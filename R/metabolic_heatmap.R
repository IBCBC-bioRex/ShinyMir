# R/metabolic_heatmap.R
# miRNA x subsystem impact matrix for metabolic heatmap

# -----------------------------------------------------------------------------
# Compute miRNA × subsystem matrix.
# df     : filtered_mr() result (must have MIRNA_NAME, NAME, SUBSYSTEM, ESS_FRAC)
# con    : DB connection (needed for total reactions per subsystem in coverage metric)
# metric : "ess_frac" | "coverage" | "n_reactions" | "fold_enrichment"
# Returns: list(mat, pmat, adjpmat) — pmat/adjpmat are NULL unless metric=="fold_enrichment"
# -----------------------------------------------------------------------------
compute_metabolic_matrix <- function(df, con, metric = "ess_frac", group_by = "subsystem") {

  if ("REACTION_NAME" %in% names(df) && !"NAME" %in% names(df))
    df <- dplyr::rename(df, NAME = REACTION_NAME)

  mirna_col <- if ("MIRNA_NAME_clean" %in% names(df)) "MIRNA_NAME_clean" else "MIRNA_NAME"
  df$mirna_ <- strip_html(df[[mirna_col]])

  group_col <- if (group_by == "reaction") "NAME" else if (group_by == "category") "CATEGORY" else "SUBSYSTEM"

  if (metric == "ess_frac" && !"ESS_FRAC" %in% names(df))
    df$ESS_FRAC <- NA_real_

  # Drop rows with missing group column
  df <- df[!is.na(df[[group_col]]) & nzchar(df[[group_col]]), ]
  if (nrow(df) == 0 || length(unique(df$mirna_)) < 1) return(NULL)

  fe_data <- NULL  # populated only for fold_enrichment

  agg <- switch(metric,

    "ess_frac" = {
      df %>%
        dplyr::group_by(mirna_, .data[[group_col]]) %>%
        dplyr::summarise(value = mean(ESS_FRAC, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(value = ifelse(is.nan(value), 0, value))
    },

    "coverage" = {
      if (group_by %in% c("subsystem", "category")) {
        subsystems_present <- unique(df$SUBSYSTEM)
        tot_sub <- dbGetQuery(con, glue_sql(
          "SELECT SUBSYSTEM, COUNT(DISTINCT REACTION_ID) AS total
             FROM {`REACTIONS_TBL`}
            WHERE SUBSYSTEM IN ({subs*})
            GROUP BY SUBSYSTEM",
          subs = subsystems_present, .con = con
        ))
        if (group_by == "category") {
          # Aggregate subsystem totals into category totals
          tot_sub$CATEGORY <- SUBSYSTEM_CATEGORY_MAP[tot_sub$SUBSYSTEM]
          tot <- tot_sub %>%
            dplyr::filter(!is.na(CATEGORY)) %>%
            dplyr::group_by(CATEGORY) %>%
            dplyr::summarise(total = sum(total, na.rm = TRUE), .groups = "drop")
          df %>%
            dplyr::group_by(mirna_, CATEGORY) %>%
            dplyr::summarise(impacted = dplyr::n_distinct(NAME), .groups = "drop") %>%
            dplyr::left_join(tot, by = "CATEGORY") %>%
            dplyr::mutate(value = impacted / total) %>%
            dplyr::select(mirna_, CATEGORY, value)
        } else {
          df %>%
            dplyr::group_by(mirna_, SUBSYSTEM) %>%
            dplyr::summarise(impacted = dplyr::n_distinct(NAME), .groups = "drop") %>%
            dplyr::left_join(tot_sub, by = "SUBSYSTEM") %>%
            dplyr::mutate(value = impacted / total) %>%
            dplyr::select(mirna_, SUBSYSTEM, value)
        }
      } else {
        reactions_present <- unique(df$NAME)
        tot <- dbGetQuery(con, glue_sql(
          "SELECT R.NAME AS NAME, COUNT(DISTINCT RG.GENE_ID) AS total
             FROM {`REACTIONS_TBL`} R
             JOIN {`REACTIONS_GENES_TBL`} RG ON R.REACTION_ID = RG.REACTION_ID
            WHERE R.NAME IN ({rxns*})
            GROUP BY R.NAME",
          rxns = reactions_present, .con = con
        ))
        df %>%
          dplyr::group_by(mirna_, NAME) %>%
          dplyr::summarise(impacted = dplyr::n_distinct(GENE_NAME), .groups = "drop") %>%
          dplyr::left_join(tot, by = "NAME") %>%
          dplyr::mutate(value = impacted / total) %>%
          dplyr::select(mirna_, NAME, value)
      }
    },

    "n_reactions" = {
      gcol <- if (group_by == "subsystem") "SUBSYSTEM" else "NAME"
      vcol <- if (group_by == "subsystem") "NAME" else "GENE_NAME"
      df %>%
        dplyr::group_by(mirna_, .data[[gcol]]) %>%
        dplyr::summarise(value = dplyr::n_distinct(.data[[vcol]]), .groups = "drop")
    },

    "fold_enrichment" = {
      # Hypergeometric per miRNA × group: universe = all distinct reactions in df
      N          <- length(unique(df$NAME))
      mirna_list <- unique(df$mirna_)

      group_sizes <- df %>%
        dplyr::group_by(.data[[group_col]]) %>%
        dplyr::summarise(K = dplyr::n_distinct(NAME), .groups = "drop")

      rows <- lapply(mirna_list, function(m) {
        sub_df <- df[df$mirna_ == m, ]
        n      <- length(unique(sub_df$NAME))
        if (n == 0) return(NULL)
        grp_counts <- sub_df %>%
          dplyr::group_by(.data[[group_col]]) %>%
          dplyr::summarise(k = dplyr::n_distinct(NAME), .groups = "drop") %>%
          dplyr::left_join(group_sizes, by = group_col) %>%
          dplyr::mutate(
            fe  = (k / n) / (K / N),
            pv  = phyper(k - 1L, K, N - K, n, lower.tail = FALSE),
            mirna_ = m
          )
        grp_counts
      })
      fe_data <- dplyr::bind_rows(rows) %>%
        dplyr::group_by(mirna_) %>%
        dplyr::mutate(adj_pv = p.adjust(pv, method = "BH")) %>%
        dplyr::ungroup()

      fe_data %>% dplyr::select(mirna_, dplyr::all_of(group_col), value = fe)
    }
  )

  # Pivot wide: rows = miRNA, cols = group_col
  wide <- tidyr::pivot_wider(agg,
    names_from  = dplyr::all_of(group_col),
    values_from = value,
    values_fill = 0
  )

  mirnas <- wide$mirna_
  mat    <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- mirnas

  # Build p-value matrices for fold_enrichment
  pmat    <- NULL
  adjpmat <- NULL
  if (metric == "fold_enrichment" && !is.null(fe_data)) {
    pivot_pmat <- function(col) {
      wide_p <- tidyr::pivot_wider(
        fe_data %>% dplyr::select(mirna_, dplyr::all_of(group_col), val = dplyr::all_of(col)),
        names_from  = dplyr::all_of(group_col),
        values_from = val,
        values_fill = NA_real_
      )
      pm <- as.matrix(wide_p[, -1, drop = FALSE])
      rownames(pm) <- wide_p$mirna_
      # align to mat dims
      common_rows <- intersect(rownames(mat), rownames(pm))
      common_cols <- intersect(colnames(mat), colnames(pm))
      pm[common_rows, common_cols, drop = FALSE]
    }
    pmat    <- pivot_pmat("pv")
    adjpmat <- pivot_pmat("adj_pv")
  }

  list(mat = mat, pmat = pmat, adjpmat = adjpmat)
}
