# =============================================================================
# robustness.R – Network robustness sweep for miRNA–Disease / miRNA–Gene
# =============================================================================

# =============================================================================
# robustness.R – Network robustness sweep for miRNA–Disease
# =============================================================================

#' Sweep evidence thresholds and compute network metrics
#'
#' Single DB query (pair-level pubmed counts), then pure-R sweep.
#' Returns data.frame with one row per threshold value.
#'
#' @param con  DBI connection
#' @param mirna_clause   SQL fragment string (already includes leading AND)
#' @param disease_clause SQL fragment string (already includes leading AND)
#' @param min_thresh integer, minimum threshold to sweep
#' @param max_thresh integer, maximum threshold to sweep
#' @param use_update logical – use updated table?
#'
#' @return data.frame: threshold, n_edges, n_mirna, n_disease, density, jaccard
compute_mda_robustness <- function(con,
                                   mirna_clause,
                                   disease_clause,
                                   min_thresh         = 1L,
                                   max_thresh         = 20L,
                                   use_update         = FALSE,
                                   disease_restricted = FALSE,
                                   human_clause       = "") {

  min_thresh <- max(1L, as.integer(min_thresh))
  max_thresh <- max(min_thresh, as.integer(max_thresh))

  table_name <- if (isTRUE(use_update)) "MIRNAS_DISEASES_ARTICLES_update"
                else                    "MIRNAS_DISEASES_ARTICLES"

  # Convert clause objects to plain character (avoids double-quoting in paste)
  mc <- as.character(mirna_clause)
  dc <- as.character(disease_clause)

  hc <- as.character(human_clause)

  # Use MIRNA_PREMATURE to collapse 3p/5p variants (mirrors main table behaviour)
  sql <- paste(
    "SELECT M.MIRNA_PREMATURE AS MIRNA_NAME, C.DISEASE,",
    "       COUNT(DISTINCT A.PUBMED_ID) AS N_PUB",
    "FROM", table_name, "A",
    "JOIN MIRNAS   M ON A.MIRNA_ID   = M.MIRNA_ID",
    "JOIN DISEASES C ON A.DISEASE_ID = C.DISEASE_ID",
    "WHERE 1=1", mc, dc, hc,
    "GROUP BY M.MIRNA_PREMATURE, C.DISEASE"
  )

  pairs <- dbGetQuery(con, sql)

  if (nrow(pairs) == 0) return(data.frame())

  thresholds <- seq(min_thresh, max_thresh)

  n_total_dis <- max(N_TOTAL_DISEASES, 2L, na.rm = TRUE)

  rows <- lapply(thresholds, function(t) {
    sub <- pairs[pairs$N_PUB >= t, ]
    n_e <- nrow(sub)
    n_m <- length(unique(sub$MIRNA_NAME))
    n_d <- length(unique(sub$DISEASE))

    # DSI per miRNA at this threshold
    dsi_mean <- NA_real_
    dsi_sd   <- NA_real_
    if (n_e > 0L && n_m > 0L) {
      if (disease_restricted && !is.null(MIRNA_GLOBAL_DIS_COUNT) &&
          nrow(MIRNA_GLOBAL_DIS_COUNT) > 0) {
        # specific disease search → use global breadth per miRNA
        gdc      <- setNames(MIRNA_GLOBAL_DIS_COUNT$N_DIS_GLOBAL, MIRNA_GLOBAL_DIS_COUNT$MIRNA_NAME)
        mir_prem <- unique(sub$MIRNA_NAME)  # already MIRNA_PREMATURE from SQL
        n_g      <- pmax(gdc[mir_prem], 1L, na.rm = TRUE)
        n_g[is.na(n_g)] <- 1L
        dsi_vals <- 1 - log2(n_g) / log2(n_total_dis)
      } else {
        n_dis_per_mirna <- tapply(sub$DISEASE, sub$MIRNA_NAME, function(x) length(unique(x)))
        dsi_vals <- 1 - log2(pmax(n_dis_per_mirna, 1L)) / log2(n_total_dis)
      }
      dsi_mean <- round(mean(dsi_vals, na.rm = TRUE), 3)
      dsi_sd   <- round(sd(dsi_vals,   na.rm = TRUE), 3)
    }

    data.frame(
      threshold = t,
      n_edges   = n_e,
      n_mirna   = n_m,
      n_disease = n_d,
      dsi_mean  = dsi_mean,
      dsi_sd    = dsi_sd,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# =============================================================================
# robustness.R – Network robustness sweep for miRNA–Gene
# =============================================================================

#' Sweep evidence thresholds for miRNA-Gene network
#'
#' @param con DBI connection
#' @param mirna_clause  SQL fragment (leading AND)
#' @param gene_clause   SQL fragment (leading AND)
#' @param proof_clause  SQL fragment (leading AND) – active evidence filter
#' @param min_thresh integer
#' @param max_thresh integer
#' @return data.frame: threshold, n_edges, n_mirna, n_gene, density
compute_mga_robustness <- function(con,
                                   mirna_clause,
                                   gene_clause,
                                   proof_clause  = "",
                                   human_clause  = "",
                                   min_thresh    = 1L,
                                   max_thresh    = 20L) {
  min_thresh <- max(1L, as.integer(min_thresh))
  max_thresh <- max(min_thresh, as.integer(max_thresh))

  mc <- as.character(mirna_clause)
  gc <- as.character(gene_clause)
  pc <- as.character(proof_clause)
  hc <- as.character(human_clause)

  proof_join <- if (nchar(trimws(pc)) > 0) "JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID" else ""
  sql <- paste(
    "SELECT M.MIRNA_NAME, G.GENE_NAME,",
    "       COUNT(DISTINCT A.PUBMED_ID) AS N_PUB",
    "FROM MIRNAS_GENES_ARTICLES A",
    "JOIN MIRNAS M ON A.MIRNA_ID = M.MIRNA_ID",
    "JOIN GENES  G ON A.GENE_ID  = G.GENE_ID",
    proof_join,
    "WHERE 1=1", mc, gc, pc, hc,
    "GROUP BY M.MIRNA_NAME, G.GENE_NAME"
  )

  pairs <- dbGetQuery(con, sql)
  if (nrow(pairs) == 0) return(data.frame())

  thresholds <- seq(min_thresh, max_thresh)

  n_total_genes <- max(N_TOTAL_GENES, 2L, na.rm = TRUE)

  rows <- lapply(thresholds, function(t) {
    sub <- pairs[pairs$N_PUB >= t, ]
    n_e <- nrow(sub)
    n_m <- length(unique(sub$MIRNA_NAME))
    n_g <- length(unique(sub$GENE_NAME))

    gsi_mean <- NA_real_
    gsi_sd   <- NA_real_
    if (n_e > 0L && n_m > 0L) {
      n_genes_per_mirna <- tapply(sub$GENE_NAME, sub$MIRNA_NAME, function(x) length(unique(x)))
      gsi_vals <- 1 - log2(pmax(n_genes_per_mirna, 1L)) / log2(n_total_genes)
      gsi_mean <- round(mean(gsi_vals, na.rm = TRUE), 3)
      gsi_sd   <- round(sd(gsi_vals,   na.rm = TRUE), 3)
    }

    data.frame(threshold = t, n_edges = n_e, n_mirna = n_m,
               n_gene = n_g, gsi_mean = gsi_mean, gsi_sd = gsi_sd,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Sweep threshold × evidence category for miRNA-Gene
#'
#' For each threshold AND each category of a chosen level (STRENGTH/TYPE/THROUGHPUT):
#' N pairs with >= 1 proof in that category AND N_pub >= threshold.
#'
#' @param con DBI connection
#' @param mirna_clause SQL fragment
#' @param gene_clause  SQL fragment
#' @param level  one of "STRENGTH", "TYPE", "THROUGHPUT"
#' @param min_thresh integer
#' @param max_thresh integer
#' @return data.frame: threshold, category, n_pairs
compute_mga_breakdown <- function(con,
                                  mirna_clause,
                                  gene_clause,
                                  proof_clause = "",
                                  level      = "TYPE",
                                  min_thresh = 1L,
                                  max_thresh = 20L) {
  min_thresh <- max(1L, as.integer(min_thresh))
  max_thresh <- max(min_thresh, as.integer(max_thresh))
  level      <- match.arg(level, c("STRENGTH", "TYPE", "THROUGHPUT"))

  mc <- as.character(mirna_clause)
  gc <- as.character(gene_clause)
  pc <- as.character(proof_clause)

  # Get all pairs with their pub count AND the level category
  col_ref <- paste0("D.", level)
  sql <- paste(
    "SELECT M.MIRNA_NAME, G.GENE_NAME,", col_ref, "AS CATEGORY,",
    "       COUNT(DISTINCT A.PUBMED_ID) AS N_PUB",
    "FROM MIRNAS_GENES_ARTICLES A",
    "JOIN MIRNAS M ON A.MIRNA_ID = M.MIRNA_ID",
    "JOIN GENES  G ON A.GENE_ID  = G.GENE_ID",
    "JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID",
    "WHERE 1=1", mc, gc, pc,
    "GROUP BY M.MIRNA_NAME, G.GENE_NAME,", col_ref
  )

  raw <- dbGetQuery(con, sql)
  if (nrow(raw) == 0) return(data.frame())

  # Label STRENGTH numeric → text
  if (level == "STRENGTH") {
    raw$CATEGORY <- dplyr::recode(as.character(raw$CATEGORY),
                                   "0" = "Low (0)",
                                   "1" = "Medium (1)",
                                   "2" = "High (2)")
  }

  categories <- sort(unique(raw$CATEGORY))
  thresholds <- seq(min_thresh, max_thresh)

  rows <- lapply(thresholds, function(t) {
    sub <- raw[raw$N_PUB >= t, ]
    lapply(categories, function(cat) {
      n_pairs <- nrow(unique(sub[sub$CATEGORY == cat, c("MIRNA_NAME", "GENE_NAME")]))
      data.frame(threshold = t, category = cat, n_pairs = n_pairs,
                 stringsAsFactors = FALSE)
    })
  })

  do.call(rbind, do.call(c, rows))
}

# =============================================================================
# Co-targeting robustness sweep
# =============================================================================
#' For each evidence threshold t, filter pairs by N_PUB >= t, recompute Jaccard
#' matrix, count co-targeting pairs above jaccard_cutoff, and track stability.
#'
#' @param pairs_df  data.frame with columns col_nodes, col_shared, N_PUB
#' @param col_nodes  column used as nodes (e.g. "MIRNA_NAME")
#' @param col_shared column used as shared entities (e.g. "DISEASE" or "GENE_NAME")
#' @param min_thresh integer
#' @param max_thresh integer
#' @param jaccard_cutoff numeric [0,1] – threshold for a pair to be "co-targeting"
#' @return list: sweep = data.frame per threshold; stability = data.frame per pair
compute_cotarget_robustness <- function(pairs_df,
                                        col_nodes     = "MIRNA_NAME",
                                        col_shared    = "DISEASE",
                                        min_thresh    = 1L,
                                        max_thresh    = 20L,
                                        jaccard_cutoff = 0.1) {
  min_thresh <- max(1L, as.integer(min_thresh))
  max_thresh <- max(min_thresh, as.integer(max_thresh))
  thresholds <- seq(min_thresh, max_thresh)

  if (!all(c(col_nodes, col_shared, "N_PUB") %in% names(pairs_df)))
    return(list(sweep = data.frame(), stability = data.frame()))

  # Track per-pair stability: key = "nodeA||nodeB" (sorted)
  pair_hits <- list()

  sweep_rows <- lapply(thresholds, function(t) {
    empty_row <- data.frame(threshold = t, n_pairs = 0L,
                            mean_jaccard = NA_real_, median_jaccard = NA_real_)
    sub <- pairs_df[pairs_df$N_PUB >= t, ]
    if (nrow(sub) == 0 || length(unique(sub[[col_nodes]])) < 2)
      return(empty_row)

    mat <- compute_jaccard_matrix(sub, col_nodes = col_nodes, col_shared = col_shared)
    if (is.null(mat))
      return(empty_row)

    mat_work <- mat; diag(mat_work) <- 0
    idx <- which(mat_work >= jaccard_cutoff & upper.tri(mat_work), arr.ind = TRUE)

    if (nrow(idx) > 0) {
      keys <- apply(idx, 1, function(r) {
        paste(sort(rownames(mat)[r]), collapse = "||")
      })
      for (k in keys) pair_hits[[k]] <<- (pair_hits[[k]] %||% 0L) + 1L
    }

    jac_vals <- mat_work[mat_work > 0 & upper.tri(mat_work)]
    data.frame(
      threshold     = t,
      n_pairs       = nrow(idx),
      mean_jaccard  = if (length(jac_vals) > 0) round(mean(jac_vals),   4) else NA_real_,
      median_jaccard = if (length(jac_vals) > 0) round(median(jac_vals), 4) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  sweep_df <- do.call(rbind, sweep_rows)

  # Baseline Jaccard + shared count (at min_thresh) for each pair
  baseline_stats <- local({
    sub <- pairs_df[pairs_df$N_PUB >= min_thresh, ]
    mat <- if (nrow(sub) >= 2 && length(unique(sub[[col_nodes]])) >= 2)
      compute_jaccard_matrix(sub, col_nodes = col_nodes, col_shared = col_shared)
    else NULL
    empty <- list(jaccard = setNames(numeric(0), character(0)),
                  n_shared = setNames(integer(0), character(0)))
    if (is.null(mat)) return(empty)
    mat_w <- mat; diag(mat_w) <- 0
    idx <- which(upper.tri(mat_w) & mat_w > 0, arr.ind = TRUE)
    if (nrow(idx) == 0) return(empty)
    keys  <- apply(idx, 1, function(r) paste(sort(rownames(mat)[r]), collapse = "||"))
    jac   <- setNames(round(mat_w[idx], 4), keys)
    sets  <- lapply(split(sub[[col_shared]], sub[[col_nodes]]), unique)
    ns    <- setNames(mapply(function(a, b) {
      sa <- sets[[a]]; sb <- sets[[b]]
      if (is.null(sa) || is.null(sb)) NA_integer_ else length(intersect(sa, sb))
    }, rownames(mat)[idx[,1]], rownames(mat)[idx[,2]], SIMPLIFY = TRUE), keys)
    list(jaccard = jac, n_shared = setNames(as.integer(ns), keys))
  })
  baseline_jaccard <- baseline_stats$jaccard
  baseline_n_shared <- baseline_stats$n_shared

  # Stability table
  n_thresh <- length(thresholds)
  stability_df <- if (length(pair_hits) > 0) {
    keys  <- names(pair_hits)
    hits  <- unlist(pair_hits)
    nodes <- strsplit(keys, "\\|\\|")
    data.frame(
      node_a           = vapply(nodes, `[[`, character(1), 1),
      node_b           = vapply(nodes, `[[`, character(1), 2),
      jaccard_baseline = baseline_jaccard[keys],
      n_shared_baseline = baseline_n_shared[keys],
      n_thresholds     = as.integer(hits),
      stability_score  = round(hits / n_thresh, 3),
      stringsAsFactors = FALSE,
      row.names = NULL
    ) |> dplyr::arrange(dplyr::desc(stability_score))
  } else {
    data.frame(node_a = character(0), node_b = character(0),
               jaccard_baseline = numeric(0),
               n_shared_baseline = integer(0),
               n_thresholds = integer(0), stability_score = numeric(0))
  }

  list(sweep = sweep_df, stability = stability_df)
}
