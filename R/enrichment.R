################################################################################

resolve_clean_col <- function(df, col) {
  clean_col <- paste0(col, "_clean")
  if (clean_col %in% names(df)) df[[clean_col]] else strip_html(df[[col]])
}

top_n_by_score <- function(df, group_col, score_col, n) {
  df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(score = max(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::slice_head(n = as.integer(n)) %>%
    dplyr::pull(.data[[group_col]])
}

# Costruisce indice invertito: gene → vettore di pathway_id che lo contengono.
# Usato da run_enrichment per saltare pathway senza overlap (ottimizzazione).
build_gene_pathway_index <- function(pathway_genes_list) {
  gp_df <- data.frame(
    gene = unlist(pathway_genes_list, use.names = FALSE),
    pid  = rep(names(pathway_genes_list), lengths(pathway_genes_list)),
    stringsAsFactors = FALSE
  )
  split(gp_df$pid, gp_df$gene)
}

# Step 1: filtra df_mga in base alle selezioni utente
filter_mga_for_ontology <- function(df_mga, genes_selected, mirnas_selected) {
  df_mga$GENE_NAME  <- resolve_clean_col(df_mga, "GENE_NAME")
  df_mga$MIRNA_NAME <- resolve_clean_col(df_mga, "MIRNA_NAME")
  df_mga <- df_mga[!duplicated(df_mga), ]

  if (length(genes_selected) > 0)
    df_mga <- df_mga[df_mga$GENE_NAME %in% unique(unlist(genes_selected)), ]

  if (length(mirnas_selected) > 0)
    df_mga <- df_mga[df_mga$MIRNA_NAME %in% unique(unlist(mirnas_selected)), ]

  df_mga
}



# =============================================================================
# Pathway annotation (no statistics) — replaces get_mirna_pathway_enrichment
# Returns: which pathways the input genes belong to (combined or per miRNA)
# =============================================================================
annotate_mirna_pathways <- function(df_mga_filtered,
                                    mode              = c("combined", "per_mirna", "per_disease"),
                                    min_pathway_genes = 2L,
                                    max_pathway_genes = 1000L,
                                    top_n             = 50L,
                                    mirna_order       = NULL,
                                    rank_by           = c("N_FOUND", "COVERAGE", "N_MIRNAS"),
                                    pathway_ids_filter   = NULL,
                                    custom_pathway_sizes = NULL,
                                    custom_genes_list    = NULL,
                                    custom_names_map     = NULL,
                                    custom_gene_index    = NULL) {
  rank_by <- match.arg(rank_by)
  mode <- match.arg(mode)

  # Use custom pathway data if provided, otherwise fall back to globals
  .pathway_sizes     <- custom_pathway_sizes %||% pathway_sizes
  .pathway_genes_list <- custom_genes_list   %||% pathway_genes_list
  .pathway_names_map  <- custom_names_map    %||% pathway_names_map
  .gene_pathway_index <- custom_gene_index   %||% gene_pathway_index

  df_mga_filtered$GENE_NAME  <- resolve_clean_col(df_mga_filtered, "GENE_NAME")
  df_mga_filtered$MIRNA_NAME <- resolve_clean_col(df_mga_filtered, "MIRNA_NAME")
  df_mga_filtered <- df_mga_filtered[!duplicated(df_mga_filtered[, c("MIRNA_NAME", "GENE_NAME")]), ]

  valid_pids <- names(.pathway_sizes[
    .pathway_sizes >= min_pathway_genes & .pathway_sizes <= max_pathway_genes
  ])
  if (!is.null(pathway_ids_filter) && length(pathway_ids_filter) > 0) {
    allowed <- as.character(pathway_ids_filter)
    valid_pids <- valid_pids[valid_pids %in% allowed]
  }

  run_annotation <- function(genes_target) {
    candidate_pids <- unique(unlist(
      .gene_pathway_index[intersect(genes_target, names(.gene_pathway_index))],
      use.names = FALSE
    ))
    candidate_pids <- intersect(candidate_pids, valid_pids)
    if (length(candidate_pids) == 0) return(NULL)

    overlaps <- lapply(.pathway_genes_list[candidate_pids],
                       function(pg) intersect(genes_target, pg))
    k_vec <- lengths(overlaps)
    keep  <- k_vec > 0
    if (!any(keep)) return(NULL)
    candidate_pids <- candidate_pids[keep]
    overlaps       <- overlaps[keep]
    k_vec          <- k_vec[keep]
    K_vec          <- as.integer(.pathway_sizes[candidate_pids])

    data.frame(
      PATHWAY_ID   = as.integer(candidate_pids),
      PATHWAY_NAME = unname(.pathway_names_map[candidate_pids]),
      K            = K_vec,
      N_FOUND      = k_vec,
      COVERAGE     = round(k_vec / pmax(K_vec, 1L), 3),
      GENES_FOUND  = vapply(overlaps, paste, collapse = ", ", FUN.VALUE = character(1L)),
      stringsAsFactors = FALSE
    )
  }

  if (mode == "combined") {
    genes_target   <- unique(df_mga_filtered$GENE_NAME)
    gene_mirna_idx <- split(df_mga_filtered$MIRNA_NAME, df_mga_filtered$GENE_NAME)

    res_df <- run_annotation(genes_target)
    if (is.null(res_df) || nrow(res_df) == 0) return(data.frame())

    genes_list <- strsplit(res_df$GENES_FOUND, ", ", fixed = TRUE)
    res_df$MIRNAS_FOUND <- vapply(genes_list, function(genes) {
      paste(unique(unlist(gene_mirna_idx[genes], use.names = FALSE)), collapse = ", ")
    }, character(1L))
    res_df$N_MIRNAS <- vapply(res_df$MIRNAS_FOUND, function(x) {
      if (nzchar(x)) length(strsplit(x, ", ", fixed = TRUE)[[1L]]) else 0L
    }, integer(1L))

    sort_col <- if (rank_by == "N_MIRNAS") "N_MIRNAS" else if (rank_by == "COVERAGE") "COVERAGE" else "N_FOUND"
    res_df <- res_df[order(-res_df[[sort_col]]), ]
    res_df[, c("PATHWAY_ID", "PATHWAY_NAME", "K", "N_FOUND", "COVERAGE",
               "GENES_FOUND", "N_MIRNAS", "MIRNAS_FOUND")]

  } else if (mode == "per_mirna") {
    mirna_order <- if (!is.null(mirna_order)) mirna_order else unique(df_mga_filtered$MIRNA_NAME)
    sets <- split(df_mga_filtered$GENE_NAME,
                  factor(df_mga_filtered$MIRNA_NAME, levels = mirna_order))

    res_list <- lapply(names(sets), function(mir) {
      res <- run_annotation(unique(sets[[mir]]))
      if (is.null(res)) return(NULL)
      res$MIRNA_NAME <- mir
      res
    })

    res_df <- dplyr::bind_rows(res_list)
    if (is.null(res_df) || nrow(res_df) == 0) return(data.frame())

    res_df[, c("PATHWAY_ID", "MIRNA_NAME", "PATHWAY_NAME", "K", "N_FOUND", "COVERAGE", "GENES_FOUND")]

  } else {
    # per_disease: requires DISEASE_GROUP column in df_mga_filtered
    if (!"DISEASE_GROUP" %in% names(df_mga_filtered))
      stop("per_disease mode requires DISEASE_GROUP column in df_mga_filtered")

    groups   <- unique(df_mga_filtered$DISEASE_GROUP)
    sort_col <- if (rank_by == "COVERAGE") "COVERAGE" else "N_FOUND"

    res_list <- lapply(groups, function(grp) {
      sub <- df_mga_filtered[df_mga_filtered$DISEASE_GROUP == grp, ]
      res <- run_annotation(unique(sub$GENE_NAME))
      if (is.null(res) || nrow(res) == 0) return(NULL)
      res$DISEASE_GROUP <- grp
      res
    })

    res_df <- dplyr::bind_rows(res_list)
    if (is.null(res_df) || nrow(res_df) == 0) return(data.frame())

    top_paths <- top_n_by_score(res_df, "PATHWAY_NAME", sort_col, top_n)

    res_df <- res_df[res_df$PATHWAY_NAME %in% top_paths, ]
    res_df[, c("PATHWAY_ID", "DISEASE_GROUP", "PATHWAY_NAME", "K", "N_FOUND", "COVERAGE", "GENES_FOUND")]
  }
}
