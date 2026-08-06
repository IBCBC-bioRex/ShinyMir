# =============================================================================
# Reactome hierarchy helper — get all descendant PATHWAY_IDs (recursive)
# =============================================================================
get_descendant_pathway_ids <- function(con, pathway_ids) {
  if (length(pathway_ids) == 0) return(integer(0))
  res <- tryCatch(
    dbGetQuery(con, glue_sql("
      WITH RECURSIVE desc_tree(id) AS (
        SELECT PATHWAY_ID FROM ONTOLOGY_PATHWAYS WHERE PATHWAY_ID IN ({ids*})
        UNION ALL
        SELECT PL.PATHWAY_ID
        FROM PATHWAY_LINK PL
        INNER JOIN desc_tree d ON PL.ANCESTOR = d.id
      )
      SELECT DISTINCT id FROM desc_tree",
      ids = as.integer(pathway_ids), .con = con
    ))$id,
    error = function(e) integer(0)
  )
  as.integer(res)
}

# =============================================================================
# Hierarchy helper — for a set of pathway_ids, return their L1 ancestor name.
# L1 = top-level nodes: appear only as ANCESTOR in PATHWAY_LINK, never as PATHWAY_ID.
# Returns named character vector: PATHWAY_ID (as char) → L1_NAME.
# =============================================================================
get_l1_ancestor_map <- function(con, pathway_ids) {
  if (length(pathway_ids) == 0) return(character(0))
  ids <- as.integer(pathway_ids)
  tryCatch(
    dbGetQuery(con, glue_sql("
      WITH RECURSIVE anc_tree(leaf, node) AS (
        SELECT PATHWAY_ID AS leaf, PATHWAY_ID AS node
          FROM ONTOLOGY_PATHWAYS WHERE PATHWAY_ID IN ({ids*})
        UNION ALL
        SELECT a.leaf, PL.ANCESTOR AS node
          FROM anc_tree a
          JOIN PATHWAY_LINK PL ON PL.PATHWAY_ID = a.node
      )
      SELECT a.leaf AS PATHWAY_ID, P.PATHWAY_NAME AS L1_NAME
        FROM anc_tree a
        JOIN ONTOLOGY_PATHWAYS P ON P.PATHWAY_ID = a.node
       WHERE a.node NOT IN (SELECT PATHWAY_ID FROM PATHWAY_LINK)",
      ids = ids, .con = con
    )) %>%
      dplyr::distinct(PATHWAY_ID, .keep_all = TRUE) %>%
      { setNames(.$L1_NAME, as.character(.$PATHWAY_ID)) },
    error = function(e) character(0)
  )
}

# =============================================================================
# Helper: build CELL_LINE_ID list from expression cell-line filters
# =============================================================================
build_cell_line_ids <- function(input, suffix, con) {
  cancer_filter <- input[[paste0("expression_cancer_filter",   suffix)]]
  primary_dis   <- input[[paste0("expression_primary_disease", suffix)]]
  disease_sub   <- input[[paste0("expression_disease_subtype", suffix)]]
  specific_cl   <- input[[paste0("expression_cell_lines",      suffix)]]

  q <- "SELECT CELL_LINE_ID FROM CELL_LINE WHERE 1=1"
  if (!is.null(cancer_filter) && cancer_filter != "all") {
    cancer_db <- if (cancer_filter == "yes") "Yes" else "No"
    q <- paste0(q, as.character(glue_sql(" AND CANCER = {cv}", cv = cancer_db, .con = con)))
  }
  if (length(primary_dis) > 0)
    q <- paste0(q, " AND ", as.character(glue_sql("PRIMARY_DISEASE IN ({pd*})", pd = primary_dis, .con = con)))
  if (length(disease_sub) > 0)
    q <- paste0(q, " AND ", as.character(glue_sql("DISEASE_SUBTYPE IN ({ds*})", ds = disease_sub, .con = con)))
  if (length(specific_cl) > 0)
    q <- paste0(q, " AND ", as.character(glue_sql("CELL_LINE_NAME IN ({cl*})", cl = specific_cl, .con = con)))

  dbGetQuery(con, q)$CELL_LINE_ID
}

expand_categories_to_subsystems <- function(cats) {
  names(SUBSYSTEM_CATEGORY_MAP)[!is.na(SUBSYSTEM_CATEGORY_MAP) & SUBSYSTEM_CATEGORY_MAP %in% cats]
}

# =============================================================================
# Helper: post-processing shared by get_filtered_mr and get_filtered_mr_by_reaction
# Adds SUBSYSTEM/GPR/HUMAN_ID/FORMULA, runs essentiality, applies
# subsystem exclusion, and appends ESSENTIALITY column when grouping by reaction.
# =============================================================================
apply_mr_postprocessing <- function(df, con, input, for_plot) {
  min_ess_frac  <- as.numeric(input$min_ess_frac_mr  %||% 0)
  only_essential <- isTRUE(input$only_essential_mr)
  # Flat (ungrouped) data needed for per-(miRNA,reaction) metrics
  is_flat <- for_plot || identical(input$group_mr %||% "None", "None")

  needs_gpr <- for_plot ||
    ("miRNA and reaction" %in% (input$group_mr %||% "None")) ||
    (is_flat && (min_ess_frac > 0 || only_essential))

  if (needs_gpr && nrow(df) > 0) {
    # The base SQL query (query_mr) already returns GPR, HUMAN_ID, SUBSYSTEM, FORMULA
    # per REACTION_ID. After dedup by RXN_KEY (=HUMAN_ID) in get_filtered_mr, each
    # row carries its own correct GPR — no re-fetch needed.
    # Only fetch columns that are genuinely absent (e.g. when called from a code path
    # that uses a different base query or an older df without these columns).
    cols_wanted <- c("SUBSYSTEM", "GPR", "HUMAN_ID", "FORMULA")
    cols_fetch  <- cols_wanted[!cols_wanted %in% colnames(df)]
    if (length(cols_fetch) > 0) {
      # Use RXN_KEY (HUMAN_ID) as the join key when available — it is unique per
      # reaction variant so no MIN() aggregation is needed.
      # Fall back to NAME-based MIN() for rows without a HUMAN_ID.
      has_rxn_key <- "RXN_KEY" %in% colnames(df)
      if (has_rxn_key && !all(is.na(df$RXN_KEY) | df$RXN_KEY == "")) {
        lookup_cols <- paste(
          c("COALESCE(HUMAN_ID, NAME) AS RXN_KEY",
            cols_fetch),
          collapse = ", "
        )
        reaction_lookup <- dbGetQuery(con,
          glue_sql(paste0("SELECT ", lookup_cols, " FROM ", REACTIONS_TBL,
                          " WHERE COALESCE(HUMAN_ID, NAME) IN ({vals*})"),
                   vals = unique(df$RXN_KEY), .con = con)
        )
        cols_to_drop <- setdiff(intersect(colnames(df), colnames(reaction_lookup)), "RXN_KEY")
        if (length(cols_to_drop) > 0)
          df <- df[, setdiff(colnames(df), cols_to_drop), drop = FALSE]
        df <- dplyr::left_join(df, reaction_lookup, by = "RXN_KEY")
      } else {
        lookup_cols <- paste(
          c("NAME AS NAME",
            sapply(cols_fetch, function(col) paste0("MIN(", col, ") AS ", col))),
          collapse = ", "
        )
        reaction_lookup <- dbGetQuery(con,
          glue_sql(paste0("SELECT ", lookup_cols, " FROM ", REACTIONS_TBL,
                          " WHERE NAME IN ({vals*}) GROUP BY NAME"),
                   vals = unique(df$NAME), .con = con)
        )
        cols_to_drop <- setdiff(intersect(colnames(df), colnames(reaction_lookup)), "NAME")
        if (length(cols_to_drop) > 0)
          df <- df[, setdiff(colnames(df), cols_to_drop), drop = FALSE]
        df <- dplyr::left_join(df, reaction_lookup, by = "NAME")
      }
    }

    # RXN_KEY may not exist if df came from a non-MR code path; derive it now.
    if (!"RXN_KEY" %in% colnames(df) && "HUMAN_ID" %in% colnames(df))
      df$RXN_KEY <- ifelse(is.na(df$HUMAN_ID) | df$HUMAN_ID == "", df$NAME, df$HUMAN_ID)

    # ESS_FRAC: compute when GPR is present and data is flat
    if (is_flat && all(c("MIRNA_NAME", "NAME") %in% colnames(df)) &&
        !"ESS_FRAC" %in% colnames(df)) {
      df <- compute_essentiality_fraction(con, df)
    }

    # IS_ESSENTIAL: Boolean GPR evaluation per (miRNA, reaction) on flat data
    if (is_flat && all(c("MIRNA_NAME", "NAME", "GPR", "GENE_NAME") %in% colnames(df)) &&
        !"IS_ESSENTIAL" %in% colnames(df)) {
      df <- compute_is_essential(df)
    }

    # Apply ESS_FRAC filter
    if (min_ess_frac > 0 && "ESS_FRAC" %in% colnames(df) && nrow(df) > 0)
      df <- df[!is.na(df$ESS_FRAC) & df$ESS_FRAC >= min_ess_frac, ]

    # Apply IS_ESSENTIAL filter
    if (only_essential && "IS_ESSENTIAL" %in% colnames(df) && nrow(df) > 0)
      df <- df[!is.na(df$IS_ESSENTIAL) & df$IS_ESSENTIAL, ]
  }

  # L1/L2 pathway filter (include-only; hierarchical: L1 expands to subsystems, intersected with L2)
  l1_cats  <- input$mr_filter_l1 %||% character(0)
  l2_subs  <- input$mr_filter_l2 %||% character(0)

  included_subs <- character(0)
  if (length(l1_cats) > 0)
    included_subs <- expand_categories_to_subsystems(l1_cats)
  if (length(l2_subs) > 0)
    included_subs <- if (length(included_subs) > 0) intersect(included_subs, l2_subs) else l2_subs

  needs_sub_filter <- length(included_subs) > 0 && nrow(df) > 0 && "NAME" %in% colnames(df)
  if (needs_sub_filter) {
    if (!"SUBSYSTEM" %in% colnames(df)) {
      sub_lookup <- dbGetQuery(con, glue_sql(
        paste0("SELECT NAME AS NAME, SUBSYSTEM FROM ", REACTIONS_TBL, " WHERE NAME IN ({vals*})"),
        vals = unique(df$NAME), .con = con))
      df <- dplyr::left_join(df, sub_lookup, by = "NAME")
    }
    df <- df[!is.na(df$SUBSYSTEM) & df$SUBSYSTEM %in% included_subs, ]
  }

  # Exclusion filter: remove specified categories/subsystems
  excl_raw <- input$mr_filter_exclude %||% character(0)
  if (length(excl_raw) > 0 && nrow(df) > 0) {
    all_cats    <- unique(na.omit(SUBSYSTEM_CATEGORY_MAP))
    excl_cats   <- excl_raw[excl_raw %in% all_cats]
    excl_subs   <- excl_raw[!excl_raw %in% all_cats]
    if (length(excl_cats) > 0)
      excl_subs <- c(excl_subs, expand_categories_to_subsystems(excl_cats))
    if (length(excl_subs) > 0) {
      if (!"SUBSYSTEM" %in% colnames(df)) {
        sub_lookup <- dbGetQuery(con, glue_sql(
          paste0("SELECT NAME AS NAME, SUBSYSTEM FROM ", REACTIONS_TBL, " WHERE NAME IN ({vals*})"),
          vals = unique(df$NAME), .con = con))
        df <- dplyr::left_join(df, sub_lookup, by = "NAME")
      }
      df <- filter_reactions_by_subsystem_exclusion(df, excl_subs)
    }
  }

  if ("miRNA and reaction" %in% (input$group_mr %||% "None") && nrow(df) > 0) {
    df$ESSENTIALITY <- if ("IS_ESSENTIAL" %in% colnames(df)) {
      dplyr::case_when(
        is.na(df$IS_ESSENTIAL)  ~ "unknown",
        df$IS_ESSENTIAL         ~ "yes",
        TRUE                    ~ "no"
      )
    } else {
      ifelse(grepl("\\bor\\b", df$GPR, ignore.case = TRUE), "no", "yes")
    }
    df <- df %>% relocate(ESSENTIALITY, .after = 2)
  }
  df
}

# =============================================================================
apply_pathway_filter_by_ids <- function(df, pathway_ids, con) {
  if (length(pathway_ids) == 0 || !"GENE_NAME" %in% names(df) || nrow(df) == 0) return(df)
  pathway_genes <- dbGetQuery(con, glue_sql(
    "SELECT DISTINCT G.GENE_NAME
       FROM ONTOLOGY_PATHWAYS_GENES PG
       JOIN GENES G ON PG.GENE_ID = G.GENE_ID
      WHERE PG.PATHWAY_ID IN ({ids*})",
    ids = as.integer(pathway_ids), .con = con
  ))$GENE_NAME
  if (length(pathway_genes) == 0) df[0, , drop = FALSE]
  else df[df$GENE_NAME %in% pathway_genes, , drop = FALSE]
}

# =============================================================================
# Expression filter helpers
# =============================================================================
get_expression_gene_names <- function(input, suffix, con) {
  mode <- input[[paste0("expression_filter_mode", suffix)]]
  if (is.null(mode) || mode == "none") return(NULL)

  if (mode == "tissue") {
    tissues <- input[[paste0("expression_tissue_select", suffix)]]
    nTPM    <- input[[paste0("nTPM_expression_filter",   suffix)]]
    # Require at least one tissue selected – nTPM alone without tissue would
    # return genes expressed in ANY tissue, which is not a meaningful filter.
    if (length(tissues) == 0) return(NULL)

    where_parts <- "1=1"
    if (length(tissues) > 0)
      where_parts <- paste0(where_parts, " AND ",
                            as.character(glue_sql("t.TISSUE_NAME IN ({tissues*})", tissues = tissues, .con = con)))
    if (!is.null(nTPM) && nTPM > 0)
      where_parts <- paste0(where_parts,
                            as.character(glue_sql(" AND h.nTPM >= {nTPM*}", nTPM = nTPM, .con = con)))

    return(dbGetQuery(con, paste0(
      "SELECT DISTINCT g.GENE_NAME
         FROM HPA_TISSUE_DATA h
         JOIN TISSUE_TYPE t ON h.TISSUE_ID = t.TISSUE_ID
         JOIN GENES g ON h.GENE_ID = g.GENE_ID
        WHERE ", where_parts
    ))$GENE_NAME)
  }

  # cell_line mode
  nTPM   <- input[[paste0("nTPM_2_expression_filter", suffix)]]
  cl_ids <- build_cell_line_ids(input, suffix, con)
  if (length(cl_ids) == 0) return(character(0))

  base_q <- paste0(
    "SELECT DISTINCT g.GENE_NAME
       FROM HPA_CELLLINE_DATA h
       JOIN GENES g ON h.GENE_ID = g.GENE_ID
      WHERE ",
    as.character(glue_sql("h.CELL_LINE_ID IN ({cl_ids*})", cl_ids = cl_ids, .con = con))
  )
  if (!is.null(nTPM) && nTPM > 0)
    base_q <- paste0(base_q,
                     as.character(glue_sql(" AND h.nTPM >= {nTPM*}", nTPM = nTPM, .con = con)))
  dbGetQuery(con, base_q)$GENE_NAME
}

# Filter pathway IDs by tissue/cell-line expression.
# Uses CTE on GENE_ID (integer, indexed) — avoids large IN(gene_names) clause.
# input/suffix: same as get_expression_gene_names; min_frac in [0,1].
# Returns integer vector of PATHWAY_IDs, integer(0) if nothing passes, NULL if filter inactive.
filter_pathways_by_expression <- function(con, input, suffix, min_frac = 0.2) {
  mode <- input[[paste0("expression_filter_mode", suffix)]] %||% "none"
  if (mode == "none") return(NULL)

  nTPM_thresh <- if (mode == "tissue")
    as.numeric(input[[paste0("nTPM_expression_filter",   suffix)]] %||% 0.1)
  else
    as.numeric(input[[paste0("nTPM_2_expression_filter", suffix)]] %||% 0.1)

  expr_subquery <- if (mode == "tissue") {
    tissues <- input[[paste0("expression_tissue_select", suffix)]]
    if (length(tissues) == 0) return(NULL)
    as.character(glue_sql(
      "SELECT DISTINCT h.GENE_ID FROM HPA_TISSUE_DATA h
         JOIN TISSUE_TYPE t ON h.TISSUE_ID = t.TISSUE_ID
        WHERE t.TISSUE_NAME IN ({tissues*}) AND h.nTPM >= {nTPM_thresh}",
      tissues = tissues, nTPM_thresh = nTPM_thresh, .con = con))
  } else {
    cl_ids <- build_cell_line_ids(input, suffix, con)
    if (length(cl_ids) == 0) return(NULL)
    as.character(glue_sql(
      "SELECT DISTINCT h.GENE_ID FROM HPA_CELLLINE_DATA h
        WHERE h.CELL_LINE_ID IN ({cl_ids*}) AND h.nTPM >= {nTPM_thresh}",
      cl_ids = cl_ids, nTPM_thresh = nTPM_thresh, .con = con))
  }

  q <- paste0(
    "WITH expressed AS (", expr_subquery, ")
     SELECT pg.PATHWAY_ID,
            COUNT(DISTINCT pg.GENE_ID)                                          AS total_genes,
            COUNT(DISTINCT CASE WHEN e.GENE_ID IS NOT NULL THEN pg.GENE_ID END) AS expressed_genes
       FROM ONTOLOGY_PATHWAYS_GENES pg
       LEFT JOIN expressed e ON pg.GENE_ID = e.GENE_ID
      GROUP BY pg.PATHWAY_ID
     HAVING expressed_genes >= 1"
  )
  df <- tryCatch(dbGetQuery(con, q), error = function(e) data.frame())
  if (nrow(df) == 0) return(data.frame(PATHWAY_ID = integer(0), expressed_genes = integer(0)))
  df[df$expressed_genes / df$total_genes >= min_frac,
     c("PATHWAY_ID", "expressed_genes"), drop = FALSE]
}

apply_expression_filter <- function(df, expr_genes) {
  if (is.null(expr_genes))   return(df)
  if (nrow(df) == 0)         return(df)
  if (!"GENE_NAME" %in% names(df)) return(df)
  if (length(expr_genes) == 0) return(df[0L, , drop = FALSE])

  gene_col   <- df$GENE_NAME[!is.na(df$GENE_NAME)]
  is_concat  <- length(gene_col) > 0 &&
                any(grepl(",", gene_col, fixed = TRUE))

  if (!is_concat) {
    # Individual gene names per row – standard filter
    return(df[!is.na(df$GENE_NAME) & df$GENE_NAME %in% expr_genes, , drop = FALSE])
  }

  # GROUP_CONCAT result: each cell is "GENE_A,GENE_B,..." – split and filter
  split_genes <- strsplit(df$GENE_NAME, ",", fixed = TRUE)
  expressed   <- lapply(split_genes, function(gs) {
    gs <- trimws(gs)
    gs[nzchar(gs) & gs %in% expr_genes]
  })
  keep <- lengths(expressed) > 0
  df   <- df[keep, , drop = FALSE]

  if (nrow(df) > 0) {
    expressed_kept <- expressed[keep]
    df$GENE_NAME   <- vapply(expressed_kept, paste, character(1L), collapse = ",")
    if ("GENE_COUNT" %in% names(df))
      df$GENE_COUNT <- as.integer(lengths(expressed_kept))
  }
  df
}

add_expression_annotation <- function(df, input, suffix, con) {
  mode <- input[[paste0("expression_filter_mode", suffix)]]
  if (is.null(mode) || mode == "none" || nrow(df) == 0) return(df)

  gene_names <- unique(df$GENE_NAME)

  if (mode == "tissue") {
    # Aggregate in R – avoids GROUP_CONCAT(ORDER BY) which is SQLite-only (breaks DuckDB)
    tissues <- input[[paste0("expression_tissue_select", suffix)]]
    q <- paste0(
      "SELECT g.GENE_NAME, t.TISSUE_NAME
       FROM HPA_TISSUE_DATA h
       JOIN TISSUE_TYPE t ON h.TISSUE_ID = t.TISSUE_ID
       JOIN GENES g ON h.GENE_ID = g.GENE_ID
       WHERE ", as.character(glue_sql("g.GENE_NAME IN ({gn*})", gn = gene_names, .con = con))
    )
    if (length(tissues) > 0)
      q <- paste0(q, " AND ", as.character(glue_sql("t.TISSUE_NAME IN ({ts*})", ts = tissues, .con = con)))
    ann_raw <- dbGetQuery(con, q)
    ann <- ann_raw %>%
      dplyr::group_by(GENE_NAME) %>%
      dplyr::summarise(TISSUES = paste(sort(unique(TISSUE_NAME)), collapse = ", "),
                       .groups = "drop")
    df <- dplyr::left_join(df, ann, by = "GENE_NAME")

  } else if (mode == "cell_line") {
    cl_ids <- build_cell_line_ids(input, suffix, con)
    if (length(cl_ids) > 0) {
      # Aggregate in R – avoids GROUP_CONCAT(ORDER BY) which is SQLite-only
      ann_raw <- dbGetQuery(con, glue_sql(
        "SELECT g.GENE_NAME, h.CELL_LINE_ID, c.PRIMARY_DISEASE
         FROM HPA_CELLLINE_DATA h
         JOIN CELL_LINE c ON h.CELL_LINE_ID = c.CELL_LINE_ID
         JOIN GENES g ON h.GENE_ID = g.GENE_ID
         WHERE h.CELL_LINE_ID IN ({cl_ids*})
           AND g.GENE_NAME IN ({gn*})",
        cl_ids = cl_ids, gn = gene_names, .con = con
      ))
      ann <- ann_raw %>%
        dplyr::group_by(GENE_NAME) %>%
        dplyr::summarise(
          N_CELL_LINES    = dplyr::n_distinct(CELL_LINE_ID),
          PRIMARY_DISEASES = paste(sort(unique(PRIMARY_DISEASE[!is.na(PRIMARY_DISEASE)])),
                                   collapse = ", "),
          .groups = "drop"
        )
      df <- dplyr::left_join(df, ann, by = "GENE_NAME")
    }
  }
  df
}

# =============================================================================
# Expression filter cascade observer
# =============================================================================
observeExpressionFilterCascade <- function(input, session, suffix, con) {
  disease_id   <- paste0("expression_primary_disease", suffix)
  subtype_id   <- paste0("expression_disease_subtype", suffix)
  celllines_id <- paste0("expression_cell_lines",      suffix)
  cancer_id    <- paste0("expression_cancer_filter",   suffix)
  reset_id     <- paste0("reset_expression",           suffix)

  # PRIMARY_DISEASE cascade on cancer filter change
  observeEvent(input[[cancer_id]], {
    cv <- input[[cancer_id]]
    q  <- "SELECT DISTINCT PRIMARY_DISEASE FROM CELL_LINE WHERE PRIMARY_DISEASE IS NOT NULL"
    if (!is.null(cv) && cv != "all") {
      cancer_db <- if (cv == "yes") "Yes" else "No"
      q <- paste0(q, as.character(glue_sql(" AND CANCER = {v}", v = cancer_db, .con = con)))
    }
    choices <- dbGetQuery(con, paste0(q, " ORDER BY PRIMARY_DISEASE"))$PRIMARY_DISEASE
    updateSelectizeInput(session, disease_id, choices = choices, selected = character(0), server = TRUE)
  }, ignoreNULL = FALSE)

  # DISEASE_SUBTYPE cascade on primary disease change
  observeEvent(input[[disease_id]], {
    cv  <- input[[cancer_id]]
    pds <- input[[disease_id]]
    q   <- "SELECT DISTINCT DISEASE_SUBTYPE FROM CELL_LINE WHERE DISEASE_SUBTYPE IS NOT NULL"
    if (!is.null(cv) && cv != "all") {
      cancer_db <- if (cv == "yes") "Yes" else "No"
      q <- paste0(q, as.character(glue_sql(" AND CANCER = {v}", v = cancer_db, .con = con)))
    }
    if (length(pds) > 0)
      q <- paste0(q, " AND ", as.character(glue_sql("PRIMARY_DISEASE IN ({pd*})", pd = pds, .con = con)))
    choices <- dbGetQuery(con, paste0(q, " ORDER BY DISEASE_SUBTYPE"))$DISEASE_SUBTYPE
    updateSelectizeInput(session, subtype_id, choices = choices, selected = character(0), server = TRUE)
  }, ignoreNULL = FALSE)

  # CELL_LINES cascade – debounced: cancer→disease→subtype changes fire this in rapid
  # succession; debounce collapses them into one DB query after 350ms of quiet.
  cell_line_trigger <- reactive({
    list(input[[cancer_id]], input[[disease_id]], input[[subtype_id]])
  }) |> debounce(350)

  observe({
    cell_line_trigger()
    cv  <- input[[cancer_id]]
    pds <- input[[disease_id]]
    ds  <- input[[subtype_id]]
    q   <- "SELECT DISTINCT CELL_LINE_NAME FROM CELL_LINE WHERE 1=1"
    if (!is.null(cv) && cv != "all") {
      cancer_db <- if (cv == "yes") "Yes" else "No"
      q <- paste0(q, as.character(glue_sql(" AND CANCER = {v}", v = cancer_db, .con = con)))
    }
    if (length(pds) > 0)
      q <- paste0(q, " AND ", as.character(glue_sql("PRIMARY_DISEASE IN ({pd*})", pd = pds, .con = con)))
    if (length(ds) > 0)
      q <- paste0(q, " AND ", as.character(glue_sql("DISEASE_SUBTYPE IN ({d*})", d = ds, .con = con)))
    choices <- dbGetQuery(con, paste0(q, " ORDER BY CELL_LINE_NAME"))$CELL_LINE_NAME
    updateSelectizeInput(session, celllines_id, choices = choices, selected = character(0), server = TRUE)
  })

  # Reset
  observeEvent(input[[reset_id]], {
    updateRadioButtons(session,  paste0("expression_filter_mode",    suffix), selected = "none")
    updateSelectizeInput(session, paste0("expression_tissue_select",  suffix), selected = character(0))
    updateRadioButtons(session,  paste0("expression_cancer_filter",   suffix), selected = "all")
    updateSelectizeInput(session, disease_id,   choices = cell_line_primary_diseases, selected = character(0), server = TRUE)
    updateSelectizeInput(session, subtype_id,   selected = character(0))
    updateSelectizeInput(session, celllines_id, selected = character(0))
  })
}

# =============================================================================
# Proof helpers (AND-mode)
# =============================================================================

# Returns SQL fragment " AND (MIRNA_ID, GENE_ID) IN (pair_subquery)"
and_proof_pair_sql <- function(selected_proofs, mirna_ids = integer(0), gene_ids = integer(0), con, alias = NULL) {
  n <- length(selected_proofs)
  if (n < 2) return(DBI::SQL(""))
  pfx <- if (!is.null(alias)) paste0(alias, ".") else ""
  wp <- as.character(glue_sql("D.PROOF_NAME IN ({proofs*})", proofs = selected_proofs, .con = con))
  if (length(mirna_ids) > 0)
    wp <- paste(wp, "AND", as.character(glue_sql("A.MIRNA_ID IN ({mids*})", mids = mirna_ids, .con = con)))
  if (length(gene_ids) > 0)
    wp <- paste(wp, "AND", as.character(glue_sql("A.GENE_ID IN ({gids*})", gids = gene_ids, .con = con)))
  DBI::SQL(paste0(
    " AND (", pfx, "MIRNA_ID, ", pfx, "GENE_ID) IN (",
    "SELECT A.MIRNA_ID, A.GENE_ID FROM MIRNAS_GENES_ARTICLES A",
    " JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID",
    " WHERE ", wp,
    " GROUP BY A.MIRNA_ID, A.GENE_ID",
    " HAVING COUNT(DISTINCT D.PROOF_NAME) >= ", n, ")"
  ))
}

# Returns df(MIRNA_NAME, GENE_NAME) for pairs that have ALL selected proofs.
and_proof_valid_pairs_names <- function(selected_proofs, mirna_ids = integer(0), gene_ids = integer(0), con) {
  n <- length(selected_proofs)
  if (n < 2) return(NULL)
  wp <- as.character(glue_sql("D.PROOF_NAME IN ({proofs*})", proofs = selected_proofs, .con = con))
  if (length(mirna_ids) > 0)
    wp <- paste(wp, "AND", as.character(glue_sql("A.MIRNA_ID IN ({mids*})", mids = mirna_ids, .con = con)))
  if (length(gene_ids) > 0)
    wp <- paste(wp, "AND", as.character(glue_sql("A.GENE_ID IN ({gids*})", gids = gene_ids, .con = con)))
  dbGetQuery(con, paste0(
    "SELECT DISTINCT B.MIRNA_NAME, C.GENE_NAME",
    " FROM MIRNAS_GENES_ARTICLES A",
    " JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID",
    " JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID",
    " JOIN GENES C ON A.GENE_ID = C.GENE_ID",
    " WHERE ", wp,
    " GROUP BY A.MIRNA_ID, A.GENE_ID, B.MIRNA_NAME, C.GENE_NAME",
    " HAVING COUNT(DISTINCT D.PROOF_NAME) >= ", n
  ))
}

# =============================================================================
# Reaction size / subsystem filters
# =============================================================================
filter_reactions_by_subsystem_exclusion <- function(df, excluded_subsystems) {
  if (length(excluded_subsystems) == 0 || !("SUBSYSTEM" %in% colnames(df))) return(df)
  df[is.na(df$SUBSYSTEM) | !df$SUBSYSTEM %in% excluded_subsystems, ]
}


get_mirna_ids <- function(vals, con, arm = NULL) {
  if (length(vals) == 0) return(integer(0))
  # If any value has an arm suffix (-3p / -5p), query on MIRNA_NAME; else MIRNA_PREMATURE
  if (any(grepl("-[35]p$", vals, perl = TRUE))) {
    dbGetQuery(con, glue_sql(
      "SELECT MIRNA_ID FROM MIRNAS WHERE MIRNA_NAME IN ({vals*})",
      vals = vals, .con = con
    ))$MIRNA_ID
  } else if (!is.null(arm) && arm == "max") {
    # Per-miRNA: pick the arm with more gene associations; fallback to both if tied/missing
    res <- dbGetQuery(con, glue_sql(
      "SELECT M.MIRNA_ID, M.MIRNA_PREMATURE, M.MIRNA_ARM,
              COUNT(DISTINCT MGA.GENE_ID) AS N_GENES
       FROM MIRNAS M
       LEFT JOIN MIRNAS_GENES_ARTICLES MGA ON M.MIRNA_ID = MGA.MIRNA_ID
       WHERE M.MIRNA_PREMATURE IN ({vals*}) AND M.MIRNA_ARM IN ('5p','3p')
       GROUP BY M.MIRNA_ID, M.MIRNA_PREMATURE, M.MIRNA_ARM",
      vals = vals, .con = con
    ))
    if (nrow(res) == 0) return(integer(0))
    # For each premature, keep the arm with max N_GENES (ties → keep both)
    best <- do.call(rbind, lapply(split(res, res$MIRNA_PREMATURE), function(g) {
      mx <- max(g$N_GENES, na.rm = TRUE)
      g[g$N_GENES == mx, , drop = FALSE]
    }))
    best$MIRNA_ID
  } else if (!is.null(arm) && arm %in% c("5p", "3p")) {
    dbGetQuery(con, glue_sql(
      "SELECT MIRNA_ID FROM MIRNAS WHERE MIRNA_PREMATURE IN ({vals*}) AND MIRNA_ARM = {arm}",
      vals = vals, arm = arm, .con = con
    ))$MIRNA_ID
  } else {
    dbGetQuery(con, glue_sql(
      "SELECT MIRNA_ID FROM MIRNAS WHERE MIRNA_PREMATURE IN ({vals*})",
      vals = vals, .con = con
    ))$MIRNA_ID
  }
}

# =============================================================================
# Disease helpers
# =============================================================================
select_disease_ancestors_fun <- function(disease_selected, filter_ancestor, con, depth = "all") {
  if (filter_ancestor && length(disease_selected) > 0) {
    query_ <- if (depth == "first") {
      glue_sql("
        SELECT DISTINCT D.DISEASE
        FROM DISEASES ROOT
        JOIN DISEASE_LINK DL ON ROOT.DISEASE_ID = DL.ANCESTOR
        JOIN DISEASES D      ON DL.DISEASE_ID   = D.DISEASE_ID
        WHERE ROOT.DISEASE IN ({vals*})
        UNION
        SELECT DISEASE FROM DISEASES WHERE DISEASE IN ({vals*})",
        vals = disease_selected, .con = con)
    } else {
      glue_sql("
        WITH RECURSIVE desc_tree(id) AS (
          SELECT DISEASE_ID FROM DISEASES WHERE DISEASE IN ({vals*})
          UNION ALL
          SELECT DL.DISEASE_ID
          FROM DISEASE_LINK DL
          INNER JOIN desc_tree d ON DL.ANCESTOR = d.id
        )
        SELECT DISTINCT DISEASE
        FROM DISEASES
        WHERE DISEASE_ID IN (SELECT id FROM desc_tree)",
        vals = disease_selected, .con = con)
    }
    unique(dbGetQuery(con, query_)$DISEASE)
  } else {
    disease_selected
  }
}

# For each descendant disease, return which originally-selected disease it belongs to.
get_disease_descendant_map <- function(disease_selected, con) {
  if (length(disease_selected) == 0) return(NULL)
  query_ <- glue_sql("
    WITH RECURSIVE desc_tree(id, root_disease) AS (
      SELECT DISEASE_ID, DISEASE
      FROM DISEASES
      WHERE DISEASE IN ({vals*})
      UNION ALL
      SELECT DL.DISEASE_ID, dt.root_disease
      FROM DISEASE_LINK DL
      INNER JOIN desc_tree dt ON DL.ANCESTOR = dt.id
    )
    SELECT
      child.DISEASE     AS DISEASE,
      MIN(dt.root_disease) AS ANCESTOR_GROUP
    FROM desc_tree dt
    JOIN DISEASES child ON dt.id = child.DISEASE_ID
    GROUP BY child.DISEASE",
    vals = disease_selected, .con = con
  )
  dbGetQuery(con, query_)
}

# =============================================================================
# miRNADisease filter
# =============================================================================
filtered_mda_fun <- function(input, botton_update, con, for_plot = FALSE, my_data, group_override = NULL) {
  if (for_plot) {
    base_query <- query_md[["query_mirna_disease"]]
  } else {
    group_sel  <- group_override %||% (input$group_mda %||% "None")
    base_query <- query_md[["query_general"]]
    if ("miRNA"             %in% group_sel) base_query <- query_md[["query_mirna"]]
    if ("disease"           %in% group_sel) base_query <- query_md[["query_disease"]]
    if ("miRNA and disease" %in% group_sel) base_query <- query_md[["query_mirna_disease"]]
  }
  mda_tbl <- if (botton_update()) "MIRNAS_DISEASES_ARTICLES_update" else "MIRNAS_DISEASES_ARTICLES"
  if (botton_update()) {
    base_query <- gsub("MIRNAS_DISEASES_ARTICLES", mda_tbl, base_query, fixed = TRUE)
  }

  # Bug fix 1+3: check SAMPLE_TYPE column exists before using it; also filter inner
  # min_assoc subquery by sample type so thresholds count only matching-type pubs.
  sample_types        <- input$sample_type_mda
  has_sample_type_col <- "SAMPLE_TYPE" %in% dbListFields(con, mda_tbl)

  if (length(sample_types) > 0 && !has_sample_type_col) {
    showNotification(
      "Sample type filter requires running the MDA_SAMPLE_TYPE build step. Filter ignored.",
      type = "warning", duration = 6
    )
  }

  inner_sample_filter <- if (length(sample_types) > 0 && has_sample_type_col) {
    safe_types <- gsub("'", "''", sample_types, fixed = TRUE)
    like_parts <- vapply(safe_types,
      function(s) sprintf("SAMPLE_TYPE LIKE '%%%s%%'", s), character(1))
    paste0("AND (", paste(like_parts, collapse = " OR "), ")")
  } else {
    ""
  }

  human_clause <- if (isTRUE(input$filter_human_mda) && isTRUE(HAS_IS_HUMAN))
    DBI::SQL("AND CAST(B.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)")
  else
    DBI::SQL("")

  min_assoc_val    <- as.integer(input$min_assoc_mda %||% 1L)
  # include human_clause inside the pair-filter subquery so the count matches the outer query
  human_inner_mda <- if (isTRUE(input$filter_human_mda) && isTRUE(HAS_IS_HUMAN))
    "AND CAST(PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)"
  else
    ""
  min_assoc_clause <- if (min_assoc_val > 1L) {
    DBI::SQL(sprintf(
      "JOIN (SELECT MIRNA_ID, DISEASE_ID FROM %s WHERE 1=1 %s %s GROUP BY MIRNA_ID, DISEASE_ID HAVING COUNT(DISTINCT PUBMED_ID) >= %d) mda_pair_filter ON A.MIRNA_ID = mda_pair_filter.MIRNA_ID AND C.DISEASE_ID = mda_pair_filter.DISEASE_ID",
      mda_tbl, inner_sample_filter, human_inner_mda, min_assoc_val
    ))
  } else {
    DBI::SQL("")
  }

  mirnas_id <- get_mirna_ids(input$mirna_select, con)
  mirna_clause <- if (length(input$mirna_select) > 0) {
    if (length(mirnas_id) == 0) DBI::SQL("AND 1=0")
    else glue_sql("AND A.MIRNA_ID IN ({vals*})", vals = mirnas_id, .con = con)
  } else if (input$filter_mirna_all_checkbox) {
    DBI::SQL("")
  } else {
    DBI::SQL("AND 1=0")
  }

  disease_list_total <- select_disease_ancestors_fun(
    disease_selected = input$disease_select,
    filter_ancestor  = input$filter_diseaseancestor_checkbox,
    con   = con,
    depth = input$disease_descendant_depth %||% "all"
  )

  # "All diseases" pool: respect DO category filter if active
  cat_sel        <- input$do_category_filter %||% ""
  all_in_scope   <- if (nzchar(cat_sel) && cat_sel %in% names(do_category_disease_map))
                      do_category_disease_map[[cat_sel]]
                    else
                      disease_choices

  disease_list_total_id <- if (length(disease_list_total) > 0)
    dbGetQuery(con,
      glue_sql("SELECT DISEASE_ID FROM DISEASES WHERE DISEASE IN ({vals*})",
               vals = disease_list_total, .con = con))$DISEASE_ID
  else integer(0)

  disease_clause <- if (length(disease_list_total) > 0) {
    glue_sql("AND C.DISEASE_ID IN ({vals*})", vals = disease_list_total_id, .con = con)
  } else if (input$filter_disease_checkbox) {
    # "use all" but scoped to category
    if (nzchar(cat_sel) && cat_sel %in% names(do_category_disease_map)) {
      scope_ids <- dbGetQuery(con,
        glue_sql("SELECT DISEASE_ID FROM DISEASES WHERE DISEASE IN ({vals*})",
                 vals = all_in_scope, .con = con))$DISEASE_ID
      if (length(scope_ids) > 0)
        glue_sql("AND C.DISEASE_ID IN ({vals*})", vals = scope_ids, .con = con)
      else
        DBI::SQL("AND 1=0")
    } else {
      DBI::SQL("")
    }
  } else {
    DBI::SQL("AND 1=0")
  }

  sample_type_clause <- if (length(sample_types) > 0 && has_sample_type_col) {
    like_clauses <- lapply(sample_types, function(s)
      glue_sql("B.SAMPLE_TYPE LIKE {pattern}", pattern = paste0("%", s, "%"), .con = con))
    DBI::SQL(paste0("AND (", paste(vapply(like_clauses, as.character, character(1)), collapse = " OR "), ")"))
  } else {
    DBI::SQL("")
  }

  query <- glue_sql(base_query,
                    mirna_clause       = mirna_clause,
                    disease_clause     = disease_clause,
                    min_assoc_clause   = min_assoc_clause,
                    min_assoc_val      = min_assoc_val,
                    sample_type_clause = sample_type_clause,
                    human_clause       = human_clause,
                    .con = con)
  df <- dbGetQuery(con, query)

  # df_common: flat miRNA-disease pair query (one row per pair, with PUBMED_COUNT)
  # used for degree filter, mirna_common logic, and validation
  df_common_query <- gsub("MIRNAS_DISEASES_ARTICLES", mda_tbl,
                          query_md[["query_mirna_disease"]], fixed = TRUE)
  df_common <- if (for_plot || identical(base_query, query_md[["query_mirna_disease"]]) ||
                   identical(base_query, gsub("MIRNAS_DISEASES_ARTICLES", mda_tbl,
                                              query_md[["query_mirna_disease"]], fixed = TRUE))) {
    df
  } else {
    dbGetQuery(con, glue_sql(df_common_query,
      mirna_clause       = mirna_clause,
      disease_clause     = disease_clause,
      min_assoc_clause   = min_assoc_clause,
      min_assoc_val      = min_assoc_val,
      sample_type_clause = sample_type_clause,
      human_clause       = human_clause,
      .con = con
    ))
  }

  # Degree filter
  min_m <- as.integer(input$min_degree_mirna_mda   %||% 1L)
  min_d <- as.integer(input$min_degree_disease_mda %||% 1L)
  if ((min_m > 1L || min_d > 1L) && nrow(df_common) > 0) {
    if (min_m > 1L) {
      deg_m   <- tapply(df_common$DISEASE,    df_common$MIRNA_NAME, function(x) length(unique(x)))
      valid_m <- names(deg_m)[deg_m >= min_m]
      df_common <- df_common[df_common$MIRNA_NAME %in% valid_m, ]
      if ("MIRNA_NAME" %in% names(df)) df <- df[df$MIRNA_NAME %in% valid_m, ]
    }
    if (min_d > 1L) {
      deg_d   <- tapply(df_common$MIRNA_NAME, df_common$DISEASE,    function(x) length(unique(x)))
      valid_d <- names(deg_d)[deg_d >= min_d]
      df_common <- df_common[df_common$DISEASE %in% valid_d, ]
      if ("DISEASE" %in% names(df)) df <- df[df$DISEASE %in% valid_d, ]
    }
  }

  mirnas_total <- unique(df_common$MIRNA_NAME)
  ndisease     <- length(unique(df_common$DISEASE))
  mirna_common <- df_common %>%
    dplyr::group_by(MIRNA_NAME) %>%
    dplyr::summarise(n_diseases = n_distinct(DISEASE), .groups = "drop") %>%
    dplyr::filter(n_diseases >= ndisease) %>%
    dplyr::pull(MIRNA_NAME)

  my_data$mirnas_total     <- mirnas_total
  my_data$mirna_common     <- mirna_common
  my_data$mirna_disease_df <- if (nrow(df_common) > 0)
    unique(df_common[, c("MIRNA_NAME", "DISEASE"), drop = FALSE])
  else
    data.frame(MIRNA_NAME = character(0), DISEASE = character(0), stringsAsFactors = FALSE)

  # DSI (miRNA-level) + PAIR_SCORE (pair-level, within filter)
  # When specific diseases are selected, df_common only contains those diseases →
  # n_dis would be artificially low → DSI artificially high.
  # In that case, fall back to MIRNA_GLOBAL_DIS_COUNT (full DB breadth per miRNA).
  # When "all diseases" scope, use df_common (respects min_assoc threshold).
  if (nrow(df_common) > 0 && !is.na(N_TOTAL_DISEASES) && N_TOTAL_DISEASES > 1) {
    cat_sel  <- input$do_category_filter %||% ""
    n_total  <- if (nzchar(cat_sel) && cat_sel %in% names(N_DISEASES_BY_CATEGORY))
                  max(N_DISEASES_BY_CATEGORY[[cat_sel]], 2L, na.rm = TRUE)
                else
                  N_TOTAL_DISEASES

    disease_restricted <- length(disease_list_total) > 0

    n_dis_ctx <- if (disease_restricted && !is.null(MIRNA_GLOBAL_DIS_COUNT) &&
                     nrow(MIRNA_GLOBAL_DIS_COUNT) > 0) {
      # specific disease selection → use global per-miRNA breadth
      # df_common$MIRNA_NAME = MIRNA_PREMATURE (from query), matches MIRNA_GLOBAL_DIS_COUNT keys
      gdc <- setNames(MIRNA_GLOBAL_DIS_COUNT$N_DIS_GLOBAL, MIRNA_GLOBAL_DIS_COUNT$MIRNA_NAME)
      pmax(gdc[df_common$MIRNA_NAME], 1L, na.rm = TRUE)
    } else {
      # all-disease scope → use filtered counts (respects min_assoc)
      pmax(
        tapply(df_common$DISEASE, df_common$MIRNA_NAME,
               function(x) length(unique(x)))[df_common$MIRNA_NAME],
        1L
      )
    }
    n_dis_ctx[is.na(n_dis_ctx)] <- 1L

    df_common$DSI <- round(1 - log2(n_dis_ctx) / log2(n_total), 3)
    has_disease_count <- "DISEASE_COUNT" %in% names(df)
    has_mirna_count   <- "MIRNA_COUNT"   %in% names(df)

    if (has_disease_count && !has_mirna_count) {
      mirna_dsi <- unique(df_common[, c("MIRNA_NAME", "DSI")])
      df <- dplyr::left_join(df, mirna_dsi, by = "MIRNA_NAME")
    } else if (has_mirna_count && !has_disease_count) {
      disease_dsi <- df_common %>%
        dplyr::group_by(DISEASE) %>%
        dplyr::summarise(DSI_MEAN = round(mean(DSI, na.rm = TRUE), 3), .groups = "drop")
      df <- dplyr::left_join(df, disease_dsi, by = "DISEASE")
    } else if (!has_disease_count && !has_mirna_count && all(c("MIRNA_NAME", "DISEASE") %in% names(df))) {
      pub_by_mirna <- tapply(df_common$PUBMED_COUNT, df_common$MIRNA_NAME, sum, na.rm = TRUE)
      df_common$PAIR_SCORE <- round(df_common$PUBMED_COUNT / pub_by_mirna[df_common$MIRNA_NAME], 3)
      pair_cols <- unique(df_common[, c("MIRNA_NAME", "DISEASE", "DSI", "PAIR_SCORE")])
      df <- dplyr::left_join(df, pair_cols, by = c("MIRNA_NAME", "DISEASE"))
    }
    # query_pubmed (both counts present): skip — no meaningful pair-level join
  }

  if (!for_plot) {
    df$MIRNA_NAME_clean <- df$MIRNA_NAME
    df$DISEASE_clean    <- df$DISEASE
    df$PUBMED_ID_clean  <- as.character(df$PUBMED_ID)
    df$MIRNA_NAME <- link_info_fun(df$MIRNA_NAME, "mirna")
    df$MIRBASE_ACC <- NULL
    df$PUBMED_ID   <- link_info_fun(df$PUBMED_ID, "pubmed")
    df$DISEASE     <- link_info_fun(df$DISEASE,   "disease")
  } else {
    df$MIRBASE_ACC <- NULL
  }
  df
}

# =============================================================================
# miRNA–Gene filter (by miRNA)
# =============================================================================
resolve_mirnas_mga <- function(input, my_data, mirnas_override = NULL) {
  if (!is.null(mirnas_override)) return(mirnas_override)
  if (identical(input$mga_mirna_source, "scratch")) {
    if (isTRUE(input$mga_use_all_scratch)) mirna_choices else input$mirna_gene_select_scratch
  } else {
    if (isTRUE(input$use_all_mirnas_from_mda))           my_data$mirnas_total
    else if (isTRUE(input$filter_mirna_common_checkbox)) my_data$mirna_common
    else                                                  input$mirna_gene_select
  }
}

get_filtered_mga <- function(con, input, my_data, mirnas_override = NULL, for_plot = FALSE, group_override = NULL) {
  selected_proofs  <- input$proofs
  proof_match_mode <- if (is.null(input$proof_match)) "any" else input$proof_match
  mirnas           <- resolve_mirnas_mga(input, my_data, mirnas_override)
  # Apply arm filter for all miRNAs coming from MDA source (premature names, no arm suffix)
  mda_arm <- if (!identical(input$mga_mirna_source, "scratch")) {
    arm_val <- input$mga_mda_arm %||% "5p"
    if (arm_val %in% c("5p", "3p", "max")) arm_val else NULL
  } else NULL
  mirnas_id <- get_mirna_ids(mirnas, con, arm = mda_arm)

  # Gene filter – parallel to miRNA filter (like miRNA-Disease tab)
  use_all_genes      <- isTRUE(input$mga_use_all_genes)
  only_metabolic     <- isTRUE(input$mga_only_metabolic_genes)
  genes              <- input$gene_search_select
  gene_ids           <- if (length(genes) > 0 && !use_all_genes && !only_metabolic) {
    dbGetQuery(con, glue_sql("SELECT GENE_ID FROM GENES WHERE GENE_NAME IN ({vals*})",
                             vals = genes, .con = con))$GENE_ID
  } else integer(0)

  gene_clause_filter <- if (use_all_genes || length(genes) == 0 && !only_metabolic) {
    DBI::SQL("")
  } else if (only_metabolic) {
    DBI::SQL(paste0("AND A.GENE_ID IN (SELECT GENE_ID FROM ", REACTIONS_GENES_TBL, ")"))
  } else if (length(gene_ids) > 0) {
    glue_sql("AND A.GENE_ID IN ({vals*})", vals = gene_ids, .con = con)
  } else {
    DBI::SQL("AND 1=0")
  }

  base_query <- if (for_plot) {
    query_mg[["query_by_gene_plot"]]
  } else {
    group_sel <- group_override %||% (input$group_mga %||% "None")
    q <- query_mg[["query_by_gene_general"]]
    if ("miRNA"          %in% group_sel) q <- query_mg[["query_by_gene_grouped_mirna"]]
    if ("gene"           %in% group_sel) q <- query_mg[["query_by_gene_grouped"]]
    if ("miRNA and gene" %in% group_sel) q <- query_mg[["query_by_gene_grouped_mirna_gene"]]
    q
  }

  # Min associations per miRNA-gene pair filter
  # query_by_gene_* queries have A = MIRNAS_GENES_ARTICLES (use A.GENE_ID)
  min_assoc_val <- as.integer(input$min_assoc_mga %||% 1L)
  min_assoc_clause <- if (min_assoc_val > 1L) {
    if (length(selected_proofs) > 0) {
      glue_sql(
        "JOIN (SELECT A2.MIRNA_ID, A2.GENE_ID FROM MIRNAS_GENES_ARTICLES A2 JOIN PROOF_MIRNA_GENE D2 ON A2.PROOF_ID = D2.PROOF_ID WHERE D2.PROOF_NAME IN ({selected_proofs*}) GROUP BY A2.MIRNA_ID, A2.GENE_ID HAVING COUNT(DISTINCT A2.PUBMED_ID) >= {min_assoc_val}) mga_pair_filter ON A.MIRNA_ID = mga_pair_filter.MIRNA_ID AND A.GENE_ID = mga_pair_filter.GENE_ID",
        selected_proofs = selected_proofs, min_assoc_val = min_assoc_val, .con = con
      )
    } else {
      DBI::SQL(sprintf(
        "JOIN (SELECT MIRNA_ID, GENE_ID FROM MIRNAS_GENES_ARTICLES GROUP BY MIRNA_ID, GENE_ID HAVING COUNT(DISTINCT PUBMED_ID) >= %d) mga_pair_filter ON A.MIRNA_ID = mga_pair_filter.MIRNA_ID AND A.GENE_ID = mga_pair_filter.GENE_ID",
        min_assoc_val
      ))
    }
  } else {
    DBI::SQL("")
  }

  mirna_clause_filter <- if (length(mirnas) > 0) {
    if (length(mirnas_id) == 0) DBI::SQL("AND 1=0")
    else glue_sql("AND A.MIRNA_ID IN ({vals*})", vals = mirnas_id, .con = con)
  } else DBI::SQL("AND 1=0")

  needs_pair_postfilter <- FALSE
  if (proof_match_mode == "all" && length(selected_proofs) > 1) {
    uses_mirna_gene_query <- base_query %in% c(
      query_mg[["query_by_gene_grouped_mirna_gene"]], query_mg[["query_by_gene_plot"]]
    )
    if (uses_mirna_gene_query) {
      needs_pair_postfilter <- TRUE
    } else {
      pair_clause <- and_proof_pair_sql(selected_proofs, mirnas_id, gene_ids, con, alias = "A")
      gene_clause_filter <- DBI::SQL(paste(as.character(gene_clause_filter),
                                           as.character(pair_clause)))
    }
  }

  proof_clause <- if (length(selected_proofs) > 0) {
    glue_sql("AND D.PROOF_NAME IN ({vals*})", vals = selected_proofs, .con = con)
  } else DBI::SQL("")

  # Human-only filter: exclude PUBMED_IDs explicitly marked IS_HUMAN = 0
  # query_mirna_gene / query_mirna_gene_plot use alias B for MIRNAS_GENES_ARTICLES,
  # all other MGA queries use alias A.
  mga_uses_b_alias <- base_query %in% c(
    query_mg[["query_mirna_gene"]], query_mg[["query_mirna_gene_plot"]]
  )
  use_human      <- isTRUE(input$filter_human_mga) && isTRUE(HAS_IS_HUMAN)
  human_clause   <- if (use_human)
    DBI::SQL("AND CAST(A.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)")
  else DBI::SQL("")
  human_clause_b <- if (use_human)
    DBI::SQL("AND CAST(B.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)")
  else DBI::SQL("")

  df <- dbGetQuery(con, glue_sql(base_query,
                                  gene_clause_filter  = gene_clause_filter,
                                  mirna_clause_filter = mirna_clause_filter,
                                  proof_clause        = proof_clause,
                                  min_assoc_clause    = min_assoc_clause,
                                  min_assoc_val       = min_assoc_val,
                                  human_clause        = human_clause,
                                  human_clause_b      = human_clause_b,
                                  .con = con))

  # L1 category filter: restricts all analysis to genes in selected top-level pathway category
  l1_sel <- as.integer(input$mga_pathway_l1[nzchar(input$mga_pathway_l1 %||% "")])
  if (length(l1_sel) > 0) {
    l1_desc_ids <- tryCatch(get_descendant_pathway_ids(con, l1_sel), error = function(e) integer(0))
    if (length(l1_desc_ids) > 0)
      df <- apply_pathway_filter_by_ids(df, l1_desc_ids, con)
  }

  # Degree filter: always computed from a flat (MIRNA × GENE) result so the
  # same valid_m / valid_g set is applied regardless of the display grouping.
  # For grouped queries GENE_NAME / MIRNA_NAME are GROUP_CONCAT strings —
  # tapply on those gives wrong counts, so we derive degrees from a separate
  # flat query when needed.
  min_m <- as.integer(input$min_degree_mirna_mga %||% 1L)
  min_g <- as.integer(input$min_degree_gene_mga  %||% 1L)

  if ((min_m > 1L || min_g > 1L) && nrow(df) > 0) {
    is_flat_mga <- base_query %in% c(
      query_mg[["query_by_gene_general"]],
      query_mg[["query_by_gene_grouped_mirna_gene"]],
      query_mg[["query_by_gene_plot"]]
    )
    flat_deg <- if (is_flat_mga) {
      df
    } else {
      dbGetQuery(con, glue_sql(
        query_mg[["query_by_gene_plot"]],
        gene_clause_filter  = gene_clause_filter,
        mirna_clause_filter = mirna_clause_filter,
        proof_clause        = proof_clause,
        min_assoc_clause    = min_assoc_clause,
        min_assoc_val       = min_assoc_val,
        .con = con
      ))
    }

    valid_m <- NULL
    if (min_m > 1L && nrow(flat_deg) > 0) {
      deg_m   <- tapply(flat_deg$GENE_NAME, flat_deg$MIRNA_NAME, function(x) length(unique(x)))
      valid_m <- names(deg_m)[deg_m >= min_m]
      flat_deg <- flat_deg[flat_deg$MIRNA_NAME %in% valid_m, ]
      if ("MIRNA_NAME" %in% names(df)) df <- df[df$MIRNA_NAME %in% valid_m, ]
    }
    if (min_g > 1L && nrow(flat_deg) > 0) {
      deg_g   <- tapply(flat_deg$MIRNA_NAME, flat_deg$GENE_NAME, function(x) length(unique(x)))
      valid_g <- names(deg_g)[deg_g >= min_g]
      if ("GENE_NAME" %in% names(df))
        df <- df[df$GENE_NAME %in% valid_g, ]
    }
  }

  if (needs_pair_postfilter) {
    valid_pairs <- and_proof_valid_pairs_names(selected_proofs, mirnas_id, gene_ids, con)
    if (is.null(valid_pairs) || nrow(valid_pairs) == 0) {
      df <- df[0, ]
    } else {
      df <- dplyr::semi_join(df, valid_pairs, by = c("MIRNA_NAME", "GENE_NAME"))
    }
  }

  if ("GENE_NAME" %in% colnames(df)) {
    df <- apply_expression_filter(df, get_expression_gene_names(input, "_mga", con))
    if (!for_plot && nrow(df) > 0) df <- add_expression_annotation(df, input, "_mga", con)
  }
  if (nrow(df) > 0) df$is_metabolic <- df$GENE_NAME %in% metabolic_genes

  if (!for_plot && nrow(df) > 0) {
    cols      <- setdiff(names(df), "is_metabolic")
    new_order <- append(cols, "is_metabolic", after = length(cols) - 1)
    df        <- df[, new_order]
    df$MIRNA_NAME_clean <- df$MIRNA_NAME
    df$GENE_NAME_clean  <- df$GENE_NAME
    df$PUBMED_ID_clean  <- as.character(df$PUBMED_ID)
    df$PUBMED_ID  <- link_info_fun(df$PUBMED_ID,  "pubmed")
    df$GENE_NAME  <- link_info_fun(df$GENE_NAME,  "gene")
    df$MIRNA_NAME <- link_info_fun(df$MIRNA_NAME, "mirna")
    df$MIRBASE_MATURE_ACC <- NULL
  }
  df
}


build_min_assoc_clause_mr <- function(n) {
  if (is.null(n) || n <= 1L) return(DBI::SQL(""))
  DBI::SQL(sprintf(
    "JOIN (SELECT MIRNA_ID AS mf_mid, GENE_ID AS mf_gid FROM MIRNAS_GENES_ARTICLES GROUP BY MIRNA_ID, GENE_ID HAVING COUNT(DISTINCT PUBMED_ID) >= %d) _mr_mf ON A.MIRNA_ID = _mr_mf.mf_mid AND A.GENE_ID = _mr_mf.mf_gid",
    n
  ))
}

# =============================================================================
# miRNA–Reaction filter (by miRNA)
# =============================================================================
get_filtered_mr <- function(con, input, for_plot, my_data, group_override = NULL) {
  mirnas <- if (isTRUE(input$use_all_mirnas_from_mga))            my_data$mirnas_from_mga
            else if (isTRUE(input$filter_rmirna_common_checkbox)) my_data$mirna_common
            else                                                   input$mirna_reaction_select

  mirnas_id           <- get_mirna_ids(mirnas, con)

  base_query <- query_mr[["query_general"]]

  mirna_clause_filter <- if (length(mirnas) > 0 && length(mirnas_id) > 0) {
    glue_sql("AND MIRNA_ID IN ({vals*})", vals = mirnas_id, .con = con)
  } else DBI::SQL("AND 1=0")

  genes <- if (isTRUE(input$filter_genes_reaction)) my_data$genes_from_mga %||% character(0)
           else input$gene_reaction_select

  genes_ids <- if (length(genes) > 0) {
    dbGetQuery(con,
      glue_sql("SELECT GENE_ID FROM GENES WHERE GENE_NAME IN ({vals*})", vals = genes, .con = con)
    )$GENE_ID
  } else integer(0)

  genes_clause_filter <- if (length(genes_ids) > 0) {
    glue_sql("AND GENE_ID IN ({vals*})", vals = genes_ids, .con = con)
  } else DBI::SQL("")

  selected_proofs_mr <- input$proofs %||% character(0)
  proof_match_mode   <- if (is.null(input$proof_match_mr)) "any" else input$proof_match_mr

  proof_clause <- if (length(selected_proofs_mr) > 0) {
    glue_sql("AND F.PROOF_NAME IN ({vals*})", vals = selected_proofs_mr, .con = con)
  } else DBI::SQL("")

  min_assoc_mr_val <- as.integer(my_data$min_assoc_mga_used %||% 1L)
  min_assoc_clause <- build_min_assoc_clause_mr(min_assoc_mr_val)

  human_clause <- if (isTRUE(input$filter_human_mr) && isTRUE(HAS_IS_HUMAN))
    DBI::SQL("AND CAST(A.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)")
  else DBI::SQL("")

  df <- dbGetQuery(con, glue_sql(base_query,
                                  mirna_clause_filter = mirna_clause_filter,
                                  genes_clause_filter = genes_clause_filter,
                                  proof_clause        = proof_clause,
                                  min_assoc_clause    = min_assoc_clause,
                                  human_clause        = human_clause,
                                  .con = con))

  # RXN_KEY = HUMAN_ID (the true unique reaction identifier in Human-GEM).
  # Falls back to NAME for reactions without a HUMAN_ID (e.g. user-uploaded models).
  # Using RXN_KEY instead of NAME as the internal key means each reaction variant
  # (same name but different GPR/compartment) is treated independently.
  if (nrow(df) > 0 && "HUMAN_ID" %in% colnames(df))
    df$RXN_KEY <- ifelse(is.na(df$HUMAN_ID) | df$HUMAN_ID == "", df$NAME, df$HUMAN_ID)

  # Deduplicate per (miRNA, gene, RXN_KEY): each HUMAN_ID is already unique in the DB
  # so duplicates here can only arise from multiple join paths.
  if (nrow(df) > 0 && all(c("MIRNA_NAME", "GENE_NAME", "RXN_KEY") %in% colnames(df)))
    df <- dplyr::distinct(df, MIRNA_NAME, GENE_NAME, RXN_KEY, .keep_all = TRUE)

  if ("GENE_NAME" %in% colnames(df)) {
    if (!for_plot && nrow(df) > 0) df <- add_expression_annotation(df, input, "_mr", con)
  }
  df <- apply_mr_postprocessing(df, con, input, for_plot)
  df
}
