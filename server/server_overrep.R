# === server_overrep.R ===
# Tab Pathway Annotation: annotation table, heatmap, sankey

observeProofsFilters(input, session, suffix = "_ovr")
# Expression context inherited from miRNA-Gene tab (_mga suffix)

# Cached reactive: computed once, reused by both info display and overrepresentation_results
ovr_expr_pathway_ids <- reactive({
  if (identical(input$ovr_pathway_source, "gmt")) return(NULL)
  filter_pathways_by_expression(
    con      = con,
    input    = input,
    suffix   = "_mga",
    min_frac = as.numeric(input$ovr_expr_min_frac %||% 20) / 100
  )
})

output$ovr_expr_context_banner <- renderUI({
  mode <- input$expression_filter_mode_mga %||% "none"
  if (mode == "none") return(NULL)
  ctx_label <- if (mode == "tissue") {
    tissues <- input$expression_tissue_select_mga %||% character(0)
    if (length(tissues) == 0) "tissue (none selected)"
    else paste(tissues, collapse = ", ")
  } else if (mode == "cell_line") {
    cls <- input$expression_cell_lines_mga %||% character(0)
    if (length(cls) == 0) "cell line (none selected)"
    else paste(cls, collapse = ", ")
  } else mode
  mode_lbl <- if (mode == "tissue") "Tissue" else "Cell line"
  tags$div(
    style = paste(
      "background:#e8f4fd; border-left:4px solid #3a8fc7;",
      "padding:8px 12px; margin-bottom:10px; border-radius:3px;",
      "font-size:12px; color:#1a5276;"
    ),
    icon("info-circle"),
    tags$strong(paste0(" Expression context active (from miRNA-Gene tab)")),
    tags$br(),
    tags$span(style = "color:#555;",
      paste0(mode_lbl, ": ", ctx_label,
             " — the two filters below apply to Reactome universe and results."))
  )
})

output$ovr_expr_filter_info <- renderUI({
  mode <- input$expression_filter_mode_mga %||% "none"
  if (mode == "none") return(NULL)
  if (identical(input$ovr_pathway_source, "gmt"))
    return(tags$p(style = "color:#b8860b; font-style:italic; margin-top:4px;",
                  icon("exclamation-triangle"),
                  " Expression filter not applied to custom GMT sets — filter manually before upload."))
  expr_df <- ovr_expr_pathway_ids()
  if (is.null(expr_df)) return(NULL)
  n_total <- tryCatch(
    dbGetQuery(con, "SELECT COUNT(DISTINCT PATHWAY_ID) AS n FROM ONTOLOGY_PATHWAYS_GENES")$n,
    error = function(e) NA_integer_
  )
  n_kept    <- nrow(expr_df)
  n_removed <- if (!is.na(n_total)) n_total - n_kept else NA_integer_
  min_frac  <- as.numeric(input$ovr_expr_min_frac %||% 20)
  tissue_lbl <- paste(
    input$expression_tissue_select_mga %||%
    input$expression_cell_lines_mga    %||% "?",
    collapse = ", ")
  tags$p(
    style = "color:#555; font-style:italic; margin-top:4px;",
    HTML(paste0(
      "<b>", n_kept, "</b> pathways pass (",
      if (!is.na(n_removed)) paste0("<b>", n_removed, "</b> removed, ") else "",
      min_frac, "% expressed, ", tissue_lbl, ")"
    ))
  )
})

# ── GMT upload: single shared reactive (gmt_data defined in server.R) ────────
observeEvent(input[["ul_ovr_gmt-reset"]], {
  gmt_data(NULL)
  shinyjs::reset("ul_ovr_gmt-file")
})

observeEvent(input[["ul_ovr_gmt-file"]], {
  f <- input[["ul_ovr_gmt-file"]]
  req(f)
  gmt <- parse_gmt_file(f$datapath)
  if (is.null(gmt)) { gmt_data(NULL); return() }
  gmt_data(gmt)
  showNotification(
    paste0("GMT loaded: ", length(gmt), " gene set(s). Click 'Run annotation' to use."),
    type = "message", duration = 5
  )
})

output$ovr_gmt_status <- render_gmt_status_ui(gmt_data)

# Build custom pathway structures from GMT (reactive, recomputes on new upload)
ovr_custom_pathway_data <- reactive({
  gmt <- gmt_data()
  if (is.null(gmt) || length(gmt) == 0) return(NULL)
  # Use gene set names as both ID and display name
  genes_list <- gmt
  sizes      <- vapply(genes_list, length, integer(1L))
  names_map  <- setNames(names(gmt), names(gmt))
  gene_index <- build_gene_pathway_index(genes_list)
  list(
    pathway_sizes      = sizes,
    pathway_genes_list = genes_list,
    pathway_names_map  = names_map,
    gene_pathway_index = gene_index
  )
})

updateSelectizeInput(session, "ovr_pathway_upload_names", choices = pathway_choices, server = TRUE)

ovr_uploaded_pathway_ids <- reactiveVal(integer(0))

upload_list_server("ul_ovr_pathway_names", session, "ovr_pathway_upload_names",
  type = NULL, valid_choices = pathway_choices)

# Convert uploaded pathway names → IDs, store separately from manual L2 selection
observe({
  names_sel <- input$ovr_pathway_upload_names
  if (length(names_sel) == 0) { ovr_uploaded_pathway_ids(integer(0)); return() }
  ids <- tryCatch(
    dbGetQuery(con, glue_sql(
      "SELECT PATHWAY_ID FROM ONTOLOGY_PATHWAYS WHERE PATHWAY_NAME IN ({vals*})",
      vals = names_sel, .con = con
    ))$PATHWAY_ID,
    error = function(e) integer(0)
  )
  ovr_uploaded_pathway_ids(as.integer(ids))
}) |> bindEvent(input$ovr_pathway_upload_names, ignoreNULL = FALSE, ignoreInit = TRUE)

# ── Hierarchy availability: defined in server.R (shared with MGA tab) ─────────

output$ovr_hierarchy_ui <- renderUI({
  if (ontology_has_hierarchy()) {
    fluidRow(
      column(5,
        h5(icon("sitemap"), "Reactome hierarchy filter (optional)"),
        selectizeInput("ovr_pathway_l1", "Filter Top-level category:",
          choices = pathway_l1_choices, selected = NULL, multiple = TRUE,
          options = list(placeholder = "All top-level pathways…")),
        selectizeInput("ovr_pathway_l3", "Search specific pathway:",
          choices = NULL, selected = NULL, multiple = TRUE,
          options = list(placeholder = "Search any pathway…", maxOptions = 500)),
        upload_list_ui("ul_ovr_pathway_names",
          tags$small(style = "color:#888;", "Or upload a .txt list of pathway names"))
      ),
      column(3, br(),
        checkboxInput("ovr_include_descendants",
                      tagList(icon("sitemap"), " Include all sub-pathways"),
                      value = TRUE),
        tags$small(style = "color:#888;",
          "When checked, enrichment includes all child pathways of the selection.")
      ),
      column(2, br(), br(),
        actionButton("reset_ovr_pathway_hierarchy", "Reset hierarchy",
                     icon = icon("undo"), class = "btn-sm")
      )
    )
  } else {
    div(
      style = "margin-bottom:10px;",
      tags$span(
        class = "label label-warning",
        style = "font-size:13px; padding:6px 10px;",
        icon("exclamation-triangle"),
        " Flat ontology loaded — hierarchy filter unavailable."
      ),
      tags$p(
        style = "color:#888; font-size:12px; margin-top:6px;",
        "Restore original ontology (Data upload tab) to re-enable hierarchy filtering."
      )
    )
  }
})

# ── L1 → L3 cascade: when L1 empty load all pathways, when L1 selected load descendants ──
observeEvent(input$ovr_pathway_l1, {
  l1_ids <- as.integer(input$ovr_pathway_l1)
  if (length(l1_ids) == 0) {
    all_pw <- tryCatch(
      dbGetQuery(con, "SELECT PATHWAY_ID, PATHWAY_NAME FROM ONTOLOGY_PATHWAYS ORDER BY PATHWAY_NAME"),
      error = function(e) data.frame(PATHWAY_ID = integer(0), PATHWAY_NAME = character(0))
    )
    updateSelectizeInput(session, "ovr_pathway_l3",
      choices  = setNames(all_pw$PATHWAY_ID, all_pw$PATHWAY_NAME),
      selected = character(0), server = TRUE,
      options  = list(placeholder = "Search any pathway…", maxOptions = 500))
    return()
  }
  desc_ids <- get_descendant_pathway_ids(con, l1_ids)
  sub_pw <- if (length(desc_ids) > 0) {
    tryCatch(
      dbGetQuery(con, glue_sql(
        "SELECT PATHWAY_ID, PATHWAY_NAME FROM ONTOLOGY_PATHWAYS
         WHERE PATHWAY_ID IN ({ids*}) ORDER BY PATHWAY_NAME",
        ids = desc_ids, .con = con
      )),
      error = function(e) data.frame(PATHWAY_ID = integer(0), PATHWAY_NAME = character(0))
    )
  } else data.frame(PATHWAY_ID = integer(0), PATHWAY_NAME = character(0))
  updateSelectizeInput(session, "ovr_pathway_l3",
    choices  = setNames(sub_pw$PATHWAY_ID, sub_pw$PATHWAY_NAME),
    selected = character(0), server = TRUE,
    options  = list(placeholder = "Search within category…", maxOptions = 500))
}, ignoreNULL = FALSE)

observeEvent(input$reset_ovr_pathway_hierarchy, {
  updateSelectizeInput(session, "ovr_pathway_l1", selected = character(0))
  updateSelectizeInput(session, "ovr_pathway_l3", selected = character(0))
  updateSelectizeInput(session, "ovr_pathway_upload_names", selected = character(0))
  ovr_uploaded_pathway_ids(integer(0))
  updateCheckboxInput(session, "ovr_include_descendants", value = TRUE)
})

# ── Core reactive: build annotation data frame ───────────────────────────────
overrepresentation_results <- eventReactive(input$run_overpresentation, {
  withProgress(message = "Running pathway annotation...", value = 0, {

    setProgress(0.1, detail = "Loading gene associations")
    req(filtered_mga())
    df_mga          <- filtered_mga_flat()
    genes_selected  <- if (isTRUE(input$filter_ontology_genes_reaction)) NULL else input$ontology_gene_select
    mirnas_selected <- if (isTRUE(input$use_ontology_all_mirnas_from_mga)) NULL else input$ontology_mirna_select

    setProgress(0.3, detail = "Filtering gene associations")
    df_mga_filtered <- filter_mga_for_ontology(
      df_mga          = df_mga,
      genes_selected  = genes_selected,
      mirnas_selected = mirnas_selected
    )
    if (nrow(df_mga_filtered) == 0) return(list(combined = data.frame(), per_mirna = data.frame()))

    my_data$ovr_mga_df <- df_mga_filtered

    use_gmt    <- identical(input$ovr_pathway_source, "gmt")
    custom_pwd <- if (use_gmt) ovr_custom_pathway_data() else NULL

    if (use_gmt) {
      shiny::validate(shiny::need(!is.null(custom_pwd),
        "No GMT file loaded. Upload a .gmt file in the Pathway filters box."))
    }

    pathway_ids_filter <- if (!use_gmt && ontology_has_hierarchy()) {
      selected_l3  <- as.integer(input$ovr_pathway_l3)
      uploaded_ids <- ovr_uploaded_pathway_ids()
      selected_l1  <- as.integer(input$ovr_pathway_l1)
      manual_ids   <- unique(c(selected_l3, uploaded_ids))
      if (length(manual_ids) > 0) {
        # L3/upload explicit: checkbox controls whether to expand to descendants
        if (isTRUE(input$ovr_include_descendants))
          get_descendant_pathway_ids(con, manual_ids)
        else
          manual_ids
      } else if (length(selected_l1) > 0) {
        # L1 only: always expand (L1 nodes have no direct enrichable pathways)
        get_descendant_pathway_ids(con, selected_l1)
      } else {
        NULL
      }
    } else {
      NULL
    }

    # ── Expression filter: keep only pathways with ≥ min_frac expressed genes ──
    expr_df <- ovr_expr_pathway_ids()  # data.frame(PATHWAY_ID, expressed_genes) or NULL
    if (!is.null(expr_df) && nrow(expr_df) > 0) {
      expr_pid_chr <- as.character(expr_df$PATHWAY_ID)
      pathway_ids_filter <- if (!is.null(pathway_ids_filter))
        intersect(as.character(pathway_ids_filter), expr_pid_chr)
      else
        expr_pid_chr
    } else {
      expr_df <- NULL
    }

    common_args <- list(
      df_mga_filtered      = df_mga_filtered,
      top_n                = as.integer(input$ontology_top_n),
      rank_by              = "COVERAGE",
      pathway_ids_filter   = pathway_ids_filter,
      custom_pathway_sizes = custom_pwd$pathway_sizes,
      custom_genes_list    = custom_pwd$pathway_genes_list,
      custom_names_map     = custom_pwd$pathway_names_map,
      custom_gene_index    = custom_pwd$gene_pathway_index
    )

    setProgress(0.55, detail = "Annotating (combined)")
    res_combined <- do.call(annotate_mirna_pathways,
                            c(common_args, list(mode = "combined")))

    setProgress(0.75, detail = "Annotating (per miRNA)")
    res_per_mirna <- do.call(annotate_mirna_pathways,
                             c(common_args, list(mode = "per_mirna")))

    # ── Add TISSUE_COVERAGE when expression filter active ────────────────────
    add_tissue_coverage <- function(res, expr_df, expressed_genes_ovr) {
      if (is.null(expr_df) || nrow(res) == 0 || !"PATHWAY_ID" %in% names(res)) return(res)
      expr_set   <- expressed_genes_ovr  # character vector of expressed gene names
      k_expr_map <- setNames(expr_df$expressed_genes, as.character(expr_df$PATHWAY_ID))
      res$TISSUE_COVERAGE <- vapply(seq_len(nrow(res)), function(i) {
        pid      <- as.character(res$PATHWAY_ID[i])
        k_expr   <- k_expr_map[pid]
        if (is.na(k_expr) || k_expr == 0) return(NA_real_)
        genes_found <- strsplit(res$GENES_FOUND[i], ", ", fixed = TRUE)[[1]]
        n_expr_found <- sum(genes_found %in% expr_set)
        round(n_expr_found / k_expr, 3)
      }, numeric(1))
      res
    }

    if (!is.null(expr_df)) {
      expressed_genes_ovr <- get_expression_gene_names(input, "_mga", con)
      res_combined  <- add_tissue_coverage(res_combined,  expr_df, expressed_genes_ovr)
      res_per_mirna <- add_tissue_coverage(res_per_mirna, expr_df, expressed_genes_ovr)
    }

    list(combined = res_combined, per_mirna = res_per_mirna)
  })
})

observe({
  res <- overrepresentation_results()

  # Attach L1 ancestor category to all result dataframes (when hierarchy available)
  .add_l1 <- function(df) {
    if (is.null(df) || nrow(df) == 0 || !"PATHWAY_ID" %in% names(df)) return(df)
    if (!ONTOLOGY_HAS_HIERARCHY) { df$L1_CATEGORY <- NA_character_; return(df) }
    l1_map <- tryCatch(get_l1_ancestor_map(con, unique(df$PATHWAY_ID)), error = function(e) character(0))
    df$L1_CATEGORY <- if (length(l1_map) > 0 && !is.null(names(l1_map)))
      l1_map[as.character(df$PATHWAY_ID)]
    else
      rep(NA_character_, nrow(df))
    df
  }

  combined  <- .add_l1(res$combined)
  per_mirna <- .add_l1(res$per_mirna)

  # If L1 filter active, drop rows whose resolved L1_CATEGORY is outside the selection.
  # Needed because a pathway can be a descendant of multiple L1 nodes; get_l1_ancestor_map
  # assigns only one primary L1, which may differ from the selected ones.
  selected_l1_ids <- as.integer(input$ovr_pathway_l1 %||% integer(0))
  if (length(selected_l1_ids) > 0 && ONTOLOGY_HAS_HIERARCHY) {
    selected_l1_names <- tryCatch(
      dbGetQuery(con, glue_sql(
        "SELECT PATHWAY_NAME FROM ONTOLOGY_PATHWAYS WHERE PATHWAY_ID IN ({ids*})",
        ids = selected_l1_ids, .con = con
      ))$PATHWAY_NAME,
      error = function(e) character(0)
    )
    if (length(selected_l1_names) > 0) {
      .filter_l1 <- function(df) {
        if (is.null(df) || nrow(df) == 0 || !"L1_CATEGORY" %in% names(df)) return(df)
        df[!is.na(df$L1_CATEGORY) & df$L1_CATEGORY %in% selected_l1_names, ]
      }
      combined  <- .filter_l1(combined)
      per_mirna <- .filter_l1(per_mirna)
    }
  }

  coverage_min <- as.numeric(input$ovr_tissue_coverage_min %||% 0) / 100
  if (coverage_min > 0) {
    .filter_coverage <- function(df) {
      if (is.null(df) || nrow(df) == 0 || !"TISSUE_COVERAGE" %in% names(df)) return(df)
      df[!is.na(df$TISSUE_COVERAGE) & df$TISSUE_COVERAGE >= coverage_min, ]
    }
    combined  <- .filter_coverage(combined)
    per_mirna <- .filter_coverage(per_mirna)
  }

  my_data$over_df_combined  <- combined
  my_data$over_df_per_mirna <- per_mirna
  my_data$over_df           <- my_data$over_df_combined
})

# ── Count bar ─────────────────────────────────────────────────────────────────
output$count_mgoa <- renderUI({
  df  <- my_data$over_df
  mga <- my_data$ovr_mga_df
  req(!is.null(df) && nrow(df) > 0)

  # build a flat df with miRNA and gene counts from ovr_mga_df
  cols <- list("pathways" = "PATHWAY_NAME")
  if (!is.null(mga) && nrow(mga) > 0) {
    mirna_col <- if ("MIRNA_NAME_clean" %in% names(mga)) "MIRNA_NAME_clean" else "MIRNA_NAME"
    gene_col  <- if ("GENE_NAME_clean"  %in% names(mga)) "GENE_NAME_clean"  else "GENE_NAME"
    combined_df <- cbind(
      df[, intersect(c("PATHWAY_NAME"), names(df)), drop = FALSE],
      data.frame(
        MIRNA_NAME = rep(mga[[mirna_col]][1], nrow(df)),
        GENE_NAME  = rep(mga[[gene_col]][1],  nrow(df))
      )
    )
    # overwrite with real counts via render_count_bar on mga
    mirna_n <- dplyr::n_distinct(strip_html(mga[[mirna_col]]))
    gene_n  <- dplyr::n_distinct(strip_html(mga[[gene_col]]))
    path_n  <- dplyr::n_distinct(df$PATHWAY_NAME)

    parts <- c(
      paste0("<b>", nrow(df),  "</b> rows"),
      paste0("<b>", path_n,   "</b> unique pathways"),
      paste0("<b>", mirna_n,  "</b> unique miRNA"),
      paste0("<b>", gene_n,   "</b> unique genes")
    )
  } else {
    path_n  <- dplyr::n_distinct(df$PATHWAY_NAME)
    parts <- c(
      paste0("<b>", nrow(df), "</b> rows"),
      paste0("<b>", path_n,  "</b> unique pathways")
    )
    if ("MIRNA_NAME" %in% names(df))
      parts <- c(parts, paste0("<b>", dplyr::n_distinct(df$MIRNA_NAME), "</b> unique miRNA"))
  }

  tags$p(style = "color:#555; font-style:italic; margin-bottom:6px;",
         HTML(paste(parts, collapse = " &nbsp;|&nbsp; ")))
})

# ── Table ─────────────────────────────────────────────────────────────────────
output$download_overrep_csv <- make_csv_download(
  "pathway_annotation", reactive(my_data$over_df), con = NULL
)

output$pathway_overrepresentation_table <- DT::renderDT({
  df <- my_data$over_df
  req(!is.null(df) && nrow(df) > 0)

  # Join pathway member genes from DB (all genes in pathway, not just found ones)
  pw_ids <- unique(df$PATHWAY_ID)
  genes_per_pathway <- tryCatch({
    if (is.null(pw_ids) || length(pw_ids) == 0)
      setNames(character(0), character(0))
    else {
      gdf <- dbGetQuery(con, glue_sql(
        "SELECT PG.PATHWAY_ID, G.GENE_NAME
           FROM ONTOLOGY_PATHWAYS_GENES PG
           JOIN GENES G ON PG.GENE_ID = G.GENE_ID
          WHERE PG.PATHWAY_ID IN ({ids*})",
        ids = pw_ids, .con = con
      ))
      if (nrow(gdf) == 0)
        setNames(character(0), character(0))
      else
        tapply(gdf$GENE_NAME, gdf$PATHWAY_ID, function(x) paste(sort(x), collapse = ", "))
    }
  }, error = function(e) setNames(character(0), character(0)))

  # Add GENES_IN_PATHWAY column; if lookup fails, skip gracefully
  tryCatch({
    pw_char <- as.character(df$PATHWAY_ID)
    gp      <- genes_per_pathway  # named character vector or character(0)
    vals    <- if (length(gp) > 0 && !is.null(names(gp)))
                 as.character(gp[pw_char])
               else
                 rep(NA_character_, nrow(df))
    df$GENES_IN_PATHWAY <- vals
  }, error = function(e) NULL)

  has_genes_col <- "GENES_IN_PATHWAY" %in% names(df)
  genes_col_idx <- if (has_genes_col) which(names(df) == "GENES_IN_PATHWAY") - 1L else NULL
  pid_col_idx   <- if ("PATHWAY_ID" %in% names(df)) which(names(df) == "PATHWAY_ID") - 1L else NULL

  genes_def <- if (has_genes_col) list(list(
    targets = genes_col_idx,
    render  = DT::JS(
      "function(data, type, row) {",
      "  if (type !== 'display' || !data || data === 'NA') return data || '';",
      "  var genes = data.split(', ');",
      "  if (genes.length <= 10) return data;",
      "  return genes.slice(0, 10).join(', ') + ' … [+' + (genes.length - 10) + ' more]';",
      "}"
    )
  )) else list()

  pid_def <- if (!is.null(pid_col_idx))
    list(list(targets = pid_col_idx, visible = FALSE))
  else list()

  col_defs <- c(genes_def, pid_def)

  DT::datatable(
    df,
    filter = "top", escape = FALSE, extensions = "Buttons",
    options = list(
      pageLength = 15, scrollX = TRUE, deferRender = TRUE,
      dom = "Bfrtip",
      buttons  = list(list(extend = "colvis", text = "Columns")),
      columnDefs = col_defs
    )
  )
}, server = TRUE)

# ── Barplot combined: pathway sorted by coverage (or grouped by L1) ───────────
output$ovr_barplot_combined <- renderPlotly({
  df <- my_data$over_df_combined
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "Run annotation to see combined results."))
  df <- head(df[order(-df$COVERAGE), ], as.integer(input$ontology_top_n %||% 20L))

  group_l1   <- isTRUE(input$ovr_barplot_group_l1)
  has_l1     <- "L1_CATEGORY" %in% names(df) && any(!is.na(df$L1_CATEGORY))
  x_label    <- "Coverage (N_found / K)"

  trunc_name <- function(x, n = 45)
    ifelse(nchar(x) > n, paste0(substr(x, 1, n - 2), "…"), x)

  if (group_l1 && has_l1) {
    df <- df %>%
      dplyr::mutate(L1_CATEGORY = dplyr::coalesce(L1_CATEGORY, "Unknown")) %>%
      dplyr::arrange(dplyr::desc(COVERAGE)) %>%
      dplyr::mutate(
        PATHWAY_SHORT = trunc_name(PATHWAY_NAME),
        tooltip = paste0(
          "<b>", PATHWAY_NAME, "</b><br>",
          "L1: ", L1_CATEGORY, "<br>",
          "Coverage: ", COVERAGE, " (", N_FOUND, "/", K, ")<br>",
          "miRNAs: ", N_MIRNAS
        )
      )
    df$PATHWAY_SHORT <- factor(df$PATHWAY_SHORT, levels = rev(unique(df$PATHWAY_SHORT)))

    l1_cats   <- unique(df$L1_CATEGORY)
    l1_colors <- setNames(
      grDevices::hcl.colors(length(l1_cats), palette = "Set2"),
      l1_cats
    )

    p <- plot_ly()
    for (cat in l1_cats) {
      sub <- df[df$L1_CATEGORY == cat, ]
      p <- add_trace(p,
        data = sub, x = ~COVERAGE, y = ~PATHWAY_SHORT,
        type = "bar", orientation = "h", name = cat,
        marker = list(color = l1_colors[[cat]]),
        text = ~tooltip, hoverinfo = "text"
      )
    }
    p <- p %>% layout(
      barmode = "overlay",
      title   = list(text = paste0("Pathway annotation — grouped by L1 — ", nrow(df), " pathways"), font = list(size = 14)),
      xaxis   = list(title = x_label, range = c(0, 1), gridcolor = "#e8e8e8"),
      yaxis   = list(title = "", tickfont = list(size = 10)),
      legend  = list(title = list(text = "L1 category"), tracegroupgap = 2),
      margin  = list(l = 300, r = 180, t = 50, b = 60),
      paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa"
    )
  } else {
    df <- df %>%
      dplyr::arrange(COVERAGE) %>%
      dplyr::mutate(
        PATHWAY_SHORT = trunc_name(PATHWAY_NAME),
        tooltip = paste0(
          "<b>", PATHWAY_NAME, "</b><br>",
          if (has_l1) paste0("L1: ", dplyr::coalesce(L1_CATEGORY, "—"), "<br>") else "",
          "Coverage: ", COVERAGE, " (", N_FOUND, "/", K, ")<br>",
          "miRNAs: ", N_MIRNAS
        )
      )
    df$PATHWAY_SHORT <- factor(df$PATHWAY_SHORT, levels = unique(df$PATHWAY_SHORT))

    p <- plot_ly(df,
      x = ~COVERAGE, y = ~PATHWAY_SHORT, type = "bar", orientation = "h",
      marker = list(
        color = ~N_FOUND, colorscale = "Viridis", showscale = TRUE,
        colorbar = list(title = "Genes found", len = 0.5, thickness = 14)
      ),
      text = ~tooltip, hoverinfo = "text"
    ) %>%
    layout(
      title  = list(text = paste0("Pathway annotation — ", nrow(df), " pathways"), font = list(size = 14)),
      xaxis  = list(title = x_label, range = c(0, 1), gridcolor = "#e8e8e8"),
      yaxis  = list(title = "", tickfont = list(size = 10)),
      margin = list(l = 300, r = 120, t = 50, b = 60),
      paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa"
    )
  }

  p %>%
    plotly_clean_config() %>%
    htmlwidgets::onRender("function(el) {
      Plotly.toImage(el, {format:'png', width:2580, height:1440}).then(function(url) {
        Shiny.setInputValue('ovr_plot_png', url, {priority:'event'});
      });
    }")
}) |> bindEvent(input$run_overpresentation, input$ovr_barplot_group_l1, ignoreNULL = TRUE)

store_png("ovr_plot_png", my_data, input)

# ── L1 summary table ──────────────────────────────────────────────────────────
output$ovr_l1_summary_table <- DT::renderDT({
  df <- my_data$over_df_combined
  shiny::validate(shiny::need(
    !is.null(df) && nrow(df) > 0 &&
      "L1_CATEGORY" %in% names(df) && any(!is.na(df$L1_CATEGORY)),
    "No L1 category data available."
  ))

  get_mirnas <- if ("MIRNAS_FOUND" %in% names(df)) {
    function(sub) unique(trimws(unlist(strsplit(
      as.character(sub$MIRNAS_FOUND[!is.na(sub$MIRNAS_FOUND)]), ","))))
  } else function(sub) character(0)

  cats <- sort(unique(df$L1_CATEGORY[!is.na(df$L1_CATEGORY)]))
  tbl <- do.call(rbind, lapply(cats, function(cat) {
    sub   <- df[!is.na(df$L1_CATEGORY) & df$L1_CATEGORY == cat, ]
    n_pw  <- dplyr::n_distinct(sub$PATHWAY_NAME)
    n_m   <- length(get_mirnas(sub))
    genes <- unique(trimws(unlist(strsplit(
      as.character(sub$GENES_FOUND[!is.na(sub$GENES_FOUND) & nzchar(sub$GENES_FOUND)]), ","))))
    n_g   <- length(genes[nzchar(genes)])
    data.frame(
      `L1 Category`      = cat,
      `Pathways enriched` = n_pw,
      `miRNAs`           = n_m,
      `Genes`            = n_g,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }))

  DT::datatable(tbl,
    rownames   = FALSE, escape = FALSE,
    extensions = "Buttons",
    options    = list(
      pageLength = 25, ordering = TRUE,
      dom     = "Bt",
      buttons = list(
        list(extend = "csv",   filename = "l1_pathway_summary",
             exportOptions = list(modifier = list(page = "all"), columns = ":all")),
        list(extend = "excel", filename = "l1_pathway_summary",
             exportOptions = list(modifier = list(page = "all"), columns = ":all"))
      )
    )
  )
}, server = FALSE) |> bindEvent(my_data$over_df_combined, ignoreNULL = TRUE)


# ── Disease annotation helper (used by MR heatmap) ───────────────────────────
# Returns named char vector: miRNA → label ("Disease A" / "Disease B" / "Shared")
# NULL if MDA not run or >2 diseases.
.ovr_disease_annotation <- reactive({
  mda <- tryCatch(filtered_mda_flat(), error = function(e) NULL)
  if (is.null(mda) || nrow(mda) == 0 || !"DISEASE" %in% names(mda)) return(NULL)
  diseases <- unique(mda$DISEASE)
  if (length(diseases) == 0 || length(diseases) > 2) return(NULL)
  if (length(diseases) == 1) {
    mirnas <- unique(gsub("-[35]p$", "", as.character(mda$MIRNA_NAME), ignore.case = TRUE))
    return(setNames(rep(diseases[1], length(mirnas)), mirnas))
  }
  to_premature <- function(x) gsub("-[35]p$", "", as.character(x), ignore.case = TRUE)
  d1 <- diseases[1]; d2 <- diseases[2]
  m1 <- unique(to_premature(mda$MIRNA_NAME[mda$DISEASE == d1]))
  m2 <- unique(to_premature(mda$MIRNA_NAME[mda$DISEASE == d2]))
  all_m <- union(m1, m2)
  label <- dplyr::case_when(
    all_m %in% m1 & all_m %in% m2 ~ "Shared",
    all_m %in% m1                  ~ d1,
    TRUE                           ~ d2
  )
  setNames(label, all_m)
})

# ── Heatmap miRNA x Pathway (per_mirna mode) ──────────────────────────────────
output$mirna_pathway_heatmap <- renderPlotly({
  df <- my_data$over_df_per_mirna
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0 && "MIRNA_NAME" %in% names(df),
                "Run annotation to see the per-miRNA heatmap."))

  value_col_hm <- "COVERAGE"

  top_n_paths <- as.integer(input$ontology_top_n %||% 20L)
  top_paths <- df %>%
    dplyr::group_by(PATHWAY_NAME) %>%
    dplyr::summarise(TOTAL = sum(.data[[value_col_hm]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(TOTAL)) %>%
    dplyr::slice_head(n = top_n_paths) %>%
    dplyr::pull(PATHWAY_NAME)

  mat_df <- df %>%
    dplyr::filter(PATHWAY_NAME %in% top_paths) %>%
    dplyr::select(MIRNA_NAME, PATHWAY_NAME, dplyr::all_of(value_col_hm)) %>%
    tidyr::pivot_wider(names_from = PATHWAY_NAME, values_from = dplyr::all_of(value_col_hm), values_fill = 0)

  # Add miRNAs present in input data but absent from any top pathway (all-zero rows)
  all_mirnas <- unique(as.character(my_data$ovr_mga_df$MIRNA_NAME))
  missing    <- setdiff(all_mirnas, mat_df$MIRNA_NAME)
  if (length(missing) > 0) {
    zero_rows        <- as.data.frame(matrix(0L, nrow = length(missing), ncol = ncol(mat_df) - 1L))
    names(zero_rows) <- names(mat_df)[-1]
    zero_rows        <- cbind(MIRNA_NAME = missing, zero_rows)
    mat_df           <- dplyr::bind_rows(mat_df, zero_rows)
  }

  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$MIRNA_NAME

  # L1 category map for columns (pathways)
  l1_col_map <- NULL
  has_l1_hm  <- "L1_CATEGORY" %in% names(df) && any(!is.na(df$L1_CATEGORY))
  if (has_l1_hm) {
    l1_col_map <- setNames(
      dplyr::coalesce(df$L1_CATEGORY, "Unknown"),
      df$PATHWAY_NAME
    )
  }

  if (!isTRUE(input$ovr_heatmap_cluster_col)) {
    if (has_l1_hm) {
      # Sort columns by L1 then by total coverage within L1
      col_l1    <- l1_col_map[colnames(mat)]
      col_score <- colSums(mat, na.rm = TRUE)
      col_ord   <- order(col_l1, -col_score)
      mat <- mat[, col_ord, drop = FALSE]
    } else {
      mat <- mat[, order(colnames(mat)), drop = FALSE]
    }
  }

  # Truncate long pathway names keeping start + end so shared-prefix names stay distinct
  trunc_mid <- function(x, max_len = 50) {
    long <- nchar(x) > max_len
    half <- (max_len - 1L) %/% 2L
    ifelse(long,
           paste0(substr(x, 1, half), "…", substr(x, nchar(x) - half + 1L, nchar(x))),
           x)
  }
  orig_colnames  <- colnames(mat)   # keep for L1 lookup after truncation
  colnames(mat)  <- trunc_mid(colnames(mat))

  transposed <- isTRUE(input$ovr_heatmap_transpose)

  # Build L1 annotation aligned to (truncated) colnames(mat) for heatmaply side strip
  l1_ann <- if (has_l1_hm && !isTRUE(input$ovr_heatmap_cluster_col)) {
    setNames(l1_col_map[orig_colnames], colnames(mat))
  } else NULL

  p <- render_mirna_target_heatmap(
    mat,
    cluster_rows          = isTRUE(input$ovr_heatmap_cluster_row),
    cluster_cols          = isTRUE(input$ovr_heatmap_cluster_col),
    x1_label              = "miRNA",
    x2_label              = "Pathway",
    row_annotation        = NULL,
    l1_pathway_annotation = l1_ann,
    zero_white            = TRUE,
    seriate               = "none",
    transpose             = transposed,
    color_palette         = input$ovr_heatmap_palette %||% "YlGnBu"
  )

  p %>% htmlwidgets::onRender("function(el) {
    Plotly.toImage(el, {format:'png', width:2580, height:1680}).then(function(url) {
      Shiny.setInputValue('ovr_hm_png', url, {priority:'event'});
    });
  }")
}) |> bindEvent(input$run_overpresentation,
               input$ovr_heatmap_cluster_row, input$ovr_heatmap_cluster_col,
               input$ovr_heatmap_transpose,   input$ovr_heatmap_palette,
               ignoreNULL = TRUE)

store_png("ovr_hm_png", my_data, input)

# (Sankey removed)
