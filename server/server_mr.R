# === server_mr.R ===
# Tab miRNA-Reaction: tabella, grafo, legenda, click-link modals

# ── L1 (category) → L2 (subsystem) cascade ───────────────────────────────────
observeEvent(input$mr_filter_l1, {
  l1_cats <- input$mr_filter_l1
  if (length(l1_cats) == 0) {
    updateSelectizeInput(session, "mr_filter_l2",
      choices = subsystem_choices, selected = character(0), server = TRUE,
      options = list(placeholder = "Select category above, or search all…", maxOptions = 500L))
  } else {
    sub_in_cat <- sort(intersect(
      names(SUBSYSTEM_CATEGORY_MAP)[!is.na(SUBSYSTEM_CATEGORY_MAP) & SUBSYSTEM_CATEGORY_MAP %in% l1_cats],
      subsystem_choices))
    updateSelectizeInput(session, "mr_filter_l2",
      choices = sub_in_cat, selected = character(0), server = TRUE,
      options = list(placeholder = "Subsystems in selected category…", maxOptions = 500L))
  }
}, ignoreNULL = FALSE, ignoreInit = TRUE)

# ── Reset pathway filter ──────────────────────────────────────────────────────
observeEvent(input$reset_mr_pathway_filter, {
  updateSelectizeInput(session, "mr_filter_l1", selected = character(0))
  updateSelectizeInput(session, "mr_filter_l2",
    choices = subsystem_choices, selected = character(0), server = TRUE,
    options = list(placeholder = "Select category above, or search all…"))
  updateSelectizeInput(session, "mr_filter_exclude",
    choices = sort(union(category_choices, subsystem_choices)), selected = character(0), server = TRUE,
    options = list(placeholder = "Exclude subsystem or category…"))
})

observeEvent(input$reset_mirna_reaction, {
  updateSelectizeInput(session, "mirna_reaction_select",         selected = character(0))
  updateCheckboxInput(session,  "filter_rmirna_common_checkbox", value = FALSE)
  updateCheckboxInput(session,  "use_all_mirnas_from_mga",       value = FALSE)
})

# Aggiorna my_data$mirnas_from_mga (nomi maturi esatti) quando filtered_mga cambia
observeEvent(filtered_mga(), {
  flat <- filtered_mga_flat()
  my_data$mirnas_from_mga <- unique(flat$MIRNA_NAME)
  my_data$genes_from_mga  <- unique(flat$GENE_NAME)
}, ignoreNULL = TRUE)

# Single DB query – one row per (MIRNA, GENE, REACTION); runs only on Search click
filtered_mr_base <- eventReactive(input$obtain_mr, {
  on.exit(removeNotification("notify_mr"))
  df <- get_filtered_mr(con = con, input = input, for_plot = TRUE, my_data, group_override = "None")
  if ("NAME" %in% names(df) && !"REACTION_NAME" %in% names(df))
    df <- dplyr::rename(df, REACTION_NAME = NAME)
  shiny::validate(shiny::need(nrow(df) > 0, "No associations found with the current filters. Try broadening your selection."))
  df
})

# Table view – instant dplyr group-by, re-runs on group_mr change without hitting DB
filtered_mr <- reactive({
  df   <- filtered_mr_base()
  mode <- input$group_mr %||% "None"
  if (!"SUBSYSTEM" %in% names(df)) {
    sub_lkp <- dbGetQuery(con, glue_sql(
      "SELECT NAME AS REACTION_NAME, MIN(SUBSYSTEM) AS SUBSYSTEM FROM REACTIONS
        WHERE NAME IN ({vals*}) GROUP BY NAME",
      vals = unique(df$REACTION_NAME), .con = con
    ))
    df <- dplyr::left_join(df, sub_lkp, by = "REACTION_NAME")
  }
  switch(mode,
    "miRNA" = df %>%
      dplyr::group_by(MIRNA_NAME) %>%
      dplyr::summarise(
        REACTION_COUNT = dplyr::n_distinct(REACTION_NAME),
        REACTION_NAME  = paste(sort(unique(REACTION_NAME)), collapse = ", "),
        GENE_COUNT     = dplyr::n_distinct(GENE_NAME),
        SUBSYSTEM      = paste(sort(unique(na.omit(SUBSYSTEM))), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(REACTION_COUNT)) %>%
      dplyr::mutate(MIRNA_NAME = link_info_fun(MIRNA_NAME, "mirna")),

    "reaction" = df %>%
      dplyr::group_by(REACTION_NAME) %>%
      dplyr::summarise(
        MIRNA_COUNT = dplyr::n_distinct(MIRNA_NAME),
        MIRNA_NAME  = paste(sort(unique(MIRNA_NAME)), collapse = ", "),
        GENE_COUNT  = dplyr::n_distinct(GENE_NAME),
        HUMAN_ID    = dplyr::first(na.omit(HUMAN_ID)),
        SUBSYSTEM   = dplyr::first(na.omit(SUBSYSTEM)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(MIRNA_COUNT)),

    "miRNA and reaction" = df %>%
      dplyr::group_by(MIRNA_NAME, REACTION_NAME) %>%
      dplyr::summarise(
        GENE_COUNT = dplyr::n_distinct(GENE_NAME),
        SUBSYSTEM  = dplyr::first(na.omit(SUBSYSTEM)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(GENE_COUNT)) %>%
      dplyr::mutate(MIRNA_NAME = link_info_fun(MIRNA_NAME, "mirna")),

    "subsystem" = {
      df %>%
        dplyr::group_by(SUBSYSTEM) %>%
        dplyr::summarise(
          REACTION_COUNT = dplyr::n_distinct(REACTION_NAME),
          MIRNA_COUNT    = dplyr::n_distinct(MIRNA_NAME),
          GENE_COUNT     = dplyr::n_distinct(GENE_NAME),
          .groups = "drop"
        ) %>%
        dplyr::filter(!is.na(SUBSYSTEM)) %>%
        dplyr::arrange(dplyr::desc(REACTION_COUNT))
    },

    "category" = {
      df %>%
        dplyr::mutate(CATEGORY = SUBSYSTEM_CATEGORY_MAP[SUBSYSTEM]) %>%
        dplyr::filter(!is.na(CATEGORY)) %>%
        dplyr::group_by(CATEGORY) %>%
        dplyr::summarise(
          SUBSYSTEM_COUNT = dplyr::n_distinct(SUBSYSTEM),
          REACTION_COUNT  = dplyr::n_distinct(REACTION_NAME),
          MIRNA_COUNT     = dplyr::n_distinct(MIRNA_NAME),
          GENE_COUNT      = dplyr::n_distinct(GENE_NAME),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(REACTION_COUNT))
    },

    # "None" – flat with HTML links
    df %>%
      dplyr::mutate(
        MIRNA_NAME = link_info_fun(MIRNA_NAME, "mirna"),
        GENE_NAME  = link_info_fun(GENE_NAME,  "gene")
      )
  )
})

mr_data_ready <- reactive({
  tryCatch({ filtered_mr_base(); TRUE }, error = function(e) FALSE)
})

output$count_mr <- renderUI({
  req(filtered_mr_base())
  render_count_bar(filtered_mr_base(), list(
    "miRNA"      = "MIRNA_NAME",
    "genes"      = "GENE_NAME",
    "reactions"  = "REACTION_NAME",
    "subsystems" = "SUBSYSTEM"
  ))
})

output$result_table_mr <- renderDT({
  req(filtered_mr())
  render_standard_dt(filtered_mr())
}, server = TRUE)

output$download_mr_csv <- make_csv_download("mirna_reaction", filtered_mr)


mra_plot_module  <- plotVisualizationServer("plot_mr_net",
  get_data_fn = function() filtered_mr_base(),
  type_network     = "REACTION_NAME",
  x1 = "MIRNA_NAME", x2 = "REACTION_NAME",
  extra_graph_args = reactive(list(
    color_essentiality  = TRUE,
    mirna_color         = "#FF6666",
    reaction_color_by   = input$mr_net_color_by %||% "subsystem"
  )),
  physics_input    = reactive(input$mr_network_layout %||% "forceAtlas2Based"),
  mode             = "network",
  external_trigger = reactive(input$obtain_mr)
)
mra_plot_data   <- mra_plot_module$data
observeEvent(mra_plot_module$net_png(), {
  png <- mra_plot_module$net_png()
  if (!is.null(png) && nzchar(png)) my_data$mr_net_png <- png
})
mra_threshold  <- mra_plot_module$filter_edge

# Shared palette — all reactions, no ESS_FRAC filtering.
mr_subsystem_palette <- reactive({
  df       <- filtered_mr_base()
  color_by <- input$mr_net_color_by %||% "subsystem"
  if (is.null(df) || nrow(df) == 0) return(character(0))
  if (color_by == "category" && !"CATEGORY" %in% colnames(df) && "SUBSYSTEM" %in% colnames(df))
    df$CATEGORY <- SUBSYSTEM_CATEGORY_MAP[df$SUBSYSTEM]
  col_col    <- if (color_by == "category" && "CATEGORY" %in% colnames(df)) "CATEGORY" else "SUBSYSTEM"
  groups_vec <- sort(unique(df[[col_col]][!is.na(df[[col_col]]) & df[[col_col]] != ""]))
  if (length(groups_vec) == 0) return(character(0))
  n <- length(groups_vec)
  pal <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Set3")[seq_len(n)]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
  setNames(pal, groups_vec)
})

# ── Click-link modals (condivisi da tutti i tab) ──────────
observeEvent(input$clicked_mirna,   { req(input$clicked_mirna);   show_mirna_modal(input$clicked_mirna,   con) })
observeEvent(input$clicked_disease, { req(input$clicked_disease); show_disease_modal(input$clicked_disease, con) })
observeEvent(input$clicked_pubmed,  { req(input$clicked_pubmed);  show_pubmed_modal(input$clicked_pubmed,  con) })
observeEvent(input$clicked_gene,    { req(input$clicked_gene);    show_gene_modal(input$clicked_gene,    con) })

# ── Legenda subsystem + essenzialità ─────────────────────
output$subsystem_legend_mr <- renderUI({
  df <- mra_plot_data()
  req(nrow(df) > 0)

  legend_blocks <- list()

  color_by <- input$mr_net_color_by %||% "subsystem"

  # derive CATEGORY if needed
  if (color_by == "category" && !"CATEGORY" %in% colnames(df) && "SUBSYSTEM" %in% colnames(df))
    df$CATEGORY <- SUBSYSTEM_CATEGORY_MAP[df$SUBSYSTEM]

  group_col   <- if (color_by == "category" && "CATEGORY" %in% colnames(df)) "CATEGORY" else "SUBSYSTEM"
  group_label <- if (group_col == "CATEGORY") "Reaction category" else "Reaction subsystem"

  if (group_col %in% colnames(df)) {
    shared_pal <- mr_subsystem_palette()
    groups_vec <- sort(unique(df[[group_col]][!is.na(df[[group_col]]) & df[[group_col]] != ""]))
    # Filter to only groups visible in the (possibly threshold-filtered) df
    groups_vec <- intersect(groups_vec, names(shared_pal))

    if (length(groups_vec) > 0) {
      pal <- shared_pal[groups_vec]

      sub_items <- mapply(function(sub, col) {
        tags$div(
          style = "display:flex; align-items:center; gap:6px; margin-bottom:6px;",
          tags$div(style = sprintf(
            "width:13px; height:13px; background:%s; border:1px solid rgba(0,0,0,0.2); border-radius:3px; flex-shrink:0;", col)),
          sub
        )
      }, groups_vec, pal, SIMPLIFY = FALSE)

      legend_blocks[["subsystem"]] <- tags$div(
        style = "margin-bottom:8px;",
        tags$p(
          style = "font-weight:600; font-size:12px; color:#555; margin:0 0 6px 0; text-transform:uppercase; letter-spacing:0.4px;",
          icon("square", style = "font-size:10px; margin-right:5px;"),
          group_label
        ),
        sub_items
      )
    }
  }

  ess_levels <- list(
    list(label = "Essential  (no 'or' in GPR – no alternative isoenzyme)", color = "#B22222"),
    list(label = "Non-essential  (alternative isoenzymes present)",          color = "#6699CC")
  )
  ess_items <- lapply(ess_levels, function(e) {
    tags$span(
      style = "display:inline-flex; align-items:center; margin-right:20px; margin-bottom:4px;",
      tags$span(style = sprintf(
        "display:inline-block; width:30px; height:4px; background-color:%s; margin-right:7px; flex-shrink:0;", e$color)),
      tags$span(e$label, style = "font-size:12px; color:#444;")
    )
  })
  legend_blocks[["essentiality"]] <- tags$div(
    tags$p(
      style = "font-weight:600; font-size:12px; color:#555; margin:0 0 6px 0; text-transform:uppercase; letter-spacing:0.4px;",
      icon("minus", style = "font-size:10px; margin-right:5px;"),
      "Edge – miRNA essentiality on that reaction"
    ),
    tags$div(style = "display:flex; flex-wrap:wrap;", ess_items)
  )

  if (length(legend_blocks) == 0) return(NULL)
  tagList(tags$div(
    style = "margin-top:10px; padding:12px 16px; background:#fafafa; border:1px solid #e0e0e0; border-radius:6px;",
    legend_blocks
  ))
}) |> bindCache(mra_plot_data(), input$mr_net_color_by)

# =============================================================================
# Metabolic heatmap
# =============================================================================

.mr_hm_trigger <- reactive({
  list(
    search   = input$obtain_mr %||% 0L,
    group_by = input$metab_group_by %||% "subsystem"
  )
})

metab_matrix <- eventReactive(.mr_hm_trigger(), {
  shiny::validate(shiny::need(
    mr_data_ready(),
    "Run a miRNA–Reaction search first (Table summary tab)."
  ))
  df <- filtered_mr_base()
  shiny::validate(shiny::need(nrow(df) > 0, "No associations in current miRNA–Reaction search."))

  metric   <- "coverage"
  group_by <- input$metab_group_by %||% "subsystem"

  # SUBSYSTEM needed for subsystem/category grouping or exclusion filter
  if (!"SUBSYSTEM" %in% colnames(df) && "REACTION_NAME" %in% colnames(df) && nrow(df) > 0) {
    sub_lookup <- dbGetQuery(con, glue_sql(
      "SELECT NAME AS REACTION_NAME, MIN(SUBSYSTEM) AS SUBSYSTEM
         FROM REACTIONS WHERE NAME IN ({vals*})
         GROUP BY NAME",
      vals = unique(df$REACTION_NAME), .con = con
    ))
    df <- dplyr::left_join(df, sub_lookup, by = "REACTION_NAME")
  }

  # Category aggregation: remap SUBSYSTEM to category labels, then treat as subsystem
  if (group_by == "category") {
    df <- df %>%
      dplyr::mutate(CATEGORY = SUBSYSTEM_CATEGORY_MAP[SUBSYSTEM]) %>%
      dplyr::filter(!is.na(CATEGORY))
    shiny::validate(shiny::need(nrow(df) > 0, "No reactions matched any biological category."))
  }

  res_cm <- compute_metabolic_matrix(df, con, metric = metric, group_by = group_by)
  shiny::validate(shiny::need(
    !is.null(res_cm) && nrow(res_cm$mat) >= 1 && ncol(res_cm$mat) >= 1,
    "No data available for selected grouping and metric."
  ))

  list(
    mat       = res_cm$mat,
    pmat      = res_cm$pmat,
    adjpmat   = res_cm$adjpmat,
    metric    = metric,
    group_by  = group_by,
    col_label = if (group_by == "category") "Category" else "Subsystem"
  )
})

output$metab_heatmap <- renderPlotly({
  res <- metab_matrix()
  req(!is.null(res))

  mat        <- res$mat
  cluster_rc <- isTRUE(input$metab_cluster)

  p <- render_mirna_target_heatmap(
    mat,
    cluster_rows   = cluster_rc,
    cluster_cols   = cluster_rc,
    x1_label       = "miRNA",
    x2_label       = res$col_label,
    row_annotation = NULL,
    zero_white     = TRUE,
    transpose      = FALSE,
    color_palette  = input$metab_heatmap_palette %||% "YlGnBu"
  )
  p
})

output$metab_heatmap_table <- DT::renderDataTable(server = TRUE, {
  res <- metab_matrix()
  req(!is.null(res))

  metric_label <- "Reaction_coverage"

  df_long <- as.data.frame(as.table(res$mat))
  colnames(df_long) <- c("miRNA", res$col_label, metric_label)
  df_long <- df_long[df_long[[metric_label]] > 0, ]
  df_long <- df_long[order(-df_long[[metric_label]]), ]
  rownames(df_long) <- NULL

  render_standard_dt(df_long)
}) |> bindEvent(.mr_hm_trigger(), ignoreNULL = TRUE)

output$download_metab_heatmap_table <- downloadHandler(
  filename = function() paste0("metabolic_heatmap_", Sys.Date(), ".csv"),
  content  = function(file) {
    res <- metab_matrix()
    req(!is.null(res))
    metric_label <- switch(res$metric,
      "ess_frac"        = "Mean_ESS_FRAC",
      "coverage"        = "Reaction_coverage",
      "n_reactions"     = "N_reactions",
      "fold_enrichment" = "Fold_enrichment",
      "Value")
    df_long <- as.data.frame(as.table(res$mat))
    colnames(df_long) <- c("miRNA", res$col_label, metric_label)
    df_long <- df_long[df_long[[metric_label]] > 0, ]
    df_long <- df_long[order(-df_long[[metric_label]]), ]
    write.csv(df_long, file, row.names = FALSE)
  }
)

# ── Subsystem/Category barplot (Overview) ─────────────────────────────────────
output$mr_barplot <- renderPlotly({
  shiny::validate(shiny::need(mr_data_ready(), "Run a miRNA–Reaction search first."))
  df <- filtered_mr_base()
  shiny::validate(shiny::need(nrow(df) > 0, "No associations to display."))

  if (!"SUBSYSTEM" %in% colnames(df)) {
    sub_lookup <- dbGetQuery(con, glue_sql(
      "SELECT NAME AS REACTION_NAME, MIN(SUBSYSTEM) AS SUBSYSTEM FROM REACTIONS
        WHERE NAME IN ({vals*}) GROUP BY NAME",
      vals = unique(df$REACTION_NAME %||% df$NAME), .con = con
    ))
    join_col <- if ("REACTION_NAME" %in% names(df)) "REACTION_NAME" else "NAME"
    df <- dplyr::left_join(df, sub_lookup, by = setNames("REACTION_NAME", join_col))
  }

  group_by_bar <- input$mr_barplot_group_by %||% "subsystem"

  # Query totals using real SUBSYSTEM names (before any category remap)
  subsystems_present <- unique(df$SUBSYSTEM[!is.na(df$SUBSYSTEM)])
  tot_sub <- if (length(subsystems_present) > 0) {
    tryCatch(
      dbGetQuery(con, glue_sql(
        "SELECT SUBSYSTEM, COUNT(DISTINCT REACTION_ID) AS total
           FROM {`REACTIONS_TBL`}
          WHERE SUBSYSTEM IN ({subs*})
          GROUP BY SUBSYSTEM",
        subs = subsystems_present, .con = con
      )),
      error = function(e) data.frame(SUBSYSTEM = character(0), total = integer(0))
    )
  } else data.frame(SUBSYSTEM = character(0), total = integer(0))

  if (group_by_bar == "category") {
    df <- df %>%
      dplyr::mutate(SUBSYSTEM = SUBSYSTEM_CATEGORY_MAP[SUBSYSTEM]) %>%
      dplyr::filter(!is.na(SUBSYSTEM))
    shiny::validate(shiny::need(nrow(df) > 0, "No reactions matched any biological category."))
    # Aggregate sub-totals to category totals
    tot_rxn <- tot_sub %>%
      dplyr::mutate(SUBSYSTEM = SUBSYSTEM_CATEGORY_MAP[SUBSYSTEM]) %>%
      dplyr::filter(!is.na(SUBSYSTEM)) %>%
      dplyr::group_by(SUBSYSTEM) %>%
      dplyr::summarise(total = sum(total, na.rm = TRUE), .groups = "drop")
  } else {
    tot_rxn <- tot_sub
  }

  agg <- df %>%
    dplyr::group_by(SUBSYSTEM) %>%
    dplyr::summarise(
      N_reactions = dplyr::n_distinct(REACTION_NAME),
      N_mirnas    = dplyr::n_distinct(MIRNA_NAME),
      .groups = "drop"
    ) %>%
    dplyr::left_join(tot_rxn, by = "SUBSYSTEM") %>%
    dplyr::mutate(
      coverage = ifelse(!is.na(total) & total > 0, N_reactions / total, NA_real_)
    ) %>%
    dplyr::arrange(dplyr::coalesce(coverage, 0)) %>%
    dplyr::mutate(
      SUBSYSTEM = factor(SUBSYSTEM, levels = SUBSYSTEM),
      tooltip   = paste0("<b>", SUBSYSTEM, "</b><br>",
                         "Coverage: ", ifelse(!is.na(coverage), scales::percent(coverage, accuracy = 0.1), "—"), "<br>",
                         "Reactions targeted: ", N_reactions, " / ", ifelse(is.na(total), "?", total), "<br>",
                         "miRNAs: ", N_mirnas)
    )

  group_label <- if (group_by_bar == "category") "Category" else "Subsystem"

  plot_ly(agg,
    x = ~coverage, y = ~SUBSYSTEM, type = "bar", orientation = "h",
    marker = list(
      color     = ~N_mirnas,
      colorscale = list(c(0, "#d4eac8"), c(1, "#2a7a3b")),
      showscale = TRUE,
      colorbar  = list(title = "N miRNAs", len = 0.5)
    ),
    text = ~tooltip, hoverinfo = "text"
  ) %>%
  layout(
    xaxis  = list(title = "Coverage (reactions targeted / total)", range = c(0, 1),
                  tickformat = ".0%", zeroline = FALSE, gridcolor = "#e8e8e8"),
    yaxis  = list(title = "", tickfont = list(size = 10)),
    paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa",
    margin = list(l = 220, r = 80, t = 30, b = 50),
    showlegend = FALSE
  ) %>%
  plotly_clean_config() %>%
  htmlwidgets::onRender("
    function(el, x) {
      Plotly.toImage(el, {format:'png', width:2580, height:1680}).then(function(url) {
        Shiny.setInputValue('mr_barplot_png', url, {priority:'event'});
      });
    }
  ")
}) |> bindEvent(input$obtain_mr, input$mr_barplot_group_by, ignoreNULL = TRUE)

observeEvent(input$mr_barplot_png, {
  png <- input$mr_barplot_png
  if (!is.null(png) && nzchar(png)) my_data$mr_barplot_png <- png
})

# =============================================================================
# Co-targeting analysis (miRNA–Reaction / subsystems)
# =============================================================================

cotarget_mr_df <- reactive({
  shiny::validate(shiny::need(
    mr_data_ready(),
    "Run a miRNA–Reaction search first (Table summary tab)."
  ))
  df <- filtered_mr_base()
  shiny::validate(shiny::need(nrow(df) > 0, "No associations in current miRNA–Reaction search."))

  df$mirna_ <- df$MIRNA_NAME

  result <- unique(data.frame(
    MIRNA_NAME = df$mirna_,
    SUBSYSTEM  = df$SUBSYSTEM,
    stringsAsFactors = FALSE
  ))
  result <- result[!is.na(result$SUBSYSTEM) & result$SUBSYSTEM != "", ]

  shiny::validate(shiny::need(nrow(result) > 0,
                "No subsystem associations found. Check that the reaction data includes subsystem info."))
  result
})

# =============================================================================
# miRNA Essentiality analysis
# =============================================================================

mr_essentiality_df <- reactive({
  shiny::validate(shiny::need(
    mr_data_ready(),
    "Run a miRNA-Reaction search first (Table tab)."
  ))
  df <- filtered_mr_base()
  if (!"ESS_FRAC"     %in% names(df)) df <- compute_essentiality_fraction(con, df)
  if (!"IS_ESSENTIAL" %in% names(df)) df <- compute_is_essential(df)

  # One row per (miRNA, reaction-key) to avoid counting gene-level duplicates.
  # RXN_KEY = HUMAN_ID (the true unique reaction ID); falls back to REACTION_NAME
  # for rows without a HUMAN_ID.
  rxn_key_col <- if ("RXN_KEY" %in% names(df)) "RXN_KEY" else "REACTION_NAME"
  pairs <- dplyr::distinct(df, MIRNA_NAME, .data[[rxn_key_col]], .keep_all = TRUE)
  pairs$mirna_ <- pairs$MIRNA_NAME

  # IS_ESSENTIAL is the GPR-boolean column (same source as table); guaranteed present above
  pairs %>%
    dplyr::group_by(MIRNA_NAME = mirna_) %>%
    dplyr::summarise(
      N_reactions   = dplyr::n_distinct(.data[[rxn_key_col]]),
      N_essential   = sum(IS_ESSENTIAL == TRUE, na.rm = TRUE),
      N_partial     = N_reactions - N_essential,
      Mean_ESS_FRAC = round(mean(ESS_FRAC, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Pct_essential = round(100 * N_essential / N_reactions, 1)) %>%
    dplyr::arrange(dplyr::desc(N_reactions))
}) |> bindEvent(input$run_mr_essentiality, ignoreNULL = TRUE)

output$mr_essentiality_plot <- renderPlotly({
  df <- mr_essentiality_df()
  shiny::validate(shiny::need(nrow(df) > 0, "No data available."))

  df$MIRNA_NAME <- factor(df$MIRNA_NAME, levels = rev(df$MIRNA_NAME))

  plot_ly(df) %>%
    add_trace(
      x = ~N_essential, y = ~MIRNA_NAME, type = "bar", orientation = "h",
      name = "Essential", marker = list(color = "#B22222"),
      text = ~paste0("Essential: ", N_essential, " (", Pct_essential, "%)"),
      hoverinfo = "text"
    ) %>%
    add_trace(
      x = ~N_partial, y = ~MIRNA_NAME, type = "bar", orientation = "h",
      name = "Partial", marker = list(color = "#6699CC"),
      text = ~paste0("Partial: ", N_partial),
      hoverinfo = "text"
    ) %>%
    layout(
      barmode = "stack",
      xaxis   = list(title = "Number of reactions", zeroline = FALSE, gridcolor = "#e8e8e8"),
      yaxis   = list(title = "", tickfont = list(size = 10)),
      legend  = list(orientation = "h", x = 0, y = -0.15),
      paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa",
      margin = list(l = 200, r = 20, t = 20, b = 60)
    ) %>%
    plotly_clean_config() %>%
    htmlwidgets::onRender("function(el) {
      Plotly.toImage(el, {format:'png', width:2580, height:1440}).then(function(url) {
        Shiny.setInputValue('mr_ess_plot_png', url, {priority:'event'});
      });
    }")
}) |> bindEvent(input$run_mr_essentiality, mr_essentiality_df(), ignoreNULL = TRUE)

observeEvent(input$mr_ess_plot_png, {
  png <- input$mr_ess_plot_png
  if (!is.null(png) && nzchar(png)) my_data$mr_ess_plot_png <- png
})

output$mr_essentiality_table <- DT::renderDataTable({
  df <- mr_essentiality_df()
  shiny::validate(shiny::need(nrow(df) > 0, "No data available."))
  DT::datatable(
    df,
    rownames  = FALSE,
    filter    = "top",
    extensions = "Buttons",
    options = list(
      pageLength = 20,
      dom = "Bfrtip",
      buttons = list(
        list(extend = "csv",   exportOptions = list(modifier = list(page = "all"), columns = ":all")),
        list(extend = "excel", exportOptions = list(modifier = list(page = "all"), columns = ":all"))
      )
    )
  )
}) |> bindEvent(input$run_mr_essentiality, mr_essentiality_df(), ignoreNULL = TRUE)

output$download_mr_essentiality_csv <- make_csv_download("mirna_essentiality", mr_essentiality_df)

# =============================================================================
# Co-targeting
# =============================================================================

cotarget_mr_mat <- reactive({
  df        <- cotarget_mr_df()
  node_type <- input$mr_cotarget_node %||% "mirna_mirna"
  MAX_DIM   <- 100L
  if (node_type == "subsystem_subsystem") {
    n <- length(unique(df$SUBSYSTEM))
    shiny::validate(shiny::need(n >= 2, "Need at least 2 subsystems. Broaden miRNA–Reaction filters."))
    shiny::validate(shiny::need(n <= MAX_DIM, paste0("Too many nodes (", n, "×", n, " matrix). Limit: ", MAX_DIM, "×", MAX_DIM, ". Raise the min-degree or reduce the selection.")))
    compute_jaccard_matrix(df, col_nodes = "SUBSYSTEM", col_shared = "MIRNA_NAME")
  } else {
    n <- length(unique(df$MIRNA_NAME))
    shiny::validate(shiny::need(n >= 2, "Need at least 2 miRNAs. Broaden miRNA–Reaction filters."))
    shiny::validate(shiny::need(n <= MAX_DIM, paste0("Too many nodes (", n, "×", n, " matrix). Limit: ", MAX_DIM, "×", MAX_DIM, ". Raise the min-degree or reduce the miRNA selection.")))
    compute_jaccard_matrix(df, col_nodes = "MIRNA_NAME", col_shared = "SUBSYSTEM")
  }
})

output$cotarget_mr_heatmap <- renderPlotly({
  mat       <- cotarget_mr_mat()
  req(!is.null(mat))
  node_type <- input$mr_cotarget_node %||% "mirna_mirna"
  node_lbl  <- if (node_type == "subsystem_subsystem") "Subsystem" else "miRNA"
  render_cotarget_heatmap(mat,
    cluster       = isTRUE(input$cotarget_mr_cluster),
    node_label    = node_lbl,
    color_palette = input$cotarget_mr_palette   %||% "YlOrRd",
    mask_zero     = isTRUE(input$cotarget_mr_mask_zero),
    mask_diag     = isTRUE(input$cotarget_mr_mask_diag)) %>%
  htmlwidgets::onRender("function(el) {
    Plotly.toImage(el, {format:'png', width:2580, height:2100}).then(function(url) {
      Shiny.setInputValue('mr_hm_raw_png', url, {priority:'event'});
    });
  }")
})

store_png_variant("mr_hm_raw_png", my_data, input,
  node_input  = "mr_cotarget_node",
  variant_map = list(mirna_mirna        = "mr_hm_mirna_png",
                     subsystem_subsystem = "mr_hm_subsystem_png"))

output$cotarget_mr_count <- renderUI({
  mat       <- cotarget_mr_mat()
  req(!is.null(mat))
  node_type <- input$mr_cotarget_node %||% "mirna_mirna"
  node_lbl  <- if (node_type == "subsystem_subsystem") "subsystems" else "miRNAs"
  tbl       <- cotarget_pairs_table(mat, threshold = 0)
  n_pairs   <- if (is.null(tbl)) 0L else nrow(tbl)
  tags$p(
    style = "color:#555; font-style:italic; margin-bottom:6px;",
    HTML(paste0(
      "<b>", nrow(mat), "</b> ", node_lbl, " &nbsp;|&nbsp; ",
      "<b>", n_pairs,  "</b> pairs"
    ))
  )
})

cotarget_mr_pairs_df <- reactive({
  mat       <- cotarget_mr_mat()
  req(!is.null(mat))
  node_type <- input$mr_cotarget_node %||% "mirna_mirna"
  if (node_type == "subsystem_subsystem") {
    cotarget_pairs_table(mat, threshold = 0,
                          df_source  = cotarget_mr_df(),
                          col_nodes  = "SUBSYSTEM",
                          col_shared = "MIRNA_NAME")
  } else {
    cotarget_pairs_table(mat, threshold = 0,
                          df_source  = cotarget_mr_df(),
                          col_nodes  = "MIRNA_NAME",
                          col_shared = "SUBSYSTEM")
  }
})

output$cotarget_mr_table <- DT::renderDataTable({
  tbl <- cotarget_mr_pairs_df()
  shiny::validate(shiny::need(!is.null(tbl), "No co-targeting pairs found."))
  render_standard_dt(tbl)
})



