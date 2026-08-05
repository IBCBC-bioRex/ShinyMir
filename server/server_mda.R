# === server_mda.R ===
# Tab miRNA-Disease: tabella, grafo, interazione con selezione

# Helper: compute 5 network metrics for a subgraph edge-list
# sub   = data.frame with MIRNA_NAME / DISEASE columns (flat, one row per edge)
# n_m   = total miRNAs in sample (fixed set size)
# n_d   = total diseases in sample (fixed or sampled set size)

# Populate DO category filter
observe({
  cats <- names(do_category_disease_map)
  choices <- c("All categories" = "", setNames(cats, cats))
  updateSelectInput(session, "do_category_filter", choices = choices, selected = "")
})

# Update disease_select choices when a DO category is selected
observe({
  cat_sel <- input$do_category_filter
  if (!is.null(cat_sel) && nzchar(cat_sel) && cat_sel %in% names(do_category_disease_map)) {
    avail <- sort(unique(do_category_disease_map[[cat_sel]]))
    updateSelectizeInput(session, "disease_select",
                         choices = avail, selected = character(0), server = FALSE)
  }
}) |> bindEvent(input$do_category_filter, ignoreNULL = TRUE, ignoreInit = TRUE)

# Shared reactive – single DB query for both legend, graph, and heatmap
mda_disease_map <- reactive({
  if (isTRUE(input$filter_diseaseancestor_checkbox) && length(input$disease_select) > 0)
    get_disease_descendant_map(input$disease_select, con)
  else NULL
})

output$mda_legend <- renderUI({
  anc_map <- mda_disease_map()
  if (!is.null(anc_map) && nrow(anc_map) > 0) {
    grp_unique <- sort(unique(na.omit(anc_map$ANCESTOR_GROUP)))
    n_grp      <- length(grp_unique)
    pal        <- GROUP_PALETTE
    grp_col    <- pal[((seq_len(n_grp) - 1L) %% length(pal)) + 1L]
    named_cols <- c(setNames("#FF6666", "miRNA"), setNames(grp_col, grp_unique))
    return(build_legend_ui(named_cols, dot_shapes = "miRNA"))
  }
  build_legend_ui(c("miRNA" = "#FF6666", "Disease" = "#6699FF"), dot_shapes = "miRNA")
})

upload_list_server("ul_mirna_select",   session, "mirna_select",   type = "mirna", valid_choices = mirna_choices)
upload_list_server("ul_disease_select", session, "disease_select",                 valid_choices = disease_choices)

# Single DB query – one row per (MIRNA, DISEASE, PUBMED_ID); runs only on Search click
filtered_mda_base <- eventReactive(input$obtain_mda, {
  on.exit(removeNotification("notify_mda"))
  df <- filtered_mda_fun(input, botton_update_mda_database, con = con,
                         for_plot = FALSE, my_data, group_override = "None")
  if (nrow(df) == 0) shiny::validate("No associations found with the current filters. Try broadening your selection.")
  df
})

# Table view – instant dplyr group-by, re-runs on group_mda change without hitting DB
filtered_mda <- reactive({
  df <- filtered_mda_base()
  mode <- input$group_mda %||% "None"
  switch(mode,
    "miRNA" = df %>%
      dplyr::group_by(MIRNA_NAME_clean) %>%
      dplyr::summarise(
        DISEASE_COUNT = dplyr::n_distinct(DISEASE_clean),
        DISEASE       = paste(sort(unique(DISEASE_clean)), collapse = ","),
        PUBMED_COUNT  = dplyr::n_distinct(PUBMED_ID_clean),
        DSI           = dplyr::first(DSI),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(DISEASE_COUNT)) %>%
      dplyr::mutate(MIRNA_NAME = link_info_fun(MIRNA_NAME_clean, "mirna")) %>%
      dplyr::select(MIRNA_NAME, DISEASE_COUNT, DISEASE, PUBMED_COUNT, DSI),
    "disease" = df %>%
      dplyr::group_by(DISEASE_clean) %>%
      dplyr::summarise(
        MIRNA_COUNT  = dplyr::n_distinct(MIRNA_NAME_clean),
        MIRNA_NAME   = paste(sort(unique(MIRNA_NAME_clean)), collapse = ","),
        PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(MIRNA_COUNT)) %>%
      dplyr::mutate(DISEASE = link_info_fun(DISEASE_clean, "disease")) %>%
      dplyr::select(DISEASE, MIRNA_COUNT, MIRNA_NAME, PUBMED_COUNT),
    "miRNA and disease" = df %>%
      dplyr::group_by(MIRNA_NAME_clean, DISEASE_clean) %>%
      dplyr::summarise(
        PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean),
        PUBMED_ID    = paste(sort(unique(PUBMED_ID_clean)), collapse = ","),
        DSI          = dplyr::first(DSI),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(PUBMED_COUNT)) %>%
      dplyr::mutate(
        MIRNA_NAME = link_info_fun(MIRNA_NAME_clean, "mirna"),
        DISEASE    = link_info_fun(DISEASE_clean,    "disease"),
        PUBMED_ID  = link_info_fun(PUBMED_ID,        "pubmed")
      ) %>%
      dplyr::select(MIRNA_NAME, DISEASE, PUBMED_COUNT, PUBMED_ID, DSI),
    df  # "None" – one row per article, HTML links already in base
  )
})

# Flat miRNA-disease pairs (one row per pair) – derived from base, no extra DB call
filtered_mda_flat <- reactive({
  df <- filtered_mda_base()
  df %>%
    dplyr::group_by(MIRNA_NAME_clean, DISEASE_clean) %>%
    dplyr::summarise(PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean), .groups = "drop") %>%
    dplyr::rename(MIRNA_NAME = MIRNA_NAME_clean, DISEASE = DISEASE_clean)
})

output$count_mda <- renderUI({
  req(filtered_mda_base())
  render_count_bar(filtered_mda_base(), list(
    "miRNA"    = "MIRNA_NAME_clean",
    "diseases" = "DISEASE_clean",
    "PubMed"   = "PUBMED_ID_clean"
  ))
})

output$result_table_mda <- renderDT({
  req(filtered_mda())
  render_standard_dt(filtered_mda())
}, server = TRUE)

output$download_mda_csv <- make_csv_download(
  "mirna_disease", filtered_mda, con = con,
  include_abstract = function() input$download_abstract_mda
)

mda_extra_graph <- reactive({
  dsi_map <- if (isTRUE(input$mda_color_nodes_by_dsi)) {
    tryCatch({
      df <- mda_dsi_df()
      if (!is.null(df) && nrow(df) > 0)
        setNames(df$DSI, df$MIRNA_NAME)
      else NULL
    }, error = function(e) NULL)
  } else NULL
  list(
    mirna_color       = "#FF6666",
    target_color      = "#6699FF",
    disease_group_map = mda_disease_map(),
    mirna_dsi_map     = dsi_map
  )
})

# Network tab module
mda_plot_module <- plotVisualizationServer("plot_md_net",
  get_data_fn      = function() filtered_mda_flat(),
  type_network     = "DISEASE",
  x1 = "MIRNA_NAME", x2 = "DISEASE",
  extra_graph_args = mda_extra_graph,
  mode             = "network",
  external_trigger = reactive(input$obtain_mda),
  physics_input    = reactive(input$mda_network_layout %||% "forceAtlas2Based")
)
observeEvent(mda_plot_module$net_png(), {
  png <- mda_plot_module$net_png()
  if (!is.null(png) && nzchar(png)) my_data$mda_net_png <- png
})


# Popola mirna_gene_select con i miRNA visibili nella tabella mda
observeEvent(input$result_table_mda_rows_all, {
  mirnas <- unique(filtered_mda_flat()$MIRNA_NAME)
  updateSelectizeInput(session, "mirna_gene_select", choices = mirnas, server = TRUE)
}, ignoreNULL = TRUE)

# =============================================================================
# Robustness analysis (miRNA–Disease)
# =============================================================================

rob_mda_data <- eventReactive(input$run_rob_mda, {
  shiny::validate(shiny::need(
    tryCatch({ filtered_mda_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA–Disease search first (Table tab)."
  ))

  withProgress(message = "Running robustness sweep...", value = 0, {
    setProgress(0.1, detail = "Building clauses")

    # Re-derive the same mirna/disease clauses used in filtered_mda_fun
    mirnas_id <- get_mirna_ids(input$mirna_select, con)
    mirna_clause <- if (length(input$mirna_select) > 0) {
      if (length(mirnas_id) == 0) DBI::SQL("AND 1=0")
      else glue_sql("AND A.MIRNA_ID IN ({vals*})", vals = mirnas_id, .con = con)
    } else if (isTRUE(input$filter_mirna_all_checkbox)) {
      DBI::SQL("")
    } else {
      DBI::SQL("AND 1=0")
    }

    disease_list <- select_disease_ancestors_fun(
      disease_selected = input$disease_select,
      filter_ancestor  = input$filter_diseaseancestor_checkbox,
      con   = con,
      depth = input$disease_descendant_depth %||% "all"
    )
    cat_sel      <- input$do_category_filter %||% ""
    all_in_scope <- if (nzchar(cat_sel) && cat_sel %in% names(do_category_disease_map))
                      do_category_disease_map[[cat_sel]]
                    else
                      disease_choices

    disease_clause <- if (length(disease_list) > 0) {
      disease_ids <- dbGetQuery(con,
        glue_sql("SELECT DISEASE_ID FROM DISEASES WHERE DISEASE IN ({vals*})",
                 vals = disease_list, .con = con))$DISEASE_ID
      if (length(disease_ids) > 0)
        glue_sql("AND C.DISEASE_ID IN ({vals*})", vals = disease_ids, .con = con)
      else
        DBI::SQL("AND 1=0")
    } else if (isTRUE(input$filter_disease_checkbox)) {
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

    human_clause_rob <- if (isTRUE(input$filter_human_mda) && isTRUE(HAS_IS_HUMAN))
      "AND CAST(A.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)"
    else
      ""

    setProgress(0.3, detail = "Sweeping thresholds")
    df <- compute_mda_robustness(
      con                = con,
      mirna_clause       = mirna_clause,
      disease_clause     = disease_clause,
      min_thresh         = input$rob_mda_min_thresh %||% 1L,
      max_thresh         = input$rob_mda_max_thresh %||% 20L,
      use_update         = botton_update_mda_database(),
      disease_restricted = length(disease_list) > 0,
      human_clause       = human_clause_rob
    )
    shiny::validate(shiny::need(nrow(df) > 0, "No associations found for this selection."))
    setProgress(1)
    df
  })
})

output$rob_mda_plot <- renderPlotly({
  df <- rob_mda_data()
  req(!is.null(df) && nrow(df) > 0)

  current_thresh <- as.integer(input$min_assoc_mda %||% 1L)

  p <- plot_ly(df, x = ~threshold) %>%
    add_lines(y = ~norm_base(n_edges),   name = "Edge retention",    line = list(color = "#2c7bb6", width = 2)) %>%
    add_lines(y = ~norm_base(n_mirna),   name = "miRNA retention",   line = list(color = "#d7191c", width = 2)) %>%
    add_lines(y = ~norm_base(n_disease), name = "Disease retention", line = list(color = "#1a9641", width = 2))

  # DSI mean ± SD on secondary Y axis
  if ("dsi_mean" %in% names(df) && any(!is.na(df$dsi_mean))) {
    df_dsi <- df[!is.na(df$dsi_mean), ]
    p <- p %>%
      add_lines(data = df_dsi, x = ~threshold, y = ~dsi_mean, name = "DSI mean",
                line = list(color = "#fdae61", width = 2, dash = "dash"),
                yaxis = "y", inherit = FALSE)
  }

  if (current_thresh >= min(df$threshold) && current_thresh <= max(df$threshold)) {
    p <- p %>% add_segments(
      x = current_thresh, xend = current_thresh, y = 0, yend = 1,
      line = list(color = "black", dash = "dash", width = 1),
      name = "Current threshold", showlegend = TRUE, inherit = FALSE
    )
  }

  p %>% layout(
    title  = list(text = "Network robustness vs evidence threshold", font = list(size = 14)),
    xaxis  = list(title = "Min evidence threshold", dtick = 1),
    yaxis  = list(title = "Fraction of baseline (normalised)", range = c(0, 1.05)),
    legend = list(orientation = "h", y = -0.2),
    hovermode = "x unified",
    paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa"
  ) %>% plotly_clean_config() %>%
  htmlwidgets::onRender("function(el) {
    Plotly.toImage(el, {format:'png', width:2580, height:1440}).then(function(url) {
      Shiny.setInputValue('mda_rob_png', url, {priority:'event'});
    });
  }")
}) |> bindEvent(input$run_rob_mda, ignoreNULL = TRUE)

store_png("mda_rob_png", my_data, input)

output$rob_mda_table <- DT::renderDataTable({
  df <- rob_mda_data()
  req(!is.null(df) && nrow(df) > 0)
  col_map <- c(
    threshold = "Threshold",
    n_edges   = "N edges",
    n_mirna   = "N miRNA",
    n_disease = "N diseases",
    dsi_mean  = "DSI mean",
    dsi_sd    = "DSI SD"
  )
  df <- df[, intersect(names(col_map), names(df)), drop = FALSE]
  names(df) <- col_map[names(df)]
  DT::datatable(df, rownames = FALSE,
                options = list(pageLength = 30, scrollX = TRUE, dom = "tip"))
}) |> bindEvent(input$run_rob_mda, ignoreNULL = TRUE)

rob_mda_cotarget_data <- eventReactive(input$run_rob_mda_cotarget, {
  shiny::validate(shiny::need(
    tryCatch({ filtered_mda_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA–Disease search first (Table tab)."
  ))
  withProgress(message = "Running co-targeting robustness...", value = 0.2, {
    df  <- filtered_mda_base()
    node_type <- input$rob_mda_cotarget_node %||% "mirna_mirna"
    if (node_type == "mirna_mirna") {
      col_nodes <- "MIRNA_NAME"; col_shared <- "DISEASE"
    } else {
      col_nodes <- "DISEASE"; col_shared <- "MIRNA_NAME"
    }
    m_col      <- if ("MIRNA_NAME_clean" %in% names(df)) "MIRNA_NAME_clean" else "MIRNA_NAME"
    d_col      <- if ("DISEASE_clean"    %in% names(df)) "DISEASE_clean"    else "DISEASE"
    obs_mirnas <- unique(strip_html(df[[m_col]]))
    obs_dis    <- unique(strip_html(df[[d_col]]))
    tbl <- if (isTRUE(botton_update_mda_database()))
      "MIRNAS_DISEASES_ARTICLES_update" else "MIRNAS_DISEASES_ARTICLES"
    human_sql_ct <- if (isTRUE(input$filter_human_mda) && isTRUE(HAS_IS_HUMAN))
      "AND CAST(A.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)"
    else
      ""
    pairs <- dbGetQuery(con, glue_sql(
      "SELECT M.MIRNA_PREMATURE AS MIRNA_NAME, C.DISEASE,
              COUNT(DISTINCT A.PUBMED_ID) AS N_PUB
       FROM {`tbl`} A
       JOIN MIRNAS   M ON A.MIRNA_ID   = M.MIRNA_ID
       JOIN DISEASES C ON A.DISEASE_ID = C.DISEASE_ID
       WHERE M.MIRNA_PREMATURE IN ({obs_mirnas*})
         AND C.DISEASE          IN ({obs_dis*})
         {human_sql_ct}
       GROUP BY M.MIRNA_PREMATURE, C.DISEASE",
      tbl = tbl, obs_mirnas = obs_mirnas, obs_dis = obs_dis,
      human_sql_ct = DBI::SQL(human_sql_ct), .con = con
    ))
    incProgress(0.6)
    compute_cotarget_robustness(
      pairs_df       = pairs,
      col_nodes      = col_nodes,
      col_shared     = col_shared,
      min_thresh     = input$rob_mda_min_thresh_ct %||% 1L,
      max_thresh     = input$rob_mda_max_thresh_ct %||% 20L,
      jaccard_cutoff = input$rob_mda_jaccard_cutoff %||% 0.1
    )
  })
})

output$rob_mda_cotarget_plot <- renderPlotly({
  res <- rob_mda_cotarget_data()
  df  <- res$sweep
  req(!is.null(df) && nrow(df) > 0)
  plot_ly(df, x = ~threshold) %>%
    add_lines(y = ~n_pairs,       name = "Co-targeting pairs", line = list(color = "#2c7bb6", width = 2)) %>%
    add_lines(y = ~mean_jaccard,  name = "Mean Jaccard",       line = list(color = "#fdae61", width = 2, dash = "dash"),
              yaxis = "y2") %>%
    layout(
      title  = list(text = "Co-targeting robustness vs evidence threshold", font = list(size = 14)),
      xaxis  = list(title = "Min evidence threshold", dtick = 1),
      yaxis  = list(title = "N co-targeting pairs", side = "left"),
      yaxis2 = list(title = "Mean Jaccard", overlaying = "y", side = "right",
                    range = c(0, 1), showgrid = FALSE),
      legend = list(orientation = "h", y = -0.2),
      hovermode = "x unified",
      paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa"
    ) %>% plotly_clean_config() %>%
    htmlwidgets::onRender("function(el) {
      Plotly.toImage(el, {format:'png', width:2580, height:1440}).then(function(url) {
        Shiny.setInputValue('mda_rob_ct_raw_png', url, {priority:'event'});
      });
    }")
}) |> bindEvent(input$run_rob_mda_cotarget, ignoreNULL = TRUE)

store_png_variant("mda_rob_ct_raw_png", my_data, input,
  node_input  = "rob_mda_cotarget_node",
  variant_map = list(mirna_mirna    = "mda_rob_ct_mirna_png",
                     disease_disease = "mda_rob_ct_disease_png"))

output$rob_mda_stability_table <- DT::renderDataTable({
  res <- rob_mda_cotarget_data()
  df  <- res$stability
  req(!is.null(df) && nrow(df) > 0)
  node_type <- isolate(input$rob_mda_cotarget_node) %||% "mirna_mirna"
  lbl_a <- if (node_type == "mirna_mirna") "miRNA A" else "Disease A"
  lbl_b <- if (node_type == "mirna_mirna") "miRNA B" else "Disease B"
  if (ncol(df) == 4) names(df) <- c(lbl_a, lbl_b, "N thresholds", "Stability score")
  DT::datatable(df, rownames = FALSE,
                options = list(pageLength = 20, scrollX = TRUE, dom = "frtip"))
}) |> bindEvent(input$run_rob_mda_cotarget, ignoreNULL = TRUE)

# =============================================================================
# Co-targeting analysis (miRNA–Disease)
# =============================================================================

cotarget_mda_df <- eventReactive(input$obtain_mda, {
  shiny::validate(shiny::need(
    tryCatch({ filtered_mda_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA–Disease search first (Table summary tab)."
  ))
  df <- filtered_mda_base()
  shiny::validate(shiny::need(nrow(df) > 0, "No associations in current miRNA–Disease search."))

  mirna_col   <- if ("MIRNA_NAME_clean" %in% names(df)) "MIRNA_NAME_clean" else "MIRNA_NAME"
  disease_col <- if ("DISEASE_clean"    %in% names(df)) "DISEASE_clean"    else "DISEASE"

  df$mirna_   <- strip_html(df[[mirna_col]])
  df$disease_ <- strip_html(df[[disease_col]])

  # Min associations per miRNA-disease pair filter
  min_assoc <- max(1L, as.integer(input$min_assoc_mda %||% 1L))
  if (min_assoc > 1L) {
    pair_counts <- df %>%
      dplyr::group_by(mirna_, disease_) %>%
      dplyr::summarise(
        n = if ("PUBMED_ID" %in% names(df))
              length(split_unique(PUBMED_ID))
            else
              dplyr::n(),
        .groups = "drop"
      )
    keep_pairs <- pair_counts[pair_counts$n >= min_assoc, c("mirna_", "disease_")]
    df <- dplyr::semi_join(df, keep_pairs, by = c("mirna_", "disease_"))
    shiny::validate(shiny::need(nrow(df) > 0,
                  "No pairs survive the minimum associations filter. Lower the threshold."))
  }

  result <- unique(data.frame(
    MIRNA_NAME = df$mirna_,
    DISEASE    = df$disease_,
    stringsAsFactors = FALSE
  ))
  shiny::validate(shiny::need(nrow(result) > 0, "No miRNA–disease pairs found."))
  result
})

cotarget_mda_mat <- reactive({
  df        <- cotarget_mda_df()
  node_type <- input$mda_cotarget_node %||% "mirna_mirna"
  MAX_DIM   <- 100L
  if (node_type == "disease_disease") {
    n <- length(unique(df$DISEASE))
    shiny::validate(shiny::need(n >= 2, "Need at least 2 diseases. Broaden miRNA–Disease filters."))
    shiny::validate(shiny::need(n <= MAX_DIM, paste0("Too many nodes (", n, "×", n, " matrix). Limit: ", MAX_DIM, "×", MAX_DIM, ". Raise the min-degree or reduce the disease selection.")))
    compute_jaccard_matrix(df, col_nodes = "DISEASE", col_shared = "MIRNA_NAME")
  } else {
    n <- length(unique(df$MIRNA_NAME))
    shiny::validate(shiny::need(n >= 2, "Need at least 2 miRNAs. Broaden miRNA–Disease filters."))
    shiny::validate(shiny::need(n <= MAX_DIM, paste0("Too many nodes (", n, "×", n, " matrix). Limit: ", MAX_DIM, "×", MAX_DIM, ". Raise the min-degree or reduce the miRNA selection.")))
    compute_jaccard_matrix(df, col_nodes = "MIRNA_NAME", col_shared = "DISEASE")
  }
})

output$cotarget_mda_heatmap <- renderPlotly({
  mat       <- cotarget_mda_mat()
  req(!is.null(mat))
  node_type <- input$mda_cotarget_node %||% "mirna_mirna"
  node_lbl  <- if (node_type == "disease_disease") "Disease" else "miRNA"
  render_cotarget_heatmap(mat,
    cluster       = isTRUE(input$cotarget_mda_cluster),
    node_label    = node_lbl,
    color_palette = input$cotarget_mda_palette   %||% "YlOrRd",
    mask_zero     = isTRUE(input$cotarget_mda_mask_zero),
    mask_diag     = isTRUE(input$cotarget_mda_mask_diag)) %>%
  htmlwidgets::onRender("function(el) {
    Plotly.toImage(el, {format:'png', width:2580, height:2100}).then(function(url) {
      Shiny.setInputValue('cotarget_mda_hm_raw_png', url, {priority:'event'});
    });
  }")
}) |> bindCache(
  cotarget_mda_mat(),
  input$mda_cotarget_node,
  input$cotarget_mda_cluster,
  input$cotarget_mda_palette,
  input$cotarget_mda_mask_zero,
  input$cotarget_mda_mask_diag
)

store_png_variant("cotarget_mda_hm_raw_png", my_data, input,
  node_input   = "mda_cotarget_node",
  variant_map  = list(mirna_mirna    = "cotarget_mda_mirna_png",
                      disease_disease = "cotarget_mda_disease_png"))

output$cotarget_mda_count <- renderUI({
  mat       <- cotarget_mda_mat()
  req(!is.null(mat))
  node_type <- input$mda_cotarget_node %||% "mirna_mirna"
  node_lbl  <- if (node_type == "disease_disease") "diseases" else "miRNAs"
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

cotarget_mda_pairs_df <- reactive({
  mat       <- cotarget_mda_mat()
  req(!is.null(mat))
  node_type <- input$mda_cotarget_node %||% "mirna_mirna"
  if (node_type == "disease_disease") {
    cotarget_pairs_table(mat, threshold = 0,
                          df_source  = cotarget_mda_df(),
                          col_nodes  = "DISEASE",
                          col_shared = "MIRNA_NAME")
  } else {
    cotarget_pairs_table(mat, threshold = 0,
                          df_source  = cotarget_mda_df(),
                          col_nodes  = "MIRNA_NAME",
                          col_shared = "DISEASE")
  }
})

output$cotarget_mda_table <- DT::renderDataTable({
  tbl <- cotarget_mda_pairs_df()
  shiny::validate(shiny::need(!is.null(tbl), "No co-targeting pairs found."))
  render_standard_dt(tbl)
})




# DSI specificity profile — used for node colouring in network view
mda_dsi_df <- reactive({
  shiny::validate(shiny::need(
    tryCatch({ filtered_mda_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA-Disease search first (Table tab)."
  ))
  df <- filtered_mda_base()
  shiny::validate(shiny::need("DSI" %in% names(df), "DSI not available. Re-run Table search."))
  unique(df[, c("MIRNA_NAME_clean", "DSI"), drop = FALSE]) %>%
    dplyr::rename(MIRNA_NAME = MIRNA_NAME_clean) %>%
    dplyr::filter(!is.na(DSI))
})
