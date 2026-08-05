# === server_mga.R ===
# Tab miRNA-Gene: tabella, grafo, propagazione scelte verso tab successivi

MGA_MIRNA_COLOR <- "#FF6666"

output$mga_legend <- renderUI({
  mirna_cols <- c("miRNA" = MGA_MIRNA_COLOR)

  gene_color <- "#6699FF"
  gene_cols  <- c("Gene" = gene_color)

  src <- input$mga_pathway_source %||% "reactome"
  has_pw_filter <- if (src == "gmt") {
    length(input$mga_gmt_set_select) > 0 && !is.null(gmt_data())
  } else {
    length(as.integer(input$mga_pathway_l2)) > 0 || length(input$mga_pathway_select) > 0
  }

  if (has_pw_filter) {
    if (src == "gmt") {
      sel_sets <- input$mga_gmt_set_select
      pw_col   <- setNames(
        GROUP_PALETTE[((seq_along(sel_sets) - 1L) %% length(GROUP_PALETTE)) + 1L], sel_sets)
    } else {
      sel_names <- if (length(input$mga_pathway_l2) > 0) {
        anchor_ids <- as.integer(input$mga_pathway_l2)
        tryCatch(
          dbGetQuery(con, glue_sql(
            "SELECT DISTINCT PATHWAY_NAME FROM ONTOLOGY_PATHWAYS WHERE PATHWAY_ID IN ({ids*}) ORDER BY PATHWAY_NAME",
            ids = anchor_ids, .con = con
          ))$PATHWAY_NAME,
          error = function(e) character(0)
        )
      } else {
        sort(input$mga_pathway_select %||% character(0))
      }
      pw_col <- setNames(
        GROUP_PALETTE[((seq_along(sel_names) - 1L) %% length(GROUP_PALETTE)) + 1L], sel_names)
    }
    named_cols <- c(mirna_cols, pw_col, c("Mixed gene sets" = "#DAA520"))
  } else {
    named_cols <- c(mirna_cols, gene_cols)
  }

  build_legend_ui(named_cols, dot_shapes = names(mirna_cols))
})

updateSelectizeInput(session, "gene_search_select", choices = total_genes, server = TRUE)

observeEvent(input$mga_pathway_source, {
  updateSelectizeInput(session, "mga_pathway_select", choices = pathway_choices, server = TRUE)
}, ignoreNULL = FALSE)

# ── GMT: uses shared gmt_data reactive (upload in Pathway Annotation tab) ────
# When GMT changes, repopulate the set selector in the MGA network panel
observeEvent(gmt_data(), {
  gmt <- gmt_data()
  if (is.null(gmt)) {
    updateSelectizeInput(session, "mga_gmt_set_select",
      choices = character(0), selected = character(0), server = TRUE)
  } else {
    updateSelectizeInput(session, "mga_gmt_set_select",
      choices = names(gmt), selected = character(0), server = TRUE)
  }
}, ignoreNULL = FALSE)

output$mga_gmt_status_banner <- renderUI({
  gmt <- gmt_data()
  if (is.null(gmt))
    tags$span(style = "color:#888;", "No GMT loaded.")
  else
    tags$span(style = "color:#1a5276;",
              icon("check-circle"), paste0(length(gmt), " gene set(s) loaded."))
})

observeEvent(input$mga_use_3p5p_scratch, {
  choices <- if (isTRUE(input$mga_use_3p5p_scratch)) mirna_choices_3p5p else mirna_choices
  updateSelectizeInput(session, "mirna_gene_select_scratch",
                       choices = choices, selected = character(0), server = TRUE)
}, ignoreNULL = FALSE)

upload_list_server("ul_mirna_gene_scratch", session, "mirna_gene_select_scratch", type = "mirna",
  valid_choices = function() if (isTRUE(input$mga_use_3p5p_scratch)) mirna_choices_3p5p else mirna_choices)
upload_list_server("ul_gene_search",        session, "gene_search_select",        type = "gene",
  valid_choices = total_genes)

observeEvent(input$reset_mga_pathway, {
  updateSelectizeInput(session, "mga_pathway_l2",     selected = character(0))
  updateSelectizeInput(session, "mga_pathway_select", selected = character(0))
  updateSelectizeInput(session, "mga_gmt_set_select", selected = character(0), server = TRUE)
})

# ── Pathway colour UI (shown in Network subtab) ───────────────────────────────
output$mga_pathway_filter_ui <- renderUI({
  tagList(
    h6(icon("sitemap"), " Colour nodes by pathway"),
    radioButtons("mga_pathway_source", NULL,
      choices  = c("Reactome" = "reactome", "Custom GMT" = "gmt"),
      selected = input$mga_pathway_source %||% "reactome", inline = TRUE),

    # ── GMT panel ────────────────────────────────────────────────────────────
    conditionalPanel(
      condition = "input.mga_pathway_source == 'gmt'",
      tags$div(
        style = "background:#e8f4fd; border-left:4px solid #3a8fc7; padding:8px 12px; margin-bottom:8px; border-radius:3px; font-size:12px; color:#1a5276;",
        icon("info-circle"),
        tags$strong(" GMT file is shared — upload in the Pathway Annotation tab."),
        tags$br(),
        uiOutput("mga_gmt_status_banner")
      ),
      selectizeInput("mga_gmt_set_select", "Gene sets to highlight:",
        choices = NULL, selected = NULL, multiple = TRUE,
        options = list(placeholder = "Load a GMT in Pathway Annotation tab first…",
                       maxOptions = 500, closeAfterSelect = FALSE)),
      actionButton("reset_mga_pathway", "Reset selection", icon = icon("undo"), class = "btn-sm")
    ),

    # ── Reactome panel ────────────────────────────────────────────────────────
    conditionalPanel(
      condition = "input.mga_pathway_source == 'reactome' || input.mga_pathway_source == null",
      if (ontology_has_hierarchy()) {
        all_pw <- tryCatch(
          dbGetQuery(con, "SELECT PATHWAY_ID, PATHWAY_NAME FROM ONTOLOGY_PATHWAYS ORDER BY PATHWAY_NAME"),
          error = function(e) data.frame(PATHWAY_ID = integer(0), PATHWAY_NAME = character(0))
        )
        tagList(
          selectizeInput("mga_pathway_l2", "Pathway:",
            choices  = setNames(all_pw$PATHWAY_ID, all_pw$PATHWAY_NAME),
            selected = NULL, multiple = TRUE,
            options  = list(placeholder = "Search all pathways…", maxOptions = 200)),
          actionButton("reset_mga_pathway", "Reset", icon = icon("undo"), class = "btn-sm")
        )
      } else {
        tagList(
          selectizeInput("mga_pathway_select", "Pathway:",
            choices = pathway_choices, selected = NULL, multiple = TRUE,
            options = list(placeholder = "All pathways (no filter)", maxOptions = 200)),
          actionButton("reset_mga_pathway", "Reset", icon = icon("undo"), class = "btn-sm")
        )
      }
    )
  )
})

observeEvent(input$reset_gene_search, {
  updateSelectizeInput(session, "gene_search_select", selected = character(0))
  updateCheckboxInput(session, "mga_use_all_genes", value = FALSE)
})


# Single DB query + GMT filter – runs only on Search click
filtered_mga_base <- eventReactive(input$obtain_mga, {
  on.exit(removeNotification("notify_mga"))
  my_data$min_assoc_mga_used <- as.integer(input$min_assoc_mga %||% 1L)
  df <- get_filtered_mga(con = con, input = input, my_data, group_override = "None")
  shiny::validate(shiny::need(nrow(df) > 0, "No associations found with the current filters. Try broadening your selection."))

  if (identical(input$mga_pathway_source, "gmt")) {
    gmt      <- gmt_data()
    sel_sets <- input$mga_gmt_set_select
    if (!is.null(gmt) && length(sel_sets) > 0) {
      gmt_genes <- unique(unlist(gmt[sel_sets], use.names = FALSE))
      df <- df[df$GENE_NAME_clean %in% gmt_genes, ]
      shiny::validate(shiny::need(nrow(df) > 0, "No associations found for the selected gene sets. Try different sets."))
    }
  }

  df
})

# Table view – instant dplyr group-by, re-runs on group_mga change without hitting DB
filtered_mga <- reactive({
  df   <- filtered_mga_base()
  mode <- input$group_mga %||% "None"
  switch(mode,
    "miRNA" = df %>%
      dplyr::group_by(MIRNA_NAME_clean) %>%
      dplyr::summarise(
        GENE_COUNT   = dplyr::n_distinct(GENE_NAME_clean),
        GENE_NAME    = paste(sort(unique(GENE_NAME_clean)), collapse = ","),
        PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(GENE_COUNT)) %>%
      dplyr::mutate(MIRNA_NAME = link_info_fun(MIRNA_NAME_clean, "mirna")) %>%
      dplyr::select(MIRNA_NAME, GENE_COUNT, GENE_NAME, PUBMED_COUNT),
    "gene" = df %>%
      dplyr::group_by(GENE_NAME_clean) %>%
      dplyr::summarise(
        MIRNA_COUNT  = dplyr::n_distinct(MIRNA_NAME_clean),
        MIRNA_NAME   = paste(sort(unique(MIRNA_NAME_clean)), collapse = ","),
        PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(MIRNA_COUNT)) %>%
      dplyr::mutate(GENE_NAME = link_info_fun(GENE_NAME_clean, "gene")) %>%
      dplyr::select(GENE_NAME, MIRNA_COUNT, MIRNA_NAME, PUBMED_COUNT),
    "miRNA and gene" = df %>%
      dplyr::group_by(MIRNA_NAME_clean, GENE_NAME_clean) %>%
      dplyr::summarise(
        PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean),
        PUBMED_ID    = paste(sort(unique(PUBMED_ID_clean)), collapse = ","),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(PUBMED_COUNT)) %>%
      dplyr::mutate(
        MIRNA_NAME = link_info_fun(MIRNA_NAME_clean, "mirna"),
        GENE_NAME  = link_info_fun(GENE_NAME_clean,  "gene"),
        PUBMED_ID  = link_info_fun(PUBMED_ID,         "pubmed")
      ) %>%
      dplyr::select(MIRNA_NAME, GENE_NAME, PUBMED_COUNT, PUBMED_ID),
    df  # "None" – one row per article, HTML links already in base
  )
})

# Flat MIRNA×GENE pairs – derived from GMT-filtered base
filtered_mga_flat <- reactive({
  filtered_mga_base() %>%
    dplyr::group_by(MIRNA_NAME_clean, GENE_NAME_clean) %>%
    dplyr::summarise(
      PUBMED_COUNT = dplyr::n_distinct(PUBMED_ID_clean),
      is_metabolic = dplyr::first(is_metabolic),
      .groups = "drop"
    ) %>%
    dplyr::rename(MIRNA_NAME = MIRNA_NAME_clean, GENE_NAME = GENE_NAME_clean)
})

output$count_mga <- renderUI({
  req(filtered_mga_base())
  df <- filtered_mga_base()

  bar <- render_count_bar(df, list(
    "miRNA"  = "MIRNA_NAME_clean",
    "genes"  = "GENE_NAME_clean",
    "PubMed" = "PUBMED_ID_clean"
  ))

  metab_line <- if ("is_metabolic" %in% names(df)) {
    gene_col <- "GENE_NAME_clean"
    genes_clean  <- df[[gene_col]]
    n_total_g    <- length(unique(genes_clean))
    n_metab_g    <- length(unique(genes_clean[df$is_metabolic == TRUE]))
    pct          <- if (n_total_g > 0) round(100 * n_metab_g / n_total_g) else 0
    tags$p(
      style = "color:#555; font-style:italic; margin-bottom:6px;",
      HTML(paste0(
        "<b>", n_metab_g, "</b> metabolic genes",
        " &nbsp;|&nbsp; ",
        "<b>", n_total_g - n_metab_g, "</b> non-metabolic",
        " &nbsp;&mdash;&nbsp; <b>", pct, "%</b> metabolic"
      ))
    )
  } else NULL

  tagList(bar, metab_line)
})

output$result_table_mga <- renderDT({
  req(filtered_mga())
  shiny::validate(shiny::need(nrow(filtered_mga()) > 0, "No associations found with the current filters. Try broadening your selection."))
  render_standard_dt(filtered_mga())
}, server = TRUE)

output$download_mga_csv <- make_csv_download(
  "mirna_gene", filtered_mga, con = con,
  include_abstract = function() input$download_abstract_mga
)

# Propaga geni e miRNA visibili verso i tab successivi
observeEvent(input$result_table_mga_rows_all, {
  flat   <- filtered_mga_flat()
  genes  <- unique(flat$GENE_NAME)
  mirnas <- unique(flat$MIRNA_NAME)
  updateSelectizeInput(session, "gene_reaction_select", choices = intersect(genes, metabolic_genes), server = TRUE)
  updateSelectizeInput(session, "mirna_reaction_select", choices = mirnas, server = TRUE)
}, ignoreNULL = TRUE)

observeEvent(input$reset_mg, {
  for (id in c("proof_type", "proof_strength", "proof_throughput", "proofs")) {
    updateSelectInput(session, inputId = id, selected = character(0))
  }
  updateRadioButtons(session, "proof_match", selected = "any")
})

observeProofsFilters(input, session, suffix = "")  # miRNA-Gene
observeExpressionFilterCascade(input, session, suffix = "_mga", con)

mga_get_data <- function() filtered_mga_flat()

mga_extra_graph <- eventReactive(input$obtain_mga, {
  pathway_labels <- NULL
  gene_color_map <- NULL

  src <- input$mga_pathway_source %||% "reactome"

  if (src == "gmt") {
    # ── GMT mode: color genes by selected gene sets ───────────────────────────
    gmt        <- gmt_data()
    sel_sets   <- input$mga_gmt_set_select
    if (!is.null(gmt) && length(sel_sets) > 0) {
      active_gmt <- gmt[sel_sets]
      pal    <- GROUP_PALETTE
      set_col <- setNames(pal[((seq_along(sel_sets) - 1L) %% length(pal)) + 1L], sel_sets)

      all_set_genes <- unique(unlist(active_gmt, use.names = FALSE))
      pw_grouped <- do.call(rbind, lapply(all_set_genes, function(g) {
        in_sets <- sel_sets[vapply(active_gmt, function(genes) g %in% genes, logical(1))]
        data.frame(GENE_NAME = g,
                   PATHWAYS  = paste(in_sets, collapse = ", "),
                   n_pw      = length(in_sets),
                   stringsAsFactors = FALSE)
      }))
      pathway_labels <- setNames(as.list(pw_grouped$PATHWAYS), pw_grouped$GENE_NAME)
      gene_color_map <- setNames(
        lapply(seq_len(nrow(pw_grouped)), function(i) {
          if (pw_grouped$n_pw[i] == 1)
            set_col[[ pw_grouped$PATHWAYS[i] ]]
          else
            "#DAA520"
        }),
        pw_grouped$GENE_NAME
      )
    }

  } else {
    # ── Reactome mode: color genes by selected pathway IDs ────────────────────
    anchor_ids <- as.integer(input$mga_pathway_l2)

    if (length(anchor_ids) > 0) {
      pw_long <- tryCatch(dbGetQuery(con, glue_sql(
        "SELECT DISTINCT G.GENE_NAME, P.PATHWAY_NAME
           FROM ONTOLOGY_PATHWAYS_GENES PG
           JOIN ONTOLOGY_PATHWAYS P ON PG.PATHWAY_ID = P.PATHWAY_ID
           JOIN GENES G             ON PG.GENE_ID    = G.GENE_ID
          WHERE PG.PATHWAY_ID IN ({ids*})",
        ids = anchor_ids, .con = con
      )), error = function(e) data.frame(GENE_NAME = character(0), PATHWAY_NAME = character(0)))

      if (nrow(pw_long) > 0) {
        pw_grouped <- pw_long |>
          dplyr::group_by(GENE_NAME) |>
          dplyr::summarise(PATHWAYS = paste(unique(PATHWAY_NAME), collapse = ", "),
                           n_pw     = dplyr::n_distinct(PATHWAY_NAME), .groups = "drop")
        pathway_labels <- setNames(as.list(pw_grouped$PATHWAYS), pw_grouped$GENE_NAME)
        all_pw_names   <- sort(unique(pw_long$PATHWAY_NAME))
        pal    <- GROUP_PALETTE
        pw_col <- setNames(pal[((seq_along(all_pw_names) - 1L) %% length(pal)) + 1L], all_pw_names)
        gene_color_map <- setNames(
          lapply(seq_len(nrow(pw_grouped)), function(i) {
            g <- pw_grouped$GENE_NAME[i]
            if (pw_grouped$n_pw[i] == 1)
              pw_col[[ pw_long$PATHWAY_NAME[pw_long$GENE_NAME == g][1] ]]
            else
              "#DAA520"
          }),
          pw_grouped$GENE_NAME
        )
      }

    } else if (length(input$mga_pathway_select) > 0) {
      selected_pw <- sort(input$mga_pathway_select)
      pw_long <- tryCatch(dbGetQuery(con, glue_sql(
        "SELECT DISTINCT G.GENE_NAME, P.PATHWAY_NAME
           FROM ONTOLOGY_PATHWAYS_GENES PG
           JOIN ONTOLOGY_PATHWAYS P ON PG.PATHWAY_ID = P.PATHWAY_ID
           JOIN GENES G             ON PG.GENE_ID    = G.GENE_ID
          WHERE P.PATHWAY_NAME IN ({vals*})",
        vals = selected_pw, .con = con
      )), error = function(e) data.frame(GENE_NAME = character(0), PATHWAY_NAME = character(0)))

      if (nrow(pw_long) > 0) {
        pw_grouped <- pw_long |>
          dplyr::group_by(GENE_NAME) |>
          dplyr::summarise(PATHWAYS = paste(unique(PATHWAY_NAME), collapse = ", "),
                           n_pw     = dplyr::n_distinct(PATHWAY_NAME), .groups = "drop")
        pathway_labels <- setNames(as.list(pw_grouped$PATHWAYS), pw_grouped$GENE_NAME)
        pal    <- GROUP_PALETTE
        pw_col <- setNames(pal[((seq_along(selected_pw) - 1L) %% length(pal)) + 1L], selected_pw)
        gene_color_map <- setNames(
          lapply(seq_len(nrow(pw_grouped)), function(i) {
            g <- pw_grouped$GENE_NAME[i]
            if (pw_grouped$n_pw[i] == 1)
              pw_col[[ pw_long$PATHWAY_NAME[pw_long$GENE_NAME == g][1] ]]
            else
              "#DAA520"
          }),
          pw_grouped$GENE_NAME
        )
      }
    }
  }

  list(
    mirna_color         = MGA_MIRNA_COLOR,
    target_color        = "#6699FF",
    gene_pathway_labels = pathway_labels,
    gene_color_map      = gene_color_map
  )
})

# Network tab module
mga_plot_module <- plotVisualizationServer("plot_mg_net",
  get_data_fn      = mga_get_data,
  type_network     = "GENE_NAME",
  x1 = "MIRNA_NAME", x2 = "GENE_NAME",
  extra_graph_args = mga_extra_graph,
  mode             = "network",
  external_trigger = reactive(input$obtain_mga),
  physics_input    = reactive(input$mga_network_layout %||% "forceAtlas2Based")
)
observeEvent(mga_plot_module$net_png(), {
  png <- mga_plot_module$net_png()
  if (!is.null(png) && nzchar(png)) my_data$mga_net_png <- png
})




# =============================================================================
# Co-targeting analysis
# =============================================================================

cotarget_gene_df <- reactive({
  shiny::validate(shiny::need(
    tryCatch({ filtered_mga_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA–Gene search first (Table summary tab)."
  ))

  df <- filtered_mga_base()
  shiny::validate(shiny::need(nrow(df) > 0, "No associations in current miRNA–Gene search."))



  # Use clean columns (strip HTML)
  mirna_col <- if ("MIRNA_NAME_clean" %in% names(df)) "MIRNA_NAME_clean" else "MIRNA_NAME"
  gene_col  <- if ("GENE_NAME_clean"  %in% names(df)) "GENE_NAME_clean"  else "GENE_NAME"

  df$mirna_ <- strip_html(df[[mirna_col]])
  df$gene_  <- strip_html(df[[gene_col]])

  # Filter by minimum associations per miRNA-gene pair
  min_assoc <- max(1L, as.integer(input$min_assoc_mga %||% 1L))
  if (min_assoc > 1L) {
    pair_counts <- df %>%
      dplyr::group_by(mirna_, gene_) %>%
      dplyr::summarise(
        n = if ("PUBMED_ID" %in% names(df))
              length(split_unique(PUBMED_ID))
            else
              dplyr::n(),
        .groups = "drop"
      )
    keep_pairs <- pair_counts[pair_counts$n >= min_assoc, c("mirna_", "gene_")]
    df <- dplyr::semi_join(df, keep_pairs, by = c("mirna_", "gene_"))
    shiny::validate(shiny::need(nrow(df) > 0,
                  "No pairs survive the minimum associations filter. Lower the threshold."))
  }

  result <- unique(data.frame(
    MIRNA_NAME = df$mirna_,
    GENE_NAME  = df$gene_,
    stringsAsFactors = FALSE
  ))

  shiny::validate(shiny::need(nrow(result) > 0, "No miRNA–gene pairs found."))
  result
})

cotarget_mat <- reactive({
  df        <- cotarget_gene_df()
  node_type <- input$mga_cotarget_node %||% "mirna_mirna"
  MAX_DIM   <- 100L
  if (node_type == "gene_gene") {
    n <- length(unique(df$GENE_NAME))
    shiny::validate(shiny::need(n >= 2, "Need at least 2 genes. Broaden miRNA–Gene filters."))
    shiny::validate(shiny::need(n <= MAX_DIM, paste0("Too many nodes (", n, "×", n, " matrix). Limit: ", MAX_DIM, "×", MAX_DIM, ". Raise the min-degree or reduce the gene selection.")))
    compute_jaccard_matrix(df, col_nodes = "GENE_NAME", col_shared = "MIRNA_NAME")
  } else {
    n <- length(unique(df$MIRNA_NAME))
    shiny::validate(shiny::need(n >= 2, "Need at least 2 miRNAs. Broaden miRNA–Gene filters."))
    shiny::validate(shiny::need(n <= MAX_DIM, paste0("Too many nodes (", n, "×", n, " matrix). Limit: ", MAX_DIM, "×", MAX_DIM, ". Raise the min-degree or reduce the miRNA selection.")))
    compute_jaccard_matrix(df, col_nodes = "MIRNA_NAME", col_shared = "GENE_NAME")
  }
})

# Named vector: MIRNA_NAME (arm-specific) → disease label or "Shared"
# Returns NULL when: clustering ON, checkbox OFF, not MDA source, or only 1 disease
cotarget_mirna_annotation <- reactive({
  if (!isTRUE(input$cotarget_show_disease_ann)) return(NULL)
  if (isTRUE(input$cotarget_cluster))           return(NULL)

  df        <- cotarget_gene_df()
  node_type <- input$mga_cotarget_node %||% "mirna_mirna"
  if (node_type != "mirna_mirna") return(NULL)

  is_mda_source <- !identical(input$mga_mirna_source, "scratch") &&
    (isTRUE(input$use_all_mirnas_from_mda) || isTRUE(input$filter_mirna_common_checkbox))
  if (!is_mda_source) return(NULL)

  md_map <- my_data$mirna_disease_df
  if (is.null(md_map) || nrow(md_map) == 0) return(NULL)

  n_diseases <- length(unique(md_map$DISEASE))
  if (n_diseases < 2) return(NULL)

  mirnas <- unique(df$MIRNA_NAME)
  prem_map <- tryCatch(
    dbGetQuery(con, glue_sql(
      "SELECT MIRNA_NAME, MIRNA_PREMATURE FROM MIRNAS WHERE MIRNA_NAME IN ({mirnas*})",
      mirnas = mirnas, .con = con
    )),
    error = function(e) NULL
  )
  if (is.null(prem_map) || nrow(prem_map) == 0) return(NULL)

  disease_order <- c(input$disease_select %||% character(0), "Shared", "Unknown")

  ann <- setNames(vapply(mirnas, function(m) {
    prem <- prem_map$MIRNA_PREMATURE[prem_map$MIRNA_NAME == m]
    dis  <- unique(md_map$DISEASE[md_map$MIRNA_NAME %in% prem])
    if (length(dis) == 0)                        return("Unknown")
    if (length(dis) >= n_diseases)               return("Shared")
    dis[1L]
  }, character(1L)), mirnas)
  # attach disease_order as attribute so render function can sort
  attr(ann, "disease_order") <- disease_order
  ann
})

output$cotarget_heatmap <- renderPlotly({
  mat       <- cotarget_mat()
  req(!is.null(mat))
  node_type <- input$mga_cotarget_node %||% "mirna_mirna"
  node_lbl  <- if (node_type == "gene_gene") "Gene" else "miRNA"
  render_cotarget_heatmap(mat,
    cluster       = isTRUE(input$cotarget_cluster),
    node_label    = node_lbl,
    color_palette = input$cotarget_palette   %||% "YlOrRd",
    mask_zero     = isTRUE(input$cotarget_mask_zero),
    mask_diag     = isTRUE(input$cotarget_mask_diag),
    annotation    = cotarget_mirna_annotation()) %>%
  htmlwidgets::onRender("function(el) {
    Plotly.toImage(el, {format:'png', width:2580, height:2100}).then(function(url) {
      Shiny.setInputValue('mga_hm_raw_png', url, {priority:'event'});
    });
  }")
}) |> bindCache(
  cotarget_mat(),
  input$mga_cotarget_node,
  input$cotarget_cluster,
  input$cotarget_palette,
  input$cotarget_mask_zero,
  input$cotarget_mask_diag,
  cotarget_mirna_annotation()
)

store_png_variant("mga_hm_raw_png", my_data, input,
  node_input  = "mga_cotarget_node",
  variant_map = list(mirna_mirna = "mga_hm_mirna_png",
                     gene_gene   = "mga_hm_gene_png"))

output$cotarget_count <- renderUI({
  mat       <- cotarget_mat()
  req(!is.null(mat))
  node_type <- input$mga_cotarget_node %||% "mirna_mirna"
  node_lbl  <- if (node_type == "gene_gene") "genes" else "miRNAs"
  tbl       <- cotarget_pairs_table(mat, threshold = 0)
  n_pairs   <- if (is.null(tbl)) 0L else nrow(tbl)
  tags$p(
    style = "color:#555; font-style:italic; margin-bottom:6px;",
    HTML(paste0(
      "<b>", nrow(mat), "</b> ", node_lbl, " &nbsp;|&nbsp; ",
      "<b>", n_pairs, "</b> pairs"
    ))
  )
})

cotarget_mga_pairs_df <- reactive({
  mat       <- cotarget_mat()
  req(!is.null(mat))
  node_type <- input$mga_cotarget_node %||% "mirna_mirna"
  if (node_type == "gene_gene") {
    cotarget_pairs_table(mat, threshold = 0,
                          df_source  = cotarget_gene_df(),
                          col_nodes  = "GENE_NAME",
                          col_shared = "MIRNA_NAME")
  } else {
    cotarget_pairs_table(mat, threshold = 0,
                          df_source  = cotarget_gene_df(),
                          col_nodes  = "MIRNA_NAME",
                          col_shared = "GENE_NAME")
  }
})

output$cotarget_table <- DT::renderDataTable({
  tbl <- cotarget_mga_pairs_df()
  shiny::validate(shiny::need(!is.null(tbl), "No co-targeting pairs found."))
  render_standard_dt(tbl)
})



# =============================================================================
# Robustness analysis (miRNA–Gene)
# =============================================================================

# Helper: build mirna + gene + proof clauses from current input state
.mga_rob_clauses <- function(input, my_data, con) {
  mirnas    <- resolve_mirnas_mga(input, my_data)
  mirnas_id <- get_mirna_ids(mirnas, con)
  use_all_mirnas <- isTRUE(input$mga_use_all_scratch) ||
    isTRUE(input$use_all_mirnas_from_mda) ||
    isTRUE(input$filter_mirna_common_checkbox)
  mirna_clause <- if (length(mirnas) > 0) {
    if (length(mirnas_id) == 0) DBI::SQL("AND 1=0")
    else glue_sql("AND A.MIRNA_ID IN ({vals*})", vals = mirnas_id, .con = con)
  } else if (use_all_mirnas) {
    DBI::SQL("")
  } else {
    DBI::SQL("AND 1=0")
  }

  use_all_genes  <- isTRUE(input$mga_use_all_genes)
  only_metabolic <- isTRUE(input$mga_only_metabolic_genes)
  genes          <- input$gene_search_select
  gene_ids       <- if (length(genes) > 0 && !use_all_genes && !only_metabolic)
    dbGetQuery(con, glue_sql("SELECT GENE_ID FROM GENES WHERE GENE_NAME IN ({vals*})",
                             vals = genes, .con = con))$GENE_ID
  else integer(0)
  gene_clause <- if (use_all_genes || (length(genes) == 0 && !only_metabolic)) DBI::SQL("")
                 else if (only_metabolic) DBI::SQL("AND A.GENE_ID IN (SELECT GENE_ID FROM REACTIONS_GENES)")
                 else if (length(gene_ids) > 0)
                   glue_sql("AND A.GENE_ID IN ({vals*})", vals = gene_ids, .con = con)
                 else DBI::SQL("AND 1=0")

  selected_proofs <- input$proofs
  proof_clause <- if (length(selected_proofs) > 0)
    glue_sql("AND D.PROOF_NAME IN ({vals*})", vals = selected_proofs, .con = con)
  else DBI::SQL("")

  human_clause <- if (isTRUE(input$filter_human_mga) && isTRUE(HAS_IS_HUMAN))
    DBI::SQL("AND CAST(A.PUBMED_ID AS VARCHAR) NOT IN (SELECT CAST(PUBMED_ID AS VARCHAR) FROM ARTICLES WHERE IS_HUMAN = 0)")
  else
    DBI::SQL("")

  list(mirna = mirna_clause, gene = gene_clause, proof = proof_clause, human = human_clause)
}

# ── Plot 1: threshold sweep ───────────────────────────────────────────────────
rob_mga_thresh_data <- eventReactive(input$run_rob_mga, {
  shiny::validate(shiny::need(
    tryCatch({ filtered_mga_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA–Gene search first (Table tab)."
  ))
  withProgress(message = "Running threshold sweep...", value = 0, {
    setProgress(0.2, detail = "Building clauses")
    cl <- .mga_rob_clauses(input, my_data, con)
    setProgress(0.5, detail = "Sweeping thresholds")
    df <- compute_mga_robustness(
      con          = con,
      mirna_clause = cl$mirna,
      gene_clause  = cl$gene,
      proof_clause = cl$proof,
      human_clause = cl$human,
      min_thresh   = input$rob_mga_min_thresh %||% 1L,
      max_thresh   = input$rob_mga_max_thresh %||% 20L
    )
    shiny::validate(shiny::need(nrow(df) > 0, "No associations found for this selection."))
    setProgress(1)
    df
  })
})

output$rob_mga_plot <- renderPlotly({
  df <- rob_mga_thresh_data()
  req(!is.null(df) && nrow(df) > 0)

  current_thresh <- as.integer(input$min_assoc_mga %||% 1L)

  p <- plot_ly(df, x = ~threshold) %>%
    add_lines(y = ~norm_base(n_edges), name = "Edge retention",  line = list(color = "#2c7bb6", width = 2)) %>%
    add_lines(y = ~norm_base(n_mirna), name = "miRNA retention", line = list(color = "#d7191c", width = 2)) %>%
    add_lines(y = ~norm_base(n_gene),  name = "Gene retention",  line = list(color = "#1a9641", width = 2))

  if (current_thresh >= min(df$threshold) && current_thresh <= max(df$threshold)) {
    p <- p %>% add_segments(
      x = current_thresh, xend = current_thresh, y = 0, yend = 1,
      line = list(color = "black", dash = "dash", width = 1),
      name = "Current threshold", showlegend = TRUE, inherit = FALSE
    )
  }
  p %>% layout(
    title     = list(text = "Network robustness vs evidence threshold", font = list(size = 14)),
    xaxis     = list(title = "Min publications per pair", dtick = 1),
    yaxis     = list(title = "Fraction of baseline (normalised)", range = c(0, 1.05)),
    legend    = list(orientation = "h", y = -0.2),
    hovermode = "x unified",
    paper_bgcolor = "#ffffff", plot_bgcolor = "#fafafa"
  ) %>% plotly_clean_config() %>%
  htmlwidgets::onRender("function(el) {
    Plotly.toImage(el, {format:'png', width:2580, height:1440}).then(function(url) {
      Shiny.setInputValue('mga_rob_png', url, {priority:'event'});
    });
  }")
}) |> bindEvent(input$run_rob_mga, ignoreNULL = TRUE)

store_png("mga_rob_png", my_data, input)

output$rob_mga_table <- DT::renderDataTable(server = TRUE, {
  df <- rob_mga_thresh_data()
  req(!is.null(df) && nrow(df) > 0)
  col_map <- c(
    threshold = "Threshold",
    n_edges   = "N edges",
    n_mirna   = "N miRNA",
    n_gene    = "N genes"
  )
  df <- df[, intersect(names(col_map), names(df)), drop = FALSE]
  names(df) <- col_map[names(df)]
  DT::datatable(df, rownames = FALSE,
                options = list(pageLength = 30, scrollX = TRUE, dom = "tip"))
}) |> bindEvent(input$run_rob_mga, ignoreNULL = TRUE)

# ── Co-targeting robustness ───────────────────────────────────────────────────
rob_mga_cotarget_data <- eventReactive(input$run_rob_mga_cotarget, {
  shiny::validate(shiny::need(
    tryCatch({ filtered_mga_base(); TRUE }, error = function(e) FALSE),
    "Run a miRNA–Gene search first (Table tab)."
  ))
  withProgress(message = "Running co-targeting robustness...", value = 0.2, {
    node_type <- input$rob_mga_cotarget_node %||% "mirna_mirna"
    if (node_type == "mirna_mirna") {
      col_nodes <- "MIRNA_NAME"; col_shared <- "GENE_NAME"
    } else {
      col_nodes <- "GENE_NAME"; col_shared <- "MIRNA_NAME"
    }
    df <- filtered_mga_base()
    req(!is.null(df) && nrow(df) > 0)
    mirna_col <- if ("MIRNA_NAME_clean" %in% names(df)) "MIRNA_NAME_clean" else "MIRNA_NAME"
    gene_col  <- if ("GENE_NAME_clean"  %in% names(df)) "GENE_NAME_clean"  else "GENE_NAME"
    # Build pairs with N_PUB per miRNA-gene pair
    pairs <- df %>%
      dplyr::mutate(MIRNA_NAME = strip_html(.data[[mirna_col]]),
                    GENE_NAME  = strip_html(.data[[gene_col]])) %>%
      dplyr::group_by(MIRNA_NAME, GENE_NAME) %>%
      dplyr::summarise(
        N_PUB = if ("PUBMED_ID" %in% names(df)) length(split_unique(PUBMED_ID)) else dplyr::n(),
        .groups = "drop"
      )
    incProgress(0.5)
    compute_cotarget_robustness(
      pairs_df       = pairs,
      col_nodes      = col_nodes,
      col_shared     = col_shared,
      min_thresh     = input$rob_mga_min_thresh_ct %||% 1L,
      max_thresh     = input$rob_mga_max_thresh_ct %||% 20L,
      jaccard_cutoff = input$rob_mga_jaccard_cutoff %||% 0.1
    )
  })
})

output$rob_mga_cotarget_plot <- renderPlotly({
  res <- rob_mga_cotarget_data()
  df  <- res$sweep
  req(!is.null(df) && nrow(df) > 0)
  plot_ly(df, x = ~threshold) %>%
    add_lines(y = ~n_pairs,      name = "Co-targeting pairs", line = list(color = "#2c7bb6", width = 2)) %>%
    add_lines(y = ~mean_jaccard, name = "Mean Jaccard",       line = list(color = "#fdae61", width = 2, dash = "dash"),
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
        Shiny.setInputValue('mga_rob_ct_raw_png', url, {priority:'event'});
      });
    }")
}) |> bindEvent(input$run_rob_mga_cotarget, ignoreNULL = TRUE)

store_png_variant("mga_rob_ct_raw_png", my_data, input,
  node_input  = "rob_mga_cotarget_node",
  variant_map = list(mirna_mirna = "mga_rob_ct_mirna_png",
                     gene_gene   = "mga_rob_ct_gene_png"))

output$rob_mga_stability_table <- DT::renderDataTable(server = TRUE, {
  res <- rob_mga_cotarget_data()
  df  <- res$stability
  req(!is.null(df) && nrow(df) > 0)
  node_type <- isolate(input$rob_mga_cotarget_node) %||% "mirna_mirna"
  lbl_a <- if (node_type == "mirna_mirna") "miRNA A" else "Gene A"
  lbl_b <- if (node_type == "mirna_mirna") "miRNA B" else "Gene B"
  if (ncol(df) == 4) names(df) <- c(lbl_a, lbl_b, "N thresholds", "Stability score")
  DT::datatable(df, rownames = FALSE,
                options = list(pageLength = 20, scrollX = TRUE, dom = "frtip"))
}) |> bindEvent(input$run_rob_mga_cotarget, ignoreNULL = TRUE)


