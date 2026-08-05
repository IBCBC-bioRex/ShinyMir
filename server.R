server <- function(input, output, session) {

  observe({

  # =========================================================
  # Reactive values condivisi tra tutti i tab
  # =========================================================
  my_data <- reactiveValues()
  my_data$disease <- dbGetQuery(con, "SELECT DISEASE_ID, DISEASE FROM DISEASES")
  my_data$genes   <- dbGetQuery(con, "SELECT GENE_ID, GENE_NAME FROM GENES")

  updateSelectizeInput(session, "mirna_select",                choices = mirna_choices,   server = TRUE)
  updateSelectizeInput(session, "disease_select",             choices = disease_choices, server = TRUE)
  updateSelectizeInput(session, "mirna_gene_select_scratch",  choices = mirna_choices,   server = TRUE)
  updateSelectizeInput(session, "mirna_reaction_select_scratch", choices = mirna_choices, server = TRUE)
  updateSelectizeInput(session, "gene_reaction_select_manual", choices = metabolic_genes, server = TRUE)

  # Documentation outputs removed: all static HTML, rendered directly in ui_helpers.R

  observeEvent(input$reset_all, {
    shinyjs::runjs("window.location.reload(true)")
  })

  # ONTOLOGY_HAS_HIERARCHY is a plain boolean pre-computed in global.R at boot
  ontology_has_hierarchy <- reactive(ONTOLOGY_HAS_HIERARCHY)

  # Shared GMT reactive — single upload point (Pathway Annotation tab), consumed by MGA + OVR
  gmt_data <- reactiveVal(NULL)

  # =========================================================
  # Mpdules server
  # =========================================================
  source("server/server_upload.R",  local = environment())
  source("server/server_toggles.R", local = environment())
  source("server/server_mda.R",     local = environment())
  source("server/server_mga.R",     local = environment())
  source("server/server_mr.R",      local = environment())
  source("server/server_overrep.R", local = environment())

  # =========================================================
  # Suspend heavy outputs when their panel is hidden
  # (default is TRUE but visNetwork/plotly inside conditionalPanel
  #  can escape Shiny's auto-detection – be explicit)
  # =========================================================
  heavy_outputs <- c(
    # Co-targeting heatmaps (networks removed)
    "cotarget_heatmap",
    "cotarget_mda_heatmap",
    "cotarget_mr_heatmap",
    # Pathway annotation
    "mirna_pathway_heatmap", "ovr_barplot_combined",
    "mr_essentiality_plot", "mr_barplot",
    # Metabolic heatmap
    "metab_heatmap",
    # Network validation
    # Main result tables
    "result_table_mda", "result_table_mga", "result_table_mr"
  )
  for (oid in heavy_outputs) {
    outputOptions(output, oid, suspendWhenHidden = TRUE)
  }

  # =========================================================
  # Generate report — two-phase: JS captures canvases, then R builds HTML
  # =========================================================
  observeEvent(input$download_report_btn, {
    mda_df      <- tryCatch(filtered_mda(),      error = function(e) NULL)
    mda_flat_df <- tryCatch(filtered_mda_flat(), error = function(e) NULL)
    mga_df      <- tryCatch(filtered_mga(),      error = function(e) NULL)
    mga_flat_df <- tryCatch(filtered_mga_flat(), error = function(e) NULL)
    mr_df       <- tryCatch(filtered_mr(),       error = function(e) NULL)
    mr_flat_df  <- tryCatch(filtered_mr_base(),  error = function(e) NULL)
    ovr_df <- tryCatch({
      df <- my_data$ovr_mga_df
      if (is.null(df) || nrow(df) == 0) NULL else df
    }, error = function(e) NULL)
    # Create temp working directory
    tmp_dir <- tempfile(pattern = "shinymir_report_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    # Helper: decode data URL → write PNG → return filename or NULL
    .write_png <- function(data_url, filename) {
      if (is.null(data_url) || !nzchar(data_url %||% "")) return(NULL)
      tryCatch({
        b64 <- sub("^data:image/[^;]+;base64,", "", data_url)
        writeBin(jsonlite::base64_dec(b64), file.path(tmp_dir, filename))
        filename
      }, error = function(e) NULL)
    }

    img_files <- list(
      # MDA
      mda_net              = .write_png(my_data$mda_net_png,              "mda_network.png"),
      mda_rob              = .write_png(my_data$mda_rob_png,              "mda_robustness.png"),
      mda_rob_ct_mirna     = .write_png(my_data$mda_rob_ct_mirna_png,    "mda_rob_ct_mirna.png"),
      mda_rob_ct_disease   = .write_png(my_data$mda_rob_ct_disease_png,  "mda_rob_ct_disease.png"),
      cotarget_mda_mirna   = .write_png(my_data$cotarget_mda_mirna_png,  "cotarget_mda_mirna.png"),
      cotarget_mda_disease = .write_png(my_data$cotarget_mda_disease_png,"cotarget_mda_disease.png"),
      # MGA
      mga_net              = .write_png(my_data$mga_net_png,             "mga_network.png"),
      mga_hm_mirna         = .write_png(my_data$mga_hm_mirna_png,       "mga_cotarget_mirna.png"),
      mga_hm_gene          = .write_png(my_data$mga_hm_gene_png,        "mga_cotarget_gene.png"),
      mga_rob              = .write_png(my_data$mga_rob_png,             "mga_robustness.png"),
      mga_rob_ct_mirna     = .write_png(my_data$mga_rob_ct_mirna_png,   "mga_rob_ct_mirna.png"),
      mga_rob_ct_gene      = .write_png(my_data$mga_rob_ct_gene_png,    "mga_rob_ct_gene.png"),
      # MR
      mr_net               = .write_png(my_data$mr_net_png,             "mr_network.png"),
      mr_barplot           = .write_png(my_data$mr_barplot_png,         "mr_barplot.png"),
      mr_hm_mirna          = .write_png(my_data$mr_hm_mirna_png,       "mr_cotarget_mirna.png"),
      mr_hm_subsystem      = .write_png(my_data$mr_hm_subsystem_png,   "mr_cotarget_subsystem.png"),
      ess_plot             = .write_png(my_data$mr_ess_plot_png,        "mr_essentiality.png"),
      # OVR
      ovr_plot             = .write_png(my_data$ovr_plot_png,           "ovr_barplot.png"),
      ovr_hm               = .write_png(my_data$ovr_hm_png,            "ovr_heatmap.png")
    )

    db_info <- list(
      filename     =  "shinyMIR.duckdb",
      n_diseases   = N_TOTAL_DISEASES,
      n_genes      = N_TOTAL_GENES,
      n_reactions  = tryCatch(dbGetQuery(con, "SELECT COUNT(DISTINCT REACTION_ID) FROM REACTIONS")[[1]], error = function(e) NA_integer_),
      n_mirnas     = length(mirna_choices)
    )
    html <- build_shinymir_report(
      mda_df = mda_df, mda_flat_df = mda_flat_df,
      mga_df = mga_df, mga_flat_df = mga_flat_df,
      mr_df  = mr_df,  mr_flat_df  = mr_flat_df,
      ovr_df = ovr_df,
      input = input, img_files = img_files, db_info = db_info
    )
    writeLines(html, file.path(tmp_dir, "report.html"), useBytes = FALSE)

    # ZIP all files in tmp_dir
    zip_path <- tempfile(fileext = ".zip")
    all_files <- list.files(tmp_dir, full.names = FALSE)
    owd <- setwd(tmp_dir)
    tryCatch(utils::zip(zip_path, files = all_files),
             finally = setwd(owd))

    raw_zip <- readBin(zip_path, "raw", file.info(zip_path)$size)
    unlink(zip_path)
    b64 <- jsonlite::base64_enc(raw_zip)
    fname <- paste0("shinymir_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    session$sendCustomMessage("trigger_zip_download", list(b64 = b64, filename = fname))
  })

  })
} # end server
