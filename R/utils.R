################################################################################
#####Funzioni per info click (disease, mirna,gene,reaction)
################################################################################


# Suppress showlegend on all plotly traces (used after heatmaply side-color injection).
suppress_sidebar_legends <- function(p) {
  p$x$data <- lapply(p$x$data, function(tr) {
    if (!is.null(tr$showlegend)) tr$showlegend <- FALSE
    tr
  })
  p
}

# Strip -3p / -5p suffix and lowercase — normalises mature → premature miRNA name
normalize_mirna_name <- function(x) tolower(gsub("-[35]p$", "", x, ignore.case = TRUE))

# Costruisce un box legenda dinamico per le network.
# items: named character vector  nome → colore hex, es. c("miRNA"="#FF6666", "Disease"="#6699FF")
# dot_shapes: nomi che usano cerchio (dot), gli altri usano quadrato (box)
build_legend_ui <- function(items, dot_shapes = "miRNA") {
  tags$div(
    style = "padding:10px; background:white; border:1px solid #ddd; border-radius:5px;",
    tags$strong("Legend"),
    tags$div(
      style = "margin-top:8px;",
      lapply(seq_along(items), function(i) {
        lbl   <- names(items)[i]
        col   <- items[i]
        is_dot <- lbl %in% dot_shapes
        tags$div(
          style = "display:flex; align-items:center; gap:6px; margin-bottom:6px;",
          tags$div(style = sprintf(
            "width:14px; height:14px; background:%s; flex-shrink:0;%s",
            col,
            if (is_dot) "border-radius:50%;" else ""
          )),
          lbl
        )
      })
    )
  )
}

customStyles <- function() {
  HTML("
  .box {
    border-radius: 10px;
    box-shadow: 0 6px 16px rgba(0,0,0,0.08);
    contain: layout style;        /* isolate reflow to this box */
  }

  .box-header {
    font-weight: 600;
  }

  h4 {
    font-weight: 600;
    margin-top: 0;
  }

  .btn {
    border-radius: 6px;
  }

  .main-header .logo {
    font-weight: 700;
    letter-spacing: 0.5px;
  }

  .vis-legend {
    background: rgba(255,255,255,0.92) !important;
    border: 1px solid #cccccc !important;
    border-radius: 5px !important;
    padding: 6px 8px !important;
  }

  /* GPU-composite layers for heavy outputs */
  .vis-network,
  .shiny-plot-output,
  .plotly {
    will-change: transform;
    transform: translateZ(0);
  }

  /* Smooth fade-in for spinner replacement */
  .shiny-bound-output {
    transition: opacity 0.15s ease;
  }
  ")
}
################################################################################
# Helper condiviso per tutti i modal informativi.
# tabella: matrix a 2 colonne (label | valore). escape=FALSE per celle con HTML.
show_info_modal <- function(title, tabella, escape = FALSE) {
  table_html <- kableExtra::kable(
    tabella,
    format      = "html",
    escape      = escape,
    table.attr  = "class='table table-bordered'"
  )
  showModal(modalDialog(
    title     = title,
    HTML(table_html),
    easyClose = TRUE,
    footer    = modalButton("Close")
  ))
}

show_disease_modal <- function(disease_name, con) {
  description <- dbGetQuery(con, glue_sql(
    "SELECT DISEASE, DOID, DESCRIPTION FROM DISEASES WHERE DISEASE = {val*}",
    val = disease_name, .con = con
  ))
  tabella <- matrix(
    c("DOID", description$DOID,
      "Description", description$DESCRIPTION),
    ncol = 2, byrow = TRUE
  )
  show_info_modal(paste("Detail about:", disease_name), tabella)
}


show_pubmed_modal <- function(pubmed_id, con) {
  description <- dbGetQuery(con, glue_sql(
    "SELECT * FROM ARTICLES WHERE PUBMED_ID = {val*}",
    val = pubmed_id, .con = con
  ))
  if (nrow(description) == 0) {
    showNotification("No article found for this PubMed ID.", type = "error")
    return()
  }
  pubmed_link <- sprintf(
    '<a href="https://pubmed.ncbi.nlm.nih.gov/%s" target="_blank">%s</a>',
    htmltools::htmlEscape(pubmed_id), htmltools::htmlEscape(pubmed_id)
  )
  tabella <- matrix(
    c("TITLE",       description$TITLE,
      "YEAR",        description$YEAR,
      "ABSTRACT",    description$ABSTRACT,
      "PUBMED Link", pubmed_link),
    ncol = 2, byrow = TRUE
  )
  show_info_modal(paste("Details for PubMed ID:", pubmed_id), tabella, escape = FALSE)
}


show_gene_modal <- function(gene_name, con) {
  description <- dbGetQuery(con, glue_sql(
    "SELECT * FROM GENES WHERE GENE_NAME = {val*}",
    val = gene_name, .con = con
  ))
  if (nrow(description) == 0) {
    showNotification("No information found for this gene.", type = "error")
    return()
  }
  tabella <- matrix(
    c("HGNC_ID",        description$hgnc_id,
      "ENSEMBL",        description$ensembl_gene_id,
      "CHROMOSOME",     description$chromosome_name,
      "START_POSITION", description$start_position,
      "END POSITION",   description$end_position,
      "STRAND",         description$strand,
      "TYPE",           description$gene_biotype,
      "DESCRIPTION",    description$description),
    ncol = 2, byrow = TRUE
  )
  show_info_modal(paste("Details for GENE:", gene_name), tabella, escape = FALSE)
}


show_mirna_modal <- function(mirna_name, con) {
  if (endsWith(mirna_name, "-3p") || endsWith(mirna_name, "-5p")) {
    description <- dbGetQuery(con, glue_sql(
      "SELECT * FROM MIRNAS WHERE MIRNA_NAME = {val*}",
      val = mirna_name, .con = con
    ))
    if (nrow(description) == 0) {
      showNotification("No information found for this miRNA", type = "error")
      return()
    }
    tabella <- matrix(
      c("SEQUENCE",     description$SEQUENCE,
        "ACCESSION_ID", description$MIRBASE_MATURE_ACC),
      ncol = 2, byrow = TRUE
    )
  } else {
    # Single query for both arms
    desc_both <- dbGetQuery(con, glue_sql(
      "SELECT * FROM MIRNAS WHERE MIRNA_NAME IN ({vals*})",
      vals = paste0(mirna_name, c("-3p", "-5p")), .con = con
    ))
    if (nrow(desc_both) == 0) {
      showNotification("No information found for this miRNA", type = "error")
      return()
    }
    description_3p <- desc_both[desc_both$MIRNA_NAME == paste0(mirna_name, "-3p"), ]
    description_5p <- desc_both[desc_both$MIRNA_NAME == paste0(mirna_name, "-5p"), ]
    tabella <- matrix(
      c("3P SEQUENCE",     if (nrow(description_3p) > 0) description_3p$SEQUENCE          else NA_character_,
        "3P ACCESSION_ID", if (nrow(description_3p) > 0) description_3p$MIRBASE_MATURE_ACC else NA_character_,
        "5P SEQUENCE",     if (nrow(description_5p) > 0) description_5p$SEQUENCE          else NA_character_,
        "5P ACCESSION_ID", if (nrow(description_5p) > 0) description_5p$MIRBASE_MATURE_ACC else NA_character_),
      ncol = 2, byrow = TRUE
    )
  }
  show_info_modal(paste("Details for miRNA:", mirna_name), tabella, escape = FALSE)
}

# ==============================================================================
# Helper condivisi – evitano duplicazione in server.R
# ==============================================================================

drop_clean_cols <- function(df) {
  clean_cols <- grep("_clean$", names(df), value = TRUE)
  df[, setdiff(names(df), clean_cols), drop = FALSE]
}

# Render DT standard usato in tutti e tre i tab principali + over-representation
render_standard_dt <- function(df) {
  df <- drop_clean_cols(df)
  # Flatten list-columns and non-standard types that crash DT filter bars
  df[] <- lapply(df, function(col) {
    if (is.list(col))
      return(vapply(col, function(x) paste(unlist(x), collapse = ", "), character(1)))
    if (inherits(col, "integer64"))
      return(as.integer(col))
    if (is.numeric(col) && !is.integer(col) && length(col) > 0 &&
        all(is.na(col) | col == floor(col)))
      return(as.integer(col))
    col
  })
  DT::datatable(
    df,
    filter     = "top",
    escape     = FALSE,
    extensions = "Buttons",
    selection  = "multiple",
    options    = list(
      pageLength  = 10,
      scrollX     = TRUE,
      deferRender = TRUE,
      dom         = "Bfrtip",
      buttons     = list(
        list(extend = "colvis", text = "Columns")
      )
    )
  )
}

# Prepare df for CSV/Excel download: strip HTML tags, drop _clean helper columns.
clean_df_for_export <- function(df) {
  df <- drop_clean_cols(df)
  df[] <- lapply(df, function(col) {
    if (is.character(col)) strip_html(col) else col
  })
  df
}

# Aggiunge una colonna ABSTRACT a df usando PUBMED_ID_clean (valore raw,
# eventualmente GROUP_CONCAT "id1,id2,...") per recuperare gli abstract da ARTICLES.
# Più ID nella stessa riga -> abstract uniti con " | ".
attach_abstracts <- function(df, con) {
  if (!"PUBMED_ID_clean" %in% names(df)) return(df)

  ids_list  <- strsplit(as.character(df$PUBMED_ID_clean), ",\\s*")
  all_ids   <- unique(trimws(unlist(ids_list)))
  all_ids   <- all_ids[nzchar(all_ids)]
  if (length(all_ids) == 0) {
    df$ABSTRACT <- NA_character_
    return(df)
  }

  abstracts <- dbGetQuery(con, glue_sql(
    "SELECT PUBMED_ID, ABSTRACT FROM ARTICLES WHERE PUBMED_ID IN ({vals*})",
    vals = all_ids, .con = con
  ))
  lookup <- setNames(abstracts$ABSTRACT, as.character(abstracts$PUBMED_ID))

  df$ABSTRACT <- vapply(ids_list, function(ids) {
    ids <- trimws(ids)
    txt <- lookup[ids]
    txt <- txt[!is.na(txt)]
    if (length(txt) == 0) NA_character_ else paste(txt, collapse = " | ")
  }, character(1))

  df
}

# Standard downloadHandler for a reactive df – produces a CSV file.
# Usage in server: output$my_dl <- make_csv_download("filename", reactive_df)
# con/include_abstract opzionali: se include_abstract() è TRUE, aggiunge colonna
# ABSTRACT recuperata da ARTICLES tramite PUBMED_ID_clean.
make_csv_download <- function(base_name, reactive_df, con = NULL, include_abstract = NULL) {
  downloadHandler(
    filename = function() paste0(base_name, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- tryCatch(reactive_df(), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) {
        showNotification("No data to download. Run a search first.", type = "warning", duration = 5)
        req(FALSE)
      }
      include_fn <- if (is.function(include_abstract)) include_abstract else function() isTRUE(include_abstract)
      if (!is.null(include_abstract) && isTRUE(include_fn()) && !is.null(con)) {
        df <- attach_abstracts(df, con)
      }
      write.csv(clean_df_for_export(df), file, row.names = FALSE)
    }
  )
}

# Standard plotly config: remove clutter, keep zoom/pan/download.
# Usage: my_plot %>% plotly_clean_config()
plotly_clean_config <- function(p) {
  plotly::config(
    p,
    displaylogo            = FALSE,
    scrollZoom             = TRUE,
    modeBarButtonsToRemove = c(
      "sendDataToCloud", "editInChartStudio",
      "select2d", "lasso2d",
      "toggleSpikelines", "hoverCompareCartesian",
      "hoverClosestCartesian"
    )
  )
}

# Barra di riepilogo conteggi sopra ogni tabella.
# cols: named list, es. list("miRNA" = "MIRNA_NAME", "diseases" = "DISEASE")
# Per PUBMED_ID gestisce le celle multi-valore separate da virgola.
render_count_bar <- function(df, cols) {
  parts <- paste0("<b>", nrow(df), "</b> rows")

  for (label in names(cols)) {
    col <- cols[[label]]
    if (!(col %in% names(df))) next

    if (col == "PUBMED_ID") {
      n <- n_distinct(
        unlist(strsplit(strip_html(paste(df[[col]], collapse = ",")), ",\\s*"))
      )
    } else {
      n <- n_distinct(strip_html(df[[col]]))
    }
    parts <- c(parts, paste0("<b>", n, "</b> unique ", label))
  }

  tags$p(
    style = "color: #555; font-style: italic; margin-bottom: 6px;",
    HTML(paste(parts, collapse = " &nbsp;|&nbsp; "))
  )
}

# Costruisce la proof_clause per glue_sql.
# alias: il prefisso della tabella PROOF_MIRNA_GENE nella query (D, F, ecc.)
build_proof_clause <- function(selected_proofs, con, alias = "D") {
  if (length(selected_proofs) > 0) {
    glue_sql(
      paste0("AND ", alias, ".PROOF_NAME IN ({vals*})"),
      vals = selected_proofs, .con = con
    )
  } else {
    DBI::SQL("")
  }
}

# Observer che aggiorna le scelte del select "proofs{suffix}" filtrando
# proof_table (definito in global.R) in base ai tre filtri avanzati:
#   proof_type{suffix}, proof_strength{suffix}, proof_throughput{suffix}.
# Usato nei tab miRNA-Gene (suffix = "" e "_gene") e miRNA-Reaction (suffix = "_mr").
observeProofsFilters <- function(input, session, suffix = "") {
  strength_id   <- paste0("proof_strength", suffix)
  type_id       <- paste0("proof_type", suffix)
  throughput_id <- paste0("proof_throughput", suffix)
  proofs_id     <- paste0("proofs", suffix)

  observe({
    strength_sel   <- input[[strength_id]]
    type_sel       <- input[[type_id]]
    throughput_sel <- input[[throughput_id]]

    proof_table_filtered <- proof_table
    if (length(type_sel) > 0) {
      proof_table_filtered <- proof_table_filtered %>% filter(TYPE %in% type_sel)
    }
    if (length(strength_sel) > 0) {
      proof_table_filtered <- proof_table_filtered %>% filter(STRENGTH %in% strength_sel)
    }
    if (length(throughput_sel) > 0) {
      proof_table_filtered <- proof_table_filtered %>% filter(THROUGHPUT %in% throughput_sel)
    }

    proofs_name       <- proof_table_filtered$PROOF_NAME
    any_filter_active <- length(type_sel) > 0 ||
                         length(strength_sel) > 0 ||
                         length(throughput_sel) > 0

    updateSelectInput(
      session, proofs_id,
      choices  = proofs_name,
      selected = if (any_filter_active) proofs_name else NULL
    )
  })
}


# =============================================================================
# Upload-list module – reusable file-upload for any selectize input
#
# UI:   upload_list_ui(id)           → compact fileInput below selectize
# Server: upload_list_server(id)     → returns reactive character vector
#
# Parent server observes the reactive and calls updateSelectizeInput.
# File format: .txt or .csv, one item per line, UTF-8.
# =============================================================================
# Internal normalizers ─────────────────────────────────────────────────────────
.norm_mirna <- function(x) {
  # 1. Fix case: mir- → miR-
  x <- gsub("(^|-)mir-", "\\1miR-", x, ignore.case = TRUE, perl = TRUE)
  # 2. Strip any non-hsa organism prefix (3 lowercase letters + hyphen)
  x <- sub("^(?!hsa-)[a-z]{3}-", "", x, perl = TRUE)
  # 3. Add hsa- if still no prefix
  ifelse(!grepl("^hsa-", x, perl = TRUE), paste0("hsa-", x), x)
}
.norm_gene <- function(x) toupper(x)


upload_list_ui <- function(id, inner_ui) {
  ns <- NS(id)
  tagList(
    inner_ui,
    tags$div(
      style = "margin-top:6px;",
      checkboxInput(ns("show_upload"),
                    tagList(icon("upload"), tags$small(" Upload from file")),
                    value = FALSE)
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == true", ns("show_upload")),
      fileInput(
        ns("file"), NULL,
        buttonLabel = tagList(icon("upload"), " Choose file"),
        placeholder = ".txt / .csv / .gmt – one item per line or GMT format",
        accept      = c(".txt", ".csv", ".gmt"),
        width       = "100%"
      ),
      uiOutput(ns("gmt_picker")),
      uiOutput(ns("preview")),
      actionButton(ns("reset"), "Clear", icon = icon("times"),
                   class = "btn-xs btn-default", style = "margin-top:2px;")
    )
  )
}

.process_upload_items <- function(raw_items, norm_fn, valid_choices,
                                  session, session_parent, target_id,
                                  stored, warn_msgs, type) {
  msgs <- character(0)

  n_raw     <- length(raw_items)
  raw_items <- unique(raw_items)
  n_dup     <- n_raw - length(raw_items)
  if (n_dup > 0)
    msgs <- c(msgs, paste0("⚠ ", n_dup, " duplicate(s) removed."))

  norm_items <- norm_fn(raw_items)
  n_norm     <- sum(norm_items != raw_items)
  if (n_norm > 0)
    msgs <- c(msgs, paste0("ℹ ", n_norm, " name(s) normalised (e.g. prefix/case fix)."))

  vc          <- if (is.function(valid_choices)) valid_choices() else valid_choices
  final_items <- norm_items

  if (!is.null(vc) && length(vc) > 0) {
    matched   <- norm_items[norm_items %in% vc]
    unmatched <- norm_items[!norm_items %in% vc]

    arm_stripped    <- character(0)
    still_unmatched <- character(0)
    if (length(unmatched) > 0 && identical(type, "mirna")) {
      has_arm  <- grepl("-[35]p$", unmatched, perl = TRUE)
      stripped <- sub("-[35]p$", "", unmatched[has_arm])
      now_ok   <- stripped[stripped %in% vc]
      arm_stripped <- now_ok
      still_unmatched <- c(unmatched[!has_arm], stripped[!stripped %in% vc])
      if (length(now_ok) > 0)
        msgs <- c(msgs, paste0(
          "⚠ ", length(now_ok),
          " item(s) had 3p/5p arm info stripped (not used in this context): ",
          paste(head(unmatched[has_arm][stripped %in% vc], 5), collapse = ", "),
          if (length(now_ok) > 5) paste0(" … and ", length(now_ok) - 5, " more") else ""
        ))
    } else {
      still_unmatched <- unmatched
    }

    final_items <- unique(c(matched, arm_stripped))

    if (length(still_unmatched) > 0)
      msgs <- c(msgs, paste0(
        "❌ ", length(still_unmatched),
        " item(s) not found in database: ",
        paste(head(still_unmatched, 5), collapse = ", "),
        if (length(still_unmatched) > 5) paste0(" … and ", length(still_unmatched) - 5, " more") else ""
      ))

    updateSelectizeInput(session_parent, target_id,
                         choices = vc, selected = final_items, server = TRUE)
  } else {
    updateSelectizeInput(session_parent, target_id,
                         selected = final_items, server = TRUE)
  }

  stored(final_items)
  warn_msgs(msgs)
}

upload_list_server <- function(id, session_parent, target_id,
                               type          = NULL,
                               valid_choices = NULL,
                               max_preview   = 15L) {
  norm_fn <- switch(type %||% "none",
    mirna = .norm_mirna,
    gene  = .norm_gene,
    function(x) x
  )

  moduleServer(id, function(input, output, session) {
    stored      <- reactiveVal(character(0))
    warn_msgs   <- reactiveVal(character(0))
    gmt_data    <- reactiveVal(NULL)   # named list: gene_set_name -> character vector of genes

    # GMT picker UI: shown only when a GMT file is loaded
    output$gmt_picker <- renderUI({
      gmt <- gmt_data()
      if (is.null(gmt) || length(gmt) == 0) return(NULL)
      tagList(
        tags$div(
          style = "margin:4px 0;",
          tags$small(tags$b(icon("layer-group"),
            paste0(" GMT file: ", length(gmt), " gene set(s) detected"))),
          selectInput(session$ns("gmt_sets"), "Select gene set(s):",
            choices  = names(gmt),
            selected = names(gmt)[1],
            multiple = TRUE,
            width    = "100%"),
          actionButton(session$ns("gmt_load"), "Load selected gene sets",
            icon = icon("check"), class = "btn-xs btn-primary",
            style = "margin-top:4px;")
        )
      )
    })

    observeEvent(input$gmt_load, {
      gmt <- gmt_data()
      req(gmt, input$gmt_sets)
      sel_sets <- input$gmt_sets[input$gmt_sets %in% names(gmt)]
      raw_items <- unique(unlist(gmt[sel_sets], use.names = FALSE))
      .process_upload_items(raw_items, norm_fn, valid_choices,
                            session, session_parent, target_id,
                            stored, warn_msgs, type)
    })

    observeEvent(input$file, {
      req(input$file)
      raw_lines <- tryCatch(
        readLines(input$file$datapath, warn = FALSE, encoding = "UTF-8"),
        error = function(e) character(0)
      )
      raw_lines <- raw_lines[nzchar(trimws(raw_lines))]

      # Detect GMT: tab-separated lines with ≥3 fields
      is_gmt <- length(raw_lines) > 0 &&
        all(vapply(head(raw_lines, 3), function(l) length(strsplit(l, "\t")[[1]]) >= 3, logical(1)))

      if (is_gmt) {
        parsed <- lapply(raw_lines, function(l) {
          fields <- strsplit(l, "\t")[[1]]
          list(name = fields[1], genes = trimws(fields[-(1:2)]))
        })
        gmt_sets <- setNames(
          lapply(parsed, `[[`, "genes"),
          vapply(parsed, `[[`, character(1), "name")
        )
        gmt_data(gmt_sets)
        return()   # wait for user to pick sets
      }

      gmt_data(NULL)
      raw_items <- raw_lines[nzchar(trimws(raw_lines))]
      .process_upload_items(raw_items, norm_fn, valid_choices,
                            session, session_parent, target_id,
                            stored, warn_msgs, type)
    })

    observeEvent(input$reset, {
      stored(character(0))
      warn_msgs(character(0))
      gmt_data(NULL)
      vc <- if (is.function(valid_choices)) valid_choices() else valid_choices
      if (!is.null(vc)) {
        updateSelectizeInput(session_parent, target_id,
                             choices = vc, selected = character(0), server = TRUE)
      } else {
        updateSelectizeInput(session_parent, target_id, selected = character(0))
      }
      shinyjs::reset("file")
      updateCheckboxInput(session, "show_upload", value = FALSE)
    })

    output$preview <- renderUI({
      items <- stored()
      msgs  <- warn_msgs()
      if (length(items) == 0 && length(msgs) == 0) return(NULL)
      n     <- length(items)
      shown <- head(items, max_preview)
      tagList(
        if (length(msgs) > 0)
          tags$div(
            style = paste0("font-size:11px; background:#fff8e1; border-left:3px solid #f39c12;",
                           " padding:4px 8px; margin-bottom:4px; border-radius:0 4px 4px 0;"),
            lapply(msgs, function(m) tags$div(m))
          ),
        if (n > 0) tags$div(
          tags$p(
            style = "font-size:12px; color:#555; margin:4px 0 2px;",
            tags$b(n), " items accepted",
            if (n > max_preview)
              tags$span(style = "color:#888; font-style:italic;",
                        paste0(" – showing first ", max_preview))
          ),
          tags$div(
            style = paste0("font-size:11px; background:#f5f5f5; border-radius:4px;",
                           " padding:4px 8px; max-height:70px; overflow-y:auto;",
                           " word-break:break-all; color:#333;"),
            paste(shown, collapse = ", ")
          )
        )
      )
    })
  })
}


link_info_fun <- function(data, type) {
  # Versione vettorializzata: accetta una singola stringa oppure un intero
  # vettore (es. tutta una colonna di un data.frame).
  # Fast-path: se nessun elemento contiene una virgola (query non raggruppate),
  # si evita lo split e si usa sprintf vettoriale – ordini di grandezza più veloce.
  # htmlEscape previene XSS nei data-attribute e nel testo visibile.
  if (length(data) == 0) return(character(0))

  if (!any(grepl(",", data, fixed = TRUE))) {
    return(sprintf(
      '<a href="#" class="%s-link" data-%s="%s">%s</a>',
      type, type,
      htmltools::htmlEscape(data, attribute = TRUE),
      htmltools::htmlEscape(data)
    ))
  }

  # Slow-path: celle contenenti ID separati da virgola (query raggruppate).
  ids_list <- strsplit(data, ",\\s*")
  vapply(ids_list, function(ids) {
    if (length(ids) == 0) return("")
    paste(
      sprintf(
        '<a href="#" class="%s-link" data-%s="%s">%s</a>',
        type, type,
        htmltools::htmlEscape(ids, attribute = TRUE),
        htmltools::htmlEscape(ids)
      ),
      collapse = ", "
    )
  }, character(1))
}

strip_html <- function(x) gsub("<.*?>", "", x)

split_unique <- function(values, sep = ",") {
  values |>
    unlist() |>
    strsplit(sep) |>
    unlist() |>
    trimws() |>
    (\(x) x[nzchar(x)])() |>
    unique()
}

preview_render_table_function <-function(preview_data){
  DT::datatable(
    preview_data,
    options = list(
      pageLength = 10,
      selection = "multiple",
      dom = 'Bfrtip',
      buttons = list(
        'csv', 'excel',
        list(extend = 'colvis', text = 'Columns')
      )
    ),
    filter = "top",
    escape = FALSE,
    extensions = "Buttons"
  )  
}

toggle_inputs <- function(condition, disable_ids = NULL, enable_ids = NULL) {
  observe({
    if (isTRUE(condition())) {
      if (!is.null(disable_ids)) lapply(disable_ids, shinyjs::disable)
      if (!is.null(enable_ids))  lapply(enable_ids, shinyjs::enable)
    } else {
      if (!is.null(disable_ids)) lapply(disable_ids, shinyjs::enable)
      if (!is.null(enable_ids))  lapply(enable_ids, shinyjs::disable)
    }
  })
}

reset_select_when <- function(condition, session, input_id) {
  observe({
    if (isTRUE(tryCatch(condition(), error = function(e) FALSE))) {
      updateSelectInput(session, input_id, selected = character(0))
    }
  })
}

norm_base <- function(x) {
  b <- x[1]
  if (is.na(b) || b == 0) return(rep(NA_real_, length(x)))
  x / b
}

store_png <- function(input_id, my_data, input, data_key = input_id) {
  observeEvent(input[[input_id]], {
    val <- input[[input_id]]
    if (!is.null(val) && nzchar(val)) my_data[[data_key]] <- val
  })
}

# Variant-aware PNG store: routes raw capture to one of two data_keys
# based on the current value of a node-type radioButton input.
# variant_map: named list e.g. list(mirna_mirna = "X_mirna_png", disease_disease = "X_disease_png")
store_png_variant <- function(input_id, my_data, input, node_input, variant_map) {
  observeEvent(input[[input_id]], {
    val <- input[[input_id]]
    if (!is.null(val) && nzchar(val)) {
      node <- isolate(input[[node_input]]) %||% names(variant_map)[1L]
      key  <- variant_map[[node]] %||% variant_map[[1L]]
      if (!is.null(key)) my_data[[key]] <- val
    }
  })
}

parse_gmt_file <- function(path) {
  raw_lines <- tryCatch(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function(e) character(0)
  )
  raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
  if (length(raw_lines) == 0) return(NULL)
  parsed <- lapply(raw_lines, function(l) {
    fields <- strsplit(l, "\t")[[1]]
    list(name = fields[1], genes = trimws(fields[-(1:2)][nzchar(trimws(fields[-(1:2)]))]))
  })
  setNames(lapply(parsed, `[[`, "genes"), vapply(parsed, `[[`, character(1), "name"))
}

render_gmt_status_ui <- function(gmt_reactive, max_preview = 8) {
  renderUI({
    gmt <- gmt_reactive()
    if (is.null(gmt)) return(tags$p(style = "color:#888; font-size:12px;",
      "No GMT loaded. Upload a .gmt file above."))
    tags$div(
      tags$span(class = "label label-success", style = "font-size:12px; padding:4px 8px;",
        icon("check"), paste0(" ", length(gmt), " gene sets loaded")),
      tags$ul(style = "margin-top:6px; font-size:12px; max-height:80px; overflow-y:auto;",
        lapply(head(names(gmt), max_preview), function(n)
          tags$li(n, tags$small(style = "color:#888;", paste0(" (", length(gmt[[n]]), " genes)")))),
        if (length(gmt) > max_preview) tags$li(style = "color:#888;",
          paste0("… and ", length(gmt) - max_preview, " more")) else NULL
      )
    )
  })
}

.network_layout_choices <- c(
  "ForceAtlas2" = "forceAtlas2Based", "BarnesHut" = "barnesHut",
  "Repulsion" = "repulsion", "Hierarchical" = "hierarchical",
  "Bipartite" = "bipartite", "Static" = "no physics"
)

.cotarget_palette_choices <- c("YlOrRd","YlGnBu","Blues","Reds","Oranges","RdPu","PuRd","Greens")
.heatmap_palette_choices <- c("YlGnBu","Blues","Greens","Reds","Oranges","RdYlBu","RdBu","PuRd","YlOrRd","Purples")

.cotarget_controls_ui <- function(prefix) {
  fluidRow(
    column(3, checkboxInput(paste0(prefix, "_cluster"),   "Cluster rows/cols",  value = FALSE)),
    column(3, checkboxInput(paste0(prefix, "_mask_zero"), "Mask zeros (white)", value = TRUE)),
    column(3, checkboxInput(paste0(prefix, "_mask_diag"), "Mask diagonal",      value = FALSE)),
    column(3, selectInput(paste0(prefix, "_palette"), "Color scale:",
                          choices = .cotarget_palette_choices,
                          selected = "YlOrRd"))
  )
}