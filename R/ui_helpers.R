# Funzione per creare un selectInput compatto con etichetta e menu a tendina sulla stessa riga
compactSelect <- function(label, inputId, choices,multiple,width) {
  tags$div(
    style = "display:flex; align-items:center; margin-bottom:6px;",
    tags$label(label, style = "width:180px; margin:0;font-weight: normal"),
    selectInput(inputId, NULL, choices,selected = NULL,multiple = multiple,width=width)
  )
}
#### UI per i filtri "Evidence" (usato nei tab miRNA-Gene e miRNA-Reaction)
# suffix: "" per by_mirna, "_gene" per by_gene, "_mr" per miRNA-Reaction
# Gli ID di input sono costruiti come: proofs{suffix}, show_advanced_evidence{suffix},
# proof_type{suffix}, proof_strength{suffix}, proof_throughput{suffix}, reset_mg{suffix}
evidenceFiltersUI <- function(suffix = "", reset_label = "Reset filters") {
  proofs_id   <- paste0("proofs", suffix)
  advanced_id <- paste0("show_advanced_evidence", suffix)
  reset_id    <- paste0("reset_mg", suffix)
  match_id <- paste0("proof_match", suffix)
  tagList(
    h5(icon("check-circle"), "Evidence filters"),
    selectInput(proofs_id, "Select evidence:", choices = proofs_all, multiple = TRUE),
    radioButtons(match_id, "Match mode:",
      choices = c("Any selected proof" = "any", "All selected proofs" = "all"),
      selected = "any", inline = TRUE),
    checkboxInput(advanced_id, "Advanced filters"),
    conditionalPanel(
      condition = paste0("input.", advanced_id),
      compactSelect("Type:",            paste0("proof_type", suffix),       proof_choices_type,       multiple = TRUE, width = "140px"),
      compactSelect("Strength:",        paste0("proof_strength", suffix),   proof_choices_strength,   multiple = TRUE, width = "140px"),
      compactSelect("High-throughput:", paste0("proof_throughput", suffix), proof_choices_Throughput, multiple = TRUE, width = "140px")
    ),
    actionButton(reset_id, reset_label, icon = icon("undo"))
  )
}

expressionFilterUI <- function(suffix = "") {
  mode_id      <- paste0("expression_filter_mode",    suffix)
  tissue_id    <- paste0("expression_tissue_select",  suffix)
  cancer_id    <- paste0("expression_cancer_filter",  suffix)
  disease_id   <- paste0("expression_primary_disease",suffix)
  subtype_id   <- paste0("expression_disease_subtype",suffix)
  celllines_id <- paste0("expression_cell_lines",     suffix)
  reset_id     <- paste0("reset_expression",          suffix)
  nTPM_id      <- paste0("nTPM_expression_filter",    suffix)
  nTPM_id2     <- paste0("nTPM_2_expression_filter",    suffix)

  tagList(
    h5(icon("flask"), "Expression context filter"),
    radioButtons(mode_id, NULL,
      choices = c("None" = "none", "Tissue" = "tissue", "Cell line" = "cell_line"),
      selected = "none", inline = TRUE),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'tissue'", mode_id),
      selectizeInput(tissue_id, "Select tissues:", choices = tissue_choices, multiple = TRUE),
      numericInput(nTPM_id, "Select minimum nTPM expression:", 0.1,
                  min=0,max=NA,step=0.1)
      
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'cell_line'", mode_id),
      radioButtons(cancer_id, "Cancer filter:",
        choices = c("All" = "all", "Cancer only" = "yes", "Non-cancer only" = "no"),
        selected = "all", inline = TRUE),
      selectizeInput(disease_id, "Primary disease:", choices = cell_line_primary_diseases, multiple = TRUE),
      selectizeInput(subtype_id, "Disease subtype:", choices = NULL, multiple = TRUE),
      selectizeInput(celllines_id, "Cell lines (optional):", choices = NULL, multiple = TRUE),
      numericInput(nTPM_id2, "Select minimum nTPM expression:", 0.1,min=0,max=NA,step=0.1)
    ),
    actionButton(reset_id, "Reset expression filter", icon = icon("undo"))
  )
}

visualizationHelpBox <- function(
    output_id,
    title = "Visualization help",
    icon_name = "info-circle",
    status = "warning",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
) {
  # Pull static content from global env if available (avoids renderUI roundtrip).
  # Falls back to htmlOutput() for any dynamic output_id.
  content <- tryCatch(
    get(output_id, envir = .GlobalEnv),
    error = function(e) htmlOutput(output_id)
  )

  box(
    title = tagList(icon(icon_name), title),
    status = status,
    solidHeader = TRUE,
    collapsible = collapsible,
    collapsed = collapsed,
    width = width,
    content
  )
}

##################################################
database_upload_tab_fixed_ids <- function(
    title,
    source_id,
    upload_id,
    update_id,
    restore_id,
    apply_id,
    help_id,
    preview_output_id,
    label,
    show_3p5p_button 
) {
  tabPanel(
    title = tagList(title),
    
    radioButtons(
      source_id,
      paste("Choose", label, "database:"),
      choices = setNames(
        c("default", "custom"),
        c(
          paste("Use default", label, "table"),
          paste("Upload custom", label, "table")
        )
      ),
      selected = "default"
    ),
    
    conditionalPanel(
      condition = sprintf("input.%s == 'custom'", source_id),
      br(),
      visualizationHelpBox(help_id, collapsible = TRUE),
      actionButton(
        update_id,
        "Update database",
        icon = icon("play"),
        class = "btn-success"
      ),
      actionButton(
        restore_id,
        "Restore original database",
        icon = icon("play"),
        class = "btn-success"
      ),
      br(),
      
      if (show_3p5p_button) {
        checkboxInput("use_3p5p_information", "Use 3p/5p information", value = FALSE)
      },
      br(),
      fileInput(
        upload_id,
        "Select table file",
        accept = c(".csv", ".tsv", ".xlsx")
      ),
    ),
    
    br(),
    
    actionButton(
      apply_id,
      "Show table",
      icon = icon("play"),
      class = "btn-success"
    ),
    
    hr(),
    
    h5(icon("table"), "Table preview"),
    DTOutput(preview_output_id)
  )
}

