# === server_toggles.R ===
# Toggle enable/disable e reset dei selectInput cross-tab

# -- miRNA-Gene --------------------------------------------
toggle_inputs(
  condition = function() length(input$mirna_gene_select) > 0,
  disable_ids = c("filter_mirna_common_checkbox", "use_all_mirnas_from_mda")
)

toggle_inputs(
  condition = function() isTRUE(input$mga_use_all_scratch),
  disable_ids = "mirna_gene_select_scratch"
)

toggle_inputs(
  condition = function() isTRUE(input$mga_use_all_genes) || isTRUE(input$mga_only_metabolic_genes),
  disable_ids = "gene_search_select"
)

toggle_inputs(
  condition = function() length(input$gene_search_select) > 0 || isTRUE(input$mga_only_metabolic_genes),
  disable_ids = "mga_use_all_genes"
)

toggle_inputs(
  condition = function() length(input$gene_search_select) > 0 || isTRUE(input$mga_use_all_genes),
  disable_ids = "mga_only_metabolic_genes"
)

toggle_inputs(
  condition = function() input$filter_mirna_common_checkbox,
  disable_ids = "mirna_gene_select"
)

toggle_inputs(
  condition = function() input$use_all_mirnas_from_mda,
  disable_ids = "mirna_gene_select"
)

toggle_inputs(
  condition = function() input$filter_genes_reaction,
  disable_ids = "gene_reaction_select"
)

toggle_inputs(
  condition = function() length(input$gene_reaction_select) > 0,
  disable_ids = "filter_genes_reaction"
)

toggle_inputs(
  condition = function() isTRUE(input$use_all_mirnas_from_mga) || isTRUE(input$filter_rmirna_common_checkbox),
  disable_ids = "mirna_reaction_select"
)

toggle_inputs(
  condition = function() length(input$mirna_reaction_select) > 0,
  disable_ids = c("filter_rmirna_common_checkbox", "use_all_mirnas_from_mga")
)


# -- miRNA-Disease -----------------------------------------
toggle_inputs(
  condition = function() length(input$mirna_select) > 0,
  disable_ids = "filter_mirna_all_checkbox"
)

toggle_inputs(
  condition = function() input$filter_mirna_all_checkbox,
  disable_ids = "mirna_select"
)

toggle_inputs(
  condition = function() length(input$disease_select) > 0,
  disable_ids = "filter_disease_checkbox"
)

toggle_inputs(
  condition = function() input$filter_disease_checkbox,
  disable_ids = "disease_select"
)

# -- Reset select ------------------------------------------
reset_select_when(
  condition = function() isTRUE(input$mga_use_all_scratch),
  session = session,
  input_id = "mirna_gene_select_scratch"
)

reset_select_when(
  condition = function() isTRUE(input$mga_use_all_genes),
  session = session,
  input_id = "gene_search_select"
)

reset_select_when(
  condition = function() isTRUE(input$filter_mirna_common_checkbox) || isTRUE(input$use_all_mirnas_from_mda),
  session = session,
  input_id = "mirna_gene_select"
)

reset_select_when(
  condition = function() isTRUE(input$use_all_mirnas_from_mga) || isTRUE(input$filter_rmirna_common_checkbox),
  session = session,
  input_id = "mirna_reaction_select"
)


reset_select_when(
  condition = function() input$filter_mirna_all_checkbox,
  session = session,
  input_id = "mirna_select"
)

# -- Reset buttons -----------------------------------------
observeEvent(input$reset_mirna, {
  updateSelectInput(session, "mirna_select", selected = character(0))
  updateCheckboxInput(session, "filter_mirna_all_checkbox", value = FALSE)
})

observeEvent(input$reset_disease, {
  updateSelectInput(session, "do_category_filter",                 selected = "")
  updateSelectizeInput(session, "disease_select", choices = disease_choices,
                       selected = character(0), server = TRUE)
  updateCheckboxInput(session, "filter_disease_checkbox",          value = TRUE)
  updateCheckboxInput(session, "filter_diseaseancestor_checkbox",  value = FALSE)
})

observeEvent(input$reset_mirna_gene, {
  updateSelectizeInput(session, "mirna_gene_select",            selected = character(0))
  updateSelectizeInput(session, "mirna_gene_select_scratch",    selected = character(0))
  updateSelectizeInput(session, "gene_search_select",           selected = character(0))
  updateCheckboxInput(session,  "mga_use_all_scratch",          value = FALSE)
  updateCheckboxInput(session,  "mga_use_3p5p_scratch",         value = FALSE)
  updateCheckboxInput(session,  "filter_mirna_common_checkbox", value = FALSE)
  updateCheckboxInput(session,  "use_all_mirnas_from_mda",      value = FALSE)
  updateCheckboxInput(session,  "mga_use_all_genes",            value = FALSE)
  updateCheckboxInput(session,  "mga_only_metabolic_genes",     value = FALSE)
})

# -- Mutual exclusion checkboxes (miRNA-Reaction) ----------
observeEvent(input$mirna_reaction_select, {
  if (length(input$mirna_reaction_select) > 0) {
    updateCheckboxInput(session, "filter_rmirna_common_checkbox", value = FALSE)
    updateCheckboxInput(session, "use_all_mirnas_from_mga",       value = FALSE)
  }
}, ignoreInit = TRUE)

observeEvent(input$filter_rmirna_common_checkbox, {
  if (isTRUE(input$filter_rmirna_common_checkbox) && isTRUE(input$use_all_mirnas_from_mga)) {
    updateCheckboxInput(session, "use_all_mirnas_from_mga", value = FALSE)
  }
}, ignoreInit = TRUE)

observeEvent(input$use_all_mirnas_from_mga, {
  if (isTRUE(input$use_all_mirnas_from_mga) && isTRUE(input$filter_rmirna_common_checkbox)) {
    updateCheckboxInput(session, "filter_rmirna_common_checkbox", value = FALSE)
  }
}, ignoreInit = TRUE)

# -- Pathway Annotation / Over-representation --------------
toggle_inputs(
  condition = function() isTRUE(input$use_ontology_all_mirnas_from_mga),
  disable_ids = "ontology_mirna_select"
)

toggle_inputs(
  condition = function() input$filter_ontology_genes_reaction,
  disable_ids = "ontology_gene_select"
)

reset_select_when(
  condition = function() isTRUE(input$use_ontology_all_mirnas_from_mga),
  session = session,
  input_id = "ontology_mirna_select"
)

reset_select_when(
  condition = function() input$filter_ontology_genes_reaction,
  session = session,
  input_id = "ontology_gene_select"
)
