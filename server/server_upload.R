# === server_upload.R ===
# Gestione upload/update/restore database: miRNA-Disease, Ontology, Metabolism

refresh_reaction_gene_map <- function(con) {
  tg <- dbGetQuery(con,
    "SELECT R.NAME AS REACTION_NAME, COUNT(DISTINCT RG.GENE_ID) AS TOTAL_GENES
       FROM REACTIONS R JOIN REACTIONS_GENES RG ON R.REACTION_ID = RG.REACTION_ID
      GROUP BY R.REACTION_ID, R.NAME")
  reaction_total_genes_map <<- setNames(as.integer(tg$TOTAL_GENES), tg$REACTION_NAME)
}

# ── miRNA-Disease ─────────────────────────────────────────
botton_update_mda_database  <- reactiveVal(FALSE)
botton_restore_mda_database <- reactiveVal(FALSE)

observeEvent(input$update_mda_database, {
  tryCatch({
    withProgress(message = "Updating miRNA-disease database...", value = 0, {

      db <- mda_preview_data()
      req(nrow(db) > 0)
      incProgress(0.2, detail = "Matching miRNA IDs")

      # ── Resolve MIRNA_ID ──────────────────────────────────────
      if (isFALSE(input$use_3p5p_information)) {
        mirna_map  <- dbGetQuery(con, "SELECT MIRNA_ID, MIRNA_PREMATURE AS MIRNA_NAME FROM MIRNAS")
        df_joined  <- dplyr::inner_join(db, mirna_map, by = "MIRNA_NAME")
      } else {
        mirna_map  <- dbGetQuery(con, "SELECT MIRNA_ID, MIRNA_PREMATURE, MIRNA_ARM FROM MIRNAS")
        df_joined  <- dplyr::inner_join(db, mirna_map, by = c("MIRNA_PREMATURE", "MIRNA_ARM"))
      }
      incProgress(0.2, detail = "Matching disease IDs")

      # ── Resolve DISEASE_ID ────────────────────────────────────
      disease_map <- dbGetQuery(con, "SELECT DISEASE_ID, DISEASE FROM DISEASES")
      df_joined   <- dplyr::inner_join(df_joined, disease_map, by = "DISEASE")
      incProgress(0.2, detail = "Writing update table")

      keep_cols <- c("PUBMED_ID", "MIRNA_ID", "DISEASE_ID")
      if ("SAMPLE_TYPE" %in% names(df_joined)) keep_cols <- c(keep_cols, "SAMPLE_TYPE")
      df_update <- dplyr::distinct(df_joined, dplyr::across(dplyr::all_of(keep_cols)))
      shiny::validate(shiny::need(
        nrow(df_update) > 0,
        "No valid associations after ID matching. Check that miRNA and disease names match the database."
      ))

      dbWriteTable(con, "MIRNAS_DISEASES_ARTICLES_update", df_update, overwrite = TRUE)
      botton_update_mda_database(TRUE)
      n_loaded <- nrow(df_update)
      incProgress(0.4)
    })
    showNotification(
      paste0("miRNA-disease database updated: ", n_loaded, " associations loaded."),
      type = "message", duration = 4
    )
  }, error = function(e) {
    showNotification(paste("Error updating miRNA-disease database:", e$message), type = "error", duration = 8)
  })
}, ignoreInit = TRUE)

observeEvent(input$restore_mda_database, {
  tryCatch({
    withProgress(message = "Restoring original miRNA-disease database...", value = 0.5, {
      botton_update_mda_database(FALSE)
      botton_restore_mda_database(TRUE)
      # Drop the update table if present so stale data is not reused
      if ("MIRNAS_DISEASES_ARTICLES_update" %in% dbListTables(con))
        dbExecute(con, "DROP TABLE MIRNAS_DISEASES_ARTICLES_update")
      incProgress(0.5)
    })
    shinyjs::reset("upload_database")
    showNotification("miRNA-disease database restored successfully.", type = "message", duration = 4)
  }, error = function(e) {
    showNotification(paste("Error restoring miRNA-disease database:", e$message), type = "error", duration = 8)
  })
}, ignoreInit = TRUE)

show_mda_upload_table <- reactiveVal(FALSE)
observeEvent(input$apply_mda_database, {
  show_mda_upload_table(!show_mda_upload_table())
})

mda_preview_data <- eventReactive(input$apply_mda_database, {
  if (input$use_3p5p_information == FALSE) {
    required_cols <- c("MIRNA_NAME", "DISEASE", "PUBMED_ID")
  } else {
    required_cols <- c("MIRNA_PREMATURE", "MIRNA_ARM", "DISEASE", "PUBMED_ID")
  }
  base_query <- query_md[["query_db_general"]]

  is_restore <- isolate(botton_restore_mda_database())
  botton_restore_mda_database(FALSE)

  preview_table_function(
    input$db_source,
    is_restore,
    input$upload_database,
    required_cols,
    base_query
  )
})

output$mda_database_preview <- renderDT({
  req(show_mda_upload_table())
  req(mda_preview_data())
  preview_render_table_function(mda_preview_data())
})

# ── Ontology ──────────────────────────────────────────────
show_ontology_upload_table <- reactiveVal(FALSE)
observeEvent(input$apply_ontology_database, {
  show_ontology_upload_table(!show_ontology_upload_table())
})

ontology_preview_data <- eventReactive(input$apply_ontology_database, {
  required_cols <- c("GENE_NAME", "PATHWAY_NAME")
  base_query    <- query_ontology[["query_db_general"]]
  is_restore    <- isolate(input$restore_ontology_database) > 0

  preview_table_function(
    input$db_ontology_source,
    is_restore,
    input$upload_ontology_database,
    required_cols,
    base_query
  )
})

output$ontology_database_preview <- renderDT({
  req(show_ontology_upload_table())
  req(ontology_preview_data())
  preview_render_table_function(ontology_preview_data())
})

observeEvent(input$update_ontology_database, {
  tryCatch({
    withProgress(message = "Updating ontology database...", value = 0, {

      db <- ontology_preview_data()
      incProgress(0.2)

      df_pathways <- db %>%
        dplyr::distinct(GENE_NAME, PATHWAY_NAME) %>%
        dplyr::distinct(PATHWAY_NAME) %>%
        dplyr::mutate(PATHWAY_ID = dplyr::row_number()) %>%
        dplyr::select(PATHWAY_ID, PATHWAY_NAME)
      incProgress(0.3)

      df_pathways_genes <- db %>%
        dplyr::distinct(GENE_NAME, PATHWAY_NAME) %>%
        dplyr::inner_join(df_pathways, by = "PATHWAY_NAME") %>%
        dplyr::inner_join(my_data$genes, by = "GENE_NAME") %>%
        dplyr::select(PATHWAY_ID, GENE_ID)
      incProgress(0.3)

      tables <- dbListTables(con)
      if (!("ONTOLOGY_PATHWAYS_GENES_origin" %in% tables)) {
        dbExecute(con, "CREATE TABLE ONTOLOGY_PATHWAYS_GENES_origin AS SELECT * FROM ONTOLOGY_PATHWAYS_GENES;")
        dbExecute(con, "CREATE TABLE ONTOLOGY_PATHWAYS_origin AS SELECT * FROM ONTOLOGY_PATHWAYS;")
      }
      dbWriteTable(con, "ONTOLOGY_PATHWAYS_GENES", df_pathways_genes, overwrite = TRUE)
      dbWriteTable(con, "ONTOLOGY_PATHWAYS", df_pathways, overwrite = TRUE)

      ensure_db_indices(con)
      refresh_ontology_globals(con)
      incProgress(0.2)
    })
    n_genes_not_mapped <- nrow(dplyr::anti_join(
      dplyr::distinct(db, GENE_NAME), my_data$genes, by = "GENE_NAME"))
    msg <- "Ontology database updated successfully."
    if (n_genes_not_mapped > 0)
      msg <- paste0(msg, " Warning: ", n_genes_not_mapped,
                    " gene name(s) not found in database and were excluded.")
    showNotification(msg, type = if (n_genes_not_mapped > 0) "warning" else "message", duration = 6)
  }, error = function(e) {
    showNotification(paste("Error updating ontology database:", e$message), type = "error", duration = 8)
  })
}, ignoreInit = TRUE)

observeEvent(input$restore_ontology_database, {
  withProgress(message = "Restoring original ontology database...", value = 0.4, {
    restore_table(con, "ONTOLOGY_PATHWAYS_GENES")
    restore_table(con, "ONTOLOGY_PATHWAYS")
    refresh_ontology_globals(con)
    incProgress(0.6)
  })
  showNotification("Ontology database restored successfully.", type = "message", duration = 4)
}, ignoreInit = TRUE)

# ── Metabolism ────────────────────────────────────────────
botton_update_metabolism_database  <- reactiveVal(FALSE)
botton_restore_metabolism_database <- reactiveVal(FALSE)

show_metabolism_upload_table <- reactiveVal(FALSE)
observeEvent(input$apply_metabolism_database, {
  show_metabolism_upload_table(!show_metabolism_upload_table())
})

metabolism_preview_data <- eventReactive(input$apply_metabolism_database, {
  # Minimum required: reaction metadata. "Genes" (comma-sep) is optional but needed
  # to rebuild gene-reaction associations when uploading a custom metabolic model.
  required_cols <- c("NAME", "REACTION_FORMULA", "GPR", "SUBSYSTEM", "COMPARTMENT")
  base_query    <- query_metabolism[["query_db_general"]]
  is_restore    <- isolate(botton_restore_metabolism_database())
  botton_restore_metabolism_database(FALSE)

  preview_table_function(
    input$db_metabolism_source,
    is_restore,
    input$upload_metabolism_database,
    required_cols,
    base_query
  )
})

output$metabolism_database_preview <- renderDT({
  req(show_metabolism_upload_table())
  req(metabolism_preview_data())
  preview_render_table_function(metabolism_preview_data())
})

observeEvent(input$update_metabolism_database, {
  tryCatch({
    withProgress(message = "Updating metabolism database...", value = 0, {
      db <- metabolism_preview_data()
      req(nrow(db) > 0)

      # Assign surrogate REACTION_ID (row position) – user does not need to provide it
      db$REACTION_ID <- seq_len(nrow(db))

      incProgress(0.15, detail = "Backing up original tables")
      tables <- dbListTables(con)
      if (!("REACTIONS_origin" %in% tables))
        dbExecute(con, "CREATE TABLE REACTIONS_origin AS SELECT * FROM REACTIONS;")
      if (!("REACTIONS_GENES_origin" %in% tables))
        dbExecute(con, "CREATE TABLE REACTIONS_GENES_origin AS SELECT * FROM REACTIONS_GENES;")

      incProgress(0.25, detail = "Writing REACTIONS")
      keep_cols <- intersect(c("REACTION_ID","NAME","REACTION_FORMULA","GPR",
                               "SUBSYSTEM","COMPARTMENT","HUMAN_ID","CONNECTION_GROUP"),
                             colnames(db))
      dbWriteTable(con, "REACTIONS", db[, keep_cols, drop = FALSE], overwrite = TRUE)

      # If "Genes" column present → rebuild REACTIONS_GENES and gene-count cache
      n_genes_not_mapped_mr <- 0L
      if ("Genes" %in% colnames(db)) {
        incProgress(0.2, detail = "Rebuilding gene–reaction associations")
        genes_map <- dbGetQuery(con, "SELECT GENE_ID, GENE_NAME FROM GENES")

        df_long <- db %>%
          dplyr::select(REACTION_ID, Genes) %>%
          tidyr::separate_rows(Genes, sep = ",\\s*") %>%
          dplyr::rename(GENE_NAME = Genes) %>%
          dplyr::filter(nchar(trimws(GENE_NAME)) > 0)

        n_genes_not_mapped_mr <- length(setdiff(unique(df_long$GENE_NAME), genes_map$GENE_NAME))

        df_rg <- df_long %>%
          dplyr::inner_join(genes_map, by = "GENE_NAME") %>%
          dplyr::select(REACTION_ID, GENE_ID) %>%
          dplyr::distinct()

        shiny::validate(shiny::need(nrow(df_rg) > 0,
          "No gene–reaction associations found. Check that gene names in 'Genes' column match the database."))

        dbWriteTable(con, "REACTIONS_GENES", df_rg, overwrite = TRUE)

        # Refresh cached TOTAL_GENES per reaction
        refresh_reaction_gene_map(con)
      }

      ensure_db_indices(con)

      incProgress(0.2, detail = "Refreshing dropdowns")
      new_subsystems     <- sort(unique(db$SUBSYSTEM[!is.na(db$SUBSYSTEM) & db$SUBSYSTEM != ""]))
      new_reaction_names <- sort(unique(db$NAME[!is.na(db$NAME)]))
      subsystem_choices     <<- new_subsystems
      reaction_name_choices <<- new_reaction_names

      botton_update_metabolism_database(TRUE)
      incProgress(0.2)
    })
    msg <- if ("Genes" %in% colnames(db)) {
      base <- "Metabolism database updated: reactions and gene associations replaced."
      if (n_genes_not_mapped_mr > 0)
        paste0(base, " Warning: ", n_genes_not_mapped_mr,
               " gene name(s) not found in database and were excluded.")
      else base
    } else {
      "Metabolism database updated: reaction metadata only (gene associations unchanged)."
    }
    showNotification(msg, type = if (n_genes_not_mapped_mr > 0) "warning" else "message", duration = 5)
  }, error = function(e) {
    showNotification(paste("Error updating metabolism database:", e$message), type = "error", duration = 8)
  })
}, ignoreInit = TRUE)

observeEvent(input$restore_metabolism_database, {
  tryCatch({
    withProgress(message = "Restoring original metabolism database...", value = 0.5, {
      botton_restore_metabolism_database(TRUE)
      restore_table(con, "REACTIONS")
      restore_table(con, "REACTIONS_GENES")
      # Refresh cached reaction→gene-count lookup (lost when table was dropped)
      refresh_reaction_gene_map(con)
      incProgress(0.5)
    })
    shinyjs::reset("upload_metabolism_database")
    showNotification("Metabolism database restored successfully.", type = "message", duration = 4)
  }, error = function(e) {
    showNotification(paste("Error restoring metabolism database:", e$message), type = "error", duration = 8)
  })
}, ignoreInit = TRUE)
