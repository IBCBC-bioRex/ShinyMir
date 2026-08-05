#####################################

exec_stmts <- function(con, stmts) {
  for (s in stmts) tryCatch(dbExecute(con, s), error = function(e) message(conditionMessage(e)))
}

# Returns TRUE if REACTIONS_v view exists in the DB connection.
reactions_view_exists <- function(con) {
  tryCatch({
    dbGetQuery(con, "SELECT 1 FROM REACTIONS_v LIMIT 0")
    TRUE
  }, error = function(e) FALSE)
}

# Creates (or re-creates after table DROP/overwrite) all performance-critical indices.
# Idempotent: IF NOT EXISTS means safe to call multiple times.
# Call at: app startup, after restore_table(), after any dbWriteTable(..., overwrite=TRUE)
# on indexed tables.
ensure_db_indices <- function(con) {
  stmts <- c(
    "CREATE INDEX IF NOT EXISTS idx_mga_mirna ON MIRNAS_GENES_ARTICLES(MIRNA_ID)",
    "CREATE INDEX IF NOT EXISTS idx_mga_gene  ON MIRNAS_GENES_ARTICLES(GENE_ID)",
    "CREATE INDEX IF NOT EXISTS idx_mga_proof ON MIRNAS_GENES_ARTICLES(PROOF_ID)",
    "CREATE INDEX IF NOT EXISTS idx_mda_mirna ON MIRNAS_DISEASES_ARTICLES(MIRNA_ID)",
    "CREATE INDEX IF NOT EXISTS idx_mda_dis   ON MIRNAS_DISEASES_ARTICLES(DISEASE_ID)",
    "CREATE INDEX IF NOT EXISTS idx_rg_gene   ON REACTIONS_GENES(GENE_ID)",
    "CREATE INDEX IF NOT EXISTS idx_rg_react  ON REACTIONS_GENES(REACTION_ID)"
  )
  exec_stmts(con, stmts)

  # ── Shadow tables + views for metabolism upload ───────────────────────────
  # REACTIONS_update / REACTIONS_GENES_update: empty by default, filled on upload.
  # REACTIONS_v / REACTIONS_GENES_v: replace-semantic views — show update table
  # when it has rows, otherwise fall back to base table.
  shadow_stmts <- list(
    "CREATE TABLE IF NOT EXISTS REACTIONS_update (
       REACTION_ID INTEGER, HUMAN_ID TEXT, NAME TEXT, FORMULA TEXT,
       GPR TEXT, SUBSYSTEM TEXT, COMPARTMENT TEXT)",
    "CREATE TABLE IF NOT EXISTS REACTIONS_GENES_update (
       REACTION_ID INTEGER, GENE_ID INTEGER)",
    "DROP VIEW IF EXISTS REACTIONS_v",
    "CREATE VIEW REACTIONS_v AS
       SELECT * FROM REACTIONS_update
       UNION ALL
       SELECT * FROM REACTIONS
       WHERE (SELECT COUNT(*) FROM REACTIONS_update) = 0",
    "DROP VIEW IF EXISTS REACTIONS_GENES_v",
    "CREATE VIEW REACTIONS_GENES_v AS
       SELECT * FROM REACTIONS_GENES_update
       UNION ALL
       SELECT * FROM REACTIONS_GENES
       WHERE (SELECT COUNT(*) FROM REACTIONS_GENES_update) = 0"
  )
  exec_stmts(con, shadow_stmts)

  invisible(NULL)
}

preview_table_function <- function(input_db,input_restore,input_upload,required_cols,query){
  if (input_db == "default" || input_restore == TRUE) {
    db_md_table <- dbGetQuery(con, query)
    return(db_md_table)
  } else {
    
    req(input_upload)
    
    ext <- tools::file_ext(input_upload$name)
    
    tmp <- switch(
      ext,
      "csv" = read.csv(
        input_upload$datapath,
        header = TRUE,
        stringsAsFactors = FALSE
      ),
      "tsv" = read.delim(
        input_upload$datapath,
        header = TRUE,
        stringsAsFactors = FALSE
      ),
      "xlsx" = read_excel(input_upload$datapath,col_names = TRUE)
    )
    if (is.null(tmp)) shiny::validate(shiny::need(FALSE, paste0("Unsupported format: ", ext, ". Use csv, tsv, or xlsx.")))
    
    missing_cols <- setdiff(required_cols, colnames(tmp))
    validate(
      need(
        length(missing_cols) == 0,
        paste0(
          "Invalid custom database. ",
          "Missing columns: ", paste(missing_cols, collapse = ", "), ". ",
          "Required columns are: ", paste(required_cols, collapse = ", "), "."
        )
      )
    )
    
    
    tmp

  }
}

restore_table <- function(con, table_name) {
  origin_name <- paste0(table_name, "_origin")
  tables <- dbListTables(con)

  if (!(origin_name %in% tables)) {
    showNotification(
      paste("No backup found for", table_name, "– nothing to restore."),
      type = "warning",
      duration = 5
    )
    return(invisible(NULL))
  }

  dbWithTransaction(con, {
    dbExecute(con, paste0("DROP TABLE IF EXISTS ", DBI::dbQuoteIdentifier(con, table_name), ";"))
    dbExecute(con, paste0("CREATE TABLE ", DBI::dbQuoteIdentifier(con, table_name),
                          " AS SELECT * FROM ", DBI::dbQuoteIdentifier(con, origin_name), ";"))
  })
  # Il backup _origin viene mantenuto per permettere ripristini multipli
  # Re-create indices lost when table was dropped
  ensure_db_indices(con)

  showNotification(
    paste(table_name, "restored successfully."),
    type = "message",
    duration = 4
  )
}
