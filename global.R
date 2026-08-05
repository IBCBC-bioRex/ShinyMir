library(shiny)
library(DBI)
library(duckdb)
library(dplyr)
library(DT)
library(tidyr)
library(htmlwidgets)
library(kableExtra)
library(glue)
library(visNetwork)
library(plotly)
library(heatmaply)
library(RColorBrewer)
library(igraph)
library(shinyjs)
library(shinycssloaders)
library(shinydashboard)
library(readxl)

# Sources
source("R/constant.R")
source("R/utils.R")
source("R/db.R")
source("R/filters.R")
source("R/graphs.R")
source("R/ui_helpers.R")
source("R/plot_module.R")
source("R/cotarget.R")
source("R/metabolic_heatmap.R")
source("R/robustness.R")
source("R/enrichment.R")   
source("R/report.R")

#######connect to db
file_db = "shinyMIR.duckdb"
if (!file.exists(file_db))
  stop("Database file not found: '", file_db, "'. Run creation_database/sqlite_to_duckdb.R first.")
con <- dbConnect(duckdb::duckdb(), file_db)

ensure_db_indices(con)
# Resolve effective table names: use views if created, else fall back to base tables
use_views           <- reactions_view_exists(con)
REACTIONS_TBL       <- if (use_views) "REACTIONS_v"       else "REACTIONS"
REACTIONS_GENES_TBL <- if (use_views) "REACTIONS_GENES_v" else "REACTIONS_GENES"
# Resolve query_mr string templates to the actual table names chosen above
query_mr <- lapply(query_mr, function(q) {
  q <- gsub("REACTIONS_GENES_v", REACTIONS_GENES_TBL, q, fixed = TRUE)
  q <- gsub("REACTIONS_v",       REACTIONS_TBL,       q, fixed = TRUE)
  q
})
HAS_IS_HUMAN <- tryCatch(
  "ARTICLES" %in% dbListTables(con) &&
    "IS_HUMAN" %in% dbListFields(con, "ARTICLES"),
  error = function(e) FALSE
)

N_TOTAL_DISEASES <- tryCatch(
  dbGetQuery(con, "SELECT COUNT(DISTINCT DISEASE_ID) FROM MIRNAS_DISEASES_ARTICLES")[[1]],
  error = function(e) NA_integer_
)
N_TOTAL_GENES <- tryCatch(
  dbGetQuery(con, "SELECT COUNT(DISTINCT GENE_ID) FROM MIRNAS_GENES_ARTICLES")[[1]],
  error = function(e) NA_integer_
)
MIRNA_GLOBAL_DIS_COUNT <- tryCatch(
  dbGetQuery(con, "
    SELECT M.MIRNA_PREMATURE AS MIRNA_NAME,
           COUNT(DISTINCT MDA.DISEASE_ID) AS N_DIS_GLOBAL
    FROM MIRNAS M
    JOIN MIRNAS_DISEASES_ARTICLES MDA ON M.MIRNA_ID = MDA.MIRNA_ID
    GROUP BY M.MIRNA_PREMATURE
  "),
  error = function(e) NULL
)
################## expression filter choices
tissue_choices <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT TISSUE_NAME FROM TISSUE_TYPE ORDER BY TISSUE_NAME")$TISSUE_NAME,
  error = function(e) character(0))
cell_line_primary_diseases <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT PRIMARY_DISEASE FROM CELL_LINE WHERE PRIMARY_DISEASE IS NOT NULL ORDER BY PRIMARY_DISEASE")$PRIMARY_DISEASE,
  error = function(e) character(0))
subsystem_choices <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT SUBSYSTEM FROM REACTIONS WHERE SUBSYSTEM IS NOT NULL AND SUBSYSTEM != '' ORDER BY SUBSYSTEM")$SUBSYSTEM,
  error = function(e) character(0))
category_choices            <- sort(unique(na.omit(SUBSYSTEM_CATEGORY_MAP)))
reaction_name_choices <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT NAME FROM REACTIONS ORDER BY NAME")$NAME,
  error = function(e) character(0))
##################
maximum_nodo_graph=500
################################################################################
# Scelte per i menu a tendina
mirna_choices <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT MIRNA_PREMATURE FROM MIRNAS ORDER BY MIRNA_PREMATURE")$MIRNA_PREMATURE,
  error = function(e) character(0))
mirna_choices_3p5p <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT MIRNA_NAME FROM MIRNAS WHERE MIRNA_NAME IS NOT NULL ORDER BY MIRNA_NAME")$MIRNA_NAME,
  error = function(e) character(0))
disease_choices <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT DISEASE FROM DISEASES ORDER BY DISEASE")$DISEASE,
  error = function(e) character(0))

# DO macro-categories: direct children of the absolute root (no ANCESTOR = root nodes,
# or first level below DOID:4). We find root by: nodes that appear as ANCESTOR but
# never as DISEASE_ID in DISEASE_LINK, then take their direct children.
do_root_categories <- tryCatch({
  res <- dbGetQuery(con, "
    WITH root AS (
      SELECT ANCESTOR AS id
      FROM DISEASE_LINK
      WHERE ANCESTOR NOT IN (SELECT DISEASE_ID FROM DISEASE_LINK)
    ),
    level1 AS (
      SELECT DL.DISEASE_ID
      FROM DISEASE_LINK DL
      JOIN root ON DL.ANCESTOR = root.id
    )
    SELECT DISTINCT D.DISEASE
    FROM DISEASES D
    JOIN level1 L ON D.DISEASE_ID = L.DISEASE_ID
    ORDER BY D.DISEASE
  ")
  res$DISEASE
}, error = function(e) character(0))

# Precompute: for each macro-category, which diseases belong to it (recursive)
do_category_disease_map <- tryCatch({
  if (length(do_root_categories) == 0) {
    list()
  } else {
    res <- dbGetQuery(con, "
      WITH RECURSIVE subtree(id, root_disease) AS (
        SELECT D.DISEASE_ID, D.DISEASE
        FROM DISEASES D
        WHERE D.DISEASE IN (SELECT ANCESTOR_D.DISEASE
          FROM DISEASE_LINK DL
          JOIN DISEASES ANCESTOR_D ON DL.DISEASE_ID = ANCESTOR_D.DISEASE_ID
          WHERE DL.ANCESTOR NOT IN (SELECT DISEASE_ID FROM DISEASE_LINK))
        UNION ALL
        SELECT DL.DISEASE_ID, st.root_disease
        FROM DISEASE_LINK DL
        JOIN subtree st ON DL.ANCESTOR = st.id
      )
      SELECT D.DISEASE, st.root_disease AS CATEGORY
      FROM subtree st
      JOIN DISEASES D ON st.id = D.DISEASE_ID
    ")
    if (nrow(res) == 0) list() else split(res$DISEASE, res$CATEGORY)
  }
}, error = function(e) list())
# Per-category DSI counts (pre-calculated at boot for each DO macro-category)
MIRNA_CATEGORY_DIS_COUNT <- tryCatch({
  if (length(do_category_disease_map) == 0) list() else
  lapply(do_category_disease_map, function(diseases) {
    if (length(diseases) == 0) return(setNames(integer(0), character(0)))
    res <- dbGetQuery(con, paste0(
      "SELECT M.MIRNA_PREMATURE AS MIRNA_NAME, COUNT(DISTINCT MDA.DISEASE_ID) AS N_DIS",
      " FROM MIRNAS M",
      " JOIN MIRNAS_DISEASES_ARTICLES MDA ON M.MIRNA_ID = MDA.MIRNA_ID",
      " JOIN DISEASES D ON MDA.DISEASE_ID = D.DISEASE_ID",
      " WHERE D.DISEASE IN (", paste(sprintf("'%s'", gsub("'", "''", diseases)), collapse = ","), ")",
      " GROUP BY M.MIRNA_PREMATURE"
    ))
    setNames(res$N_DIS, res$MIRNA_NAME)
  })
}, error = function(e) list())

N_DISEASES_BY_CATEGORY <- tryCatch(
  sapply(do_category_disease_map, length),
  error = function(e) integer(0)
)

# Full arm-specific → premature lookup (static, replaces per-request DB queries)
MIRNA_ARM_TO_PREMATURE <- tryCatch({
  m <- dbGetQuery(con, "SELECT MIRNA_NAME, MIRNA_PREMATURE FROM MIRNAS WHERE MIRNA_NAME IS NOT NULL")
  setNames(m$MIRNA_PREMATURE, m$MIRNA_NAME)
}, error = function(e) character(0))

# Pre-computed DSI reference means (one per category + global) — constant across session
MIRNA_GLOBAL_REF_MEAN <- tryCatch({
  if (is.null(MIRNA_GLOBAL_DIS_COUNT) || is.na(N_TOTAL_DISEASES) || N_TOTAL_DISEASES <= 1) return(NA_real_)
  round(mean(1 - log2(pmax(MIRNA_GLOBAL_DIS_COUNT$N_DIS_GLOBAL, 1)) / log2(N_TOTAL_DISEASES)), 3)
}, error = function(e) NA_real_)

MIRNA_CATEGORY_REF_MEAN <- tryCatch({
  if (length(MIRNA_CATEGORY_DIS_COUNT) == 0) return(numeric(0))
  sapply(names(MIRNA_CATEGORY_DIS_COUNT), function(cat) {
    v       <- MIRNA_CATEGORY_DIS_COUNT[[cat]]
    n_total <- max(N_DISEASES_BY_CATEGORY[[cat]], 2L, na.rm = TRUE)
    if (length(v) == 0 || is.na(n_total) || n_total <= 1) return(NA_real_)
    round(mean(1 - log2(pmax(v, 1)) / log2(n_total)), 3)
  })
}, error = function(e) numeric(0))

# Boolean: PATHWAY_LINK table exists and has rows (Reactome hierarchy available)
ONTOLOGY_HAS_HIERARCHY <- tryCatch(
  "PATHWAY_LINK" %in% dbListTables(con) &&
    dbGetQuery(con, "SELECT COUNT(*) FROM PATHWAY_LINK")[[1]] > 0,
  error = function(e) FALSE
)

# Helper: return list(gdc, n_total) for DSI computation given current category selection.
# cat_sel: value of input$do_category_filter (empty string = no category selected).
resolve_dsi_context <- function(cat_sel) {
  if (nzchar(cat_sel) &&
      cat_sel %in% names(MIRNA_CATEGORY_DIS_COUNT) &&
      length(MIRNA_CATEGORY_DIS_COUNT[[cat_sel]]) > 0) {
    list(
      gdc     = MIRNA_CATEGORY_DIS_COUNT[[cat_sel]],
      n_total = max(N_DISEASES_BY_CATEGORY[[cat_sel]], 2L, na.rm = TRUE),
      cat     = cat_sel
    )
  } else if (!is.null(MIRNA_GLOBAL_DIS_COUNT) && nrow(MIRNA_GLOBAL_DIS_COUNT) > 0) {
    list(
      gdc     = setNames(MIRNA_GLOBAL_DIS_COUNT$N_DIS_GLOBAL, MIRNA_GLOBAL_DIS_COUNT$MIRNA_NAME),
      n_total = N_TOTAL_DISEASES,
      cat     = ""
    )
  } else {
    list(gdc = NULL, n_total = N_TOTAL_DISEASES, cat = "")
  }
}

metabolic_genes <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT GENE_NAME FROM GENES A JOIN REACTIONS_GENES B ON A.GENE_ID=B.GENE_ID")$GENE_NAME,
  error = function(e) character(0))
total_genes <- tryCatch(
  dbGetQuery(con, "SELECT DISTINCT GENE_NAME FROM GENES")$GENE_NAME,
  error = function(e) character(0))
pathway_choices <- tryCatch(
  dbGetQuery(con, "SELECT PATHWAY_NAME FROM ONTOLOGY_PATHWAYS ORDER BY PATHWAY_NAME")$PATHWAY_NAME,
  error = function(e) character(0))

# ── Ontology enrichment globals (used by R/enrichment.R) ────────────────────
df_all_ontology <- tryCatch(
  dbGetQuery(con, "
    SELECT DISTINCT C.GENE_NAME, B.PATHWAY_NAME, B.PATHWAY_ID
    FROM ONTOLOGY_PATHWAYS_GENES A
    JOIN ONTOLOGY_PATHWAYS B ON A.PATHWAY_ID = B.PATHWAY_ID
    JOIN GENES C ON A.GENE_ID = C.GENE_ID
  "),
  error = function(e) data.frame(GENE_NAME = character(0), PATHWAY_NAME = character(0), PATHWAY_ID = integer(0))
)
N_tot <- tryCatch(
  dbGetQuery(con, "SELECT COUNT(DISTINCT GENE_ID) FROM GENES")[[1]],
  error = function(e) 0L
)
if (nrow(df_all_ontology) > 0) {
  pathway_sizes      <- tapply(df_all_ontology$GENE_NAME, df_all_ontology$PATHWAY_ID, function(x) length(unique(x)))
  pathway_genes_list <- lapply(split(df_all_ontology$GENE_NAME, df_all_ontology$PATHWAY_ID), unique)
  gene_pathway_index <- build_gene_pathway_index(pathway_genes_list)
  pathway_names_map  <- setNames(df_all_ontology$PATHWAY_NAME[!duplicated(df_all_ontology$PATHWAY_ID)],
                                 df_all_ontology$PATHWAY_ID[!duplicated(df_all_ontology$PATHWAY_ID)])
} else {
  pathway_sizes      <- integer(0)
  pathway_genes_list <- list()
  gene_pathway_index <- list()
  pathway_names_map  <- character(0)
}

pathway_l1_df <- tryCatch({
  dbGetQuery(con, "
    SELECT P.PATHWAY_ID, P.PATHWAY_NAME
    FROM ONTOLOGY_PATHWAYS P
    WHERE P.PATHWAY_ID NOT IN (SELECT PATHWAY_ID FROM PATHWAY_LINK)
    ORDER BY P.PATHWAY_NAME
  ")
}, error = function(e) data.frame(PATHWAY_ID = integer(0), PATHWAY_NAME = character(0)))
pathway_l1_choices <- setNames(pathway_l1_df$PATHWAY_ID, pathway_l1_df$PATHWAY_NAME)

proof_table <- tryCatch(
  dbGetQuery(con, "SELECT * FROM PROOF_MIRNA_GENE ORDER BY PROOF_NAME"),
  error = function(e) data.frame(PROOF_ID = integer(0), PROOF_NAME = character(0),
                                  TYPE = character(0), STRENGTH = integer(0), THROUGHPUT = character(0)))
proofs_all           <- sort(unique(proof_table$PROOF_NAME))
proof_choices_strength    <- sort(unique(proof_table$STRENGTH))
proof_choices_type        <- sort(unique(proof_table$TYPE))
proof_choices_Throughput  <- sort(unique(proof_table$THROUGHPUT))

#######################################################
group_mdp_choices<-c("None","miRNA","disease","miRNA and disease")
group_mga_choices<-c("None","miRNA","gene","miRNA and gene")
group_mr_choices <- c("None","miRNA","reaction","miRNA and reaction","subsystem","category")


# ── Cached TOTAL_GENES per reaction – used by compute_essentiality_fraction to avoid
# repeated DB queries. Refreshed on metabolism restore/upload.
reaction_total_genes_map <- tryCatch(local({
  tg <- dbGetQuery(con,
    "SELECT R.NAME AS NAME, COUNT(DISTINCT RG.GENE_ID) AS TOTAL_GENES
       FROM REACTIONS R
       JOIN REACTIONS_GENES RG ON R.REACTION_ID = RG.REACTION_ID
      GROUP BY R.REACTION_ID, R.NAME")
  setNames(as.integer(tg$TOTAL_GENES), tg$NAME)
}), error = function(e) setNames(integer(0), character(0)))

# Close DB connection when app process exits (app-level, not per-session)
shiny::onStop(function() tryCatch(dbDisconnect(con), error = function(e) NULL))


