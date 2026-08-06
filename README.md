This repository contains the data and scripts related to "ShinyMIR to investigate miRNA-target associations".

# ShinyMIR

An interactive R Shiny application for exploring miRNA regulatory networks across disease, gene, metabolic reaction, and pathway dimensions.

---

## Overview

ShinyMIR integrates curated miRNA–disease and miRNA–gene associations with metabolic network data and Reactome pathway annotations. It provides interactive tables, network visualisations, co-targeting analyses, robustness diagnostics, and an exportable HTML report — all from a single DuckDB database.

---

## Analysis modules

### 1. miRNA – Disease
Search curated miRNA–disease associations. Filter by miRNA, disease, Disease Ontology (DO) macro-category, minimum publications, sample type, and human-only studies. Results include an interactive table (groupable by miRNA, disease, or both), a force-directed network, co-targeting heatmaps (miRNA×miRNA and disease×disease), and network robustness/validation diagnostics.

### 2. miRNA – Gene
Explore validated miRNA–gene regulatory relationships. miRNAs can be entered manually or propagated from the miRNA–Disease search (with arm selection: 5p, 3p, both, or most associated). Evidence filters cover type, strength, throughput, and match mode (any / all). An expression filter based on Human Protein Atlas nTPM values restricts genes to those expressed above a threshold in a chosen tissue or cell line. The network view supports node colouring by Reactome pathway or by a custom GMT gene-set file. Co-targeting section scores miRNA pairs (Jaccard over shared genes) or gene pairs (Jaccard over shared miRNAs). Network Analysis tab runs threshold-sweep robustness and co-targeting robustness diagnostics.

### 3. Pathway Annotation
Annotates the miRNA–gene result set against Reactome or a custom GMT file (MSigDB, Hallmarks, KEGG, …). Outputs coverage per pathway (fraction of pathway genes targeted), number of implicated miRNAs, an overview bar chart (top N, optionally grouped by L1 macro-category), and a miRNA×pathway heatmap. When an expression filter is active in the miRNA–Gene tab, the same tissue/cell-line context is inherited and used to restrict the pathway universe.

### 4. miRNA – Reaction
Links miRNAs to metabolic reactions via shared target genes. miRNAs and genes are propagated from the miRNA–Gene search. Results can be filtered by metabolic category (L1) and subsystem (L2) and grouped by miRNA, reaction, subsystem, or category. Includes a coverage bar chart, a **Subsystem Analysis** tab (metabolic heatmap, pathway enrichment, gene ranking), and co-targeting analysis (miRNA×miRNA and subsystem×subsystem).

---

## Report

Clicking **Generate report** captures all currently rendered plots (network, heatmaps, robustness curves, barplots) as PNG images and builds a self-contained HTML report covering all four modules. The report is bundled with the images into a `.zip` file and downloaded automatically.

> Note: plots are captured at the moment of rendering. Visit each subtab (Network, Co-targeting, Network Analysis) before generating the report to ensure all images are included.

---

## Database

The application reads from a single DuckDB database file (`shinyMIR.duckdb`) placed in the project root. 

### Key tables

| Table | Content |
|---|---|
| `DISEASES` | Disease catalogue |
| `DISEASE_LINK` | DO ancestor–descendant hierarchy |
| `MIRNAS` | miRNA catalogue (premature + arm-specific names) |
| `GENES` | Gene catalogue |
| `ARTICLES` | PubMed article metadata (title, abstract, year) |
| `MIRNAS_DISEASES_ARTICLES` | miRNA–disease associations with PubMed IDs |
| `MIRNAS_GENES_ARTICLES` | miRNA–gene associations with evidence |
| `PROOF_MIRNA_GENE` | Evidence metadata (type, strength, throughput) |
| `REACTIONS` / `REACTIONS_GENES` | Metabolic reactions and gene–reaction links |
| `ONTOLOGY_PATHWAYS` / `ONTOLOGY_PATHWAYS_GENES` | Reactome pathway hierarchy and gene membership |
| `HPA_TISSUE_DATA` / `HPA_CELLLINE_DATA` | Human Protein Atlas expression data (nTPM) |

## Requirements

R ≥ 4.2. Open the project in RStudio (`ShinyMir.Rproj`) to ensure the working directory is set correctly.

### Option A — renv (recommended for reproducibility)

Locks exact package versions so the environment is identical across machines and time.

```r
install.packages("renv")
renv::restore()
```

### Option B — manual install

Installs the latest available versions (simpler, but may break with future package updates).

```r
source("install_packages.R")
```

---

## Running locally

```r
shiny::runApp(".")
```

---

## Project structure

```
ShinyMir/
├── global.R                    # Package loading, DB connection, global variables
├── ui.R                        # Dashboard UI definition
├── server.R                    # Server entry point, sources all server modules
├── server/
│   ├── server_mda.R            # miRNA–Disease logic
│   ├── server_mga.R            # miRNA–Gene logic
│   ├── server_mr.R             # miRNA–Reaction logic
│   ├── server_overrep.R        # Pathway Annotation logic
│   ├── server_upload.R         # Data upload handling
│   └── server_toggles.R        # UI toggle/disable helpers
├── R/
│   ├── constant.R              # SQL query templates
│   ├── filters.R               # DB query functions (get_filtered_*)
│   ├── graphs.R                # Network rendering (visNetwork + igraph)
│   ├── cotarget.R              # Jaccard matrix, co-targeting heatmaps
│   ├── metabolic_heatmap.R     # miRNA×subsystem coverage matrix
│   ├── robustness.R            # Threshold sweep, co-targeting robustness
│   ├── enrichment.R            # Pathway over-representation (hypergeometric)
│   ├── report.R                # HTML report builder
│   ├── plot_module.R           # Reusable plotVisualizationServer module
│   ├── db.R                    # DB utilities, index creation, upload handlers
│   ├── ui_helpers.R            # Shared UI components
│   └── utils.R                 # Shared helpers, upload_list_server module
```
