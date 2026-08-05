# =============================================================================
# ShinyMIR — HTML report builder
# =============================================================================

.report_css <- '
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:system-ui,-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;
     font-size:14px;line-height:1.6;color:#222;background:#f0f2f5;
     padding:32px 16px 60px}
.report{max-width:940px;margin:0 auto;background:#fff;
        border:1px solid #dde3ec;border-radius:4px;overflow:hidden}
.report-header{background:#1c3557;color:#fff;padding:28px 40px;
               display:flex;justify-content:space-between;align-items:flex-end;gap:24px}
.report-header h1{font-size:20px;font-weight:700;color:#fff}
.report-header .sub{font-size:12px;color:#a8c4de;margin-top:3px}
.report-header .meta{text-align:right;font-size:12px;color:#a8c4de;
                     line-height:1.9;white-space:nowrap}
.report-header .meta strong{color:#d4e8f7}
.toc{padding:16px 40px;background:#f7f9fc;border-bottom:1px solid #dde3ec}
.toc h2{font-size:10px;font-weight:700;text-transform:uppercase;
         letter-spacing:.8px;color:#6b7a8d;margin-bottom:8px}
.toc ul{list-style:none;display:flex;flex-wrap:wrap;gap:4px 18px}
.toc li{font-size:13px}
.toc a{color:#2878b5;text-decoration:none}
.section{padding:32px 40px;border-bottom:1px solid #e8ecf2}
.section:last-of-type{border-bottom:none}
.sec-head{display:flex;align-items:center;gap:12px;margin-bottom:20px;
          padding-bottom:10px;border-bottom:2px solid #e8ecf2}
.sec-badge{display:inline-flex;align-items:center;justify-content:center;
           width:30px;height:30px;border-radius:5px;color:#fff;
           font-size:12px;font-weight:700;flex-shrink:0}
.sec-head h2{font-size:16px;font-weight:700;color:#1c3557}
.sec-ts{margin-left:auto;font-size:11px;color:#6b7a8d;
        background:#f0f4fb;padding:2px 9px;border-radius:10px;white-space:nowrap}
.filter-box{background:#f7f9fc;border:1px solid #dde3ec;border-left:3px solid #2878b5;
            border-radius:4px;padding:12px 16px;margin-bottom:16px}
.filter-box h3{font-size:10px;font-weight:700;text-transform:uppercase;
               letter-spacing:.7px;color:#6b7a8d;margin-bottom:8px}
.frow{display:flex;align-items:baseline;gap:8px;margin-bottom:5px;font-size:12.5px}
.frow:last-child{margin-bottom:0}
.flabel{color:#6b7a8d;width:170px;flex-shrink:0;font-size:11.5px}
.fval{display:flex;flex-wrap:wrap;gap:3px}
.badge{display:inline-block;padding:2px 7px;border-radius:9px;
       font-size:11px;font-weight:500}
.b-blue{background:#e0eefa;color:#1a5a9a}
.b-green{background:#e2f5ea;color:#1a7a3a}
.b-grey{background:#eaedf2;color:#4a5568}
.b-orange{background:#fef3e2;color:#a05a00}
.count-bar{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:16px}
.chip{display:flex;align-items:center;gap:6px;padding:5px 13px;
      border-radius:18px;border:1px solid #dde3ec;background:#f7f9fc;font-size:13px}
.chip .n{font-size:16px;font-weight:700;color:#1c3557;line-height:1}
.chip .l{color:#6b7a8d;font-size:11.5px}
.sublabel{font-size:10px;font-weight:700;text-transform:uppercase;
          letter-spacing:.7px;color:#8a9bb5;margin:16px 0 6px}
.twrap{overflow-x:auto;border:1px solid #dde3ec;border-radius:4px;margin-bottom:6px}
table{width:100%;border-collapse:collapse;font-size:12px}
thead th{background:#1c3557;color:#fff;text-align:left;padding:8px 11px;
         font-weight:600;white-space:nowrap}
tbody tr:nth-child(even){background:#f4f7fb}
tbody td{padding:6px 11px;border-bottom:1px solid #eaecf0;
         max-width:280px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.tnote{font-size:11px;color:#6b7a8d;font-style:italic;margin-bottom:12px}
.no-data{font-size:13px;color:#8a9bb5;font-style:italic;
         padding:18px;background:#f7f9fc;border-radius:4px;
         border:1px dashed #c8d0de;text-align:center}
.fig-wrap{margin-bottom:8px}
.fig-wrap img{max-width:100%;border:1px solid #dde3ec;border-radius:4px;display:block}
.fig-cap{font-size:11.5px;color:#6b7a8d;font-style:italic;margin-top:4px}
.report-footer{background:#f7f9fc;border-top:1px solid #dde3ec;
               padding:14px 40px;font-size:11px;color:#8a9bb5;
               display:flex;justify-content:space-between}
.report-footer strong{color:#4a5568}
.net-stats{display:flex;flex-wrap:wrap;gap:16px;padding:8px 14px;
           background:#f0f4fb;border:1px solid #dde3ec;border-radius:4px;
           font-size:12px;color:#4a5568;margin-top:6px;margin-bottom:14px}
.net-stats strong{color:#1c3557}
.db-info-box{display:flex;flex-wrap:wrap;gap:8px 24px;padding:10px 16px;
             background:#f0f4fb;border:1px solid #dde3ec;border-radius:4px;
             font-size:12px;color:#4a5568;margin-bottom:4px}
.db-info-box strong{color:#1c3557}
.plot-desc{font-size:12.5px;color:#4a5568;line-height:1.6;
           background:#f7f9fc;border-left:3px solid #c8d0de;
           border-radius:0 4px 4px 0;padding:8px 12px;margin:10px 0 6px}
.subsec{margin:20px 0 4px;padding-top:14px;border-top:1px dashed #dde3ec}
.subsec-title{font-size:11px;font-weight:700;text-transform:uppercase;
              letter-spacing:.8px;color:#2878b5;margin-bottom:10px}
.variant-label{font-size:10.5px;font-weight:600;color:#5a6a7d;
               background:#eef2fb;padding:2px 8px;border-radius:8px;
               display:inline-block;margin:8px 0 4px}
</style>
'

# ── Internal helpers ─────────────────────────────────────────────────────────

.strip_html <- function(x) gsub("<[^>]+>", "", as.character(x))

.badge_html <- function(vals, cls = "b-grey") {
  vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
  if (length(vals) == 0)
    return('<span class="badge b-grey">All / none</span>')
  paste(sprintf('<span class="badge %s">%s</span>',
                cls, htmltools::htmlEscape(as.character(vals))),
        collapse = " ")
}

.filter_row <- function(label, vals, cls = "b-grey") {
  sprintf('<div class="frow"><span class="flabel">%s</span><div class="fval">%s</div></div>',
          label, .badge_html(vals, cls))
}

.chip <- function(n, label) {
  sprintf('<div class="chip"><span class="n">%s</span><span class="l">%s</span></div>',
          format(n, big.mark = ","), label)
}


.table_html <- function(df, max_rows = 10) {
  if (is.null(df) || nrow(df) == 0)
    return('<p class="no-data">No data.</p>')
  df_show <- head(df, max_rows)
  df_show[] <- lapply(df_show, function(col) .strip_html(col))
  ths  <- paste(sprintf("<th>%s</th>", names(df_show)), collapse = "")
  rows <- apply(df_show, 1, function(r) {
    tds <- paste(sprintf('<td title="%s">%s</td>',
                         htmltools::htmlEscape(as.character(r)),
                         htmltools::htmlEscape(substr(as.character(r), 1, 55))),
                 collapse = "")
    sprintf("<tr>%s</tr>", tds)
  })
  note <- if (nrow(df) > max_rows)
    sprintf('<p class="tnote">Showing %d of %s rows</p>',
            max_rows, format(nrow(df), big.mark = ","))
  else
    sprintf('<p class="tnote">%s rows total</p>', format(nrow(df), big.mark = ","))
  paste0('<div class="twrap"><table><thead><tr>', ths,
         '</tr></thead><tbody>', paste(rows, collapse = ""),
         '</tbody></table></div>', note)
}

.section_wrap <- function(id, n, title, color, ts, filters, counts, body) {
  sprintf('
<div class="section" id="%s">
  <div class="sec-head">
    <div class="sec-badge" style="background:%s;">%d</div>
    <h2>%s</h2>
    <span class="sec-ts">%s</span>
  </div>
  %s
  <div class="count-bar">%s</div>
  %s
</div>', id, color, n, title, ts, filters, counts, body)
}

# ── Static network via igraph → base64 PNG ────────────────────────────────────

# Embed a relative PNG filename as <img> (for ZIP export).
.captured_img_html <- function(filename, caption) {
  if (is.null(filename) || !nzchar(filename %||% "")) return("")
  sprintf('<div class="fig-wrap">
    <p class="sublabel" style="margin-top:0;">%s</p>
    <img src="%s" alt="%s" style="max-width:100%%;border:1px solid #dde3ec;border-radius:6px;display:block;">
  </div>', caption, htmltools::htmlEscape(filename), htmltools::htmlEscape(caption))
}

.plot_block <- function(img_files, key, desc, caption = "") {
  if (is.null(img_files[[key]])) return("")
  paste0('<p class="plot-desc">', desc, '</p>', .captured_img_html(img_files[[key]], caption))
}

.subsec <- function(title, ...) {
  body <- paste(..., sep = "\n")
  if (!nzchar(trimws(body))) return("")
  sprintf('<div class="subsec"><div class="subsec-title">%s</div>%s</div>', title, body)
}

.variant_label <- function(label) {
  sprintf('<span class="variant-label">%s</span>', htmltools::htmlEscape(label))
}

# Build paired variant block: two optional images under a common subsection
.variant_pair_block <- function(img_files, key_a, key_b, label_a, label_b, desc_a, desc_b) {
  a <- if (!is.null(img_files[[key_a]])) paste0(.variant_label(label_a), .plot_block(img_files, key_a, desc_a, label_a)) else ""
  b <- if (!is.null(img_files[[key_b]])) paste0(.variant_label(label_b), .plot_block(img_files, key_b, desc_b, label_b)) else ""
  paste0(a, b)
}

# Derive co-targeting pair table from a flat df (MIRNA_NAME + second col)
.cotarget_pairs_html <- function(flat_df, col_nodes, col_shared, label) {
  if (is.null(flat_df) || !all(c(col_nodes, col_shared) %in% names(flat_df))) return("")
  mat <- tryCatch(
    compute_jaccard_matrix(flat_df, col_nodes = col_nodes, col_shared = col_shared),
    error = function(e) NULL
  )
  if (is.null(mat)) return("")
  tbl <- tryCatch(
    cotarget_pairs_table(mat, threshold = 0,
                         df_source = flat_df, col_nodes = col_nodes, col_shared = col_shared),
    error = function(e) NULL
  )
  if (is.null(tbl) || nrow(tbl) == 0) return("")
  paste0('<div class="sublabel">', label, '</div>', .table_html(tbl, max_rows = 10))
}

.section_empty <- function(id, n, title, color) {
  sprintf('<div class="section" id="%s"><div class="sec-head"><div class="sec-badge" style="background:%s;">%d</div><h2>%s</h2></div><p class="no-data">No search performed in this session.</p></div>',
          id, color, n, title)
}

# ── Helpers ──────────────────────────────────────────────────────────────────

# Resolve pathway L1 IDs → names using global pathway_l1_choices (named vector: name→id)
.l1_names <- function(ids) {
  ids <- suppressWarnings(as.integer(ids))
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return(character(0))
  # pathway_l1_choices: setNames(id_vec, name_vec) — values are IDs
  nm <- names(pathway_l1_choices)[pathway_l1_choices %in% ids]
  if (length(nm) == 0) as.character(ids) else nm
}

# Resolve pathway L2 IDs → names from pathway_names_map (names_map: id→name)
.l2_names <- function(ids) {
  ids <- suppressWarnings(as.integer(ids))
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return(character(0))
  nm <- pathway_names_map[as.character(ids)]
  nm[is.na(nm)] <- as.character(ids[is.na(nm)])
  unname(nm)
}

# ── Section: MDA ─────────────────────────────────────────────────────────────

.section_mda <- function(df, input, img_files = list(), ts, flat_df = NULL) {
  if (is.null(df)) return(.section_empty("mda", 1, "miRNA - Disease", "#2878b5"))

  use_all_m   <- isTRUE(input$filter_mirna_all_checkbox)
  use_all_d   <- isTRUE(input$filter_disease_checkbox)
  mirnas      <- input$mirna_select       %||% character(0)
  diseases    <- input$disease_select     %||% character(0)
  do_cat      <- input$do_category_filter %||% ""
  ancestors   <- isTRUE(input$filter_diseaseancestor_checkbox)
  min_pub     <- input$min_assoc_mda      %||% 1L
  sample_type <- input$sample_type_mda   %||% character(0)
  human_only  <- isTRUE(input$filter_human_mda)
  grouping    <- input$group_mda          %||% "None"

  filters <- paste0(
    '<div class="filter-box" style="border-left-color:#2878b5;"><h3>Filters applied</h3>',
    .filter_row("miRNAs",
      if (use_all_m) "All miRNAs" else if (length(mirnas) == 0) "None selected" else mirnas,
      cls = "b-blue"),
    .filter_row("Diseases",
      if (use_all_d) "All diseases" else if (length(diseases) == 0) "None selected" else diseases,
      cls = "b-blue"),
    if (nzchar(do_cat)) .filter_row("Disease category (DO)", do_cat, cls = "b-orange") else "",
    if (ancestors) .filter_row("Include disease descendants", "Yes", cls = "b-grey") else "",
    .filter_row("Min. publications",  as.character(min_pub), cls = "b-grey"),
    if (length(sample_type) > 0) .filter_row("Sample type", sample_type, cls = "b-orange") else "",
    .filter_row("Human studies only", if (human_only) "Yes" else "No", cls = "b-grey"),
    .filter_row("Table grouping",     grouping, cls = "b-grey"),
    '</div>')

  n_mirna  <- if ("MIRNA_NAME" %in% names(df)) dplyr::n_distinct(.strip_html(df$MIRNA_NAME)) else "—"
  n_disease <- if ("DISEASE"   %in% names(df)) dplyr::n_distinct(.strip_html(df$DISEASE))    else "—"
  n_pairs  <- if (!is.null(flat_df)) nrow(flat_df) else "—"
  counts   <- paste(.chip(n_mirna, "miRNAs"), .chip(n_disease, "diseases"), .chip(n_pairs, "pairs"))

  # ── Table-derived summaries ──
  flat <- flat_df
  mirna_col   <- if (!is.null(flat) && "MIRNA_NAME" %in% names(flat)) "MIRNA_NAME" else NULL
  disease_col <- if (!is.null(flat) && "DISEASE"    %in% names(flat)) "DISEASE"    else NULL

  grp_by_mirna <- if (!is.null(flat) && !is.null(mirna_col) && !is.null(disease_col))
    tryCatch(flat %>%
      dplyr::group_by(MIRNA = .data[[mirna_col]]) %>%
      dplyr::summarise(
        N_diseases   = dplyr::n_distinct(.data[[disease_col]]),
        Diseases     = paste(sort(unique(.data[[disease_col]])), collapse = ", "),
        .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_diseases)),
      error = function(e) NULL)
  else NULL

  grp_by_disease <- if (!is.null(flat) && !is.null(mirna_col) && !is.null(disease_col))
    tryCatch(flat %>%
      dplyr::group_by(Disease = .data[[disease_col]]) %>%
      dplyr::summarise(
        N_miRNAs = dplyr::n_distinct(.data[[mirna_col]]),
        miRNAs   = paste(sort(unique(.data[[mirna_col]])), collapse = ", "),
        .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_miRNAs)),
      error = function(e) NULL)
  else NULL

  # ── Subsections ──
  tab_results <- .subsec("Results",
    '<div class="sublabel">Main result table</div>',
    .table_html(df),
    if (!is.null(grp_by_mirna)) paste0('<div class="sublabel">Grouped by miRNA</div>', .table_html(grp_by_mirna)) else "",
    if (!is.null(grp_by_disease)) paste0('<div class="sublabel">Grouped by disease</div>', .table_html(grp_by_disease)) else ""
  )

  tab_network <- .subsec("Network",
    if (!is.null(img_files[["mda_net"]]))
      paste0('<p class="plot-desc">Network connecting miRNAs to diseases. Edge thickness represents number of supporting publications.</p>',
             .captured_img_html(img_files[["mda_net"]], "miRNA-Disease network"))
    else
      '<p class="no-data">Network image not captured — open the Network tab before generating the report.</p>'
  )

  ct_pairs_mirna   <- .cotarget_pairs_html(flat, "MIRNA_NAME", "DISEASE",  "Co-targeting pairs (miRNA-miRNA)")
  ct_pairs_disease <- .cotarget_pairs_html(flat, "DISEASE",    "MIRNA_NAME", "Co-targeting pairs (disease-disease)")
  ct_hm_mirna      <- .plot_block(img_files, "cotarget_mda_mirna",
    "miRNA-miRNA co-targeting heatmap: rows/columns are miRNAs; cell colour = Jaccard similarity based on shared disease targets. High Jaccard = miRNAs target the same diseases.",
    "Co-targeting heatmap (miRNA-miRNA)")
  ct_hm_disease    <- .plot_block(img_files, "cotarget_mda_disease",
    "Disease-disease co-targeting heatmap: rows/columns are diseases; cell colour = Jaccard similarity based on shared regulating miRNAs. High Jaccard = diseases share many regulating miRNAs.",
    "Co-targeting heatmap (disease-disease)")

  tab_cotarget <- .subsec("Co-targeting",
    if (nzchar(ct_pairs_mirna))   paste0(.variant_label("miRNA-miRNA"), ct_pairs_mirna)   else "",
    if (nzchar(ct_hm_mirna))      paste0(.variant_label("miRNA-miRNA"), ct_hm_mirna)      else "",
    if (nzchar(ct_pairs_disease)) paste0(.variant_label("disease-disease"), ct_pairs_disease) else "",
    if (nzchar(ct_hm_disease))    paste0(.variant_label("disease-disease"), ct_hm_disease) else "",
    if (!nzchar(ct_pairs_mirna) && !nzchar(ct_hm_mirna) && !nzchar(ct_pairs_disease) && !nzchar(ct_hm_disease))
      '<p class="no-data">Co-targeting images not captured — open the Co-targeting tab and switch between miRNA-miRNA and disease-disease views before generating the report.</p>'
    else ""
  )

  rob_img     <- .plot_block(img_files, "mda_rob",
    "Robustness curve: network edge count as minimum evidence threshold increases. Steep drop = many low-confidence associations.",
    "Network robustness")
  rob_ct_mirna  <- .plot_block(img_files, "mda_rob_ct_mirna",
    "Co-targeting robustness (miRNA-miRNA): number of miRNA pairs that co-target the same diseases across evidence thresholds.",
    "Co-targeting robustness (miRNA-miRNA)")
  rob_ct_disease <- .plot_block(img_files, "mda_rob_ct_disease",
    "Co-targeting robustness (disease-disease): number of disease pairs co-regulated by the same miRNAs across evidence thresholds.",
    "Co-targeting robustness (disease-disease)")

  tab_robustness <- .subsec("Network Analysis",
    rob_img,
    if (nzchar(rob_ct_mirna))   paste0(.variant_label("miRNA-miRNA"),    rob_ct_mirna)   else "",
    if (nzchar(rob_ct_disease)) paste0(.variant_label("disease-disease"), rob_ct_disease) else "",
    if (!nzchar(rob_img) && !nzchar(rob_ct_mirna) && !nzchar(rob_ct_disease))
      '<p class="no-data">Robustness images not captured — run the Network Analysis in the miRNA-Disease tab before generating the report.</p>'
    else ""
  )

  body <- paste0(tab_results, tab_network, tab_cotarget, tab_robustness)

  .section_wrap("mda", 1, "miRNA - Disease", "#2878b5", ts, filters, counts, body)
}

# ── Section: MGA ─────────────────────────────────────────────────────────────

.section_mga <- function(df, input, img_files = list(), ts, flat_df = NULL) {
  if (is.null(df)) return(.section_empty("mga", 2, "miRNA - Gene", "#27876a"))

  src        <- input$mga_mirna_source          %||% "disease"
  mirnas_sc  <- input$mirna_gene_select_scratch %||% character(0)
  mirnas_mda <- input$mirna_gene_select         %||% character(0)
  use_all_sc <- isTRUE(input$mga_use_all_scratch)
  arm        <- input$mga_mda_arm               %||% "5p/3p"
  proofs     <- input$proofs                    %||% character(0)
  min_pub    <- input$min_assoc_mga             %||% 1L
  only_met   <- isTRUE(input$mga_only_metabolic_genes)
  grouping   <- input$group_mga                 %||% "None"
  l1_ids     <- input$mga_pathway_l1            %||% character(0)
  l2_ids     <- input$mga_pathway_l2            %||% character(0)
  l1_names   <- tryCatch(.l1_names(l1_ids), error = function(e) as.character(l1_ids))
  l2_names   <- tryCatch(.l2_names(l2_ids), error = function(e) as.character(l2_ids))

  mirna_label <- switch(src,
    "scratch" = if (use_all_sc) "All miRNAs (scratch)"
                else if (length(mirnas_sc) > 0) mirnas_sc else "None selected",
    "disease" = if (length(mirnas_mda) > 0) mirnas_mda
                else paste0("From MDA (arm: ", arm, ")"),
    paste0("Source: ", src)
  )

  proof_match <- input$proof_match     %||% "any"
  proof_type  <- input$proof_type      %||% character(0)
  proof_str   <- input$proof_strength  %||% character(0)
  proof_tp    <- input$proof_throughput %||% character(0)

  filters <- paste0(
    '<div class="filter-box" style="border-left-color:#27876a;"><h3>Filters applied</h3>',
    .filter_row("miRNA source", switch(src, "scratch"="Manual input", "disease"="From miRNA-Disease", src), cls = "b-grey"),
    .filter_row("miRNAs",       mirna_label, cls = "b-blue"),
    if (src == "disease") .filter_row("miRNA arm", arm, cls = "b-grey") else "",
    .filter_row("Evidence types", if (length(proofs) == 0) "All" else proofs, cls = "b-green"),
    .filter_row("Evidence match mode",
      if (proof_match == "all") "All selected proofs required" else "Any selected proof", cls = "b-grey"),
    if (length(proof_type) > 0) .filter_row("Evidence — type",       proof_type, cls = "b-green") else "",
    if (length(proof_str)  > 0) .filter_row("Evidence — strength",   proof_str,  cls = "b-green") else "",
    if (length(proof_tp)   > 0) .filter_row("Evidence — throughput", proof_tp,   cls = "b-green") else "",
    .filter_row("Min. publications",    as.character(min_pub), cls = "b-grey"),
    .filter_row("Only metabolic genes", if (only_met) "Yes" else "No", cls = "b-grey"),
    if (length(l1_names) > 0) .filter_row("Pathway L1 (category)", l1_names, cls = "b-orange") else
      .filter_row("Pathway L1 (category)", "All", cls = "b-grey"),
    if (length(l2_names) > 0) .filter_row("Pathway L2 (specific)",  l2_names, cls = "b-orange") else
      .filter_row("Pathway L2 (specific)", "All", cls = "b-grey"),
    .filter_row("Table grouping", grouping, cls = "b-grey"),
    '</div>')

  flat <- flat_df
  n_mirna <- if ("MIRNA_NAME" %in% names(df)) dplyr::n_distinct(.strip_html(df$MIRNA_NAME)) else "—"
  n_gene  <- if ("GENE_NAME"  %in% names(df)) dplyr::n_distinct(.strip_html(df$GENE_NAME))  else "—"
  n_pairs <- if (!is.null(flat)) nrow(flat) else "—"
  counts  <- paste(.chip(n_mirna, "miRNAs"), .chip(n_gene, "genes"), .chip(n_pairs, "pairs"))

  # Summary tables from flat
  mirna_col <- if (!is.null(flat) && "MIRNA_NAME" %in% names(flat)) "MIRNA_NAME" else NULL
  gene_col  <- if (!is.null(flat) && "GENE_NAME"  %in% names(flat)) "GENE_NAME"  else NULL

  grp_by_mirna <- if (!is.null(flat) && !is.null(mirna_col) && !is.null(gene_col))
    tryCatch(flat %>%
      dplyr::group_by(miRNA = .data[[mirna_col]]) %>%
      dplyr::summarise(N_genes = dplyr::n_distinct(.data[[gene_col]]),
                       Genes   = paste(sort(unique(.data[[gene_col]])), collapse = ", "),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_genes)),
      error = function(e) NULL)
  else NULL

  grp_by_gene <- if (!is.null(flat) && !is.null(mirna_col) && !is.null(gene_col))
    tryCatch(flat %>%
      dplyr::group_by(Gene = .data[[gene_col]]) %>%
      dplyr::summarise(N_miRNAs = dplyr::n_distinct(.data[[mirna_col]]),
                       miRNAs   = paste(sort(unique(.data[[mirna_col]])), collapse = ", "),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_miRNAs)),
      error = function(e) NULL)
  else NULL

  metab_tbl <- if (!is.null(flat) && "is_metabolic" %in% names(flat) &&
                   !is.null(mirna_col) && !is.null(gene_col))
    tryCatch(flat %>%
      dplyr::filter(is_metabolic == TRUE) %>%
      dplyr::group_by(miRNA = .data[[mirna_col]]) %>%
      dplyr::summarise(N_metabolic_genes = dplyr::n_distinct(.data[[gene_col]]),
                       Genes = paste(sort(unique(.data[[gene_col]])), collapse = ", "),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_metabolic_genes)),
      error = function(e) NULL)
  else NULL

  tab_results <- .subsec("Results",
    '<div class="sublabel">Main result table</div>', .table_html(df),
    if (!is.null(grp_by_mirna)) paste0('<div class="sublabel">Grouped by miRNA</div>', .table_html(grp_by_mirna)) else "",
    if (!is.null(grp_by_gene))  paste0('<div class="sublabel">Grouped by gene</div>',  .table_html(grp_by_gene))  else "",
    if (!is.null(metab_tbl))    paste0('<div class="sublabel">Metabolic genes per miRNA</div>', .table_html(metab_tbl)) else ""
  )

  tab_network <- .subsec("Network",
    if (!is.null(img_files[["mga_net"]]))
      paste0('<p class="plot-desc">Network linking miRNAs to their target genes. Highly connected hub genes are candidate key regulators.</p>',
             .captured_img_html(img_files[["mga_net"]], "miRNA-Gene network"))
    else
      '<p class="no-data">Network image not captured — open the Network tab before generating the report.</p>'
  )

  ct_pairs_mirna <- .cotarget_pairs_html(flat, "MIRNA_NAME", "GENE_NAME",  "Co-targeting pairs (miRNA-miRNA)")
  ct_pairs_gene  <- .cotarget_pairs_html(flat, "GENE_NAME",  "MIRNA_NAME", "Co-targeting pairs (gene-gene)")
  ct_hm_mirna    <- .plot_block(img_files, "mga_hm_mirna",
    "miRNA-miRNA co-targeting heatmap: Jaccard similarity based on shared gene targets. Clustered miRNAs share overlapping target sets.",
    "Co-targeting heatmap (miRNA-miRNA)")
  ct_hm_gene     <- .plot_block(img_files, "mga_hm_gene",
    "Gene-gene co-targeting heatmap: Jaccard similarity based on shared regulating miRNAs. Genes in the same cluster are co-regulated.",
    "Co-targeting heatmap (gene-gene)")

  tab_cotarget <- .subsec("Co-targeting",
    if (nzchar(ct_pairs_mirna)) paste0(.variant_label("miRNA-miRNA"), ct_pairs_mirna) else "",
    if (nzchar(ct_hm_mirna))    paste0(.variant_label("miRNA-miRNA"), ct_hm_mirna)    else "",
    if (nzchar(ct_pairs_gene))  paste0(.variant_label("gene-gene"),   ct_pairs_gene)  else "",
    if (nzchar(ct_hm_gene))     paste0(.variant_label("gene-gene"),   ct_hm_gene)     else "",
    if (!nzchar(ct_pairs_mirna) && !nzchar(ct_hm_mirna) && !nzchar(ct_pairs_gene) && !nzchar(ct_hm_gene))
      '<p class="no-data">Co-targeting images not captured — open the Co-targeting tab and switch between variants before generating the report.</p>'
    else ""
  )

  rob_img    <- .plot_block(img_files, "mga_rob",
    "Network robustness: edge/node retention across evidence thresholds. Steep drop = many low-confidence edges.",
    "Network robustness")
  rob_ct_mirna  <- .plot_block(img_files, "mga_rob_ct_mirna",
    "Co-targeting robustness (miRNA-miRNA): co-targeting pair stability as evidence threshold increases.",
    "Co-targeting robustness (miRNA-miRNA)")
  rob_ct_gene   <- .plot_block(img_files, "mga_rob_ct_gene",
    "Co-targeting robustness (gene-gene): gene co-regulation stability across evidence thresholds.",
    "Co-targeting robustness (gene-gene)")

  tab_robustness <- .subsec("Network Analysis",
    rob_img,
    if (nzchar(rob_ct_mirna)) paste0(.variant_label("miRNA-miRNA"), rob_ct_mirna) else "",
    if (nzchar(rob_ct_gene))  paste0(.variant_label("gene-gene"),   rob_ct_gene)  else "",
    if (!nzchar(rob_img) && !nzchar(rob_ct_mirna) && !nzchar(rob_ct_gene))
      '<p class="no-data">Robustness images not captured — run Network Analysis in the miRNA-Gene tab before generating the report.</p>'
    else ""
  )

  body <- paste0(tab_results, tab_network, tab_cotarget, tab_robustness)

  .section_wrap("mga", 2, "miRNA - Gene", "#27876a", ts, filters, counts, body)
}

# ── Section: Pathway Annotation ───────────────────────────────────────────────

.section_ovr <- function(df, input, img_files = list(), ts) {
  if (is.null(df)) return(.section_empty("pathway", 3, "Pathway Annotation", "#b05a00"))
  
  pw_source         <- input$ovr_pathway_source      %||% "reactome"
  top_n             <- input$ontology_top_n           %||% 20L
  proofs            <- input$proofs_ovr               %||% character(0)
  ovr_l1_ids        <- input$ovr_pathway_l1           %||% character(0)
  ovr_l3_ids        <- input$ovr_pathway_l3           %||% character(0)
  ovr_descendants   <- isTRUE(input$ovr_include_descendants)
  ovr_group_l1      <- isTRUE(input$ovr_barplot_group_l1)
  ovr_l1_names      <- tryCatch(.l1_names(ovr_l1_ids), error = function(e) as.character(ovr_l1_ids))
  ovr_l3_names      <- tryCatch(.l2_names(ovr_l3_ids), error = function(e) as.character(ovr_l3_ids))

  filters <- paste0(
    '<div class="filter-box" style="border-left-color:#b05a00;"><h3>Filters applied</h3>',
    .filter_row("Evidence types",
                if (length(proofs) == 0) "All" else proofs, cls = "b-green"),
    .filter_row("Pathway database", switch(pw_source, "reactome"="Reactome", "gmt"="GMT file", pw_source), cls = "b-grey"),
    .filter_row("Top N pathways",   as.character(top_n), cls = "b-grey"),
    if (length(ovr_l1_names) > 0) .filter_row("L1 category filter", ovr_l1_names, cls = "b-orange") else
      .filter_row("L1 category filter", "All", cls = "b-grey"),
    if (length(ovr_l3_names) > 0) .filter_row("L3 specific pathway filter", ovr_l3_names, cls = "b-orange") else "",
    if (length(ovr_l1_ids) > 0 || length(ovr_l3_ids) > 0)
      .filter_row("Include sub-pathways", if (ovr_descendants) "Yes" else "No", cls = "b-grey") else "",
    .filter_row("Overview — group by L1", if (ovr_group_l1) "Yes" else "No", cls = "b-grey"),
    local({
      expr_mode <- input$expression_filter_mode_ovr %||% "none"
      if (expr_mode == "none" || identical(pw_source, "gmt")) return("")
      tissue_lbl <- paste(
        input$expression_tissue_select_ovr %||%
        input$expression_cell_lines_ovr    %||% "—",
        collapse = ", ")
      min_frac <- as.numeric(input$ovr_expr_min_frac %||% 20)
      nTPM_val <- if (expr_mode == "tissue")
        as.numeric(input$nTPM_expression_filter_ovr   %||% 0.1)
      else
        as.numeric(input$nTPM_2_expression_filter_ovr %||% 0.1)
      .filter_row("Expression filter",
        paste0(if (expr_mode == "tissue") "Tissue" else "Cell line", ": ", tissue_lbl,
               " | nTPM ≥ ", nTPM_val,
               " | ≥ ", min_frac, "% pathway genes expressed"),
        cls = "b-blue")
    }),
    '</div>')
  
  n_path  <- if ("PATHWAY_NAME" %in% names(df)) dplyr::n_distinct(df$PATHWAY_NAME) else "—"
  n_genes <- if ("N_FOUND"      %in% names(df)) sum(df$N_FOUND, na.rm = TRUE) else "—"
  counts  <- paste(.chip(n_path, "pathways"), .chip(n_genes, "gene hits total"))
  
  top_df <- tryCatch({
    df %>%
      dplyr::arrange(dplyr::desc(COVERAGE)) %>%
      dplyr::select(dplyr::any_of(c("PATHWAY_NAME","COVERAGE","N_FOUND","K","N_MIRNAS","MIRNA_NAME"))) %>%
      dplyr::slice_head(n = as.integer(top_n))
  }, error = function(e) head(df, 10L))
  
  ovr_img <- .plot_block(img_files, "ovr_plot", "Horizontal bar chart of the top enriched pathways ranked by coverage (fraction of pathway genes targeted by the selected miRNAs). Bars are coloured by pathway L1 category. Longer bars indicate a higher proportion of pathway members under miRNA regulation.", "Pathway annotation chart")

  ovr_hm_img <- .plot_block(img_files, "ovr_hm", "Heatmap showing miRNA-by-pathway coverage. Each cell represents the fraction of a pathway's genes targeted by a given miRNA. Row and column clustering reveals groups of miRNAs that co-regulate the same pathways, and pathways with broad vs. narrow miRNA regulation.", "Pathway annotation heatmap")

  l1_summary_html <- tryCatch({
    if ("L1_CATEGORY" %in% names(df) && any(!is.na(df$L1_CATEGORY)) &&
        "GENES_FOUND" %in% names(df)) {
      # miRNA source: per-mirna mode has MIRNA_NAME col; combined has MIRNAS_FOUND
      get_mirnas <- if ("MIRNA_NAME" %in% names(df)) {
        function(sub) unique(as.character(sub$MIRNA_NAME))
      } else if ("MIRNAS_FOUND" %in% names(df)) {
        function(sub) unique(unlist(strsplit(as.character(sub$MIRNAS_FOUND[!is.na(sub$MIRNAS_FOUND)]), ",\\s*")))
      } else {
        function(sub) character(0)
      }
      cats <- sort(unique(df$L1_CATEGORY[!is.na(df$L1_CATEGORY)]))
      rows_html <- paste(vapply(cats, function(cat) {
        sub  <- df[!is.na(df$L1_CATEGORY) & df$L1_CATEGORY == cat, ]
        n_m  <- length(get_mirnas(sub))
        genes <- unique(unlist(strsplit(
          as.character(sub$GENES_FOUND[!is.na(sub$GENES_FOUND) & nzchar(sub$GENES_FOUND)]),
          ",\\s*")))
        genes <- genes[nzchar(genes)]
        n_g  <- length(genes)
        n_pw <- dplyr::n_distinct(sub$PATHWAY_NAME)
        sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
                htmltools::htmlEscape(cat),
                format(n_pw, big.mark = ","),
                format(n_m, big.mark = ","),
                format(n_g, big.mark = ","))
      }, character(1)), collapse = "")
      paste0(
        '<div class="sublabel">Summary by L1 pathway category</div>',
        '<div class="twrap"><table><thead><tr>',
        '<th>L1 Category</th><th>Pathways enriched</th><th>miRNAs</th><th>Genes</th>',
        '</tr></thead><tbody>', rows_html, '</tbody></table></div>',
        '<p class="tnote">miRNA and gene counts are distinct across all enriched pathways in that category.</p>')
    } else ""
  }, error = function(e) "")

  body <- paste0(
    '<div class="sublabel">Top pathways by coverage</div>', .table_html(top_df),
    ovr_img, l1_summary_html, ovr_hm_img)

  .section_wrap("pathway", 3, "Pathway Annotation", "#b05a00",
                ts, filters, counts, body)
}


# ── Section: MR ──────────────────────────────────────────────────────────────

.section_mr <- function(df, input, img_files = list(), ts, flat_df = NULL) {
  if (is.null(df)) return(.section_empty("mr", 4, "miRNA - Reaction", "#7a3fbf"))

  mirnas_sel <- input$mirna_reaction_select         %||% character(0)
  genes_sel  <- input$gene_reaction_select          %||% character(0)
  filter_g   <- isTRUE(input$filter_genes_reaction)
  l1         <- input$mr_filter_l1           %||% character(0)
  l2         <- input$mr_filter_l2           %||% character(0)
  excl       <- input$mr_filter_exclude      %||% character(0)
  grouping   <- input$group_mr               %||% "None"
  bar_group  <- input$mr_barplot_group_by    %||% "subsystem"
  hm_group   <- input$metab_group_by         %||% "subsystem"
  hm_cluster <- isTRUE(input$metab_cluster)
  hm_palette <- input$metab_heatmap_palette  %||% "YlGnBu"

  mirna_label <- if (isTRUE(input$use_all_mirnas_from_mga)) "All miRNAs from miRNA-Gene"
                 else if (isTRUE(input$filter_rmirna_common_checkbox)) "Common miRNAs"
                 else if (length(mirnas_sel) > 0) mirnas_sel else "From miRNA-Gene"
  gene_label  <- if (filter_g) "All genes from miRNA-Gene table"
                 else if (length(genes_sel) > 0) genes_sel else "From miRNA-Gene"

  filters <- paste0(
    '<div class="filter-box" style="border-left-color:#7a3fbf;"><h3>Filters applied</h3>',
    .filter_row("miRNAs",     mirna_label, cls = "b-blue"),
    .filter_row("Genes",      gene_label,  cls = "b-blue"),
    if (length(l1) > 0) .filter_row("Category (L1)",  l1, cls = "b-orange") else .filter_row("Category (L1)", "All", cls = "b-grey"),
    if (length(l2) > 0) .filter_row("Subsystem (L2)", l2, cls = "b-orange") else .filter_row("Subsystem (L2)", "All", cls = "b-grey"),
    if (length(excl) > 0) .filter_row("Excluded subsystems", excl, cls = "b-orange") else "",
    .filter_row("Table grouping",         grouping,                           cls = "b-grey"),
    .filter_row("Overview — group by",    tools::toTitleCase(bar_group),      cls = "b-grey"),
    .filter_row("Heatmap — group by",     tools::toTitleCase(hm_group),       cls = "b-grey"),
    .filter_row("Heatmap — clustering",   if (hm_cluster) "Yes" else "No",   cls = "b-grey"),
    .filter_row("Heatmap — colour scale", hm_palette,                         cls = "b-grey"),
    '</div>')

  flat <- flat_df
  react_col <- if ("REACTION_NAME" %in% names(df)) "REACTION_NAME"
               else if ("NAME" %in% names(df)) "NAME" else NULL
  n_mirna <- if ("MIRNA_NAME" %in% names(df)) dplyr::n_distinct(.strip_html(df$MIRNA_NAME))   else "—"
  n_react <- if (!is.null(react_col))          dplyr::n_distinct(df[[react_col]])              else "—"
  n_sub   <- if ("SUBSYSTEM"  %in% names(df)) dplyr::n_distinct(df$SUBSYSTEM, na.rm = TRUE)   else "—"
  counts  <- paste(.chip(n_mirna, "miRNAs"), .chip(n_react, "reactions"), .chip(n_sub, "subsystems"))

  # Summary tables from flat
  grp_by_mirna <- if (!is.null(flat) && "MIRNA_NAME" %in% names(flat) && "SUBSYSTEM" %in% names(flat))
    tryCatch(flat %>%
      dplyr::group_by(miRNA = MIRNA_NAME) %>%
      dplyr::summarise(N_subsystems = dplyr::n_distinct(SUBSYSTEM, na.rm = TRUE),
                       N_reactions  = dplyr::n_distinct(REACTION_NAME %||% ".", na.rm = TRUE),
                       Subsystems   = paste(sort(unique(SUBSYSTEM[!is.na(SUBSYSTEM)])), collapse = ", "),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_subsystems)),
      error = function(e) NULL)
  else NULL

  grp_by_subsystem <- if (!is.null(flat) && "MIRNA_NAME" %in% names(flat) && "SUBSYSTEM" %in% names(flat))
    tryCatch(flat %>%
      dplyr::filter(!is.na(SUBSYSTEM) & nzchar(SUBSYSTEM)) %>%
      dplyr::group_by(Subsystem = SUBSYSTEM) %>%
      dplyr::summarise(N_miRNAs    = dplyr::n_distinct(MIRNA_NAME),
                       N_reactions = if ("REACTION_NAME" %in% names(flat)) dplyr::n_distinct(REACTION_NAME, na.rm = TRUE) else NA_integer_,
                       miRNAs      = paste(sort(unique(MIRNA_NAME)), collapse = ", "),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(N_miRNAs)),
      error = function(e) NULL)
  else NULL

  tab_results <- .subsec("Results",
    '<div class="sublabel">Main result table</div>', .table_html(df),
    if (!is.null(grp_by_mirna))     paste0('<div class="sublabel">Grouped by miRNA</div>',     .table_html(grp_by_mirna))     else "",
    if (!is.null(grp_by_subsystem)) paste0('<div class="sublabel">Grouped by subsystem</div>', .table_html(grp_by_subsystem)) else ""
  )

  tab_network <- .subsec("Network",
    if (!is.null(img_files[["mr_net"]]))
      paste0('<p class="plot-desc">Metabolic reaction network: miRNAs connected to metabolic reactions via shared target genes. Reaction nodes coloured by subsystem; essentiality score available on hover.</p>',
             .captured_img_html(img_files[["mr_net"]], "miRNA-Reaction network"))
    else
      '<p class="no-data">Network image not captured — open the Network tab before generating the report.</p>'
  )

  tab_subsystem <- .subsec("Subsystem Analysis",
    .plot_block(img_files, "mr_barplot",
      paste0("Overview bar chart: miRNA reaction coverage per metabolic subsystem/category (group by: ", tools::toTitleCase(bar_group), "). Longer bar = higher fraction of subsystem reactions regulated."),
      "Coverage overview"),
    ""
  )

  # Co-targeting pair tables from flat (miRNA×SUBSYSTEM)
  ct_flat <- if (!is.null(flat) && all(c("MIRNA_NAME", "SUBSYSTEM") %in% names(flat)))
    flat[!is.na(flat$SUBSYSTEM) & nzchar(flat$SUBSYSTEM), c("MIRNA_NAME", "SUBSYSTEM")]
  else NULL

  ct_pairs_mirna    <- .cotarget_pairs_html(ct_flat, "MIRNA_NAME", "SUBSYSTEM", "Co-targeting pairs (miRNA-miRNA, shared subsystems)")
  ct_pairs_subsys   <- .cotarget_pairs_html(ct_flat, "SUBSYSTEM",  "MIRNA_NAME", "Co-targeting pairs (subsystem-subsystem, shared miRNAs)")
  ct_hm_mirna       <- .plot_block(img_files, "mr_hm_mirna",
    "miRNA-miRNA co-targeting heatmap: Jaccard similarity based on shared metabolic subsystem targets.",
    "Co-targeting heatmap (miRNA-miRNA)")
  ct_hm_subsystem   <- .plot_block(img_files, "mr_hm_subsystem",
    "Subsystem-subsystem co-targeting heatmap: Jaccard similarity based on shared regulating miRNAs.",
    "Co-targeting heatmap (subsystem-subsystem)")

  tab_cotarget <- .subsec("Co-targeting",
    if (nzchar(ct_pairs_mirna))   paste0(.variant_label("miRNA-miRNA"),         ct_pairs_mirna)   else "",
    if (nzchar(ct_hm_mirna))      paste0(.variant_label("miRNA-miRNA"),         ct_hm_mirna)      else "",
    if (nzchar(ct_pairs_subsys))  paste0(.variant_label("subsystem-subsystem"), ct_pairs_subsys)  else "",
    if (nzchar(ct_hm_subsystem))  paste0(.variant_label("subsystem-subsystem"), ct_hm_subsystem)  else "",
    if (!nzchar(ct_pairs_mirna) && !nzchar(ct_hm_mirna) && !nzchar(ct_pairs_subsys) && !nzchar(ct_hm_subsystem))
      '<p class="no-data">Co-targeting images not captured — open the Co-targeting tab and switch between variants before generating the report.</p>'
    else ""
  )

  tab_essentiality <- .subsec("Essentiality",
    .plot_block(img_files, "ess_plot",
      "Essentiality profile: essential vs. partial reactions per miRNA. Higher essential fraction = miRNA preferentially regulates metabolically critical reactions.",
      "Essentiality plot")
  )

  body <- paste0(tab_results, tab_network, tab_subsystem, tab_cotarget, tab_essentiality)

  .section_wrap("mr", 4, "miRNA - Reaction", "#7a3fbf", ts, filters, counts, body)
}


# ── Main entry point ─────────────────────────────────────────────────────────

build_shinymir_report <- function(mda_df, mga_df, mr_df, ovr_df,
                                   input, img_files = list(),
                                   db_info = NULL, mr_enr = NULL,
                                   mda_flat_df = NULL,
                                   mga_flat_df = NULL,
                                   mr_flat_df  = NULL) {
  ts    <- format(Sys.time(), "%Y-%m-%d %H:%M")
  n_sec <- sum(!is.null(mda_df), !is.null(mga_df),
               !is.null(mr_df),  !is.null(ovr_df))

  db_block <- if (!is.null(db_info)) {
    fmt <- function(x) if (is.na(x) || is.null(x)) "—" else format(x, big.mark = ",")
    sprintf(
      '<div class="db-info-box">
        <span><strong>Database:</strong> %s</span>
        <span><strong>miRNAs:</strong> %s</span>
        <span><strong>Diseases:</strong> %s</span>
        <span><strong>Genes:</strong> %s</span>
        <span><strong>Reactions:</strong> %s</span>
      </div>',
      htmltools::htmlEscape(db_info$filename %||% "—"),
      fmt(db_info$n_mirnas),
      fmt(db_info$n_diseases),
      fmt(db_info$n_genes),
      fmt(db_info$n_reactions)
    )
  } else ""

  paste0(
    '<!DOCTYPE html><html lang="en"><head>',
    '<meta charset="UTF-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>ShinyMIR - Analysis Report</title>',
    .report_css,
    '</head><body><div class="report">',

    sprintf('
<div class="report-header">
  <div><h1>ShinyMIR - Analysis Report</h1>
    <div class="sub">Multi-omics miRNA association analysis</div></div>
  <div class="meta">
    <div><strong>Generated:</strong> %s</div>
    <div><strong>Sections with data:</strong> %d / 4</div>
  </div>
</div>', ts, n_sec),

    if (nzchar(db_block)) paste0('<div style="padding:12px 40px 0;">', db_block, '</div>') else "",

    '<div class="toc"><h2>Contents</h2><ul>
      <li><a href="#mda">1. miRNA-Disease</a></li>
      <li><a href="#mga">2. miRNA-Gene</a></li>
      <li><a href="#pathway">3. Pathway Annotation</a></li>
      <li><a href="#mr">4. miRNA-Reaction</a></li>
    </ul></div>',

    .section_mda(mda_df, input, img_files, ts, flat_df = mda_flat_df),
    .section_mga(mga_df, input, img_files, ts, flat_df = mga_flat_df),
    .section_ovr(ovr_df, input, img_files, ts),
    .section_mr(mr_df,  input, img_files, ts, flat_df = mr_flat_df),
    

    sprintf('<div class="report-footer">
      <span>Generated by <strong>ShinyMIR</strong></span>
      <span>%s</span></div>', ts),

    '</div></body></html>'
  )
}
