mirna_target_graph_fun <- function(df, maximum_nodo_graph = 200,
                                   physics_type = "forceAtlas2Based",
                                   filter_edge,
                                   type_network = "DISEASE",
                                   color_essentiality = FALSE,
                                   show_notifications = TRUE,
                                   mirna_color        = "#FF6666",
                                   target_color       = "#6699FF",
                                   metabolic_color    = "#66CC66",
                                   nonmetabolic_color = "#6699FF",
                                   disease_group_map   = NULL,
                                   gene_pathway_labels = NULL,
                                   gene_color_map      = NULL,
                                   mirna_dsi_map       = NULL,
                                   reaction_color_by   = "subsystem",
                                   subsystem_palette   = NULL,
                                   sel_input_id        = NULL,
                                   node_click_input_id = NULL) {
  # normalize: REACTION_NAME column → NAME so the NAME branch below matches
  if (type_network == "REACTION_NAME") {
    if ("REACTION_NAME" %in% names(df)) df <- dplyr::rename(df, NAME = REACTION_NAME)
    type_network <- "NAME"
  }

  use_per_node_disease_color <- FALSE  # may be set TRUE in the DISEASE branch below

  # DSI → per-node miRNA color (red=promiscuous/DSI≈0, blue=specific/DSI≈1)
  dsi_pal <- colorRampPalette(c("#d73027", "#fee090", "#4575b4"))
  mirna_node_colors <- function(names) {
    if (is.null(mirna_dsi_map) || length(mirna_dsi_map) == 0)
      return(rep(mirna_color, length(names)))
    vals <- mirna_dsi_map[names]
    vals[is.na(vals)] <- 0.5
    dsi_pal(101)[pmin(pmax(round(vals * 100L), 0L), 100L) + 1L]
  }
  mirna_node_titles <- function(names) {
    base <- paste("miRNA:", names)
    if (is.null(mirna_dsi_map) || length(mirna_dsi_map) == 0) return(base)
    vals <- mirna_dsi_map[names]
    ifelse(is.na(vals), base, paste0(base, "<br>DSI: ", round(vals, 3)))
  }

  clean_col <- paste0(type_network, "_clean")
  if (clean_col %in% names(df)) {
    df[[type_network]] <- df[[clean_col]]
  }
  
  if ("PUBMED_ID" %in% colnames(df)){
    df$PUBMED_ID <- NULL
  }
  if (nrow(df) == 0) return(NULL)
  
  
  nodes_mirna <- unique(df$MIRNA_NAME)
  
  if (type_network == "NAME") {
    # ── Colora i nodi reazione per SUBSYSTEM o CATEGORY ───────
    if (!"HUMAN_ID"         %in% colnames(df)) df$HUMAN_ID         <- df$NAME
    # derive CATEGORY from SUBSYSTEM via global map
    if (!"CATEGORY" %in% colnames(df) && "SUBSYSTEM" %in% colnames(df))
      df$CATEGORY <- SUBSYSTEM_CATEGORY_MAP[df$SUBSYSTEM]

    color_col <- if (reaction_color_by == "category" && "CATEGORY" %in% colnames(df) &&
                     any(!is.na(df$CATEGORY))) "CATEGORY" else "SUBSYSTEM"

    optional_cols <- intersect(c("FORMULA", "GPR"), colnames(df))
    reaction_subsystem <- df %>%
      dplyr::select(dplyr::all_of(c("NAME", "SUBSYSTEM", "CATEGORY", "HUMAN_ID", optional_cols))) %>%
      dplyr::distinct(NAME, .keep_all = TRUE)

    reactions  <- unique(df[[type_network]])
    groups_vec <- unique(reaction_subsystem[[color_col]])
    groups_vec <- sort(groups_vec[!is.na(groups_vec) & groups_vec != ""])

    if (!is.null(subsystem_palette) && length(subsystem_palette) > 0) {
      # Use pre-built palette from server (guarantees legend ↔ network consistency)
      subsystem_color <- subsystem_palette
      # Extend for any group not in shared palette (fallback grey)
      missing_groups <- setdiff(groups_vec, names(subsystem_color))
      if (length(missing_groups) > 0)
        subsystem_color[missing_groups] <- "#AAAAAA"
    } else {
      n_sub <- length(groups_vec)
      if (n_sub == 0) {
        sub_palette <- character(0)
      } else if (n_sub <= 8) {
        sub_palette <- RColorBrewer::brewer.pal(max(3, n_sub), "Set2")[seq_len(n_sub)]
      } else if (n_sub <= 12) {
        sub_palette <- RColorBrewer::brewer.pal(n_sub, "Set3")[seq_len(n_sub)]
      } else {
        sub_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_sub)
      }
      subsystem_color <- setNames(sub_palette, groups_vec)
    }

    reaction_nodes <- reaction_subsystem %>%
      dplyr::filter(NAME %in% reactions) %>%
      dplyr::mutate(
        subsystem_label = dplyr::if_else(
          !is.na(.data[[color_col]]) & .data[[color_col]] != "",
          .data[[color_col]], "Unknown"
        ),
        node_color = dplyr::if_else(
          subsystem_label != "Unknown",
          subsystem_color[subsystem_label],
          "#AAAAAA"
        )
      )
    
    # Label mostrato sul nodo: HUMAN_ID (più corto); fallback a NAME se mancante
    reaction_label <- dplyr::coalesce(
      dplyr::na_if(as.character(reaction_nodes$HUMAN_ID), ""),
      reaction_nodes$NAME
    )
    
    nodes <- data.frame(
      id    = c(nodes_mirna, reaction_nodes$NAME),
      label = c(nodes_mirna, reaction_label),
      group = c(rep("miRNA", length(nodes_mirna)),
                paste0("sub_", reaction_nodes$subsystem_label)),
      title = c(mirna_node_titles(nodes_mirna),
                paste0("<b>", reaction_label, "</b><br>",
                       "Name: ",      reaction_nodes$NAME, "<br>",
                       "Subsystem: ", reaction_nodes$subsystem_label, "<br>",
                       "Formula: ",   dplyr::coalesce(reaction_nodes$FORMULA, ""), "<br>",
                       "GPR: ",       dplyr::coalesce(reaction_nodes$GPR, ""))),
      shape = c(rep("dot",  length(nodes_mirna)),
                rep("box",  nrow(reaction_nodes))),
      color = c(mirna_node_colors(nodes_mirna),
                reaction_nodes$node_color),
      size  = 20,
      stringsAsFactors = FALSE
    )
    
  } else if (type_network != "GENE_NAME") {
    disease_nodes <- unique(df[[type_network]])

    # ── Per-group colouring when disease_group_map is supplied ──────────────
    use_per_node_disease_color <- !is.null(disease_group_map) &&
                                  is.data.frame(disease_group_map) &&
                                  nrow(disease_group_map) > 0

    if (use_per_node_disease_color) {
      grp_unique <- sort(unique(na.omit(disease_group_map$ANCESTOR_GROUP)))
      n_grp      <- length(grp_unique)
      grp_col    <- setNames(
        GROUP_PALETTE[((seq_len(n_grp) - 1L) %% length(GROUP_PALETTE)) + 1L],
        grp_unique
      )

      # O(1) lookup via named vector – avoids match() inside sapply (was O(n*m))
      anc_lookup      <- setNames(disease_group_map$ANCESTOR_GROUP, disease_group_map$DISEASE)
      disease_col_vec <- vapply(disease_nodes, function(d) {
        grp <- anc_lookup[[d]]
        if (!is.null(grp) && !is.na(grp) && grp %in% names(grp_col)) grp_col[[grp]] else target_color
      }, character(1))

      nodes <- data.frame(
        id    = c(nodes_mirna, disease_nodes),
        label = c(nodes_mirna, disease_nodes),
        group = c(rep("miRNA", length(nodes_mirna)), rep(type_network, length(disease_nodes))),
        title = c(mirna_node_titles(nodes_mirna), paste(type_network, ":", disease_nodes)),
        shape = c(rep("dot", length(nodes_mirna)), rep("box", length(disease_nodes))),
        color = c(mirna_node_colors(nodes_mirna), disease_col_vec),
        size  = 20,
        stringsAsFactors = FALSE
      )
    } else {
      nodes <- data.frame(
        id    = c(nodes_mirna, disease_nodes),
        label = c(nodes_mirna, disease_nodes),
        group = c(rep("miRNA", length(nodes_mirna)), rep(type_network, length(disease_nodes))),
        title = c(mirna_node_titles(nodes_mirna), paste(type_network, ":", disease_nodes)),
        shape = c(rep("dot", length(nodes_mirna)), rep("box", length(disease_nodes))),
        color = c(mirna_node_colors(nodes_mirna), rep(target_color, length(disease_nodes))),
        size  = 20,
        stringsAsFactors = FALSE
      )
    }
  } else {
    
    nodes_gene_metabolic <- df %>%
      filter(is_metabolic == TRUE) %>%
      pull({{ type_network }}) %>%  # oppure pull(!!sym(type_network))
      unique()
    nodes_gene_nonmetabolic <- df %>%
      filter(is_metabolic == FALSE) %>%
      pull({{ type_network }}) %>%  # oppure pull(!!sym(type_network))
      unique()
    
    # ── Build gene tooltips (pathway info only when filter active) ──────────
    gene_tooltip <- function(gene_vec) {
      if (!is.null(gene_pathway_labels) && length(gene_pathway_labels) > 0) {
        sapply(gene_vec, function(g) {
          pw <- gene_pathway_labels[[g]]
          if (!is.null(pw) && nzchar(pw))
            paste0("<b>", g, "</b><br>Pathways: ", pw)
          else
            paste(type_network, ":", g)
        }, USE.NAMES = FALSE)
      } else {
        paste(type_network, ":", gene_vec)
      }
    }

    title_mirna    <- mirna_node_titles(nodes_mirna)
    title_metabolic    <- if (length(nodes_gene_metabolic)    > 0) gene_tooltip(nodes_gene_metabolic)    else character(0)
    title_nonmetabolic <- if (length(nodes_gene_nonmetabolic) > 0) gene_tooltip(nodes_gene_nonmetabolic) else character(0)

    all_genes <- c(nodes_gene_metabolic, nodes_gene_nonmetabolic)
    nodes <- data.frame(
      id = c(nodes_mirna, all_genes),
      label = c(nodes_mirna, all_genes),
      group = c(rep("miRNA", length(nodes_mirna)),
                rep("metabolic",     length(nodes_gene_metabolic)),
                rep("non_metabolic", length(nodes_gene_nonmetabolic))),
      title = c(title_mirna, title_metabolic, title_nonmetabolic),
      shape = c(rep("dot", length(nodes_mirna)),
                rep("box", length(all_genes))),
      size  = rep(20, length(nodes_mirna) + length(all_genes)),
      color = c(mirna_node_colors(nodes_mirna),
                rep(target_color, length(all_genes))),
      stringsAsFactors = FALSE
    )

    # ── Per-gene pathway colouring (overrides fixed gene color) ──────────────
    if (!is.null(gene_color_map) && length(gene_color_map) > 0) {
      gene_cols <- sapply(all_genes, function(g) {
        col <- gene_color_map[[g]]
        if (!is.null(col) && !is.na(col)) col else target_color
      }, USE.NAMES = FALSE)
      nodes$color <- c(mirna_node_colors(nodes_mirna), gene_cols)
    }

  }
  if (type_network!="NAME" && type_network!="GENE_NAME"){
    groups <- list(
      miRNA = list(color = list(background = mirna_color,  border = "#B22222")),
      target = list(color = list(background = target_color, border = "#1E90FF"))
    )
  }else if (type_network == "GENE_NAME") {
    groups <- list(
      miRNA         = list(color = list(background = mirna_color,        border = "#B22222")),
      metabolic     = list(color = list(background = metabolic_color,    border = "#228B22")),
      non_metabolic = list(color = list(background = nonmetabolic_color, border = "#6699FF"))
    )
  }
  if (type_network == "NAME") {
    if (!"GENES_CONTROLLED" %in% colnames(df)) df$GENES_CONTROLLED <- NA_integer_
    if (!"TOTAL_GENES"      %in% colnames(df)) df$TOTAL_GENES      <- NA_integer_
    if (!"ESS_FRAC"         %in% colnames(df)) df$ESS_FRAC         <- NA_real_

    edges_filtered <- df %>%
      group_by(from = MIRNA_NAME, to = !!sym(type_network)) %>%
      summarise(
        ess_frac         = suppressWarnings(max(ESS_FRAC,         na.rm = TRUE)),
        genes_controlled = suppressWarnings(max(GENES_CONTROLLED, na.rm = TRUE)),
        total_genes      = suppressWarnings(max(TOTAL_GENES,      na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      dplyr::mutate(
        ess_frac         = ifelse(is.finite(ess_frac),         ess_frac,         NA_real_),
        genes_controlled = ifelse(is.finite(genes_controlled), genes_controlled, NA_integer_),
        total_genes      = ifelse(is.finite(total_genes),      total_genes,      NA_integer_)
      ) %>%
      dplyr::filter(is.na(ess_frac) | ess_frac >= 0) %>%  # keep all; ess_frac filter is display-only
      dplyr::mutate(
        width = 1 + 7 * dplyr::coalesce(ess_frac, 0),
        title = sprintf(
          "Controlled genes: %s / %s<br>Fraction: %s",
          ifelse(is.na(genes_controlled), "?", as.character(genes_controlled)),
          ifelse(is.na(total_genes),      "?", as.character(total_genes)),
          ifelse(is.na(ess_frac),         "?", sprintf("%.2f", ess_frac))
        )
      )
  } else {
    df$COUNTER <- as.numeric(df$PUBMED_COUNT)
    df$PUBMED_COUNT = NULL
    
    edges_filtered <- df %>%
      group_by(from = MIRNA_NAME, to = !!sym(type_network)) %>%
      summarise(
        width = sum(COUNTER, na.rm = TRUE),
        title = paste0("Weight: ", width),
        .groups = 'drop'
      ) %>%
      filter(width >= filter_edge)
  }
  
  # ── Essenzialità: colora gli archi per coppia (miRNA, reazione) ──────
  if (type_network == "NAME" && color_essentiality &&
      ("IS_ESSENTIAL" %in% colnames(df) || "GPR" %in% colnames(df))) {
    edge_essentiality <- if ("IS_ESSENTIAL" %in% colnames(df)) {
      # Use pre-computed IS_ESSENTIAL (Boolean GPR evaluation)
      df %>%
        dplyr::select(MIRNA_NAME, NAME, IS_ESSENTIAL) %>%
        dplyr::distinct() %>%
        dplyr::group_by(MIRNA_NAME, NAME) %>%
        dplyr::summarise(
          is_ess     = any(IS_ESSENTIAL == TRUE, na.rm = TRUE),
          .groups    = "drop"
        ) %>%
        dplyr::mutate(
          edge_color = dplyr::if_else(is_ess, "#B22222", "#6699CC"),
          edge_label = dplyr::if_else(is_ess, "Essential", "Non-essential")
        ) %>%
        dplyr::select(MIRNA_NAME, NAME, edge_color, edge_label)
    } else {
      # Fallback: GPR binary check (no 'or' = essential)
      df %>%
        dplyr::select(MIRNA_NAME, NAME, GPR) %>%
        dplyr::distinct() %>%
        dplyr::group_by(MIRNA_NAME, NAME) %>%
        dplyr::summarise(
          is_ess  = any(!grepl("\\bor\\b", GPR, ignore.case = TRUE)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          edge_color = dplyr::if_else(is_ess, "#B22222", "#6699CC"),
          edge_label = dplyr::if_else(is_ess, "Essential", "Non-essential")
        ) %>%
        dplyr::select(MIRNA_NAME, NAME, edge_color, edge_label)
    }
    edges_filtered <- edges_filtered %>%
      dplyr::left_join(edge_essentiality,
                       by = c("from" = "MIRNA_NAME", "to" = "NAME")) %>%
      dplyr::mutate(
        color = dplyr::coalesce(edge_color, "#A9A9A9"),
        title = paste0(title, "<br>Essentiality: ", dplyr::coalesce(edge_label, "Unknown"))
      ) %>%
      dplyr::select(-edge_color, -edge_label)
  }

  # Rimuovi i nodi non connessi
  nodes_filtered <- nodes %>%
    filter(id %in% c(edges_filtered$from, edges_filtered$to))
  
  # Numero effettivo di nodi del grafo (era `length(nodes_filtered)`, che però su
  # un data.frame restituisce il numero di colonne, non di righe – bug storico).
  total_nodes <- nrow(nodes_filtered)
  
  # Sopra la soglia `maximum_nodo_graph` (500) il rendering animato con la fisica
  # lato browser diventa molto lento: si passa a un layout statico pre-calcolato
  # in R con igraph, che visualizza quasi istantaneamente anche grafi grandi.
  use_igraph_layout <- total_nodes > maximum_nodo_graph

  
  if (use_igraph_layout) {
    if (show_notifications) {
      showNotification(paste0("Large network (", total_nodes, " nodes): ..."),
                       type = "message", duration = 4)
    }
  }
  
  if (type_network!="NAME"){
    edges_filtered <- edges_filtered %>%
      dplyr::mutate(width = {
        w <- log2(width + 1)
        mn <- min(w, na.rm = TRUE)
        mx <- max(w, na.rm = TRUE)
        if (mx > mn) 1 + 7 * (w - mn) / (mx - mn) else rep(2, length(w))
      })
  }
  if (type_network == "NAME") {
    # I colori sono già nella colonna `color` dei nodi – visNetwork li usa direttamente
    net <- visNetwork(nodes_filtered, edges_filtered)
  } else if (type_network != "GENE_NAME") {
    net <- visNetwork(nodes_filtered, edges_filtered) %>%
      visGroups(groupname = "miRNA", color = groups$miRNA$color)
    if (!use_per_node_disease_color) {
      # Single colour for all disease nodes (default behaviour)
      net <- net %>% visGroups(groupname = type_network, color = groups$target$color)
    }
    # When use_per_node_disease_color = TRUE, per-node colours are in the `color`
    # column of nodes_filtered – visNetwork uses them directly, no visGroups needed.
  } else {
    # Gene nodes: colors set per-node in nodes$color column (fixed or pathway-based).
    # No visGroups for gene groups — group-level color would override per-node color.
    net <- visNetwork(nodes_filtered, edges_filtered) %>%
      visGroups(groupname = "miRNA", color = groups$miRNA$color)
  }
  net = net %>%
    visNodes(font = list(size = 18, color = "#343434")) %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
               nodesIdSelection = TRUE,
               selectedBy = "group") %>%
    visInteraction(navigationButtons = TRUE, dragNodes = TRUE, zoomView = FALSE,
                   multiselect = TRUE) %>%
    visLayout(randomSeed = 42)
  
  # Se l'essenzialità è attiva, i colori degli archi vengono dalla colonna `color`
  # del dataframe edges_filtered; altrimenti applichiamo il grigio di default.
  # NB: non impostiamo `width` qui per non sovrascrivere i pesi per-arco
  # (per NAME la width rappresenta la essentiality fraction).
  if (!(type_network == "NAME" && color_essentiality)) {
    net <- net %>%
      visEdges(color = list(color = "#A9A9A9", highlight = "#000000"))
  } else {
    net <- net %>%
      visEdges(smooth = list(enabled = TRUE, type = "continuous"))
  }
  
  if (use_igraph_layout || physics_type %in% c("no physics", "bipartite")) {
    if (physics_type == "bipartite") {
      # Bipartite: miRNA left column, targets right column, sorted by degree
      is_mirna <- nodes_filtered$group == "miRNA"
      mirna_ids  <- nodes_filtered$id[is_mirna]
      target_ids <- nodes_filtered$id[!is_mirna]
      n_mirna  <- length(mirna_ids)
      n_target <- length(target_ids)
      y_mirna  <- if (n_mirna  > 1) seq(1, n_mirna)  else 1
      y_target <- if (n_target > 1) seq(1, n_target) else 1
      coords <- data.frame(
        id = c(mirna_ids, target_ids),
        x  = c(rep(1,  n_mirna), rep(3, n_target)),
        y  = c(y_mirna * (max(n_mirna, n_target) / n_mirna),
               y_target)
      )
      nodes_filtered <- merge(nodes_filtered, coords, by = "id", all.x = TRUE)
      net <- visNetwork(nodes_filtered, edges_filtered) %>%
        visPhysics(enabled = FALSE) %>%
        visLayout(randomSeed = 42)
    } else {
      net <- net %>%
        visIgraphLayout(layout = "layout_with_fr") %>%
        visPhysics(enabled = FALSE)
    }
  } else if (physics_type == "hierarchical") {
    net <- net %>%
      visHierarchicalLayout(direction = "LR", sortMethod = "directed") %>%
      visPhysics(enabled = FALSE)
  } else {
    net <- net %>%
      visPhysics(
        solver        = physics_type,
        stabilization = list(
          enabled          = TRUE,
          iterations       = 100,
          updateInterval   = 50,
          onlyDynamicEdges = FALSE,
          fit              = TRUE
        )
      ) %>%
      visEvents(stabilizationIterationsDone =
        "function() { this.setOptions({ physics: { enabled: false } }); }")
  }
  
  # ── Node selection events for delete-node feature ──────────────────────────
  if (!is.null(sel_input_id)) {
    net <- net %>% visEvents(
      selectNode   = sprintf(
        "function(p){ Shiny.setInputValue('%s', p.nodes, {priority:'event'}); }",
        sel_input_id),
      deselectNode = sprintf(
        "function(p){ Shiny.setInputValue('%s', [],      {priority:'event'}); }",
        sel_input_id)
    )
  }

  # ── Node click → info modal ─────────────────────────────────────────────────
  if (!is.null(node_click_input_id)) {
    net <- net %>% visEvents(
      click = sprintf(
        "function(p){
           if(p.nodes.length > 0){
             var nd = this.body.data.nodes.get(p.nodes[0]);
             Shiny.setInputValue('%s', {id: p.nodes[0], group: nd ? nd.group : null}, {priority:'event'});
           }
         }",
        node_click_input_id
      )
    )
  }

  return(net)
}


# =============================================================================
# GPR Boolean evaluation — IS_ESSENTIAL
# =============================================================================
# Evaluate whether a miRNA can fully disrupt a reaction given its GPR logic.
# targeted_genes: character vector of genes silenced by this miRNA.
# Strategy: substitute each gene with TRUE (active) or FALSE (silenced),
# replace 'and'→'&' and 'or'→'|', then eval() the Boolean expression.
# Returns TRUE  = reaction is disrupted (miRNA has essential control).
# Returns FALSE = reaction can still proceed via untargeted genes.
# Returns NA    = GPR missing / parse error.
evaluate_gpr_disruption <- function(gpr, targeted_genes) {
  if (is.null(gpr) || is.na(gpr) || !nzchar(trimws(gpr))) return(NA)

  # Extract gene tokens (word chars, not 'and'/'or')
  tokens <- regmatches(gpr, gregexpr("[A-Za-z][A-Za-z0-9_.:-]*", gpr))[[1]]
  gene_tokens <- unique(tokens[!tolower(tokens) %in% c("and", "or")])
  if (length(gene_tokens) == 0) return(NA)

  expr <- gpr
  for (g in gene_tokens) {
    val  <- if (g %in% targeted_genes) "FALSE" else "TRUE"
    expr <- gsub(paste0("\\b", g, "\\b"), val, expr, perl = TRUE)
  }
  expr <- gsub("\\band\\b", " & ", expr, ignore.case = TRUE, perl = TRUE)
  expr <- gsub("\\bor\\b",  " | ", expr, ignore.case = TRUE, perl = TRUE)

  result <- tryCatch(eval(parse(text = expr)), error = function(e) NA)
  if (is.na(result)) return(NA)
  !isTRUE(result)   # disrupted = expression FALSE = IS_ESSENTIAL TRUE
}

# Adds IS_ESSENTIAL column to df (one value per (MIRNA_NAME, NAME) pair).
# Requires flat rows with GENE_NAME, NAME, MIRNA_NAME, GPR columns.
compute_is_essential <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  needed <- c("MIRNA_NAME", "NAME", "GPR", "GENE_NAME")
  if (!all(needed %in% colnames(df))) return(df)

  pairs <- df %>%
    dplyr::group_by(MIRNA_NAME, NAME) %>%
    dplyr::summarise(
      gpr_val  = dplyr::first(GPR),
      targeted = list(unique(GENE_NAME)),
      .groups  = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(IS_ESSENTIAL = evaluate_gpr_disruption(gpr_val, targeted)) %>%
    dplyr::ungroup() %>%
    dplyr::select(MIRNA_NAME, NAME, IS_ESSENTIAL)

  dplyr::left_join(df, pairs, by = c("MIRNA_NAME", "NAME"))
}

# =============================================================================
# ── Essentiality fraction per coppia (miRNA, reazione) ──────────────────────
# Aggiunge al df le colonne GENES_CONTROLLED, TOTAL_GENES, ESS_FRAC.
# GENES_CONTROLLED: geni che il miRNA targetta E che catalizzano la reazione.
# TOTAL_GENES:      geni totali della reazione nel GPR (indipendente dal miRNA).
# ESS_FRAC:         GENES_CONTROLLED / TOTAL_GENES ∈ (0, 1].
compute_essentiality_fraction <- function(con, df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  # --- Numeratore ---------------------------------------------------------
  if (!"GENES_CONTROLLED" %in% colnames(df)) {
    if ("GENE_COUNT" %in% colnames(df)) {
      # df già aggregato per (miRNA, reazione), con GENE_COUNT pronto
      df$GENES_CONTROLLED <- as.integer(df$GENE_COUNT)
    } else if ("GENE_NAME" %in% colnames(df)) {
      # Row-level o GROUP_CONCAT: conta geni distinti per coppia
      df <- df %>%
        dplyr::group_by(MIRNA_NAME, NAME) %>%
        dplyr::mutate(GENES_CONTROLLED = length(split_unique(GENE_NAME))) %>%
        dplyr::ungroup()
    } else {
      # Nessuna info sui geni nel df → query mirata
      gc <- dbGetQuery(con, glue_sql(
        "SELECT B.MIRNA_NAME, E.NAME AS NAME,
                COUNT(DISTINCT C.GENE_ID) AS GENES_CONTROLLED
           FROM MIRNAS_GENES_ARTICLES A
           JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
           JOIN GENES C ON A.GENE_ID = C.GENE_ID
           JOIN {`REACTIONS_GENES_TBL`} D ON C.GENE_ID = D.GENE_ID
           JOIN {`REACTIONS_TBL`} E ON D.REACTION_ID = E.REACTION_ID
          WHERE B.MIRNA_NAME IN ({mirnas*}) AND E.NAME IN ({reactions*})
          GROUP BY B.MIRNA_ID, B.MIRNA_NAME, E.REACTION_ID, E.NAME",
        mirnas    = unique(df$MIRNA_NAME),
        reactions = unique(df$NAME), .con = con))
      df <- dplyr::left_join(df, gc, by = c("MIRNA_NAME", "NAME"))
    }
  }
  
  # --- Denominatore -------------------------------------------------------
  # Use pre-loaded reaction_total_genes_map (global.R) to avoid per-call DB query.
  # Falls back to DB query only if a reaction name is not in the cache (e.g. after
  # metabolism upload before cache refresh).
  if (!"TOTAL_GENES" %in% colnames(df)) {
    rxns        <- unique(df$NAME)
    cached_vals <- reaction_total_genes_map[rxns]
    if (any(is.na(cached_vals))) {
      # Some reactions not in cache – fall back to DB for missing ones
      missing_rxns <- rxns[is.na(cached_vals)]
      tg_missing   <- dbGetQuery(con, glue_sql(
        "SELECT E.NAME AS NAME, COUNT(DISTINCT RG.GENE_ID) AS TOTAL_GENES
           FROM {`REACTIONS_TBL`} E JOIN {`REACTIONS_GENES_TBL`} RG ON E.REACTION_ID = RG.REACTION_ID
          WHERE E.NAME IN ({vals*}) GROUP BY E.REACTION_ID, E.NAME",
        vals = missing_rxns, .con = con))
      extra <- setNames(as.integer(tg_missing$TOTAL_GENES), tg_missing$NAME)
      cached_vals[names(extra)] <- extra
    }
    df$TOTAL_GENES <- as.integer(cached_vals[df$NAME])
  }
  
  # --- Frazione -----------------------------------------------------------
  if (all(c("GENES_CONTROLLED", "TOTAL_GENES") %in% colnames(df))) {
    df$ESS_FRAC <- ifelse(
      is.na(df$TOTAL_GENES) | df$TOTAL_GENES == 0,
      NA_real_,
      df$GENES_CONTROLLED / df$TOTAL_GENES
    )
  }
  
  df
}

safe_brewer_pal <- function(name, n = 9) {
  tryCatch(RColorBrewer::brewer.pal(n, name),
           error = function(e) RColorBrewer::brewer.pal(n, "YlGnBu"))
}

render_mirna_target_heatmap <- function(mat, cluster_rows, cluster_cols,
                                        x1_label, x2_label,
                                        col_annotation = NULL,
                                        row_annotation = NULL,
                                        l1_pathway_annotation = NULL,
                                        seriate = "mean",
                                        zero_white = FALSE,
                                        transpose = FALSE,
                                        color_palette = "YlGnBu") {
  # seriate = "mean"  → O(n log n), good default for interactive use
  # seriate = "OLO"   → O(n³), better quality but slow on large matrices (>50 nodes)
  # seriate = "none"  → no reordering
  if (is.null(mat)) return(NULL)

  if (transpose) {
    mat          <- t(mat)
    tmp          <- x1_label;     x1_label     <- x2_label;     x2_label     <- tmp
    tmp          <- cluster_rows; cluster_rows <- cluster_cols; cluster_cols <- tmp

    # col_annotation is a data.frame(DISEASE, ANCESTOR_GROUP); convert to named vector for row slot
    # row_annotation is a named vector (miRNA→label); convert to data.frame for col slot
    new_row_ann <- if (!is.null(col_annotation) && is.data.frame(col_annotation) &&
                       "DISEASE" %in% names(col_annotation) && "ANCESTOR_GROUP" %in% names(col_annotation)) {
      setNames(col_annotation$ANCESTOR_GROUP, col_annotation$DISEASE)
    } else {
      col_annotation
    }
    new_col_ann <- if (!is.null(row_annotation) && is.character(row_annotation) && !is.null(names(row_annotation))) {
      data.frame(DISEASE        = names(row_annotation),
                 ANCESTOR_GROUP = unname(row_annotation),
                 stringsAsFactors = FALSE)
    } else {
      row_annotation
    }
    col_annotation <- new_col_ann
    row_annotation <- new_row_ann
  }

  n_row <- nrow(mat)
  n_col <- ncol(mat)

  dendrogram      <- "none"
  show_dendrogram <- c(FALSE, FALSE)
  if (cluster_rows & cluster_cols) {
    dendrogram <- "both"; show_dendrogram <- c(TRUE, TRUE)
  } else if (cluster_rows) {
    dendrogram <- "row";    show_dendrogram <- c(TRUE, FALSE)
  } else if (cluster_cols) {
    dendrogram <- "column"; show_dendrogram <- c(FALSE, TRUE)
  }

  # ── Optional ancestor-group annotation bar above columns ─────────────────
  # col_annotation: data.frame(DISEASE, ANCESTOR_GROUP) from get_disease_descendant_map
  col_side <- NULL
  if (!is.null(col_annotation) && is.data.frame(col_annotation) &&
      nrow(col_annotation) > 0 && n_col > 0) {
    grp_vec <- col_annotation$ANCESTOR_GROUP[
      match(colnames(mat), col_annotation$DISEASE)
    ]
    grp_vec[is.na(grp_vec)] <- "Other"
    col_side <- data.frame(
      Ancestor = factor(grp_vec, levels = unique(grp_vec)),
      row.names = colnames(mat),
      stringsAsFactors = FALSE
    )
  }

  # ── Optional disease-membership annotation bar beside rows ───────────────
  # row_annotation: named character vector miRNA_NAME → disease label
  row_side <- NULL
  if (!is.null(row_annotation) && length(row_annotation) > 0 && n_row > 0) {
    labels  <- row_annotation[match(rownames(mat), names(row_annotation))]
    labels[is.na(labels)] <- "Unknown"
    row_side <- data.frame(
      Disease = factor(labels, levels = unique(na.omit(labels))),
      row.names = rownames(mat),
      stringsAsFactors = FALSE
    )
  }

  # L1 category annotation strip — uses heatmaply col_side / row_side (reliable inside subplots)
  if (!is.null(l1_pathway_annotation) && length(l1_pathway_annotation) > 0) {
    l1_idx <- if (transpose) rownames(mat) else colnames(mat)
    l1_vec <- l1_pathway_annotation[l1_idx]
    l1_vec[is.na(l1_vec)] <- "Unknown"
    l1_df <- data.frame(
      `L1 Category` = factor(l1_vec, levels = unique(na.omit(l1_vec))),
      row.names = l1_idx,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (transpose) {
      row_side <- if (is.null(row_side)) l1_df else cbind(row_side, l1_df)
    } else {
      col_side <- if (is.null(col_side)) l1_df else cbind(col_side, l1_df)
    }
  }

  safe_pal <- safe_brewer_pal(color_palette)
  pal <- if (zero_white)
    colorRampPalette(c("#FFFFFF", safe_pal[-1]))(255)
  else
    colorRampPalette(safe_pal)(255)

  heatmaply::heatmaply(
    mat,
    dendrogram      = dendrogram,
    show_dendrogram = show_dendrogram,
    colors          = pal,
    xlab            = x2_label,
    ylab            = x1_label,
    col_side_colors = col_side,
    row_side_colors = row_side,
    seriate         = if (cluster_rows || cluster_cols) seriate else "none",
    fontsize_row    = 9,
    fontsize_col    = 9
  ) %>%
    layout(
      margin     = list(l = 150, r = 80, t = if (!is.null(col_side)) 120 else 60,
                        b = max(150L, 6L * max(nchar(colnames(mat)), 0L))),
      showlegend = !is.null(col_side) || !is.null(row_side)
    ) %>%
    plotly_clean_config()
}


prepare_heatmap_data <- function(df, filter_count, x1, x2) {
  if (nrow(df) == 0) return(NULL)
  
  # Assicurati che x1 e x2 siano colonne del data frame
  if (!(x1 %in% colnames(df)) || !(x2 %in% colnames(df))) {
    stop("x1 or x2 not found in data frame columns.")
  }
  
  # Ordine dei pesi per la cella:
  #   1. ESS_FRAC  (frazione di essenzialità, per grafo miRNA-Reaction)
  #   2. PUBMED_COUNT (somma articoli, per disease/gene)
  #   3. n() (fallback)
  count_df <- if ("ESS_FRAC" %in% colnames(df)) {
    df %>%
      group_by(!!sym(x1), !!sym(x2)) %>%
      summarise(
        count = suppressWarnings(max(ESS_FRAC, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(count = ifelse(is.finite(count), count, 0)) %>%
      filter(count >= filter_count)
  } else {
    df %>%
      group_by(!!sym(x1), !!sym(x2)) %>%
      summarise(
        count = if ("PUBMED_COUNT" %in% colnames(df)) sum(PUBMED_COUNT, na.rm = TRUE) else n(),
        .groups = "drop"
      ) %>%
      filter(count >= filter_count)
  }
  
  # Pivot → matrix (avoid redundant tibble→data.frame→matrix conversions)
  wide <- tidyr::pivot_wider(count_df,
                             names_from  = all_of(x2),
                             values_from = "count",
                             values_fill = 0)
  mat         <- as.matrix(wide[, -1])
  rownames(mat) <- wide[[x1]]
  mat
}



