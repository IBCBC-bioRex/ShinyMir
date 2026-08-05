# R/cotarget.R
# Co-targeting analysis: Jaccard similarity between nodes based on shared entities

# -----------------------------------------------------------------------------
# Compute Jaccard similarity matrix.
# col_nodes:  column of "A" entities that become matrix rows/cols (e.g. "MIRNA_NAME")
# col_shared: column of "B" entities used as shared targets     (e.g. "GENE_NAME")
# Backward-compatible defaults kept for existing call sites.
# -----------------------------------------------------------------------------
compute_jaccard_matrix <- function(df, col_nodes = NULL, col_shared = NULL) {
  if (is.null(col_nodes))  col_nodes  <- if ("MIRNA_PREMATURE" %in% names(df)) "MIRNA_PREMATURE" else "MIRNA_NAME"
  if (is.null(col_shared)) col_shared <- "GENE_NAME"

  pairs <- unique(data.frame(
    node   = df[[col_nodes]],
    shared = df[[col_shared]],
    stringsAsFactors = FALSE
  ))
  pairs <- pairs[!is.na(pairs$node) & !is.na(pairs$shared), ]
  nodes <- sort(unique(pairs$node))
  if (length(nodes) < 2) return(NULL)

  # Binary shared x node matrix
  bin_mat     <- table(pairs$shared, pairs$node)
  bin_mat     <- (bin_mat > 0L) * 1L

  # Intersection counts via cross-product (node x node)
  intersection <- crossprod(bin_mat)
  sizes        <- diag(intersection)

  # Jaccard = |A ∩ B| / (|A| + |B| - |A ∩ B|)
  union_mat <- outer(sizes, sizes, "+") - intersection
  jaccard   <- intersection / union_mat
  diag(jaccard) <- 1

  jaccard
}


# -----------------------------------------------------------------------------
# Long-format table of co-targeting pairs above threshold.
# col_nodes, col_shared: same semantics as compute_jaccard_matrix.
# df_source: optional original data frame to count shared entities per pair.
# (Legacy param df_genes is an alias for df_source.)
# -----------------------------------------------------------------------------
cotarget_pairs_table <- function(mat, threshold = 0.1,
                                  df_source  = NULL,
                                  col_nodes  = "MIRNA_NAME",
                                  col_shared = "GENE_NAME",
                                  # legacy alias
                                  df_genes   = NULL) {
  if (is.null(mat)) return(NULL)

  # Legacy compat
  if (is.null(df_source) && !is.null(df_genes)) {
    df_source  <- df_genes
    col_nodes  <- if ("MIRNA_PREMATURE" %in% names(df_genes)) "MIRNA_PREMATURE" else "MIRNA_NAME"
    col_shared <- "GENE_NAME"
  }

  mat_work <- mat
  diag(mat_work) <- 0

  idx <- which(mat_work >= threshold & upper.tri(mat_work), arr.ind = TRUE)
  if (nrow(idx) == 0) return(NULL)

  node_names <- rownames(mat)
  label_a    <- paste0(col_nodes, "_A")
  label_b    <- paste0(col_nodes, "_B")

  result <- data.frame(
    A       = node_names[idx[, 1]],
    B       = node_names[idx[, 2]],
    Jaccard = round(mat_work[idx], 3),
    stringsAsFactors = FALSE
  )
  colnames(result)[1:2] <- c(label_a, label_b)

  if (!is.null(df_source) &&
      col_nodes  %in% names(df_source) &&
      col_shared %in% names(df_source)) {
    # Pre-deduplicate sets once – intersect on unique elements is faster
    shared_sets     <- lapply(
      split(df_source[[col_shared]], df_source[[col_nodes]]),
      unique
    )
    shared_col_name <- paste0("Common_", col_shared)
    result[[shared_col_name]] <- mapply(function(a, b) {
      sa <- shared_sets[[a]]; sb <- shared_sets[[b]]
      if (is.null(sa) || is.null(sb)) NA_integer_
      else length(intersect(sa, sb))
    }, result[[label_a]], result[[label_b]], SIMPLIFY = TRUE)
  }

  result[order(-result$Jaccard), ]
}


# -----------------------------------------------------------------------------
# Render co-targeting heatmaply (Jaccard matrix).
# cluster: logical – cluster rows and columns.
# node_label: axis label (e.g. "miRNA", "Disease", "Gene", "Subsystem").
# -----------------------------------------------------------------------------
render_cotarget_heatmap <- function(mat, cluster = FALSE, node_label = "miRNA",
                                    color_palette = "YlOrRd", mask_zero = TRUE,
                                    mask_diag = TRUE, annotation = NULL) {
  if (is.null(mat) || nrow(mat) < 2) return(NULL)

  if (cluster) {
    dist_mat <- as.dist(1 - mat)
    hc       <- hclust(dist_mat, method = "average")
    dend     <- as.dendrogram(hc)
  } else {
    dend <- FALSE
  }

  mat_plot <- mat
  if (mask_diag)  diag(mat_plot) <- NA
  if (mask_zero)  mat_plot[mat_plot == 0] <- NA

  safe_pal <- safe_brewer_pal(color_palette)

  # When annotation provided (no clustering): reorder rows/cols by disease then alpha
  has_ann <- !is.null(annotation) && length(annotation) > 0 &&
    all(rownames(mat_plot) %in% names(annotation))

  if (has_ann) {
    disease_order <- attr(annotation, "disease_order") %||% sort(unique(annotation))
    ann_labels    <- annotation[rownames(mat_plot)]
    rank_dis      <- match(ann_labels, disease_order)
    rank_dis[is.na(rank_dis)] <- length(disease_order) + 1L
    ord           <- order(rank_dis, rownames(mat_plot))
    mat_plot      <- mat_plot[ord, ord, drop = FALSE]
    ann_labels    <- ann_labels[ord]
  }

  make_ann_df <- function(labels, nms) data.frame(Disease = labels, row.names = nms, stringsAsFactors = FALSE)

  rsc <- if (has_ann) make_ann_df(ann_labels, rownames(mat_plot)) else NULL
  csc <- if (has_ann) make_ann_df(ann_labels, colnames(mat_plot)) else NULL

  p <- heatmaply::heatmaply(
    mat_plot,
    Rowv            = dend,
    Colv            = dend,
    colors          = colorRampPalette(safe_pal)(255),
    xlab            = node_label,
    ylab            = node_label,
    na.value        = "white",
    key.title       = "Jaccard",
    dendrogram      = if (cluster) "both" else "none",
    show_dendrogram = c(cluster, cluster),
    row_side_colors = rsc,
    col_side_colors = csc
  )

  if (has_ann) p <- suppress_sidebar_legends(p)

  p %>% plotly_clean_config()
}
