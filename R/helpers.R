# ============================================================
# Gene Expression Explorer - Pure Helper + Plot Functions
# Sourced by app.R; tested in tests/testthat/.
# Assumes these packages are attached by the caller:
#   SummarizedExperiment, ggplot2, dplyr, ggrepel
# ============================================================

# --- Helper Functions ---

#' Sanitize column names (remove Factor, dots, whitespace)
sanitize_name <- function(x) {
  x %>%
    gsub("Factor", "", ., fixed = TRUE) %>%
    gsub("\\.", "", .) %>%
    trimws()
}

#' Generate a default color palette for N conditions
generate_palette <- function(n) {
  if (n <= 8) {
    cols <- c("#808080", "#D4E6F1", "#7F9FB6", "#337996",
              "#085778", "#023E57", "#A3BE8C", "#BF616A")
    return(cols[seq_len(n)])
  }
  hcl.colors(n, palette = "Dark 3")
}

#' Auto-detect DE columns in rowData
detect_de_columns <- function(rd_cols) {
  fc_candidates <- c("log2Ratio", "log2FoldChange", "logFC", "log2ratio")
  fdr_candidates <- c("fdr", "padj", "FDR", "adj.P.Val", "p_adjusted")
  pval_candidates <- c("pValue", "pvalue", "PValue", "P.Value")

  fc_col <- intersect(fc_candidates, rd_cols)[1]
  fdr_col <- intersect(fdr_candidates, rd_cols)[1]
  pval_col <- intersect(pval_candidates, rd_cols)[1]

  list(fc = fc_col, fdr = fdr_col, pval = pval_col)
}

#' Auto-detect gene symbol column in rowData
detect_gene_col <- function(rd_cols) {
  candidates <- c("gene_name", "gene_symbol", "Symbol", "SYMBOL",
                   "external_gene_name", "hgnc_symbol")
  match <- intersect(candidates, rd_cols)[1]
  if (is.na(match)) NULL else match
}

#' Get categorical colData columns
get_categorical_cols <- function(se) {
  md <- as.data.frame(colData(se), check.names = FALSE)
  colnames(md) <- sanitize_name(colnames(md))
  cat_cols <- names(md)[sapply(md, function(x) {
    is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 10)
  })]
  cat_cols
}

# --- Plot Functions ---

#' Create barplot for gene expression
create_barplot <- function(
    se, gene, assay_name, group_col, second_group = NULL,
    group_levels = NULL, group_colors = NULL,
    ctrl_alpha = 0.5, gene_symbol_col = NULL,
    point_color = "black", jitter_width = 0.12,
    dodge_width = 0.8, bar_width = 0.7, err_width = 0.2
) {
  stopifnot(gene %in% rownames(se))
  stopifnot(assay_name %in% assayNames(se))

  md <- as.data.frame(colData(se), check.names = FALSE)
  colnames(md) <- sanitize_name(colnames(md))
  group_col <- sanitize_name(group_col)
  stopifnot(group_col %in% colnames(md))

  rd <- as.data.frame(rowData(se), check.names = FALSE)

  # Gene symbol for title
  gene_symbol <- gene
  if (!is.null(gene_symbol_col) && gene_symbol_col %in% colnames(rd)) {
    tmp <- rd[gene, gene_symbol_col]
    if (!is.na(tmp) && nchar(tmp) > 0) gene_symbol <- tmp
  }

  df <- md %>%
    dplyr::mutate(
      sample = colnames(se),
      value = as.numeric(assay(se, assay_name)[gene, ])
    ) %>%
    dplyr::transmute(
      sample = sample,
      value = value,
      value_log2 = log2(value + 1),
      Group = .data[[group_col]]
    )

  # Add second group if specified
  has_second <- !is.null(second_group) && second_group != "None"
  if (has_second) {
    second_group <- sanitize_name(second_group)
    if (second_group %in% colnames(md)) {
      df$SecondGroup <- as.factor(md[[second_group]])
    } else {
      has_second <- FALSE
    }
  }

  # Filter/order by selected group levels
  if (!is.null(group_levels) && length(group_levels) > 0) {
    df <- df %>%
      dplyr::filter(Group %in% group_levels) %>%
      dplyr::mutate(Group = factor(Group, levels = group_levels))
  } else {
    df <- df %>% dplyr::mutate(Group = as.factor(Group))
  }

  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data for selected groups") +
        theme_void()
    )
  }

  # Summarise
  if (has_second) {
    df_sum <- df %>%
      group_by(Group, SecondGroup) %>%
      summarise(mean = mean(value_log2, na.rm = TRUE),
                sd = sd(value_log2, na.rm = TRUE), .groups = "drop")
  } else {
    df_sum <- df %>%
      group_by(Group) %>%
      summarise(mean = mean(value_log2, na.rm = TRUE),
                sd = sd(value_log2, na.rm = TRUE), .groups = "drop")
  }

  pd <- position_dodge(width = dodge_width)

  p <- ggplot()

  if (has_second) {
    # Dodged bars by second group with alpha
    second_levels <- levels(df$SecondGroup)
    alpha_values <- setNames(
      c(ctrl_alpha, rep(1, length(second_levels) - 1)),
      second_levels
    )

    p <- p +
      geom_col(data = df_sum,
               aes(x = Group, y = mean, fill = Group, alpha = SecondGroup),
               width = bar_width, position = pd) +
      geom_errorbar(data = df_sum,
                    aes(x = Group, ymin = mean - sd, ymax = mean + sd,
                        group = SecondGroup),
                    width = err_width, position = pd) +
      geom_point(data = df,
                 aes(x = Group, y = value_log2, group = SecondGroup),
                 position = position_jitterdodge(
                   jitter.width = jitter_width, dodge.width = dodge_width),
                 size = 2, alpha = 0.9, color = point_color) +
      scale_alpha_manual(values = alpha_values, guide = "none")
  } else {
    p <- p +
      geom_col(data = df_sum,
               aes(x = Group, y = mean, fill = Group),
               width = bar_width) +
      geom_errorbar(data = df_sum,
                    aes(x = Group, ymin = mean - sd, ymax = mean + sd),
                    width = err_width) +
      geom_point(data = df,
                 aes(x = Group, y = value_log2),
                 position = position_jitter(width = jitter_width),
                 size = 2, alpha = 0.9, color = point_color)
  }

  if (!is.null(group_colors)) {
    p <- p + scale_fill_manual(values = group_colors, drop = FALSE)
  }

  p + theme_bw(base_size = 14) +
    labs(title = gene_symbol, x = NULL, y = "log2(normalized counts + 1)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    guides(fill = guide_legend(title = sanitize_name(group_col)))
}

#' Build volcano data frame from rowData (with gene_id, log2FC, fdr, symbol, Status)
build_volcano_df <- function(se, fc_col, fdr_col, gene_symbol_col = NULL,
                             fc_thresh = 1, fdr_thresh = 0.05) {
  rd <- as.data.frame(rowData(se), check.names = FALSE)
  if (is.null(fc_col) || is.null(fdr_col)) return(NULL)

  df <- data.frame(
    gene_id = rownames(rd),
    log2FC = as.numeric(rd[[fc_col]]),
    fdr = as.numeric(rd[[fdr_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(gene_symbol_col) && gene_symbol_col %in% colnames(rd)) {
    df$symbol <- rd[[gene_symbol_col]]
    df$symbol[is.na(df$symbol)] <- df$gene_id[is.na(df$symbol)]
  } else {
    df$symbol <- df$gene_id
  }
  df <- df %>% filter(!is.na(log2FC), !is.na(fdr), fdr > 0)
  df$neg_log10_fdr <- -log10(df$fdr)
  df$Status <- case_when(
    df$fdr < fdr_thresh & df$log2FC > fc_thresh ~ "UP",
    df$fdr < fdr_thresh & df$log2FC < -fc_thresh ~ "DOWN",
    TRUE ~ "NS"
  )
  df$Status <- factor(df$Status, levels = c("UP", "DOWN", "NS"))
  df
}

#' Create volcano plot from a prepared data frame
create_volcano <- function(df, fc_col = "log2FC", fdr_col = "FDR",
                           fc_thresh = 1, fdr_thresh = 0.05, n_labels = 10) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "DE columns not detected in rowData") +
             theme_void())
  }

  df_sig <- df %>%
    filter(Status != "NS") %>%
    arrange(fdr) %>%
    head(n_labels)

  status_colors <- c(UP = "#E74C3C", DOWN = "#3498DB", NS = "#CCCCCC")

  ggplot(df, aes(x = log2FC, y = neg_log10_fdr, color = Status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = status_colors) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(fdr_thresh),
               linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = df_sig, aes(label = symbol),
      size = 3, max.overlaps = 20, color = "black"
    ) +
    theme_bw(base_size = 14) +
    labs(x = paste0("log2 Fold Change (", fc_col, ")"),
         y = paste0("-log10 FDR (", fdr_col, ")"),
         title = "Volcano Plot") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

#' Create PCA plot
create_pca <- function(
    se, assay_name, group_col, second_group = NULL,
    group_colors = NULL, ntop = 500,
    pc_x = 1, pc_y = 2, show_ellipses = FALSE
) {
  mat <- log2(assay(se, assay_name) + 1)

  # Select top variable genes
  rv <- apply(mat, 1, var, na.rm = TRUE)
  ntop <- min(ntop, nrow(mat))
  select_genes <- order(rv, decreasing = TRUE)[seq_len(ntop)]
  mat_top <- mat[select_genes, ]

  # PCA
  pca <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)
  pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  md <- as.data.frame(colData(se), check.names = FALSE)
  colnames(md) <- sanitize_name(colnames(md))
  group_col <- sanitize_name(group_col)

  pca_df <- data.frame(
    Sample = colnames(se),
    PCx = pca$x[, pc_x],
    PCy = pca$x[, pc_y],
    Group = md[[group_col]],
    stringsAsFactors = FALSE
  )

  has_second <- !is.null(second_group) && second_group != "None"
  if (has_second) {
    second_group <- sanitize_name(second_group)
    if (second_group %in% colnames(md)) {
      pca_df$SecondGroup <- as.factor(md[[second_group]])
    } else {
      has_second <- FALSE
    }
  }

  p <- ggplot(pca_df, aes(x = PCx, y = PCy, color = Group))

  if (has_second) {
    p <- p + geom_point(aes(shape = SecondGroup), size = 4, alpha = 0.85)
  } else {
    p <- p + geom_point(size = 4, alpha = 0.85)
  }

  if (show_ellipses && length(unique(pca_df$Group)) > 2) {
    p <- p + stat_ellipse(aes(group = Group), type = "norm", level = 0.95,
                          alpha = 0.4)
  }

  if (!is.null(group_colors)) {
    p <- p + scale_color_manual(values = group_colors)
  }

  p + theme_bw(base_size = 14) +
    labs(
      x = paste0("PC", pc_x, " (", pct_var[pc_x], "%)"),
      y = paste0("PC", pc_y, " (", pct_var[pc_y], "%)"),
      title = "PCA Plot",
      color = group_col
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
