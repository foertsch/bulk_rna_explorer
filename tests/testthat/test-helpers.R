# Tests for pure helper + plot functions in R/helpers.R

# --- sanitize_name --------------------------------------------------------

test_that("sanitize_name strips 'Factor', dots, and trims whitespace", {
  expect_equal(sanitize_name("Condition [Factor]"), "Condition []")
  expect_equal(sanitize_name("Dox [Factor]"), "Dox []")
  expect_equal(sanitize_name("gene.name"), "genename")
  expect_equal(sanitize_name("  col  "), "col")
  expect_equal(sanitize_name("abc"), "abc")
})

test_that("sanitize_name handles a vector", {
  expect_equal(
    sanitize_name(c("A [Factor]", "B.", " C ")),
    c("A []", "B", "C")
  )
})

# --- generate_palette -----------------------------------------------------

test_that("generate_palette returns n colors", {
  expect_length(generate_palette(1), 1)
  expect_length(generate_palette(5), 5)
  expect_length(generate_palette(8), 8)
  expect_length(generate_palette(15), 15)
})

test_that("generate_palette uses the scanpy vega palette for small n", {
  expect_equal(generate_palette(3), c("#1f77b4", "#ff7f0e", "#279e68"))
})

test_that("generate_palette returns valid hex colors", {
  pal <- generate_palette(10)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}", pal)))
})

# --- detect_de_columns ----------------------------------------------------

test_that("detect_de_columns finds FGCZ conventions", {
  res <- detect_de_columns(c("gene_name", "log2Ratio", "fdr", "pValue"))
  expect_equal(res$fc, "log2Ratio")
  expect_equal(res$fdr, "fdr")
  expect_equal(res$pval, "pValue")
})

test_that("detect_de_columns finds standard DESeq2 conventions", {
  res <- detect_de_columns(c("log2FoldChange", "padj", "pvalue"))
  expect_equal(res$fc, "log2FoldChange")
  expect_equal(res$fdr, "padj")
  expect_equal(res$pval, "pvalue")
})

test_that("detect_de_columns returns NA when no match", {
  res <- detect_de_columns(c("foo", "bar"))
  expect_true(is.na(res$fc))
  expect_true(is.na(res$fdr))
  expect_true(is.na(res$pval))
})

# --- detect_gene_col ------------------------------------------------------

test_that("detect_gene_col finds common conventions", {
  expect_equal(detect_gene_col(c("gene_id", "gene_name")), "gene_name")
  expect_equal(detect_gene_col(c("Symbol", "x", "y")), "Symbol")
  expect_equal(detect_gene_col(c("hgnc_symbol", "foo")), "hgnc_symbol")
})

test_that("detect_gene_col returns NULL when no match", {
  expect_null(detect_gene_col(c("foo", "bar")))
  expect_null(detect_gene_col(character(0)))
})

# --- get_categorical_cols -------------------------------------------------

test_that("get_categorical_cols picks up factor and character columns", {
  se <- make_mock_se()
  cols <- get_categorical_cols(se)
  expect_true("Condition" %in% cols)
  expect_true("Dox" %in% cols)
  expect_true("Sample" %in% cols)
})

test_that("get_categorical_cols treats few-unique numerics as categorical", {
  se <- make_mock_se(n_samples = 6)
  # numeric_cov has 6 unique values (<= 10) so it's classified categorical
  cols <- get_categorical_cols(se)
  expect_true("numeric_cov" %in% cols)
})

test_that("colData/rowData readers emit no '...' ignored warnings", {
  # Guards the removal of the silently-ignored `check.names = FALSE` arg from
  # as.data.frame() on SummarizedExperiment DataFrames.
  se <- make_mock_se()
  expect_no_warning(get_categorical_cols(se))
  expect_no_warning(
    create_barplot(se, "ENSG00001", "xNorm", "Condition",
                   gene_symbol_col = "gene_name")
  )
  expect_no_warning(build_volcano_df(se, "log2Ratio", "fdr", "gene_name"))
})

# --- build_volcano_df -----------------------------------------------------

test_that("build_volcano_df produces expected columns", {
  se <- make_mock_se()
  df <- build_volcano_df(se, fc_col = "log2Ratio", fdr_col = "fdr",
                          gene_symbol_col = "gene_name")
  expect_s3_class(df, "data.frame")
  expect_true(all(c("gene_id", "log2FC", "fdr", "symbol",
                    "neg_log10_fdr", "Status") %in% colnames(df)))
  expect_s3_class(df$Status, "factor")
  expect_equal(levels(df$Status), c("UP", "DOWN", "NS"))
})

test_that("build_volcano_df classifies UP/DOWN correctly", {
  se <- make_mock_se()
  df <- build_volcano_df(se, fc_col = "log2Ratio", fdr_col = "fdr",
                          gene_symbol_col = "gene_name",
                          fc_thresh = 1, fdr_thresh = 0.05)
  # Gene 1: log2 = 3, fdr = 0.001 -> UP
  # Gene 2: log2 = -3, fdr = 0.001 -> DOWN
  # Gene 3: log2 = 2.5, fdr = 0.01 -> UP
  row1 <- df[df$gene_id == "ENSG00001", ]
  row2 <- df[df$gene_id == "ENSG00002", ]
  row3 <- df[df$gene_id == "ENSG00003", ]
  expect_equal(as.character(row1$Status), "UP")
  expect_equal(as.character(row2$Status), "DOWN")
  expect_equal(as.character(row3$Status), "UP")
})

test_that("build_volcano_df returns NULL when DE cols are missing", {
  se <- make_mock_se()
  expect_null(build_volcano_df(se, fc_col = NULL, fdr_col = "fdr"))
  expect_null(build_volcano_df(se, fc_col = "log2Ratio", fdr_col = NULL))
})

test_that("build_volcano_df uses gene_id as fallback when symbol col missing", {
  se <- make_mock_se()
  df <- build_volcano_df(se, fc_col = "log2Ratio", fdr_col = "fdr",
                          gene_symbol_col = NULL)
  expect_equal(df$symbol, df$gene_id)
})

# --- create_barplot -------------------------------------------------------

test_that("create_barplot returns a ggplot (no second group)", {
  se <- make_mock_se()
  p <- create_barplot(
    se = se, gene = "ENSG00001",
    assay_name = "xNorm", group_col = "Condition",
    gene_symbol_col = "gene_name"
  )
  expect_s3_class(p, "ggplot")
})

test_that("create_barplot returns a ggplot (with second group)", {
  se <- make_mock_se()
  p <- create_barplot(
    se = se, gene = "ENSG00001",
    assay_name = "xNorm", group_col = "Condition",
    second_group = "Dox", gene_symbol_col = "gene_name"
  )
  expect_s3_class(p, "ggplot")
})

test_that("create_barplot fails on invalid gene or assay", {
  se <- make_mock_se()
  expect_error(create_barplot(se, "NOT_A_GENE", "xNorm", "Condition"))
  expect_error(create_barplot(se, "ENSG00001", "bad_assay", "Condition"))
})

test_that("create_barplot y-axis label reflects the log_transform toggle", {
  se <- make_mock_se()
  p_log <- create_barplot(se, "ENSG00001", "xNorm", "Condition",
                          log_transform = TRUE)
  p_raw <- create_barplot(se, "ENSG00001", "xNorm", "Condition",
                          log_transform = FALSE)
  expect_match(p_log$labels$y, "^log2\\(xNorm")
  expect_equal(p_raw$labels$y, "xNorm")
})

test_that("create_barplot tolerates FGCZ '[Factor]' colData names", {
  se <- make_mock_se()
  # Rename Condition -> "Condition [Factor]" as FGCZ Sushi output does
  cd <- colData(se)
  colnames(cd)[colnames(cd) == "Condition"] <- "Condition [Factor]"
  colData(se) <- cd
  # get_categorical_cols cleans the name to "Condition"; that cleaned name is
  # what the app feeds back to create_barplot, which re-derives the match.
  cats <- get_categorical_cols(se)
  expect_true("Condition" %in% cats)
  p <- create_barplot(se, "ENSG00001", "xNorm", "Condition",
                      gene_symbol_col = "gene_name")
  expect_s3_class(p, "ggplot")
})

# --- create_volcano -------------------------------------------------------

test_that("create_volcano returns a ggplot from prepared df", {
  se <- make_mock_se()
  df <- build_volcano_df(se, "log2Ratio", "fdr", "gene_name")
  p <- create_volcano(df)
  expect_s3_class(p, "ggplot")
})

test_that("create_volcano returns placeholder for NULL/empty df", {
  expect_s3_class(create_volcano(NULL), "ggplot")
  expect_s3_class(create_volcano(data.frame()), "ggplot")
})

# --- create_pca -----------------------------------------------------------

test_that("create_pca returns a ggplot", {
  se <- make_mock_se()
  p <- create_pca(se, assay_name = "xNorm", group_col = "Condition")
  expect_s3_class(p, "ggplot")
})

test_that("create_pca respects ntop larger than gene count", {
  se <- make_mock_se(n_genes = 20)
  p <- create_pca(se, assay_name = "xNorm", group_col = "Condition",
                  ntop = 1000)
  expect_s3_class(p, "ggplot")
})

test_that("create_pca with ellipses does not error", {
  se <- make_mock_se()
  p <- create_pca(se, assay_name = "xNorm", group_col = "Condition",
                  show_ellipses = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("create_pca with second group (shape aesthetic)", {
  se <- make_mock_se()
  p <- create_pca(se, assay_name = "xNorm", group_col = "Condition",
                  second_group = "Dox")
  expect_s3_class(p, "ggplot")
})

# --- compute_pca / plot_pca -----------------------------------------------

test_that("compute_pca returns projection, pct_var, and sample names", {
  se <- make_mock_se()
  pca <- compute_pca(se, "xNorm", ntop = 100)
  expect_true(all(c("x", "pct_var", "samples") %in% names(pca)))
  expect_equal(nrow(pca$x), ncol(se))
  expect_equal(pca$samples, colnames(se))
})

test_that("plot_pca builds a ggplot from a compute_pca result", {
  se <- make_mock_se()
  pca <- compute_pca(se, "xNorm")
  expect_s3_class(plot_pca(pca, se, group_col = "Condition"), "ggplot")
  expect_s3_class(
    plot_pca(pca, se, group_col = "Condition", second_group = "Dox"),
    "ggplot"
  )
})

test_that("plot_pca clamps axis choices beyond the available PCs", {
  # 3 samples -> at most 3 PCs; asking for PC4 must not error
  se <- make_mock_se(n_genes = 30, n_samples = 3)
  pca <- compute_pca(se, "xNorm", ntop = 30)
  expect_lte(ncol(pca$x), 3)
  expect_s3_class(
    plot_pca(pca, se, group_col = "Condition", pc_x = 1, pc_y = 4),
    "ggplot"
  )
})

test_that("compute_pca log_transform toggle changes the projection", {
  se <- make_mock_se()
  pca_log <- compute_pca(se, "xNorm", ntop = 50, log_transform = TRUE)
  pca_raw <- compute_pca(se, "xNorm", ntop = 50, log_transform = FALSE)
  expect_equal(dim(pca_log$x), dim(pca_raw$x))
  expect_false(isTRUE(all.equal(pca_log$x, pca_raw$x)))
})

# --- save_plot_file -------------------------------------------------------

test_that("save_plot_file writes valid PNG/PDF to a Shiny-style temp path", {
  p <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  # Shiny hands content() tempfile(fileext = file_ext(name)); for "plot.png"
  # that is tempfile(fileext = "png") -> ".../fileXXXXpng", a path with no real
  # extension. ggsave() then needs `device` explicitly -- this guards the
  # download bug where the buttons produced no file.
  png_path <- tempfile(fileext = tools::file_ext("plot.png"))
  save_plot_file(p, png_path, "png", width = 6, height = 4, dpi = 150)
  expect_true(file.exists(png_path) && file.info(png_path)$size > 0)
  expect_identical(readBin(png_path, "raw", 4L),
                   as.raw(c(0x89, 0x50, 0x4e, 0x47)))  # PNG magic bytes

  pdf_path <- tempfile(fileext = tools::file_ext("plot.pdf"))
  save_plot_file(p, pdf_path, "pdf", width = 6, height = 4)
  expect_true(file.exists(pdf_path) && file.info(pdf_path)$size > 0)
  expect_identical(rawToChar(readBin(pdf_path, "raw", 5L)), "%PDF-")
})

# --- classify_de_status ---------------------------------------------------

test_that("classify_de_status labels UP/DOWN/NS by thresholds", {
  status <- classify_de_status(
    log2fc = c(2, -2, 0.5, 2),
    fdr = c(0.01, 0.01, 0.01, 0.5),
    fc_thresh = 1, fdr_thresh = 0.05
  )
  expect_equal(status, c("UP", "DOWN", "NS", "NS"))
})

# --- none_to_null ---------------------------------------------------------

test_that("none_to_null maps sentinels to NULL and passes values through", {
  expect_null(none_to_null("(none)"))
  expect_null(none_to_null(NULL))
  expect_null(none_to_null("None", "None"))
  expect_equal(none_to_null("log2FC"), "log2FC")
  expect_equal(none_to_null("None"), "None")
})

# --- select_cols ----------------------------------------------------------

test_that("select_cols keeps requested columns and drops sentinels/absent", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  expect_equal(colnames(select_cols(df, c("a", "c"))), c("a", "c"))
  expect_equal(colnames(select_cols(df, c("a", "(none)", "zzz"))), "a")
})

test_that("select_cols falls back to first n columns when none usable", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9, d = 1:3)
  expect_equal(colnames(select_cols(df, "(none)", fallback_n = 2)), c("a", "b"))
})
