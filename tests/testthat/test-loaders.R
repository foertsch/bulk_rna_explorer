# Tests for input loading + format conversion in R/helpers.R

# --- detect_file_type -----------------------------------------------------

test_that("detect_file_type maps extensions (case-insensitive, .gz aware)", {
  expect_equal(detect_file_type("x.rds"), "rds")
  expect_equal(detect_file_type("x.RDS"), "rds")
  expect_equal(detect_file_type("x.qs2"), "qs2")
  expect_equal(detect_file_type("x.qs"), "qs")
  expect_equal(detect_file_type("x.csv"), "delim")
  expect_equal(detect_file_type("x.tsv"), "delim")
  expect_equal(detect_file_type("x.txt"), "delim")
  expect_equal(detect_file_type("counts.csv.gz"), "delim")
  expect_equal(detect_file_type("x.parquet"), "parquet")
  expect_equal(detect_file_type("x.feather"), "feather")
  expect_equal(detect_file_type("x.bam"), "unknown")
})

# --- sniff_delim ----------------------------------------------------------

test_that("sniff_delim detects comma, tab, and semicolon", {
  csv <- tempfile(fileext = ".csv")
  writeLines(c("a,b,c", "1,2,3"), csv)
  expect_equal(sniff_delim(csv), ",")

  tsv <- tempfile(fileext = ".tsv")
  writeLines(c("a\tb\tc", "1\t2\t3"), tsv)
  expect_equal(sniff_delim(tsv), "\t")

  scsv <- tempfile(fileext = ".csv")
  writeLines(c("a;b;c", "1;2;3"), scsv)
  expect_equal(sniff_delim(scsv), ";")
})

# --- classify_object ------------------------------------------------------

test_that("classify_object identifies known shapes", {
  expect_equal(classify_object(make_mock_se()), "se")
  expect_equal(classify_object(matrix(1:4, nrow = 2)), "matrix")
  expect_equal(classify_object(data.frame(a = 1)), "data.frame")
  dge <- structure(list(counts = matrix(1:4, nrow = 2)), class = "DGEList")
  expect_equal(classify_object(dge), "dgelist")
  expect_equal(classify_object("hello"), "unknown")
})

# --- classify_table -------------------------------------------------------

test_that("classify_table distinguishes counts from DE results", {
  counts <- data.frame(gene_id = c("g1", "g2"), s1 = c(1, 2), s2 = c(3, 4))
  expect_equal(classify_table(counts), "counts")

  de <- data.frame(gene_name = c("A", "B"), log2FoldChange = c(1, -1),
                   padj = c(0.01, 0.2))
  expect_equal(classify_table(de), "de_results")
})

# --- matrix_to_se ---------------------------------------------------------

test_that("matrix_to_se builds a counts assay with sample colData", {
  m <- matrix(1:6, nrow = 3,
              dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  se <- matrix_to_se(m)
  expect_true(is(se, "SummarizedExperiment"))
  expect_equal(assayNames(se), "counts")
  expect_equal(dim(se), c(3, 2))
  expect_equal(se$Sample, c("s1", "s2"))
})

test_that("matrix_to_se synthesizes names when dimnames are missing", {
  se <- matrix_to_se(matrix(1:6, nrow = 3))
  expect_equal(rownames(se), c("gene_1", "gene_2", "gene_3"))
  expect_equal(colnames(se), c("sample_1", "sample_2"))
})

# --- dgelist_to_se --------------------------------------------------------

test_that("dgelist_to_se carries samples and genes", {
  dge <- structure(list(
    counts = matrix(1:6, nrow = 3,
                    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))),
    samples = data.frame(group = c("A", "B"), row.names = c("s1", "s2")),
    genes = data.frame(gene_name = c("X1", "X2", "X3"))
  ), class = "DGEList")
  se <- dgelist_to_se(dge)
  expect_equal(dim(se), c(3, 2))
  expect_true("group" %in% colnames(colData(se)))
  expect_true("gene_name" %in% colnames(rowData(se)))
})

# --- df_counts_to_se ------------------------------------------------------

test_that("df_counts_to_se splits annotation columns from the counts assay", {
  df <- data.frame(gene_id = c("g1", "g2"), s1 = c(1, 2), s2 = c(3, 4),
                   stringsAsFactors = FALSE)
  se <- df_counts_to_se(df)
  expect_equal(assayNames(se), "counts")
  expect_equal(dim(se), c(2, 2))
  expect_equal(rownames(se), c("g1", "g2"))
  expect_true("gene_id" %in% colnames(rowData(se)))
})

test_that("df_counts_to_se errors when there are no numeric columns", {
  df <- data.frame(a = c("x", "y"), b = c("p", "q"), stringsAsFactors = FALSE)
  expect_error(df_counts_to_se(df), "numeric")
})

# --- df_de_to_se ----------------------------------------------------------

test_that("df_de_to_se keeps DE columns in rowData and has no assay", {
  de <- data.frame(gene_name = c("A", "B"), log2FoldChange = c(1, -2),
                   padj = c(0.01, 0.2), stringsAsFactors = FALSE)
  se <- df_de_to_se(de)
  expect_equal(length(assayNames(se)), 0)
  expect_true(all(c("gene_name", "log2FoldChange", "padj") %in%
                    colnames(rowData(se))))
  expect_equal(nrow(se), 2)
  expect_equal(ncol(se), 0)
})

# --- coerce_to_se ---------------------------------------------------------

test_that("coerce_to_se dispatches and passes SE through unchanged", {
  se <- make_mock_se()
  expect_identical(coerce_to_se(se), se)
  expect_true(is(coerce_to_se(matrix(1:4, nrow = 2)), "SummarizedExperiment"))
  expect_error(coerce_to_se("not convertible"))
})

# --- load_dataset ---------------------------------------------------------

test_that("load_dataset reads a CSV counts matrix", {
  df <- data.frame(gene_id = c("g1", "g2"), s1 = c(1, 2), s2 = c(3, 4),
                   stringsAsFactors = FALSE)
  tmp <- tempfile(fileext = ".csv")
  write.csv(df, tmp, row.names = FALSE)
  se <- load_dataset(tmp)
  expect_true(is(se, "SummarizedExperiment"))
  expect_equal(dim(se), c(2, 2))
  expect_true("counts" %in% assayNames(se))
})

test_that("load_dataset reads a TSV counts matrix", {
  df <- data.frame(gene_id = c("g1", "g2"), s1 = c(1, 2), s2 = c(3, 4),
                   stringsAsFactors = FALSE)
  tmp <- tempfile(fileext = ".tsv")
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  se <- load_dataset(tmp)
  expect_equal(dim(se), c(2, 2))
  expect_true("counts" %in% assayNames(se))
})

test_that("load_dataset reads a DE-results CSV into rowData", {
  de <- data.frame(gene_name = c("A", "B"), log2FoldChange = c(3, -3),
                   padj = c(0.001, 0.002), stringsAsFactors = FALSE)
  tmp <- tempfile(fileext = ".csv")
  write.csv(de, tmp, row.names = FALSE)
  se <- load_dataset(tmp)
  expect_equal(length(assayNames(se)), 0)
  expect_true("log2FoldChange" %in% colnames(rowData(se)))
  expect_equal(nrow(se), 2)
})

test_that("load_dataset round-trips an RDS matrix and SE", {
  m <- matrix(1:6, nrow = 3,
              dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  tmp_m <- tempfile(fileext = ".rds")
  saveRDS(m, tmp_m)
  se_m <- load_dataset(tmp_m)
  expect_true(is(se_m, "SummarizedExperiment"))
  expect_equal(dim(se_m), c(3, 2))

  se_mock <- make_mock_se()
  tmp_se <- tempfile(fileext = ".rds")
  saveRDS(se_mock, tmp_se)
  expect_equal(dim(load_dataset(tmp_se)), dim(se_mock))
})

test_that("load_dataset errors on missing files and unknown types", {
  expect_error(load_dataset(tempfile(fileext = ".rds")), "does not exist")
  bad <- tempfile(fileext = ".xyz")
  writeLines("nonsense", bad)
  expect_error(load_dataset(bad), "Unsupported")
})

test_that("load_dataset reads qs2 files", {
  skip_if_not_installed("qs2")
  se_mock <- make_mock_se()
  tmp <- tempfile(fileext = ".qs2")
  qs2::qs_save(se_mock, tmp)
  expect_equal(dim(load_dataset(tmp)), dim(se_mock))
})

test_that("load_dataset reads qs files", {
  skip_if_not_installed("qs")
  se_mock <- make_mock_se()
  tmp <- tempfile(fileext = ".qs")
  qs::qsave(se_mock, tmp)
  expect_equal(dim(load_dataset(tmp)), dim(se_mock))
})

test_that("load_dataset reads parquet tables", {
  skip_if_not_installed("arrow")
  df <- data.frame(gene_id = c("g1", "g2"), s1 = c(1, 2), s2 = c(3, 4),
                   stringsAsFactors = FALSE)
  tmp <- tempfile(fileext = ".parquet")
  arrow::write_parquet(df, tmp)
  se <- load_dataset(tmp)
  expect_true("counts" %in% assayNames(se))
  expect_equal(dim(se), c(2, 2))
})
