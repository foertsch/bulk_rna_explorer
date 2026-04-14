# testthat auto-sources files named helper-*.R before running tests.
# Builds a tiny synthetic SummarizedExperiment mimicking FGCZ DESeq2 output.

make_mock_se <- function(n_genes = 20, n_samples = 12, seed = 42) {
  set.seed(seed)

  sample_ids <- paste0("sample_", seq_len(n_samples))
  gene_ids <- paste0("ENSG", sprintf("%05d", seq_len(n_genes)))

  # counts + xNorm assays
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  xNorm <- counts * runif(length(counts), 0.8, 1.2)

  # colData with Condition and Dox factors (matches v1 app's expected schema)
  conditions <- rep(c("cond_A", "cond_B", "cond_C"),
                    length.out = n_samples)
  dox <- rep(c("ctrl", "Dox"), length.out = n_samples)

  col_data <- DataFrame(
    Sample = sample_ids,
    Condition = factor(conditions, levels = c("cond_A", "cond_B", "cond_C")),
    Dox = factor(dox, levels = c("ctrl", "Dox")),
    numeric_cov = rnorm(n_samples)
  )

  # rowData with gene_name + DE statistics
  row_data <- DataFrame(
    gene_id = gene_ids,
    gene_name = paste0("SYM", seq_len(n_genes)),
    log2Ratio = rnorm(n_genes, sd = 2),
    fdr = runif(n_genes, 0, 1),
    pValue = runif(n_genes, 0, 1)
  )
  # Guarantee a few hits in each direction for volcano tests
  row_data$log2Ratio[1:3] <- c(3, -3, 2.5)
  row_data$fdr[1:3] <- c(0.001, 0.001, 0.01)

  SummarizedExperiment(
    assays = list(counts = counts, xNorm = xNorm),
    rowData = row_data,
    colData = col_data
  )
}
