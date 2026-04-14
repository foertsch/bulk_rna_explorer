# Generates two small toy contrasts as .rds files for CI + demo use.
# Uses FGCZ-style colData naming ("Condition [Factor]", "Dox [Factor]") so
# the app's sanitize_name() path is exercised. Run from the app root:
#   Rscript scripts/make_toy_data.R

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
})

make_toy_se <- function(seed, up_idx, down_idx,
                        n_genes = 200, n_samples = 12) {
  set.seed(seed)

  sample_ids <- paste0("s", sprintf("%02d", seq_len(n_samples)))
  gene_ids <- paste0("ENSG", sprintf("%05d", seq_len(n_genes)))

  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 80),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  xNorm <- counts * runif(length(counts), 0.8, 1.2)

  conditions <- rep(c("cond_A", "cond_B", "cond_C"), length.out = n_samples)
  dox <- rep(c("ctrl", "Dox"), length.out = n_samples)

  cd <- DataFrame(Sample = sample_ids, check.names = FALSE)
  cd$`Condition [Factor]` <- factor(conditions,
                                     levels = c("cond_A", "cond_B", "cond_C"))
  cd$`Dox [Factor]` <- factor(dox, levels = c("ctrl", "Dox"))
  cd$numeric_cov <- rnorm(n_samples)

  rd <- DataFrame(
    gene_id = gene_ids,
    gene_name = paste0("SYM", seq_len(n_genes)),
    log2Ratio = rnorm(n_genes, sd = 1.2),
    fdr = runif(n_genes, 0.05, 1),
    pValue = runif(n_genes, 0.05, 1)
  )
  # Seed deterministic UP / DOWN hits so each contrast looks different
  rd$log2Ratio[up_idx] <- runif(length(up_idx), 1.5, 4)
  rd$fdr[up_idx] <- runif(length(up_idx), 1e-5, 0.01)
  rd$log2Ratio[down_idx] <- runif(length(down_idx), -4, -1.5)
  rd$fdr[down_idx] <- runif(length(down_idx), 1e-5, 0.01)
  rd$pValue <- rd$fdr * runif(n_genes, 0.5, 1)

  SummarizedExperiment(
    assays = list(counts = counts, xNorm = xNorm),
    rowData = rd,
    colData = cd
  )
}

out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

se_A <- make_toy_se(seed = 1, up_idx = 1:15, down_idx = 16:30)
se_B <- make_toy_se(seed = 2, up_idx = 40:60, down_idx = 61:75)

saveRDS(se_A, file.path(out_dir, "toy_contrast_A.rds"))
saveRDS(se_B, file.path(out_dir, "toy_contrast_B.rds"))

for (f in c("toy_contrast_A.rds", "toy_contrast_B.rds")) {
  size_kb <- round(file.info(file.path(out_dir, f))$size / 1024, 1)
  message(sprintf("wrote data/%s (%s KB)", f, size_kb))
}
