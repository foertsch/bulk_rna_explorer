#!/usr/bin/env Rscript
# Convert a supported input file into an .rds SummarizedExperiment.
#
# Usage:
#   Rscript scripts/convert_to_se.R <input> [output.rds]
#
# Supported inputs (autodetected):
#   .rds .qs2 .qs            -> SummarizedExperiment, DESeqDataSet, DGEList,
#                               matrix, or data.frame held inside the container
#   .csv .tsv .txt(.gz) .tab -> counts matrix or DE-results table (delimiter
#                               and table type are sniffed automatically)
#   .parquet .feather        -> counts matrix or DE-results table
#
# Run from the repo root so R/helpers.R resolves. If no output path is given,
# the input path is reused with an .rds extension.

suppressPackageStartupMessages(library(SummarizedExperiment))
source(file.path("R", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/convert_to_se.R <input> [output.rds]",
       call. = FALSE)
}

input <- args[[1]]
output <- if (length(args) >= 2) {
  args[[2]]
} else {
  paste0(tools::file_path_sans_ext(input), ".rds")
}

se <- load_dataset(input)
saveRDS(se, output)

assays_txt <- if (length(assayNames(se)) > 0) {
  paste(assayNames(se), collapse = ", ")
} else {
  "(none)"
}

message(sprintf("Converted %s -> %s", input, output))
message(sprintf("  class : %s", paste(class(se), collapse = "/")))
message(sprintf("  dim   : %d genes x %d samples", nrow(se), ncol(se)))
message(sprintf("  assays: %s", assays_txt))
