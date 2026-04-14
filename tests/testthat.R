# Runner: source with Rscript from the app root directory.
#   cd gene_expression_app && Rscript tests/testthat.R
library(testthat)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(ggrepel)

source(file.path("R", "helpers.R"))
test_dir(file.path("tests", "testthat"), reporter = "summary")
