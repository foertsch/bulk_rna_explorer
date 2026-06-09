# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Input autodetection + format conversion: `load_dataset()` reads `.rds` /
  `.qs2` / `.qs` containers (SummarizedExperiment, DESeqDataSet, DGEList,
  matrix, data.frame) and `.csv` / `.tsv` / `.txt` / `.parquet` / `.feather`
  tables (counts matrix or DE-results, with delimiter and table-type sniffing)
  and converts each to a `SummarizedExperiment`. Optional readers (`qs2`,
  `qs`, `arrow`) are required lazily
- `scripts/convert_to_se.R`: CLI to batch-convert any supported file to `.rds`
- `docs/AGENT_GUIDE.md`: agent-facing guide to the data contract, supported
  formats, and conversion API
- `tests/testthat/test-loaders.R`: unit tests for detection, classification,
  coercion, and `load_dataset()` round-trips (optional formats skipped when
  the package is absent)
- CI: `lintr` job with project-specific `.lintr` config (line length 120,
  object-usage / object-name / pipe-consistency / indentation / brace
  linters disabled to suit the Shiny + Bioconductor codebase)
- `CHANGELOG.md` following Keep a Changelog format
- README: macOS / Linux launch instructions

### Changed
- App now accepts non-RDS uploads and auto-loads supported formats from
  `data/`; the lazy loader routes through `load_dataset()` instead of
  `readRDS()`
- CI: extended `testthat` job to a matrix covering Ubuntu and macOS
  (previously Ubuntu only). Windows Shiny construction smoke remains a
  separate job.
- README: replaced two em-dashes with comma and period for consistency

## [0.1.0] - 2026-04-14

Initial release.

### Added
- R Shiny app for interactive exploration of bulk RNA-seq DE results on
  `SummarizedExperiment` objects
- Multi-dataset loading and contrast switching from `.rds` files
- Auto-detected column mapping for FGCZ Sushi (`log2Ratio`, `fdr`) and
  standard Bioconductor (`log2FoldChange`, `padj`) conventions
- Expression barplots: mean, SD, jittered points, optional second grouping
  for dodged display
- Volcano plots with adjustable thresholds, click-to-plot, top-gene labels
- PCA plots with selectable axes, optional 95% ellipses, second-group
  shape encoding
- Filterable gene-search table with click-to-plot
- Per-group dynamic color pickers
- PNG and PDF download for all plots
- Windows desktop launcher (`launch_windows.bat`) with custom icon
- 20 unit tests (testthat) on pure helper functions
- GitHub Actions CI: `testthat` on Ubuntu, Shiny construction smoke on
  Windows
- Two toy datasets shipped for demo and CI use

### Fixed
- CI Windows smoke replaced port-polling launch with a construction-only
  check after Start-Process buffering and Bioconductor dependency setup
  caused intermittent timeouts (see PR #1, 2026-05-20)
