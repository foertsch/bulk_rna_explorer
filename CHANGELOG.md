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
- README: "Use with an AI coding agent" section with a copy-paste prompt,
  after validating that a cold agent can install and launch the app from the
  docs on a clean machine
- CI: `cleanroom-linux` job that installs the full CRAN + Bioconductor stack
  from zero in a bare `rocker/r-ver` container, launches via `run.R`, and
  probes the running app over HTTP, the only job that exercises the
  self-installer cold-start path end to end
- README / AGENT_GUIDE: Linux system-library prerequisites for from-source
  installs (`libuv1-dev` and the curl/ssl/xml2 + ragg/textshaping `-dev` set)
- Expression / PCA: "Assay is already log-scale" toggle that skips the
  `log2(x + 1)` transform for assays already on a log scale (vst, rlog,
  logCPM), avoiding a double transform. The barplot y-axis label now reflects
  the selected assay and whether it was log-transformed

### Changed
- App now accepts non-RDS uploads and auto-loads supported formats from
  `data/`; the lazy loader routes through `load_dataset()` instead of
  `readRDS()`
- `data/` auto-discovery now also picks up compressed delimited files
  (`.gz` / `.bz2` / `.xz`), matching the upload accept list and
  `detect_file_type()`
- PCA and Expression tabs now show an explanatory message for datasets with no
  expression assay (e.g. a DE-results-only table) instead of rendering blank
- CI: extended `testthat` job to a matrix covering Ubuntu and macOS
  (previously Ubuntu only). Windows Shiny construction smoke remains a
  separate job.
- README: replaced two em-dashes with comma and period for consistency
- README: document installing dev dependencies (`testthat`, `lintr`) before
  running the test/lint commands; they are not auto-installed like the app's
  runtime packages

### Fixed
- Plot downloads (PNG/PDF for the barplot, volcano, and PCA) produced no file:
  Shiny hands `downloadHandler`'s `content()` a temp path whose extension is
  appended with no separator (`.../fileXXXXpng`), so `ggsave()` could not infer
  the format and errored. The handlers now pass `device` explicitly
- PCA crashed with "subscript out of bounds" when a selected PC axis (e.g. PC4)
  exceeded the number of components a small dataset yields. The axis dropdowns
  are now capped to the available PCs and `plot_pca()` clamps out-of-range axes
- Removed ineffective `check.names = FALSE` arguments from `as.data.frame()`
  calls on `SummarizedExperiment` colData/rowData; the `DataFrame` method
  silently ignored them while emitting an "arguments in '...' ignored" warning
  on every render
- Uploaded or discovered datasets that share a basename (e.g. `results.csv` and
  `results.rds`) no longer silently overwrite each other; labels are
  de-duplicated
- A transient error when switching datasets (the previous gene briefly not
  present in the new dataset) is now guarded instead of flashing a raw error
- Cold-start launch on a machine without `shiny`: `shiny::runApp(...)` and
  `launch_windows.bat` both required `shiny` to be installed already, because
  the in-app auto-installer only runs once `app.R` is sourced. Added a `run.R`
  launcher that bootstraps `shiny`, then hands off to the app's installer for
  the rest; README and the Windows launcher now invoke `Rscript run.R`
- Cold-start install on a bare Linux machine: `fs` (>= 2.1.0) links system
  `libuv` and fails to compile without `libuv1-dev`, which cascades through
  `sass` -> `bslib` -> `shiny` so the app never launches. Found via a bare
  `rocker/r-ver` Docker run. Documented the required Linux `-dev` libraries
  and pinned them in the new `cleanroom-linux` CI job

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
