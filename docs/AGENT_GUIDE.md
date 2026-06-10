# Agent guide: using bulk_rna_explorer with your data

This file is written for AI coding agents (and humans) who need to get a user's
bulk RNA-seq data into the Gene Expression Explorer Shiny app, or convert it
programmatically. It documents the **data contract**, the **supported input
formats**, and the **conversion API**.

> Scope: this is a *usage* guide for working with the app and its data. It is
> intentionally separate from any developer-facing `AGENTS.md` / `CLAUDE.md`
> (project conventions, contribution rules). If you are modifying the codebase,
> read those instead.

## What the app does

`app.R` is a single-file Shiny app that explores bulk RNA-seq differential
expression results stored as a [`SummarizedExperiment`][se]. It renders four
tabs: PCA, Volcano, Gene Search (a filterable DE table), and an Expression
barplot. Pure plotting/loading logic lives in `R/helpers.R` and is unit-tested
in `tests/testthat/`.

[se]: https://bioconductor.org/packages/SummarizedExperiment/

## Launch

```bash
# From the repo root. Installs any missing packages on first run, then serves.
Rscript run.R
```

`run.R` bootstraps `shiny` and hands off to the app's installer for the rest, so
it works from a clean machine where `shiny::runApp(...)` would not (that needs
`shiny` already installed). Set `LAUNCH_BROWSER=false` for headless / CI use.
From an R session that already has `shiny`, `shiny::runApp("/path/to/bulk_rna_explorer")`
also works.

**On Linux, install system libraries first.** Packages compile from source there
(macOS/Windows get binaries), and a few link system libraries a minimal box
lacks, most importantly `libuv` for `fs` (its absence cascades to `sass` ->
`bslib` -> `shiny`). On Debian / Ubuntu:

```bash
sudo apt-get install -y libuv1-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
  libpng-dev libjpeg-dev libtiff5-dev libfontconfig1-dev libfreetype6-dev \
  libharfbuzz-dev libfribidi-dev
```

Files in `data/` auto-load (paths only; objects are read lazily when selected).
Users can also upload files from the panel at the top of the app.

## The data contract

Internally, every dataset is a `SummarizedExperiment` (SE). Which tabs work
depends on what the SE contains:

| Tab            | Requires                                                        |
|----------------|----------------------------------------------------------------|
| PCA            | At least one numeric `assay` (genes x samples)                  |
| Expression     | An `assay`, plus a categorical `colData` column to group by     |
| Volcano        | `rowData` columns for log2 fold-change **and** FDR              |
| Gene Search    | `rowData` (DE columns optional; shown/filterable when present)  |

Column mapping is done in the UI, but sensible defaults are auto-detected
(`R/helpers.R`):

- **Assay** (`detect`: prefers `xNorm`, else the first assay)
- **Group column** — any categorical `colData` column; prefers `Condition`.
  Trailing `[Factor]` suffixes are tolerated (`sanitize_name`).
- **Gene symbols** (`detect_gene_col`): `gene_name`, `gene_symbol`, `Symbol`,
  `SYMBOL`, `external_gene_name`, `hgnc_symbol`
- **DE columns** (`detect_de_columns`):
  - log2FC: `log2Ratio`, `log2FoldChange`, `logFC`, `log2ratio`
  - FDR: `fdr`, `padj`, `FDR`, `adj.P.Val`, `p_adjusted`
  - p-value: `pValue`, `pvalue`, `PValue`, `P.Value`

These cover FGCZ Sushi DESeq2/EdgeR output and standard Bioconductor DE results.
If your data uses different names, rename the columns or pick them in the UI.

The Expression and PCA tabs `log2(x + 1)`-transform the chosen assay by default.
If the assay is already on a log scale (vst, rlog, logCPM), tick **"Assay is
already log-scale"** in the Column Mapping panel (or pass `log_transform = FALSE`
to `create_barplot()` / `compute_pca()`) to skip it.

## Supported input formats (autodetection)

`load_dataset(path)` in `R/helpers.R` is the single entry point. It detects the
file type by extension, reads it, and returns an SE (or errors with an
actionable message). Mapping:

| Extension(s)                       | Read as           | Coerced to SE by                |
|------------------------------------|-------------------|---------------------------------|
| `.rds`                             | R object          | see object rules below          |
| `.qs2`                             | R object (`qs2`)  | see object rules below          |
| `.qs`                              | R object (`qs`)   | see object rules below          |
| `.csv` `.tsv` `.txt` `.tab` (`.gz`)| table (delim sniff)| table rules below              |
| `.parquet` `.feather`              | table (`arrow`)   | table rules below               |

**Object rules** (for `.rds` / `.qs2` / `.qs` contents):

- `SummarizedExperiment` / `RangedSummarizedExperiment` / `DESeqDataSet`
  (all inherit from SE) -> passed through unchanged.
- `DGEList` (edgeR) -> `$counts` becomes the `counts` assay, `$samples` ->
  `colData`, `$genes` -> `rowData`. (edgeR need not be installed.)
- `matrix` -> single `counts` assay; sample/gene names synthesized if absent.
- `data.frame` -> classified as a counts matrix or DE table (see below).

**Table rules** (delimited text, parquet, feather, or a `data.frame`):

- If the columns include recognizable DE fields (log2FC / FDR / p-value), the
  table is treated as **DE results**: rows are genes, the whole table becomes
  `rowData`, and there is **no assay**. Volcano and Gene Search work; PCA and
  Expression do not (a DE table carries no per-sample expression).
- Otherwise a mostly-numeric table is treated as a **counts matrix**:
  non-numeric columns (e.g. `gene_id`, `gene_name`) become `rowData` (the first
  supplies gene IDs); numeric columns become the `counts` assay; `colData` gets
  a `Sample` column from the column names.

> Note on grouping: a bare counts matrix has no sample metadata, so grouped
> plots (Expression, colored PCA) have nothing to group by. To get grouping,
> start from an SE / DESeqDataSet / DGEList that already carries `colData`, or
> add a categorical `colData` column before saving.

Optional readers (`qs2`, `qs`, `arrow`) are required lazily and only when a file
of that type is opened. If one is missing, the error tells you what to install.

## Conversion API

### In R

```r
source("R/helpers.R")          # defines load_dataset() and friends
se <- load_dataset("counts.csv")
se <- load_dataset("results.parquet")
se <- load_dataset("dataset.qs2")
```

Lower-level helpers (all pure, all tested) are available if you need them:
`detect_file_type()`, `sniff_delim()`, `read_any()`, `classify_object()`,
`classify_table()`, `coerce_to_se()`.

### Command line (batch conversion to .rds)

```bash
# Run from the repo root so R/helpers.R resolves.
Rscript scripts/convert_to_se.R <input> [output.rds]
```

If no output path is given, the input path is reused with an `.rds` extension.
The script prints the detected class, dimensions, and assay names. Drop the
resulting `.rds` into `data/` and it will auto-load on the next launch.

## Repo layout

| Path                          | What                                              |
|-------------------------------|---------------------------------------------------|
| `app.R`                       | Shiny UI + server (single file)                   |
| `R/helpers.R`                 | Pure helpers: detection, loaders, plot builders   |
| `scripts/convert_to_se.R`     | CLI: any supported file -> `.rds` SE              |
| `scripts/make_toy_data.R`     | Regenerates the bundled toy datasets              |
| `data/`                       | Auto-loaded datasets (toy contrasts ship here)    |
| `tests/testthat/`             | Unit tests (`test-helpers.R`, `test-loaders.R`)   |

## Running the tests

```bash
Rscript tests/testthat.R
```

Tests for optional formats (`qs2`, `qs`, `parquet`/`feather`) skip cleanly when
the relevant package is not installed.
