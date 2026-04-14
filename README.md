# bulk_rna_explorer

Interactive Shiny app for exploring bulk RNA-seq DESeq2/EdgeR results. Supports multiple contrasts, dynamic column mapping, expression barplots, volcano plots, and PCA plots.

## Features

- **Multi-dataset support**: Load multiple RDS files on demand, switch between contrasts
- **Dynamic column mapping**: Auto-detects assays, grouping columns, DE statistics, and gene symbols — works with any `SummarizedExperiment`
- **Expression barplots**: Mean + SD + jittered points, grouped by any `colData` column, optional second grouping
- **Volcano plots**: log2FC vs -log10(FDR) with adjustable thresholds, click-to-plot, top-gene labels
- **PCA plots**: Top-N variable genes, selectable PC axes, colored by group, shaped by optional second group
- **Gene search table**: Filterable DE statistics with click-to-plot
- **Dynamic color pickers**: One per detected group level, auto-generated palette
- **Download**: PNG and PDF for all plots

## Quick start

```r
shiny::runApp("/path/to/bulk_rna_explorer")
```

Upload `.rds` files via the browser to load contrasts on demand.

## Data format

Input `.rds` files must contain a **`SummarizedExperiment`** object. Column mapping is done in the UI, but these defaults are auto-detected:

- **Assays**: Prefers `xNorm` if present, otherwise the first available
- **Group column**: Prefers `Condition` (any categorical `colData` column works; `[Factor]` suffixes are tolerated)
- **Second group**: Prefers `Dox` (optional, for dodged barplots)
- **Gene symbols**: Detects `gene_name`, `gene_symbol`, `Symbol`, `SYMBOL`, `external_gene_name`, `hgnc_symbol`
- **DE columns**: `log2Ratio` / `log2FoldChange` for FC; `fdr` / `padj` for adjusted p-value; `pValue` / `pvalue` for raw p-value

Compatible with FGCZ Sushi DESeq2/EdgeR output and standard Bioconductor DE results.

## Requirements

Auto-installed on first run.

- **CRAN:** shiny, shinythemes, ggplot2, dplyr, colourpicker, shinyWidgets, DT, ggrepel
- **Bioconductor:** SummarizedExperiment

## Running the app

```r
# R console
shiny::runApp("/path/to/bulk_rna_explorer")
```

Or open `app.R` in RStudio and click **Run App**.

## Development

```bash
Rscript tests/testthat.R
```

## License

MIT — see [LICENSE](LICENSE).
