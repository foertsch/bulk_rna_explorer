# ============================================================
# Gene Expression Explorer v2 - Shiny App
# Author: Arion Foertsch @ FGCZ
# Supports multiple DESeq2 datasets, dynamic conditions,
# volcano plots, and PCA plots.
# ============================================================

# --- Dependency Installation ---
# Only install what's actually missing. In particular, do NOT pull BiocManager
# from CRAN when all bioc_packages are already present (this matters for CI:
# pak installs SummarizedExperiment directly but not BiocManager, and the
# unconditional CRAN call delayed the Shiny Windows smoke past its timeout).
install_if_missing <- function(packages, bioc_packages = NULL) {
  missing_cran <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  for (pkg in missing_cran) {
    message("Installing ", pkg, "...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!is.null(bioc_packages)) {
    missing_bioc <- bioc_packages[
      !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_bioc) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      for (pkg in missing_bioc) {
        message("Installing ", pkg, " from Bioconductor...")
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      }
    }
  }
}

install_if_missing(
  packages = c("shiny", "ggplot2", "dplyr", "colourpicker",
               "shinyWidgets", "DT", "shinythemes", "ggrepel"),
  bioc_packages = c("SummarizedExperiment")
)

# --- Load Libraries ---
library(shiny)
library(shinythemes)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(colourpicker)
library(shinyWidgets)
library(DT)
library(ggrepel)

# --- Load helper + plot functions ---
# shiny::runApp() sets getwd() to the app directory, so R/helpers.R resolves.
# For other entry points (e.g. sourcing app.R directly), fall back to file.path
# relative to the running script.
source(file.path("R", "helpers.R"))

# --- UI ---

# Try to load FGCZ header
fgcz_header <- tryCatch({
  header_path <- system.file("templates/fgcz_header.html",
                             package = "ezRun", lib.loc = .libPaths())
  if (file.exists(header_path)) includeHTML(header_path) else NULL
}, error = function(e) NULL)

ui <- fluidPage(
  theme = shinytheme("paper"),

  if (!is.null(fgcz_header)) fgcz_header,

  tags$div(
    style = "margin-bottom: 20px; padding: 15px; background-color: #f8f9fa; border-bottom: 2px solid #dee2e6;",
    tags$h2("Gene Expression Explorer", style = "margin: 0 0 5px 0; color: #333;"),
    tags$p(
      style = "margin: 0; color: #666; font-size: 14px;",
      "Author: Arion Foertsch @ Functional Genomics Center Zurich"
    )
  ),

  # --- Data Loading Panel (always visible) ---
  wellPanel(
    fluidRow(
      column(5,
        fileInput("data_files", "Load RDS file(s):",
                  accept = c(".rds", ".RDS"), multiple = TRUE)
      ),
      column(5,
        selectInput("active_dataset", "Active dataset:", choices = NULL)
      ),
      column(2,
        br(),
        actionButton("remove_dataset", "Remove", class = "btn-sm btn-danger",
                     style = "margin-top: 8px;")
      )
    ),
    helpText("Files in data/ folder auto-load; uploads are added on demand. Datasets are read lazily when selected.")
  ),

  # --- Main App (shown when data is loaded) ---
  conditionalPanel(
    condition = "output.data_loaded",

    # --- Column Mapping Panel ---
    wellPanel(
      style = "background-color: #f0f4f7;",
      h4("Column Mapping"),
      fluidRow(
        column(3, selectInput("assay_name", "Assay:", choices = NULL)),
        column(3, selectInput("group_col", "Group column:", choices = NULL)),
        column(3, selectInput("second_group", "2nd group (optional):", choices = NULL)),
        column(3, selectInput("gene_symbol_col", "Gene symbol column:", choices = NULL))
      ),
      fluidRow(
        column(4, selectInput("fc_col", "log2FC column:", choices = NULL)),
        column(4, selectInput("fdr_col", "FDR column:", choices = NULL)),
        column(4, selectInput("pval_col", "p-value column:", choices = NULL))
      )
    ),

    sidebarLayout(
      sidebarPanel(

        # --- Expression Plot controls ---
        conditionalPanel(
          condition = "input.main_tabs == 'Expression Plot'",
          selectizeInput("gene", "Quick Gene Select:", choices = NULL,
                         options = list(placeholder = "Type gene symbol...",
                                        maxOptions = 100)),
          hr(),
          h4("Filter Groups"),
          uiOutput("group_checkboxes_ui"),
          actionButton("select_all", "Select All", class = "btn-sm"),
          actionButton("select_none", "Select None", class = "btn-sm"),
          hr(),
          h4("Group Colors"),
          uiOutput("color_pickers_ui"),
          br(),
          sliderInput("ctrl_alpha", "2nd group first-level transparency:",
                      min = 0.2, max = 1, value = 0.5, step = 0.1),
          actionButton("reset_colors", "Reset Colors", class = "btn-sm btn-secondary"),
          hr(),
          h4("Download Plot"),
          fluidRow(
            column(6, downloadButton("download_barplot_png", "PNG", class = "btn-block")),
            column(6, downloadButton("download_barplot_pdf", "PDF", class = "btn-block"))
          )
        ),

        # --- Volcano Plot controls ---
        conditionalPanel(
          condition = "input.main_tabs == 'Volcano Plot'",
          h4("Volcano Options"),
          sliderInput("volcano_fc", "log2FC threshold:",
                      min = 0, max = 5, value = 1, step = 0.25),
          numericInput("volcano_fdr", "FDR threshold:",
                       value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput("volcano_nlabels", "Top N gene labels:",
                       value = 10, min = 0, max = 50, step = 1),
          hr(),
          h4("Download Plot"),
          fluidRow(
            column(6, downloadButton("download_volcano_png", "PNG", class = "btn-block")),
            column(6, downloadButton("download_volcano_pdf", "PDF", class = "btn-block"))
          )
        ),

        # --- PCA Plot controls ---
        conditionalPanel(
          condition = "input.main_tabs == 'PCA Plot'",
          h4("PCA Options"),
          numericInput("pca_ntop", "Top N variable genes:",
                       value = 500, min = 10, max = 10000, step = 50),
          fluidRow(
            column(6, selectInput("pca_pc_x", "X axis:", choices = paste0("PC", 1:4),
                                  selected = "PC1")),
            column(6, selectInput("pca_pc_y", "Y axis:", choices = paste0("PC", 1:4),
                                  selected = "PC2"))
          ),
          checkboxInput("pca_ellipses", "Show group ellipses (95%)",
                        value = FALSE),
          hr(),
          h4("Download Plot"),
          fluidRow(
            column(6, downloadButton("download_pca_png", "PNG", class = "btn-block")),
            column(6, downloadButton("download_pca_pdf", "PDF", class = "btn-block"))
          )
        ),

        # --- Gene Search controls ---
        conditionalPanel(
          condition = "input.main_tabs == 'Gene Search'",
          h4("Filter DE Results"),
          numericInput("fdr_threshold", "Max FDR:", value = 0.05,
                       min = 0, max = 1, step = 0.01),
          numericInput("log2fc_threshold", "Min |log2FC|:", value = 0,
                       min = 0, max = 10, step = 0.1),
          selectInput("de_filter", "Show:",
                      choices = c("All genes" = "all", "Significant only" = "sig",
                                  "Upregulated" = "up", "Downregulated" = "down"),
                      selected = "all")
        ),

        width = 3
      ),

      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          type = "tabs",

          tabPanel(
            "Expression Plot",
            br(),
            plotOutput("barplot", height = "500px"),
            hr(),
            h4("Selected Gene Information"),
            DTOutput("gene_info")
          ),

          tabPanel(
            "Gene Search",
            br(),
            helpText("Click a row to select that gene for plotting."),
            DTOutput("gene_search_table")
          ),

          tabPanel(
            "Volcano Plot",
            br(),
            helpText("Click a point to jump to that gene's expression plot."),
            plotOutput("volcano", height = "600px",
                       click = "volcano_click")
          ),

          tabPanel(
            "PCA Plot",
            br(),
            plotOutput("pca", height = "600px")
          )
        ),
        width = 9
      )
    )
  )
)

# --- Server ---
server <- function(input, output, session) {

  # ---- Dataset registry: paths + cache ----
  # dataset_paths: named list of file paths (NOT yet loaded)
  # dataset_cache: named list of loaded SE objects (lazy)
  dataset_paths <- reactiveVal(list())
  dataset_cache <- reactiveVal(list())

  # Discover RDS files in data/ folder on startup (paths only, no loading)
  observe({
    data_dirs <- c(file.path(getwd(), "data"), "data")
    for (data_dir in data_dirs) {
      if (dir.exists(data_dir)) {
        rds_files <- list.files(data_dir, pattern = "\\.rds$",
                                full.names = TRUE, ignore.case = TRUE)
        if (length(rds_files) > 0) {
          paths <- setNames(as.list(rds_files),
                            tools::file_path_sans_ext(basename(rds_files)))
          dataset_paths(paths)
          updateSelectInput(session, "active_dataset", choices = names(paths),
                            selected = names(paths)[1])
          break
        }
      }
    }
  }, priority = 100)

  # Handle file uploads (register paths; validate on first read)
  observeEvent(input$data_files, {
    req(input$data_files)
    paths <- dataset_paths()
    for (i in seq_len(nrow(input$data_files))) {
      file_info <- input$data_files[i, ]
      label <- tools::file_path_sans_ext(file_info$name)
      paths[[label]] <- file_info$datapath
    }
    dataset_paths(paths)
    updateSelectInput(session, "active_dataset", choices = names(paths),
                      selected = names(paths)[length(paths)])
  })

  # Remove dataset
  observeEvent(input$remove_dataset, {
    req(input$active_dataset)
    paths <- dataset_paths()
    cache <- dataset_cache()
    paths[[input$active_dataset]] <- NULL
    cache[[input$active_dataset]] <- NULL
    dataset_paths(paths)
    dataset_cache(cache)
    new_sel <- if (length(paths) > 0) names(paths)[1] else character(0)
    updateSelectInput(session, "active_dataset",
                      choices = if (length(paths) > 0) names(paths) else character(0),
                      selected = new_sel)
  })

  output$data_loaded <- reactive({
    length(dataset_paths()) > 0 && !is.null(input$active_dataset) &&
      input$active_dataset %in% names(dataset_paths())
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  # Active SE object (lazy-loaded with cache + validation)
  se <- reactive({
    req(input$active_dataset)
    paths <- dataset_paths()
    cache <- dataset_cache()
    req(input$active_dataset %in% names(paths))

    name <- input$active_dataset
    if (!is.null(cache[[name]])) return(cache[[name]])

    path <- paths[[name]]
    data <- withProgress(
      message = paste("Loading", name, "..."), value = 0.5,
      {
        tryCatch(readRDS(path), error = function(e) e)
      }
    )

    if (inherits(data, "error")) {
      showNotification(paste("Failed to read", name, ":", conditionMessage(data)),
                       type = "error", duration = 10)
      req(FALSE)
    }
    if (!is(data, "SummarizedExperiment")) {
      showNotification(
        paste0("'", name, "' is not a SummarizedExperiment (got: ",
               paste(class(data), collapse = ", "), ")"),
        type = "error", duration = 10
      )
      # Remove invalid entry from paths
      paths[[name]] <- NULL
      dataset_paths(paths)
      updateSelectInput(session, "active_dataset",
                        choices = if (length(paths) > 0) names(paths) else character(0),
                        selected = if (length(paths) > 0) names(paths)[1] else character(0))
      req(FALSE)
    }

    cache[[name]] <- data
    dataset_cache(cache)
    data
  })

  # ---- Column mapping: update choices when dataset changes ----
  observeEvent(se(), {
    se_obj <- se()
    rd_cols <- colnames(as.data.frame(rowData(se_obj), check.names = FALSE))
    cat_cols <- get_categorical_cols(se_obj)
    assay_choices <- assayNames(se_obj)

    # Assay
    assay_sel <- if ("xNorm" %in% assay_choices) "xNorm" else assay_choices[1]
    updateSelectInput(session, "assay_name", choices = assay_choices,
                      selected = assay_sel)

    # Group column
    group_sel <- if ("Condition" %in% cat_cols) "Condition" else cat_cols[1]
    updateSelectInput(session, "group_col", choices = cat_cols,
                      selected = group_sel)

    # Second group
    second_choices <- c("None", cat_cols)
    second_sel <- if ("Dox" %in% cat_cols) "Dox" else "None"
    updateSelectInput(session, "second_group", choices = second_choices,
                      selected = second_sel)

    # Gene symbol column
    gene_col <- detect_gene_col(rd_cols)
    gene_choices <- c("(rownames)", rd_cols)
    gene_sel <- if (!is.null(gene_col)) gene_col else "(rownames)"
    updateSelectInput(session, "gene_symbol_col", choices = gene_choices,
                      selected = gene_sel)

    # DE columns
    de <- detect_de_columns(rd_cols)
    rd_choices_fc <- c("(none)", rd_cols)
    rd_choices_fdr <- c("(none)", rd_cols)
    rd_choices_pval <- c("(none)", rd_cols)

    updateSelectInput(session, "fc_col", choices = rd_choices_fc,
                      selected = if (!is.na(de$fc)) de$fc else "(none)")
    updateSelectInput(session, "fdr_col", choices = rd_choices_fdr,
                      selected = if (!is.na(de$fdr)) de$fdr else "(none)")
    updateSelectInput(session, "pval_col", choices = rd_choices_pval,
                      selected = if (!is.na(de$pval)) de$pval else "(none)")
  })

  # ---- Dynamic group levels ----
  group_levels <- reactive({
    req(se(), input$group_col)
    md <- as.data.frame(colData(se()), check.names = FALSE)
    colnames(md) <- sanitize_name(colnames(md))
    col <- sanitize_name(input$group_col)
    req(col %in% colnames(md))
    unique(as.character(md[[col]]))
  })

  # ---- Dynamic UI: group checkboxes ----
  output$group_checkboxes_ui <- renderUI({
    lvls <- group_levels()
    checkboxGroupInput("conditions", NULL, choices = lvls, selected = lvls)
  })

  # Select all / none
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(session, "conditions", selected = group_levels())
  })
  observeEvent(input$select_none, {
    updateCheckboxGroupInput(session, "conditions", selected = character(0))
  })

  # ---- Dynamic UI: color pickers ----
  output$color_pickers_ui <- renderUI({
    lvls <- group_levels()
    palette <- generate_palette(length(lvls))
    tagList(
      lapply(seq_along(lvls), function(i) {
        colourInput(
          inputId = paste0("color_", i),
          label = lvls[i],
          value = palette[i],
          showColour = "both"
        )
      })
    )
  })

  # Reset colors
  observeEvent(input$reset_colors, {
    lvls <- group_levels()
    palette <- generate_palette(length(lvls))
    for (i in seq_along(lvls)) {
      updateColourInput(session, paste0("color_", i), value = palette[i])
    }
    updateSliderInput(session, "ctrl_alpha", value = 0.5)
  })

  # Reactive: current color vector
  current_colors <- reactive({
    lvls <- group_levels()
    colors <- sapply(seq_along(lvls), function(i) {
      val <- input[[paste0("color_", i)]]
      if (is.null(val)) generate_palette(length(lvls))[i] else val
    })
    setNames(colors, lvls)
  })

  # ---- Gene symbol helper ----
  gene_sym_col <- reactive({
    if (is.null(input$gene_symbol_col) || input$gene_symbol_col == "(rownames)") {
      NULL
    } else {
      input$gene_symbol_col
    }
  })

  # ---- Populate gene choices ----
  observe({
    req(se())
    rd <- as.data.frame(rowData(se()), check.names = FALSE)
    gene_ids <- rownames(rd)

    sym_col <- gene_sym_col()
    if (!is.null(sym_col) && sym_col %in% colnames(rd)) {
      gene_symbols <- rd[[sym_col]]
      gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]
    } else {
      gene_symbols <- gene_ids
    }
    choices <- setNames(gene_ids, gene_symbols)
    choices <- choices[order(names(choices))]
    updateSelectizeInput(session, "gene", choices = choices, server = TRUE)
  })

  # ---- Expression barplot ----
  current_barplot <- reactive({
    req(input$gene, input$assay_name, input$group_col)
    req(length(input$conditions) > 0)

    second <- if (!is.null(input$second_group) && input$second_group != "None") {
      input$second_group
    } else {
      NULL
    }

    create_barplot(
      se = se(),
      gene = input$gene,
      assay_name = input$assay_name,
      group_col = input$group_col,
      second_group = second,
      group_levels = input$conditions,
      group_colors = current_colors(),
      ctrl_alpha = input$ctrl_alpha,
      gene_symbol_col = gene_sym_col()
    )
  })

  output$barplot <- renderPlot({ current_barplot() })

  # ---- Gene info table ----
  output$gene_info <- renderDT({
    req(input$gene)
    rd <- as.data.frame(rowData(se()), check.names = FALSE)
    info <- rd[input$gene, , drop = FALSE]

    # Show relevant columns
    sym_col <- gene_sym_col()
    display_cols <- c(sym_col, input$fc_col, input$fdr_col, input$pval_col)
    display_cols <- display_cols[!is.null(display_cols) & display_cols != "(none)"]
    display_cols <- intersect(display_cols, colnames(info))

    if (length(display_cols) == 0) {
      info_display <- info[, seq_len(min(5, ncol(info))), drop = FALSE]
    } else {
      info_display <- info[, display_cols, drop = FALSE]
    }
    rownames(info_display) <- NULL

    num_cols <- names(info_display)[sapply(info_display, is.numeric)]

    dt <- datatable(info_display,
                    options = list(dom = "t", paging = FALSE, searching = FALSE),
                    rownames = FALSE)
    if (length(num_cols) > 0) dt <- formatSignif(dt, columns = num_cols, digits = 4)
    dt
  })

  # ---- Gene search table ----
  gene_table_data <- reactive({
    req(se())
    rd <- as.data.frame(rowData(se()), check.names = FALSE)

    fc <- if (input$fc_col != "(none)") input$fc_col else NULL
    fdr <- if (input$fdr_col != "(none)") input$fdr_col else NULL
    sym_col <- gene_sym_col()

    display_cols <- c(sym_col, fc, fdr, input$pval_col)
    display_cols <- display_cols[!is.null(display_cols) & display_cols != "(none)"]
    display_cols <- intersect(display_cols, colnames(rd))

    if (length(display_cols) == 0) {
      df <- rd[, seq_len(min(5, ncol(rd))), drop = FALSE]
    } else {
      df <- rd[, display_cols, drop = FALSE]
    }

    # Add status if FC and FDR columns exist
    if (!is.null(fc) && !is.null(fdr) && fc %in% colnames(df) && fdr %in% colnames(df)) {
      df$Status <- case_when(
        df[[fdr]] < input$fdr_threshold & df[[fc]] > input$log2fc_threshold ~ "UP",
        df[[fdr]] < input$fdr_threshold & df[[fc]] < -input$log2fc_threshold ~ "DOWN",
        TRUE ~ "NS"
      )

      if (input$de_filter == "sig") df <- df[df$Status != "NS", ]
      else if (input$de_filter == "up") df <- df[df$Status == "UP", ]
      else if (input$de_filter == "down") df <- df[df$Status == "DOWN", ]
    }

    df$ensembl_id <- rownames(df)
    rownames(df) <- NULL
    df
  })

  output$gene_search_table <- renderDT({
    df <- gene_table_data()
    hidden_col <- which(colnames(df) == "ensembl_id") - 1

    dt <- datatable(
      df, selection = "single", filter = "top",
      options = list(
        pageLength = 15, scrollX = TRUE,
        columnDefs = list(list(visible = FALSE, targets = hidden_col))
      ),
      rownames = FALSE
    )

    # Format numeric columns
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) > 0) dt <- formatSignif(dt, columns = num_cols, digits = 3)

    # Style status column if present
    if ("Status" %in% colnames(df)) {
      dt <- formatStyle(dt, "Status",
                        backgroundColor = styleEqual(
                          c("UP", "DOWN", "NS"),
                          c("#FFB2B2", "#99d8ff", "#f0f0f0")),
                        fontWeight = styleEqual(c("UP", "DOWN"), c("bold", "bold")))
    }
    dt
  })

  # Click row in gene search -> select gene and switch tab
  observeEvent(input$gene_search_table_rows_selected, {
    req(input$gene_search_table_rows_selected)
    df <- gene_table_data()
    selected_gene <- df$ensembl_id[input$gene_search_table_rows_selected]
    updateSelectizeInput(session, "gene", selected = selected_gene)
    updateTabsetPanel(session, "main_tabs", selected = "Expression Plot")
  })

  # ---- Volcano plot ----
  volcano_df <- reactive({
    req(se())
    fc <- if (!is.null(input$fc_col) && input$fc_col != "(none)") input$fc_col else NULL
    fdr <- if (!is.null(input$fdr_col) && input$fdr_col != "(none)") input$fdr_col else NULL
    build_volcano_df(
      se = se(), fc_col = fc, fdr_col = fdr,
      gene_symbol_col = gene_sym_col(),
      fc_thresh = input$volcano_fc,
      fdr_thresh = input$volcano_fdr
    )
  })

  current_volcano <- reactive({
    fc <- if (!is.null(input$fc_col) && input$fc_col != "(none)") input$fc_col else "log2FC"
    fdr <- if (!is.null(input$fdr_col) && input$fdr_col != "(none)") input$fdr_col else "FDR"
    create_volcano(
      df = volcano_df(),
      fc_col = fc, fdr_col = fdr,
      fc_thresh = input$volcano_fc,
      fdr_thresh = input$volcano_fdr,
      n_labels = input$volcano_nlabels
    )
  })

  output$volcano <- renderPlot({ current_volcano() })

  # Click on volcano -> select gene and jump to Expression Plot
  observeEvent(input$volcano_click, {
    df <- volcano_df()
    req(!is.null(df), nrow(df) > 0)
    clicked <- nearPoints(df, input$volcano_click,
                          xvar = "log2FC", yvar = "neg_log10_fdr",
                          maxpoints = 1, threshold = 10)
    if (nrow(clicked) > 0) {
      updateSelectizeInput(session, "gene", selected = clicked$gene_id[1])
      updateTabsetPanel(session, "main_tabs", selected = "Expression Plot")
    }
  })

  # ---- PCA plot ----
  current_pca <- reactive({
    req(se(), input$assay_name, input$group_col)

    second <- if (!is.null(input$second_group) && input$second_group != "None") {
      input$second_group
    } else {
      NULL
    }

    create_pca(
      se = se(),
      assay_name = input$assay_name,
      group_col = input$group_col,
      second_group = second,
      group_colors = current_colors(),
      ntop = input$pca_ntop,
      pc_x = as.integer(gsub("PC", "", input$pca_pc_x)),
      pc_y = as.integer(gsub("PC", "", input$pca_pc_y)),
      show_ellipses = isTRUE(input$pca_ellipses)
    )
  })

  output$pca <- renderPlot({ current_pca() })

  # ---- Download handlers ----

  # Barplot
  output$download_barplot_png <- downloadHandler(
    filename = function() {
      sym <- tryCatch({
        rd <- as.data.frame(rowData(se()), check.names = FALSE)
        s <- gene_sym_col()
        if (!is.null(s)) rd[input$gene, s] else input$gene
      }, error = function(e) input$gene)
      if (is.na(sym)) sym <- input$gene
      paste0(sym, "_barplot.png")
    },
    content = function(file) {
      ggsave(file, plot = current_barplot(), width = 10, height = 6, dpi = 300, bg = "white")
    }
  )

  output$download_barplot_pdf <- downloadHandler(
    filename = function() {
      sym <- tryCatch({
        rd <- as.data.frame(rowData(se()), check.names = FALSE)
        s <- gene_sym_col()
        if (!is.null(s)) rd[input$gene, s] else input$gene
      }, error = function(e) input$gene)
      if (is.na(sym)) sym <- input$gene
      paste0(sym, "_barplot.pdf")
    },
    content = function(file) {
      ggsave(file, plot = current_barplot(), width = 10, height = 6)
    }
  )

  # Volcano
  output$download_volcano_png <- downloadHandler(
    filename = function() paste0(input$active_dataset, "_volcano.png"),
    content = function(file) {
      ggsave(file, plot = current_volcano(), width = 10, height = 8, dpi = 300, bg = "white")
    }
  )

  output$download_volcano_pdf <- downloadHandler(
    filename = function() paste0(input$active_dataset, "_volcano.pdf"),
    content = function(file) {
      ggsave(file, plot = current_volcano(), width = 10, height = 8)
    }
  )

  # PCA
  output$download_pca_png <- downloadHandler(
    filename = function() paste0(input$active_dataset, "_pca.png"),
    content = function(file) {
      ggsave(file, plot = current_pca(), width = 10, height = 8, dpi = 300, bg = "white")
    }
  )

  output$download_pca_pdf <- downloadHandler(
    filename = function() paste0(input$active_dataset, "_pca.pdf"),
    content = function(file) {
      ggsave(file, plot = current_pca(), width = 10, height = 8)
    }
  )
}

# --- Run App ---
shinyApp(ui = ui, server = server)
