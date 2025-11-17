library(shiny)
library(tidyverse)

#-----------------------------
# Load and preprocess data
#-----------------------------
params  <- read_csv("IMPC_parameters_with_categories.csv")
cleaned <- read_csv("IMPC_cleaned3.csv")

data_merged <- cleaned %>%
  inner_join(params, by = c("parameter_id" = "parameterId")) %>%
  mutate(pvalue = as.numeric(pvalue)) %>%
  filter(!is.na(pvalue))

# Look-up vectors for UI
all_categories <- sort(unique(params$category))
all_genes      <- sort(unique(data_merged$gene_symbol))
all_parameters <- data_merged %>%
  distinct(parameter_id) %>%
  arrange(parameter_id) %>%
  pull(parameter_id)

#-----------------------------
# UI
#-----------------------------
ui <- fluidPage(
  titlePanel("Knockout Mouse Phenotype Dashboard"),
  sidebarLayout(
    sidebarPanel(
      # Used by all tabs
      selectInput(
        inputId  = "categorySelect",
        label    = "Select phenotype categories:",
        choices  = all_categories,
        selected = all_categories,
        multiple = TRUE
      ),
      
      # Controls for Gene View
      conditionalPanel(
        condition = "input.tabs == 'Gene View'",
        selectInput(
          inputId = "geneSelect",
          label   = "Select knockout gene(s):",
          choices = all_genes,
          multiple = TRUE
        )
      ),
      
      # Controls for Phenotype View
      conditionalPanel(
        condition = "input.tabs == 'Phenotype View'",
        tagList(
          selectInput(
            inputId  = "paramSelect",
            label    = "Select phenotype (parameter):",
            choices  = all_parameters,
            multiple = FALSE
          ),
          sliderInput(
            inputId = "topGenes",
            label   = "Number of top genes to show:",
            min     = 5,
            max     = 50,
            value   = 20,
            step    = 1
          )
        )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Gene View",
                 plotOutput("genePlot")
        ),
        tabPanel("Phenotype View",
                 plotOutput("phenoPlot")
        ),
        tabPanel("Gene Clusters",
                 plotOutput("clusterPlot")
        )
      )
    )
  )
)

#-----------------------------
# Server
#-----------------------------
server <- function(input, output, session) {
  
  # Filter by selected categories – shared across views
  filtered_data <- reactive({
    req(input$categorySelect)
    data_merged %>%
      filter(category %in% input$categorySelect)
  })
  
  # Keep gene choices in sync with category filter
  observeEvent(input$categorySelect, {
    genes_available <- filtered_data() %>%
      distinct(gene_symbol) %>%
      arrange(gene_symbol) %>%
      pull(gene_symbol)
    
    choices_list   <- setNames(as.list(genes_available), genes_available)
    selected_genes <- intersect(input$geneSelect, genes_available)
    
    updateSelectInput(
      session, "geneSelect",
      choices  = choices_list,
      selected = selected_genes
    )
  })
  
  sig_line <- -log10(0.05)  # p = 0.05 threshold
  
  #----------------------------------
  # (1) GENE VIEW – per-gene phenotype summary
  #----------------------------------
  output$genePlot <- renderPlot({
    req(input$geneSelect)
    
    df_gene <- filtered_data() %>%
      filter(gene_symbol %in% input$geneSelect)
    
    validate(need(nrow(df_gene) > 0,
                  "No data for the selected gene(s) and category(ies)."))
    
    df_gene_summary <- df_gene %>%
      group_by(gene_symbol, category) %>%
      summarise(min_p = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        min_p = ifelse(min_p <= 0 | is.na(min_p), NA_real_, min_p),
        logp  = -log10(min_p)
      )
    
    df_gene_summary$category <- factor(
      df_gene_summary$category,
      levels = sort(unique(df_gene_summary$category))
    )
    
    ggplot(df_gene_summary, aes(x = category, y = logp)) +
      geom_col(width = 0.7, aes(fill = category), show.legend = FALSE) +
      geom_hline(yintercept = sig_line, linetype = "dashed") +
      facet_wrap(~ gene_symbol, scales = "free_y", ncol = 4) +
      labs(
        x = "Phenotype category",
        y = "-log10(p-value)",
        title = "Phenotype effects per knockout gene"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        strip.text  = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5)
      )
  })
  
  #----------------------------------
  # (2) PHENOTYPE VIEW – single phenotype vs genes
  #----------------------------------
  output$phenoPlot <- renderPlot({
    req(input$paramSelect)
    
    df_pheno <- filtered_data() %>%
      filter(parameter_id == input$paramSelect)
    
    validate(need(nrow(df_pheno) > 0,
                  "No data for the selected phenotype in the chosen categories."))
    
    df_single_gene <- df_pheno %>%
      group_by(gene_symbol) %>%
      summarise(min_p = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        min_p = ifelse(min_p <= 0 | is.na(min_p), NA_real_, min_p),
        logp  = -log10(min_p)
      ) %>%
      filter(!is.na(logp)) %>%
      arrange(desc(logp))
    
    validate(need(nrow(df_single_gene) > 0,
                  "No valid p-values for this phenotype."))
    
    # Show top N genes (highest -log10(p))
    n_show <- min(input$topGenes, nrow(df_single_gene))
    df_plot <- df_single_gene %>%
      slice_max(order_by = logp, n = n_show, with_ties = FALSE)
    
    ggplot(df_plot,
           aes(x = logp,
               y = reorder(gene_symbol, logp))) +
      geom_col(fill = "steelblue") +
      geom_vline(xintercept = sig_line, linetype = "dashed") +
      labs(
        x = "-log10(p-value)",
        y = "Knockout gene",
        title = paste("Knockout effects for phenotype:", input$paramSelect)
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  })
  
  #----------------------------------
  # (3) GENE CLUSTERS – clustered heatmap
  #----------------------------------
  output$clusterPlot <- renderPlot({
    df_heat <- filtered_data()
    validate(need(nrow(df_heat) > 0,
                  "No data for the selected category(ies)."))
    
    df_heat_summary <- df_heat %>%
      group_by(category, gene_symbol) %>%
      summarise(min_p = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        min_p = ifelse(min_p <= 0 | is.na(min_p), NA_real_, min_p),
        logp  = -log10(min_p)
      )
    
    # Wide matrix for clustering
    mat <- df_heat_summary %>%
      select(gene_symbol, category, logp) %>%
      pivot_wider(names_from = category, values_from = logp, values_fill = 0)
    
    gene_mat <- as.matrix(mat[,-1])
    rownames(gene_mat) <- mat$gene_symbol
    
    # Hierarchical clustering on genes
    hc_genes   <- hclust(dist(gene_mat))
    gene_order <- rownames(gene_mat)[hc_genes$order]
    
    df_heat_summary$gene_symbol <- factor(
      df_heat_summary$gene_symbol,
      levels = gene_order
    )
    
    df_heat_summary$category <- factor(
      df_heat_summary$category,
      levels = sort(unique(df_heat_summary$category))
    )
    
    ggplot(df_heat_summary,
           aes(x = gene_symbol, y = category, fill = logp)) +
      geom_tile(color = "white") +
      scale_fill_gradient(
        low = "white", high = "steelblue", na.value = "grey90",
        name = "-log10(p-value)"
      ) +
      labs(
        x = "Knockout gene (clustered)",
        y = "Phenotype category",
        title = "Clusters of genes with similar phenotype profiles"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5)
      )
  })
}

#-----------------------------
# Run the app
#-----------------------------
shinyApp(ui = ui, server = server)
