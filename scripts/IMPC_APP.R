library(shiny)
library(tidyverse)
library(DT)
library(plotly)
library(viridis)


setwd("")


# 0. Custom Theme

theme_light_md3 <- function(base_size = 13){
  theme_minimal(base_size = base_size) +
    theme(
      plot.background     = element_rect(fill = "white", color = NA),
      panel.grid.major.y  = element_line(color = "#E5E7EB"),
      panel.grid.major.x  = element_blank(),
      axis.text           = element_text(color = "grey20"),
      axis.title          = element_text(color = "grey20", face = "bold")
    )
}


# 1. Load and prepare data


# Parameter metadata
param_meta <- readr::read_csv(
  "IMPC_parameter_groupings.csv",
  show_col_types = FALSE
) %>%
  distinct(parameterId, .keep_all = TRUE) %>%
  rename(category = groupings)

# Phenotype scores
pheno_raw <- readr::read_csv(
  "IMPC_cleaned3.csv",
  show_col_types = FALSE
)

# Coursework genes (optional default selection)
coursework_genes <- c("Trim25", "Cpz", "Pibf1", "Sulf1")

# Safe -log10(p)
safe_log10 <- function(p) {
  p2 <- ifelse(is.na(p) | p <= 0, 1e-16, p)
  -log10(p2)
}

# Helper to resolve "All ..." selections
resolve_selection <- function(selected, all_choices, all_label) {
  if (is.null(selected) || length(selected) == 0 || all_label %in% selected) {
    all_choices
  } else {
    setdiff(selected, all_label)
  }
}

alpha    <- 0.05
sig_line <- -log10(alpha)

# Merge data and basic cleaning
data_all <- pheno_raw %>%
  mutate(
    mouse_strain = dplyr::if_else(is.na(mouse_strain), "Unknown/NA", mouse_strain),
    pvalue       = as.numeric(pvalue)
  ) %>%
  left_join(param_meta, by = c("parameter_id" = "parameterId")) %>%
  mutate(
    category = dplyr::if_else(is.na(category), "Uncategorized", category)
  ) %>%
  filter(
    !is.na(pvalue),
    !is.na(gene_symbol),
    !is.na(parameter_id)
  )


# Create choices for dropdowns

all_genes      <- sort(unique(data_all$gene_symbol))
all_categories <- sort(unique(data_all$category))
all_strains    <- sort(unique(data_all$mouse_strain))
all_parameters <- sort(unique(data_all$parameter_id))

# Parameter name lookup (for titles)
param_map <- data_all %>%
  distinct(parameter_name, parameter_id) %>%
  arrange(parameter_name)

# Labels for "All"
ALL_STRAINS_LABEL <- "All strains"
ALL_PARAMS_LABEL  <- "All phenotypes"

# 2. UI

ui <- fluidPage(
  titlePanel("IMPC Knockout Mouse Phenotype Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Global filters
      selectInput(
        "strain",
        "Mouse strain:",
        choices  = c(ALL_STRAINS_LABEL, all_strains),
        selected = ALL_STRAINS_LABEL,
        multiple = TRUE
      ),
      
      selectInput(
        "phenogroup",
        "Phenotype group:",
        choices  = c("All", all_categories),
        selected = "All"
      ),
      
      tags$hr(),
      
      # Gene controls – only on Gene View
      conditionalPanel(
        condition = "input.main_tabs == 'Gene View'",
        
        selectInput(
          "gene",
          "Select Single Gene:",
          choices  = all_genes,
          selected = coursework_genes[1],
          multiple = FALSE
        )
      ),
      
      # Phenotype controls – only on Phenotype View
      conditionalPanel(
        condition = "input.main_tabs == 'Phenotype View'",
        
        h4("Phenotype view options"),
        
        selectInput(
          "parameter",
          "Select Phenotype (ID):",
          choices  = c("All phenotypes" = ALL_PARAMS_LABEL, all_parameters),
          selected = ALL_PARAMS_LABEL,
          multiple = TRUE
        ),
        
        sliderInput(
          "max_genes_param",
          "Number of genes to show:",
          min    = 10,
          max    = 60,
          value  = 30,
          step   = 5
        )
      ),
      
      # Cluster controls – for Gene Clusters & PCA View
      conditionalPanel(
        condition = "input.main_tabs == 'Gene Clusters' || input.main_tabs == 'PCA View'",
        
        sliderInput(
          "k_max",
          "Max clusters for elbow:",
          min    = 3,
          max    = 12,
          value  = 8,
          step   = 1
        ),
        
        sliderInput(
          "k_clusters",
          "Number of clusters (k):",
          min    = 2,
          max    = 8,
          value  = 3,
          step   = 1
        )
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "Gene View",
          br(),
          p("Visualising full phenotypic profile. Dashed black line indicates p = 0.05."),
          plotlyOutput("gene_bar", height = "600px"),
          br(),
          h4("Raw data for selected gene"),
          DTOutput("gene_table")
        ),
        
        tabPanel(
          "Phenotype View",
          br(),
          p("Choose one or more phenotypes to see which genes show the strongest effects."),
          plotOutput("param_bar", height = "400px"),
          br(),
          h4("Raw data for these phenotypes"),
          DTOutput("param_table")
        ),
        
        tabPanel(
          "Gene Clusters",
          br(),
          p("Interactive heatmap of genes × phenotype categories."),
          plotlyOutput("cluster_plot", height = "500px")
        ),
        
        tabPanel(
          "PCA View",
          br(),
          p("Interactive PCA of genes based on −log10(p) across parameters."),
          plotlyOutput("pca_plot", height = "500px"),
          br(),
          h4("Elbow plot for k-means clustering"),
          plotOutput("elbow_plot", height = "300px"),
          br(),
          h4("Top phenotype categories per cluster"),
          DTOutput("cluster_category_table")
        )
      )
    )
  )
)


# 3. Server


server <- function(input, output, session) {
  
 #Shared filtered dataset
  
  data_filtered <- reactive({
    req(input$strain, input$phenogroup)
    
    sel_strains <- resolve_selection(input$strain, all_strains, ALL_STRAINS_LABEL)
    
    df <- data_all %>%
      filter(mouse_strain %in% sel_strains)
    
    if (!is.null(input$phenogroup) && input$phenogroup != "All") {
      df <- df %>% filter(category == input$phenogroup)
    }
    
    validate(
      need(nrow(df) > 0, "No data for these filters. Try different strains or phenotype group.")
    )
    
    df
  })
  
#Gene View 
  
  gene_summary <- reactive({
    req(input$gene)
    
    df <- data_filtered() %>%
      filter(gene_symbol == input$gene)
    
    validate(
      need(nrow(df) > 0, "Selected gene has no data for the current filters.")
    )
    
    df_sum <- df %>%
      group_by(gene_symbol, parameter_id, parameter_name, category) %>%
      summarise(min_p = min(pvalue), .groups = "drop") %>%
      mutate(logp = safe_log10(min_p))
    
    df_sum %>% arrange(desc(logp))
  })
  
  output$gene_bar <- renderPlotly({
    df <- gene_summary()
    
    p <- ggplot(df, aes(
      x     = category,
      y     = logp,
      color = category,
      text  = paste(
        "Gene:", gene_symbol,
        "<br>Phenotype:", parameter_name,
        "<br>P-value:", signif(min_p, 3),
        "<br>Category:", category
      )
    )) +
      geom_jitter(width = 0.25, size = 3, alpha = 0.8) +
      geom_hline(yintercept = sig_line, linetype = "dashed", color = "black", size = 0.5) +
      annotate("text", x = 1, y = sig_line + 0.2, label = "p = 0.05",
               color = "black", size = 3.5, hjust = 0) +
      scale_color_viridis(discrete = TRUE, option = "D") +
      labs(
        title = paste("Phenotypic Profile:", input$gene),
        x     = "Phenotype Category",
        y     = "Significance Score (-log10 p-value)",
        color = "Category"
      ) +
      theme_light_md3() +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1),
        legend.position  = "none"
      )
    
    ggplotly(p, tooltip = "text")
  })
  
  output$gene_table <- renderDT({
    req(input$gene)
    
    df <- data_filtered() %>%
      filter(gene_symbol == input$gene)
    
    df <- df %>%
      mutate(pvalue = if_else(pvalue == 0, 1e-16, pvalue)) %>%
      select(
        gene_symbol, mouse_strain,
        parameter_id, parameter_name,
        category, pvalue
      ) %>%
      arrange(parameter_id, pvalue)
    
    datatable(
      df,
      options  = list(pageLength = 10),
      rownames = FALSE
    ) %>%
      formatSignif("pvalue", digits = 4)
  })
  
  #---------------- Phenotype View ----------------
  
  param_summary <- reactive({
    req(input$parameter)
    sel_params <- resolve_selection(input$parameter, all_parameters, ALL_PARAMS_LABEL)
    
    df <- data_filtered() %>%
      filter(parameter_id %in% sel_params)
    
    validate(
      need(nrow(df) > 0, "No data for these phenotypes under the current filters.")
    )
    
    df_gene <- df %>%
      group_by(gene_symbol) %>%
      summarise(min_p = min(pvalue), .groups = "drop")
    
    df_gene %>%
      mutate(logp = safe_log10(min_p)) %>%
      arrange(desc(logp)) %>%
      slice_head(n = input$max_genes_param)
  })
  
  output$param_bar <- renderPlot({
    df <- param_summary()
    sel_params <- resolve_selection(input$parameter, all_parameters, ALL_PARAMS_LABEL)
    
    if (length(sel_params) == length(all_parameters)) {
      title_txt <- "Genes affecting any phenotype (all parameters selected)"
    } else if (length(sel_params) == 1) {
      p_name <- param_map$parameter_name[param_map$parameter_id == sel_params]
      if (length(p_name) == 0) p_name <- "Unknown"
      title_txt <- paste0(p_name, " (", sel_params, ")")
    } else {
      title_txt <- paste("Genes affecting any of", length(sel_params), "selected phenotypes")
    }
    
    ggplot(df, aes(
      x = reorder(gene_symbol, logp),
      y = logp
    )) +
      geom_col(fill = "steelblue") +
      geom_hline(yintercept = sig_line, linetype = "dashed", color = "black", size = 0.8) +
      annotate("text", x = 1, y = sig_line + 0.2, label = "p = 0.05",
               color = "black", size = 4, hjust = 0) +
      coord_flip() +
      labs(
        x     = "Gene",
        y     = "-log10(p-value)",
        title = title_txt
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold")
      )
  })
  
  output$param_table <- renderDT({
    req(input$parameter)
    sel_params <- resolve_selection(input$parameter, all_parameters, ALL_PARAMS_LABEL)
    
    df <- data_filtered() %>%
      filter(parameter_id %in% sel_params)
    
    df <- df %>%
      mutate(pvalue = if_else(pvalue == 0, 1e-16, pvalue)) %>%
      select(
        gene_symbol, mouse_strain,
        parameter_id, parameter_name,
        category, pvalue
      ) %>%
      arrange(parameter_id, pvalue)
    
    datatable(
      df,
      options  = list(pageLength = 10),
      rownames = FALSE
    ) %>%
      formatSignif("pvalue", digits = 4)
  })
  
  #---------------- Gene × parameter matrix for clustering/PCA ----------------
  
  gene_param_data <- reactive({
    df <- data_filtered()
    
    df_mat <- df %>%
      group_by(gene_symbol, parameter_id, category) %>%
      summarise(min_p = min(pvalue), .groups = "drop") %>%
      mutate(
        logp       = safe_log10(min_p),
        # truncate non-significant to 0 to emphasise significant patterns
        logp_trunc = if_else(min_p <= alpha, logp, 0)
      )
    
    # Keep only genes with at least one significant phenotype
    sig_genes <- df_mat %>%
      group_by(gene_symbol) %>%
      summarise(any_sig = any(min_p <= alpha), .groups = "drop") %>%
      filter(any_sig) %>%
      pull(gene_symbol)
    
    if (length(sig_genes) == 0) {
      return(list(
        mat    = matrix(nrow = 0, ncol = 0),
        genes  = character(0),
        df_mat = df_mat[0, ]
      ))
    }
    
    df_mat_use <- df_mat %>% filter(gene_symbol %in% sig_genes)
    
    mat <- df_mat_use %>%
      select(gene_symbol, parameter_id, logp_trunc) %>%
      tidyr::pivot_wider(
        names_from  = parameter_id,
        values_from = logp_trunc,
        values_fill = 0
      )
    
    mat_num <- as.matrix(mat[, -1, drop = FALSE])
    rownames(mat_num) <- mat$gene_symbol
    
    list(
      mat    = mat_num,
      genes  = mat$gene_symbol,
      df_mat = df_mat_use
    )
  })
  
  # Elbow method: WSS vs k 
  
  elbow_data <- reactive({
    gp      <- gene_param_data()
    mat_num <- gp$mat
    n_genes <- nrow(mat_num)
    
    validate(
      need(n_genes >= 2, "Need at least two genes with significant phenotypes for elbow method.")
    )
    
    req(input$k_max)
    
    max_k <- min(input$k_max, n_genes - 1)  # k cannot exceed n_genes - 1
    
    wss_vals <- numeric(max_k)
    
    # k = 1: total variance
    wss_vals[1] <- sum(stats::var(mat_num) * (nrow(mat_num) - 1))
    
    if (max_k >= 2) {
      for (k in 2:max_k) {
        set.seed(123)
        km <- kmeans(mat_num, centers = k, nstart = 10)
        wss_vals[k] <- km$tot.withinss
      }
    }
    
    tibble(
      k   = 1:max_k,
      wss = wss_vals
    )
  })
  
  output$elbow_plot <- renderPlot({
    df <- elbow_data()
    
    ggplot(df, aes(x = k, y = wss)) +
      geom_line() +
      geom_point() +
      labs(
        x     = "Number of clusters (k)",
        y     = "Total within-cluster sum of squares",
        title = "Elbow plot for k-means clustering"
      ) +
      theme_minimal(base_size = 13)
  })
  
#Cluster assignment (k-means)
  
  cluster_info <- reactive({
    gp      <- gene_param_data()
    mat_num <- gp$mat
    genes   <- gp$genes
    n_genes <- nrow(mat_num)
    
    if (n_genes == 0) {
      return(tibble(
        gene_symbol = character(0),
        cluster     = factor()
      ))
    }
    
    req(input$k_clusters, input$k_max)
    
    # Clamp k to valid range
    k <- input$k_clusters
    k <- max(2, min(k, n_genes))
    k <- min(k, input$k_max)
    
    set.seed(123)
    km <- kmeans(mat_num, centers = k, nstart = 25)
    
    tibble(
      gene_symbol = genes,
      cluster     = factor(
        km$cluster,
        labels = paste("Cluster", seq_len(k))
      )
    )
  })
  
# Gene Clusters heatmap 
  
  output$cluster_plot <- renderPlotly({
    df <- data_filtered()
    
    df_sum <- df %>%
      group_by(gene_symbol, category) %>%
      summarise(min_p = min(pvalue), .groups = "drop") %>%
      mutate(logp = safe_log10(min_p))
    
    mat <- df_sum %>%
      tidyr::pivot_wider(
        names_from  = category,
        values_from = logp,
        values_fill = 0
      )
    
    gene_mat <- as.matrix(mat[, -1, drop = FALSE])
    rownames(gene_mat) <- mat$gene_symbol
    
    validate(
      need(nrow(gene_mat) >= 2, "Need at least two genes for clustering.")
    )
    
    hc <- hclust(dist(gene_mat))
    gene_order <- unique(rownames(gene_mat)[hc$order])
    df_sum$gene_symbol <- factor(df_sum$gene_symbol, levels = gene_order)
    
    plot_ly(
      data  = df_sum,
      x     = ~gene_symbol,
      y     = ~category,
      z     = ~logp,
      type  = "heatmap",
      colors = "Blues",
      text  = ~paste0(
        "Gene: ", gene_symbol,
        "<br>Category: ", category,
        "<br>-log10(p): ", round(logp, 2)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis = list(title = "Gene (clustered)", showticklabels = FALSE),
        yaxis = list(title = "Phenotype category"),
        title = "Clusters of genes with similar phenotype profiles"
      )
  })
  
  # PCA
  
  output$pca_plot <- renderPlotly({
    gp      <- gene_param_data()
    mat_num <- gp$mat
    genes   <- gp$genes
    n_genes <- nrow(mat_num)
    
    validate(
      need(n_genes >= 2, "Need at least two genes with significant phenotypes for PCA.")
    )
    
    clusters <- cluster_info()
    
    pca_res <- prcomp(mat_num, scale. = TRUE)
    scores  <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
    scores$gene_symbol <- genes
    scores <- dplyr::left_join(scores, clusters, by = "gene_symbol")
    
    var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
    pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
    pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
    
    scores$tooltip <- paste0(
      "Gene: ", scores$gene_symbol,
      "<br>Cluster: ", scores$cluster,
      "<br>PC1: ", round(scores$PC1, 2),
      "<br>PC2: ", round(scores$PC2, 2)
    )
    
    plot_ly(
      data  = scores,
      x     = ~PC1,
      y     = ~PC2,
      type  = "scatter",
      mode  = "markers",
      color = ~cluster,
      text  = ~tooltip,
      hoverinfo = "text"
    ) %>%
      layout(
        title = "PCA of genes based on phenotype score profiles",
        xaxis = list(title = pc1_lab),
        yaxis = list(title = pc2_lab)
      )
  })
  
  #Cluster summary table 
  output$cluster_category_table <- renderDT({
    gp      <- gene_param_data()
    df_mat  <- gp$df_mat
    mat_num <- gp$mat
    n_genes <- nrow(mat_num)
    
    validate(
      need(n_genes >= 2, "Need at least two genes with significant phenotypes to summarise clusters.")
    )
    
    clusters <- cluster_info()
    validate(
      need(nrow(clusters) > 0, "No clusters to summarise.")
    )
    
    df_cat <- df_mat %>%
      select(gene_symbol, category, logp_trunc) %>%
      left_join(clusters, by = "gene_symbol") %>%
      group_by(cluster, category) %>%
      summarise(mean_logp = mean(logp_trunc), .groups = "drop") %>%
      filter(mean_logp > 0) %>%
      arrange(cluster, desc(mean_logp))
    
    validate(
      need(nrow(df_cat) > 0, "No significant categories to summarise for these clusters.")
    )
    
    top_n <- 3
    table_data <- df_cat %>%
      group_by(cluster) %>%
      slice_max(mean_logp, n = top_n, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(mean_logp = round(mean_logp, 2)) %>%
      arrange(cluster, desc(mean_logp))
    
    datatable(
      table_data,
      options  = list(
        pageLength = nrow(table_data),
        searching  = FALSE,
        dom        = "t"
      ),
      rownames = FALSE
    )
  })
}


#4 Run the app


shinyApp(ui = ui, server = server)
