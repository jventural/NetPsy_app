# app.R

# ---- 1. Librerías ----
library(shiny)
library(shinydashboard)
library(shinyWidgets)    # progressBar
library(shinyjs)         # enable/disable botones
library(readxl)
library(ThesiStats)
library(PsyMetricTools)
library(dplyr)
library(InterconectaR)
library(psych)
library(huge)
library(corrplot)
library(matrixcalc)
library(networktools)
library(qgraph)
library(bootnet)
library(mgm)
library(tidyverse)       # ggplot2, tidyr, stringr…
library(forcats)         # fct_reorder
library(tibble)          # enframe()
library(parallel)        # detectCores()
library(cowplot)         # ggdraw, draw_image, plot_grid

# ---- 2. Detectar núcleos físicos disponibles ----
max_cores <- detectCores(logical = FALSE)
if (is.na(max_cores) || max_cores < 1) max_cores <- 1

# ---- 3. Función centrality_plots2_fixed MEJORADA ----
centrality_plots2_fixed <- function(qgraph_obj,
                                    network,
                                    groups = NULL,
                                    measure0 = "ExpectedInfluence",
                                    measure1 = NULL,
                                    color_palette = c("#F8766D", "#00BFC4"),
                                    use_abbrev = TRUE) {
  library(qgraph); library(dplyr); library(networktools)
  library(tibble); library(ggplot2); library(tidyr)
  library(stringr); library(forcats); library(purrr)
  
  # 1) Centralidad estándar
  ctbl <- centralityTable(network)
  cents_expect <- ctbl %>%
    filter(measure == !!measure0) %>%
    select(full_name = node, value) %>%
    arrange(full_name) %>%
    rename(!!measure0 := value) %>%
    mutate(
      Abrev = full_name %>%
        str_split(" ") %>% map_chr(~ .x[[1]]) %>%
        substr(1,3) %>% str_to_title()
    )
  
  # 2) Puente (opcional)
  if (!is.null(measure1)) {
    b_obj <- bridge(qgraph_obj,
                    communities    = groups,
                    useCommunities = "all",
                    normalize      = FALSE)
    qgraph_labels <- qgraph_obj$graphAttributes$Nodes$labels
    
    bridge_data <- tibble(
      qgraph_label = qgraph_labels,
      raw_bridge   = as.numeric(b_obj[[measure1]])
    ) %>%
      mutate(!!measure1 := as.numeric(scale(raw_bridge))) %>%
      select(-raw_bridge)
    
    full_names <- rownames(network$graph)
    if (is.null(full_names)) full_names <- colnames(network$graph)
    name_mapping <- tibble(
      full_name    = full_names,
      qgraph_label = qgraph_labels
    )
    
    bridge_data <- bridge_data %>%
      left_join(name_mapping, by = "qgraph_label") %>%
      mutate(
        Abrev = full_name %>%
          str_split(" ") %>% map_chr(~ .x[[1]]) %>%
          substr(1,3) %>% str_to_title()
      )
    
    cents2 <- inner_join(
      cents_expect,
      bridge_data %>% select(full_name, !!sym(measure1)),
      by = "full_name"
    )
    if (nrow(cents2) == 0) {
      if (nrow(cents_expect) == nrow(bridge_data)) {
        cents2 <- bind_cols(
          cents_expect,
          bridge_data %>% select(!!sym(measure1))
        )
      } else stop("No se pueden unir centralidad y puente")
    }
    
    # Preparamos los datos para ggplot
    y_var     <- if (use_abbrev) quo(Abrev) else quo(full_name)
    cents_long <- cents2 %>%
      pivot_longer(
        cols      = c(!!sym(measure0), !!sym(measure1)),
        names_to  = "Measure",
        values_to = "Value"
      )
    
    # Aquí mapeamos correctamente métricas → colores
    pal    <- setNames(color_palette, c(measure0, measure1))
    breaks <- c(measure0, measure1)
    labels <- c(measure0, measure1)
    
    Figura <- ggplot(cents_long,
                     aes(x = Value,
                         y = fct_reorder(!!y_var, Value),
                         color = Measure, group = Measure
                     )) +
      geom_point(size = 3) +
      geom_line(size = 0.5) +
      theme_minimal() +
      labs(x = "z-score", y = "Nodos", color = "Métrica") +
      scale_color_manual(
        values = pal,
        breaks = breaks,
        labels = labels
      ) +
      theme(
        axis.text.y  = element_text(size = 12),
        axis.text.x  = element_text(size = 12),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
      )
    
    class(Figura) <- c("silent_gg", class(Figura))
    assign("print.silent_gg",
           function(x, ...) suppressWarnings(NextMethod()),
           envir = .GlobalEnv)
    
    return(list(
      table = cents2 %>% arrange(desc(!!sym(measure0))),
      plot  = Figura
    ))
  }
  
  # 3) Solo measure0
  cents_single <- cents_expect %>%
    select(full_name, Abrev, !!sym(measure0)) %>%
    arrange(desc(!!sym(measure0))) %>%
    rename(Value = !!sym(measure0))
  
  pal_single <- setNames(color_palette[1], measure0)
  
  Figura <- ggplot(cents_single,
                   aes(x = Value,
                       y = fct_reorder(Abrev, Value),
                       group = 1)) +
    geom_line(color = pal_single, size = 0.5) +
    geom_point(color = pal_single, size = 3) +
    theme_minimal() +
    labs(x = "z-score", y = "Nodos", color = "Métrica") +
    scale_color_manual(
      values = pal_single,
      breaks = measure0,
      labels = measure0
    ) +
    theme(
      axis.text.y  = element_text(size = 12),
      axis.text.x  = element_text(size = 12),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  class(Figura) <- c("silent_gg", class(Figura))
  assign("print.silent_gg",
         function(x, ...) suppressWarnings(NextMethod()),
         envir = .GlobalEnv)
  
  return(list(table = cents_single, plot = Figura))
}

# ---- 3.5 Función auxiliar para crear figura combinada ----
create_combined_plot <- function(qgraph_obj, cent_plot, method = "base") {
  if (method == "base") {
    return(function() {
      oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
      layout(matrix(c(1,2), 1, 2), widths = c(2, 1))
      par(mar = c(2, 2, 2, 2)); plot(qgraph_obj)
      par(mar = c(5, 8, 4, 2)); plot.new()
      text(0.5, 0.5, "Centralidad\n(Ver pestaña Centralidad)",
           cex = 1.2, col = "darkblue")
    })
  } else if (method == "cowplot") {
    tryCatch({
      tmp <- tempfile(fileext = ".png")
      png(tmp, width = 800, height = 600, res = 120)
      plot(qgraph_obj); dev.off()
      qgraph_gg <- cowplot::ggdraw() +
        cowplot::draw_image(tmp) +
        theme(plot.margin = margin(0, 0, 0, 0))
      combined <- cowplot::plot_grid(
        qgraph_gg, cent_plot,
        ncol       = 2,
        rel_widths = c(2, 1),
        labels     = c("A", "B")
      )
      unlink(tmp)
      return(combined)
    }, error = function(e) {
      return(NULL)
    })
  }
}

# ---- 4. UI ----
ui <- dashboardPage(
  dashboardHeader(title = "Análisis de Redes Psicológicas"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Configuración", tabName = "config", icon = icon("sliders")),
      menuItem("Datos", tabName = "data", icon = icon("table")),
      menuItem("Descriptivos", tabName = "desc", icon = icon("chart-bar")),
      menuItem("Estimación de la red", tabName = "report", icon = icon("chart-pie")),
      menuItem("Estabilidad y Precisión", tabName = "boot", icon = icon("refresh")),
      menuItem("Comparación Grupos", tabName = "groups", icon = icon("users"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(
      tabItem("config",
              fluidRow(
                box(width = 6, fileInput("file", "Cargar datos (.xlsx/.csv)", accept = c(".xlsx",".csv"))),
                box(width = 6, selectInput("vars", "Seleccionar variables", choices = NULL, multiple = TRUE)),
                box(width = 6, textInput("group_names", "Nombres de grupos", "Dependencia,Miedo Soledad,Bienestar"),
                    textInput("group_values", "Tamaños de grupos", "3,1,1")),
                box(width = 6, selectInput("measure1", "Medida puente (opcional)",
                                           choices = c("Ninguna" = "None",
                                                       "Bridge Strength" = "bridgeStrength",
                                                       "Bridge Betweenness" = "bridgeBetweenness",
                                                       "Bridge Closeness" = "bridgeCloseness",
                                                       "Bridge Expected Influence" = "Bridge Expected Influence (1-step)"),
                                           selected = "None")),
                box(width = 12, align = "center", actionButton("run", "Ejecutar análisis", icon = icon("play")))
              )
      ),
      tabItem("data",
              fluidRow(box(width = 12, status = "primary", solidHeader = TRUE, tableOutput("df_vars")))
      ),
      tabItem("desc",
              fluidRow(
                tabBox(width = 12, title = "Descriptivos",
                       tabPanel("Estadísticos", tableOutput("desc_table")),
                       tabPanel("Goldbricker", verbatimTextOutput("goldbr")),
                       tabPanel("Edges", tableOutput("edge_summary")),
                       tabPanel("Density", verbatimTextOutput("density_report"))
                )
              )
      ),
      tabItem("report",
              fluidRow(
                tabBox(width = 12, id = "network_tabs", title = "Estimación de la Red",
                       tabPanel("Grafo",
                                box(width = 12, status = "info", solidHeader = TRUE, title = "Visualización de la Red",
                                    plotOutput("network_plot", height = "600px"))
                       ),
                       tabPanel("Centralidad",
                                box(width = 12, status = "info", solidHeader = TRUE, title = "Centralidad de Nodos",
                                    tableOutput("cent_table"),
                                    plotOutput("cent_plot", height = "400px"))
                       ),
                       tabPanel("Figura Combinada",
                                box(width = 12, status = "success", solidHeader = TRUE, title = "Red + Centralidad",
                                    plotOutput("combined_plot", height = "600px", width = "100%"),
                                    verbatimTextOutput("combined_plot_status")),
                                box(width = 12, align = "center",
                                    downloadButton("downloadReport", "Descargar Figura Combinada", class = "btn-success"))
                       )
                )
              )
      ),
      tabItem("boot",
              fluidRow(
                tabBox(width = 12, id = "stability_tabs", title = "Evaluación de Estabilidad y Precisión",
                       tabPanel("Case Dropping Bootstrap",
                                fluidRow(
                                  box(width = 6, status = "primary", solidHeader = TRUE, title = "Configuración Case Dropping",
                                      numericInput("case_boot","Bootstraps:",10,1,1000,1),
                                      numericInput("case_nCores","Núcleos:",4,1,max_cores,1),
                                      selectInput("case_statistics","Estadísticos:",
                                                  choices = list("Strength"="strength",
                                                                 "EI"="expectedInfluence",
                                                                 "Bridge"="bridgeStrength",
                                                                 "Strength+EI"="strength_ei",
                                                                 "Strength+Bridge"="strength_bridge",
                                                                 "EI+Bridge"="ei_bridge",
                                                                 "Todos"="all"),
                                                  selected="strength"),
                                      actionButton("run_case_boot","Ejecutar", icon=icon("play"), class="btn-primary")
                                  ),
                                  box(width = 6, status = "info", solidHeader = TRUE, title = "Estado",
                                      shinyWidgets::progressBar(id="case_progress", value=0, total=100, display_pct=TRUE),
                                      verbatimTextOutput("case_status"),
                                      verbatimTextOutput("case_time")
                                  )
                                ),
                                fluidRow(
                                  box(width=6, status="success", solidHeader=TRUE, title="Resultados", verbatimTextOutput("case_results")),
                                  box(width=6, status="success", solidHeader=TRUE, title="Gráfico", plotOutput("case_plot", height="400px"))
                                )
                       ),
                       tabPanel("Non-Parametric Bootstrap",
                                fluidRow(
                                  box(width=6, status="warning", solidHeader=TRUE, title="Configuración Non-Parametric",
                                      numericInput("nonparam_boot","Bootstraps:",10,1,1000,1),
                                      numericInput("nonparam_nCores","Núcleos:",4,1,max_cores,1),
                                      selectInput("nonparam_statistics","Estadísticos:",
                                                  choices=list("Todos"="all",
                                                               "Strength"="strength",
                                                               "EI"="expectedInfluence",
                                                               "Bridge"="bridgeStrength"),
                                                  selected="all"),
                                      actionButton("run_nonparam_boot","Ejecutar", icon=icon("play"), class="btn-warning")
                                  ),
                                  box(width=6, status="info", solidHeader=TRUE, title="Estado",
                                      shinyWidgets::progressBar(id="nonparam_progress", value=0, total=100, display_pct=TRUE),
                                      verbatimTextOutput("nonparam_status"),
                                      verbatimTextOutput("nonparam_time")
                                  )
                                ),
                                fluidRow(
                                  box(width=12, status="success", solidHeader=TRUE, title="Gráfico", plotOutput("nonparam_plot", height="500px"))
                                )
                       ),
                       tabPanel("Figura Combinada",
                                fluidRow(
                                  box(width=6, status="success", solidHeader=TRUE, title="Configuración Figura Combinada",
                                      selectInput("combined_statistics", "Estadísticos para combinar:",
                                                  choices=list("Strength"="strength",
                                                               "Expected Influence"="expectedInfluence",
                                                               "Bridge Strength"="bridgeStrength"),
                                                  selected="strength"),
                                      actionButton("generate_combined", "Generar Figura Combinada",
                                                   icon=icon("chart-line"), class="btn-success"),
                                      br(), br(),
                                      downloadButton("download_stability", "Descargar Figura", class="btn-success")
                                  ),
                                  box(width=6, status="info", solidHeader=TRUE, title="Información",
                                      p("Esta pestaña combina los resultados de Case Dropping y Non-Parametric Bootstrap."),
                                      p("Primero ejecute ambos tipos de bootstrap en las pestañas anteriores."),
                                      verbatimTextOutput("combined_status")
                                  )
                                ),
                                fluidRow(
                                  box(width=12, status="success", solidHeader=TRUE, title="Figura Combinada: Estabilidad y Precisión",
                                      plotOutput("combined_stability_plot", height="600px"))
                                )
                       )
                )
              )
      ),
      tabItem("groups",
              fluidRow(
                box(width=6, status="primary", solidHeader=TRUE, title="Configuración",
                    selectInput("group_var","Variable:", choices=NULL),
                    textInput("group1_name","Grupo 1:","Grupo 1"),
                    textInput("group2_name","Grupo 2:","Grupo 2"),
                    textInput("group_colors","Colores:","#FF5733,#33FFCE"),
                    actionButton("run_groups","Analizar", icon=icon("users"), class="btn-success")
                ),
                box(width=6, status="info", solidHeader=TRUE, title="Estado",
                    shinyWidgets::progressBar(id="groups_progress", value=0, total=100, display_pct=TRUE),
                    verbatimTextOutput("groups_status")
                )
              ),
              fluidRow(
                tabBox(width=12, id="groups_tabs", title="Grupos",
                       tabPanel("Redes",       plotOutput("networks_by_group_plot", height="600px")),
                       tabPanel("Centralidad", plotOutput("centrality_group_plot", height="500px")),
                       tabPanel("Bridge",      plotOutput("bridge_group_plot", height="500px")),
                       tabPanel("Combinada",
                                plotOutput("combined_groups_plot", height="700px"),
                                br(),
                                downloadButton("download_groups","Descargar", class="btn-success")
                       )
                )
              )
      )
    )
  )
)

# ---- 5. Server ----
server <- function(input, output, session) {
  rv <- reactiveValues(df_full = NULL, network_obj = NULL, groups = NULL,
                       errorMod = NULL, qgraph_obj = NULL,
                       caseDroppingBoot = NULL, nonParametricBoot = NULL,
                       combined_stability = NULL, combined_plot_final = NULL,
                       case_running = FALSE, nonparam_running = FALSE,
                       case_start_time = NULL, nonparam_start_time = NULL,
                       networks_groups = NULL, plot_centralidad_group = NULL,
                       bridge_plot_group = NULL, combined_groups_final = NULL,
                       errores_groups = NULL, pie_values = NULL, groups_running = FALSE)
  
  # Carga datos
  observeEvent(input$file, {
    df <- if (tools::file_ext(input$file$name)=="xlsx")
      read_excel(input$file$datapath)%>% slice(1:1000)
    else
      read.csv(input$file$datapath)%>% slice(1:1000)
    rv$df_full <- df
    updateSelectInput(session, "vars",
                      choices = names(df),
                      selected = names(df)[1:min(5,ncol(df))])
    fvs <- df %>% select_if(~ is.factor(.)||is.character(.)||
                              (is.numeric(.)&&length(unique(.))<=10)) %>% names()
    updateSelectInput(session, "group_var", choices = fvs)
  })
  
  # Datos filtrados
  df_vars <- eventReactive(input$run, {
    req(rv$df_full, input$vars)
    rv$df_full[, input$vars, drop = FALSE]
  })
  output$df_vars <- renderTable(df_vars())
  
  # Descriptivos y goldbricker
  output$desc_table   <- renderTable(psych::describe(df_vars()))
  output$goldbr       <- renderPrint(
    networktools::goldbricker(df_vars(), p=0.05,
                              method="hittner2003",
                              threshold=0.25,
                              corMin=0.50,
                              progressbar=FALSE)
  )
  output$edge_summary <- renderTable({
    req(rv$network_obj)
    InterconectaR::get_edge_weights_summary(rv$network_obj)
  })
  output$density_report <- renderPrint({
    req(rv$network_obj)
    InterconectaR::Density_report(rv$network_obj$graph)
  })
  
  # Estimar red y qgraph
  observeEvent(input$run, {
    req(df_vars())
    nm   <- strsplit(input$group_names, ",")[[1]]
    vals <- as.numeric(strsplit(input$group_values, ",")[[1]])
    rv$groups      <- InterconectaR::structure_groups(nm, vals)
    rv$network_obj <- estimateNetwork(df_vars(),
                                      default   = "ggmModSelect",
                                      stepwise  = TRUE,
                                      corMethod = "spearman")
    n <- ncol(df_vars())
    rv$errorMod   <- mgm_error_metrics(data=df_vars(),
                                       type = rep("g", n),
                                       level = rep(1, n))$errorCon$R2
    rv$qgraph_obj <- qgraph(rv$network_obj$graph,
                            groups      = rv$groups,
                            curveAll    = 2,
                            vsize       = 12,
                            esize       = 12,
                            palette     = "ggplot2",
                            layout      = "spring",
                            edge.labels = TRUE,
                            pie         = rv$errorMod,
                            layoutScale = c(0.8,0.8),
                            legend.cex  = 0.5,
                            labels      = abbreviate(names(df_vars()),3))
  })
  
  output$network_plot <- renderPlot({ req(rv$qgraph_obj); plot(rv$qgraph_obj) })
  
  # Centralidad
  cent <- eventReactive(input$run, {
    req(rv$qgraph_obj, rv$network_obj, rv$groups)
    m1 <- if (input$measure1=="None") NULL else input$measure1
    centrality_plots2_fixed(
      qgraph_obj    = rv$qgraph_obj,
      network       = rv$network_obj,
      groups        = rv$groups,
      measure0      = "ExpectedInfluence",
      measure1      = m1,
      color_palette = c("#F8766D","#00BFC4"),
      use_abbrev    = TRUE
    )
  })
  output$cent_table <- renderTable(cent()$table)
  output$cent_plot  <- renderPlot(print(cent()$plot))
  
  # Case Dropping Bootstrap
  observeEvent(input$run_case_boot, {
    req(rv$network_obj, rv$groups)
    if (rv$case_running) return()
    rv$case_running    <- TRUE
    rv$case_start_time <- Sys.time()
    updateProgressBar(session, "case_progress", value = 0)
    shinyjs::disable("run_case_boot")
    withProgress(message = "Ejecutando Case Dropping...", value = 0, {
      tryCatch({
        stats <- switch(input$case_statistics,
                        "strength"="strength",
                        "expectedInfluence"="expectedInfluence",
                        "bridgeStrength"="bridgeStrength",
                        "strength_ei"=c("strength","expectedInfluence"),
                        "strength_bridge"=c("strength","bridgeStrength"),
                        "ei_bridge"=c("expectedInfluence","bridgeStrength"),
                        "all"="all")
        cores <- min(input$case_nCores, max_cores)
        rv$caseDroppingBoot <- bootnet(
          rv$network_obj,
          boot       = input$case_boot,
          type       = "case",
          nCores     = cores,
          statistics = stats,
          communities= rv$groups,
          verbose    = FALSE
        )
        rv$case_running <- FALSE
        updateProgressBar(session, "case_progress", value = 100)
        shinyjs::enable("run_case_boot")
        showNotification("Case Dropping completado!", type = "message")
      }, error = function(e) {
        rv$case_running <- FALSE
        updateProgressBar(session, "case_progress", value = 0)
        shinyjs::enable("run_case_boot")
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  output$case_status  <- renderText({
    if (rv$case_running) "Case Dropping en progreso..."
    else if (!is.null(rv$caseDroppingBoot)) "Case Dropping completado"
    else "Presiona 'Ejecutar Case Dropping'"
  })
  output$case_time    <- renderText({
    req(rv$case_start_time)
    secs <- as.numeric(difftime(Sys.time(), rv$case_start_time, units="secs"))
    if (secs>60) paste0("Tiempo: ",round(secs/60,1)," min")
    else paste0("Tiempo: ",round(secs,0)," seg")
  })
  output$case_results <- renderPrint({
    req(rv$caseDroppingBoot)
    InterconectaR::filter_correlation_stability(rv$caseDroppingBoot)
  })
  output$case_plot    <- renderPlot({ req(rv$caseDroppingBoot); plot(rv$caseDroppingBoot, plot="area") })
  
  # Non-Parametric Bootstrap
  observeEvent(input$run_nonparam_boot, {
    req(rv$network_obj, rv$groups)
    if (rv$nonparam_running) return()
    rv$nonparam_running    <- TRUE
    rv$nonparam_start_time <- Sys.time()
    updateProgressBar(session, "nonparam_progress", value = 0)
    shinyjs::disable("run_nonparam_boot")
    withProgress(message = "Ejecutando Non-Parametric...", value = 0, {
      tryCatch({
        cores <- min(input$nonparam_nCores, max_cores)
        rv$nonParametricBoot <- bootnet(
          rv$network_obj,
          boot       = input$nonparam_boot,
          type       = "nonparametric",
          nCores     = cores,
          statistics = input$nonparam_statistics,
          communities= rv$groups,
          verbose    = FALSE
        )
        rv$nonparam_running <- FALSE
        updateProgressBar(session, "nonparam_progress", value = 100)
        shinyjs::enable("run_nonparam_boot")
        showNotification("Non-Parametric completado!", type = "message")
      }, error = function(e) {
        rv$nonparam_running <- FALSE
        updateProgressBar(session, "nonparam_progress", value = 0)
        shinyjs::enable("run_nonparam_boot")
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  output$nonparam_status <- renderText({
    if (rv$nonparam_running) "Non-Parametric en progreso..."
    else if (!is.null(rv$nonParametricBoot)) "Non-Parametric completado"
    else "Presiona 'Ejecutar Non-Parametric'"
  })
  output$nonparam_time   <- renderText({
    req(rv$nonparam_start_time)
    secs <- as.numeric(difftime(Sys.time(), rv$nonparam_start_time, units="secs"))
    if (secs>60) paste0("Tiempo: ",round(secs/60,1)," min")
    else paste0("Tiempo: ",round(secs,0)," seg")
  })
  output$nonparam_plot   <- renderPlot({ req(rv$nonParametricBoot); plot(rv$nonParametricBoot) })
  
  # Generar figura combinada (Estabilidad y Precisión)
  observeEvent(input$generate_combined, {
    req(rv$caseDroppingBoot, rv$nonParametricBoot)
    tryCatch({
      rv$combined_stability <- InterconectaR::plot_centrality_stability(
        rv$caseDroppingBoot,
        rv$nonParametricBoot,
        statistics = input$combined_statistics
      )
      showNotification("Figura combinada generada!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error al generar figura combinada:", e$message), type = "error")
    })
  })
  output$combined_status <- renderText({
    if (is.null(rv$caseDroppingBoot)) {
      "⚠️ Ejecute primero Case Dropping Bootstrap"
    } else if (is.null(rv$nonParametricBoot)) {
      "⚠️ Ejecute primero Non-Parametric Bootstrap"
    } else if (is.null(rv$combined_stability)) {
      "✅ Listo para generar figura combinada"
    } else {
      "✅ Figura combinada generada"
    }
  })
  output$combined_stability_plot <- renderPlot({ req(rv$combined_stability); print(rv$combined_stability) })
  output$download_stability <- downloadHandler(
    filename = function() paste0("Figura_2_Estabilidad_", Sys.Date(), ".jpg"),
    content = function(file) {
      req(rv$combined_stability)
      ggsave(file, rv$combined_stability, width = 9, height = 5, dpi = 600)
    }
  )
  
  # Figura Combinada en "Estimación de la red"
  output$combined_plot <- renderPlot({
    req(rv$network_obj, rv$groups, rv$errorMod, rv$qgraph_obj)
    tryCatch({
      m1        <- if (input$measure1=="None") NULL else input$measure1
      cent_data <- centrality_plots2_fixed(
        qgraph_obj    = rv$qgraph_obj,
        network       = rv$network_obj,
        groups        = rv$groups,
        measure0      = "ExpectedInfluence",
        measure1      = m1,
        color_palette = c("#F8766D","#00BFC4"),
        use_abbrev    = TRUE
      )
      p <- InterconectaR::combine_graphs_centrality2(
        Figura1_Derecha   = cent_data$plot,
        network           = rv$network_obj,
        groups            = rv$groups,
        error_Model       = rv$errorMod,
        ncol              = 2,
        widths            = c(0.5,0.25),
        dpi               = 300,
        legend.cex        = 0.35,
        abbreviate_labels = TRUE
      )
      print(p)
    }, error = function(e) {
      plot.new()
      text(0.5,0.5,paste("Error al generar figura combinada:\n",e$message),
           cex=1.2,col="red",adj=c(0.5,0.5))
    })
  },
  height = 600,
  width  = function() min(session$clientData$output_combined_plot_width, 1000)
  )
  
  
  output$downloadReport <- downloadHandler(
    filename = function() paste0("Figura_1_Final_", Sys.Date(), ".png"),
    content = function(file) {
      req(rv$network_obj, rv$groups, rv$errorMod, rv$qgraph_obj)
      cent_data <- centrality_plots2_fixed(
        qgraph_obj    = rv$qgraph_obj,
        network       = rv$network_obj,
        groups        = rv$groups,
        measure0      = "ExpectedInfluence",
        measure1      = if (input$measure1=="None") NULL else input$measure1,
        color_palette = c("#F8766D", "#00BFC4"),
        use_abbrev    = TRUE
      )
      png(file, width=1500, height=650, res=150)
      oldpar <- par(no.readonly=TRUE); on.exit({ dev.off(); par(oldpar) })
      layout(matrix(c(1,2),1,2), widths=c(2,1))
      par(mar=c(2,2,2,1)); plot(rv$qgraph_obj)
      par(mar=c(5,8,4,2))
      data_plot <- cent_data$table[order(cent_data$table$ExpectedInfluence), ]
      n <- nrow(data_plot)
      plot(1, type="n",
           xlim=range(data_plot$ExpectedInfluence, na.rm=TRUE)*c(1.1,1.1),
           ylim=c(0.5, n+0.5),
           xlab="z-score", ylab="", axes=FALSE)
      axis(1); axis(2, at=1:n, labels=data_plot$Abrev, las=1, cex.axis=0.9)
      points(data_plot$ExpectedInfluence, 1:n, pch=19, col="#F8766D", cex=1.5)
      if (!is.null(cent_data$table[[ input$measure1 ]])) {
        bridge_vals <- data_plot[[ input$measure1 ]]
        points(bridge_vals, 1:n, pch=19, col="#00BFC4", cex=1.5)
        for (i in 1:n) {
          lines(c(data_plot$ExpectedInfluence[i], bridge_vals[i]), c(i,i), col="gray60", lwd=1)
        }
        legend("topright", legend=c("Expected Influence", input$measure1),
               col=c("#F8766D","#00BFC4"), pch=19, bty="n", cex=0.9)
      }
      abline(h=1:n, col="gray90", lty=3); grid()
    }
  )
  
  # Análisis por grupos
  observeEvent(input$run_groups, {
    req(rv$df_full, input$group_var, input$vars)
    rv$groups_running <- TRUE
    updateProgressBar(session, "groups_progress", value = 0)
    withProgress(message="Analizando por grupos...", value=0, {
      tryCatch({
        df_for_groups <- rv$df_full %>%
          select(all_of(c(input$vars, input$group_var))) %>% 
          na.omit() %>%
          rename(.group = !!input$group_var)
        
        incProgress(0.1, detail="Estimando redes...")
        rv$networks_groups <- InterconectaR::estimate_networks_by_group(
          data      = df_for_groups,
          group_var = ".group",
          columns   = input$vars,
          default   = "ggmModSelect",
          stepwise  = TRUE,
          corMethod = "spearman",
          abbreviate_vars = TRUE,
          abbr_minlength   = 3
        )
        incProgress(0.3, detail="Centralidad por grupo...")
        colors <- strsplit(input$group_colors,",")[[1]] %>% trimws()
        if (length(colors)<2) colors <- c("#FF5733","#33FFCE")
        group_names <- c(input$group1_name,input$group2_name)
        rv$plot_centralidad_group <- InterconectaR::plot_centrality_by_group(
          rv$networks_groups,
          replacements = group_names,
          measure_spec = c("ExpectedInfluence"),
          color_palette = colors
        )
        incProgress(0.5, detail="Errores por grupo...")
        rv$errores_groups <- InterconectaR::mgm_errors_groups(
          data   = df_for_groups,
          type   = rep("g", length(input$vars)),
          level  = rep(1, length(input$vars)),
          group  = .group,
          columns= input$vars
        )
        rv$pie_values <- lapply(rv$errores_groups, `[[`, "R2")
        incProgress(0.7, detail="Redes por grupo...")
        rv$combined_plot_groups <- InterconectaR::plot_networks_by_group(
          res               = 300,
          networks_by_group = rv$networks_groups,
          groups            = rv$groups,
          pie               = rv$pie_values,
          legend.cex        = 0.8
        )
        incProgress(0.8, detail="Bridge por grupo...")
        rv$bridge_plot_group <- InterconectaR::centrality_bridge_plot(
          networks_groups = rv$networks_groups,
          group_names     = group_names,
          measure         = "Bridge Expected Influence (1-step)",
          color_palette   = colors
        )
        incProgress(0.9, detail="Combinando grupos...")
        rv$combined_groups_final <- InterconectaR::combine_groupBy(
          red_group              = rv$combined_plot_groups,
          plot_centralidad_group = rv$plot_centralidad_group,
          bridge_plot_group      = rv$bridge_plot_group$plot,
          width_a                = 12,
          width_bc               = 4.5,
          show_plot              = FALSE
        )
        incProgress(1, detail="¡Listo!")
        rv$groups_running <- FALSE
        updateProgressBar(session, "groups_progress", value = 100)
        showNotification("Análisis por grupos completado!", type="message")
      }, error = function(e) {
        rv$groups_running <- FALSE
        updateProgressBar(session, "groups_progress", value = 0)
        showNotification(paste("Error:",e$message), type="error")
      })
    })
  })
  output$groups_status <- renderText({
    if (rv$groups_running) "Analizando..."
    else if (!is.null(rv$networks_groups)) "Completado"
    else "Presiona 'Analizar'"
  })
  output$networks_by_group_plot   <- renderPlot({ req(rv$combined_plot_groups);   print(rv$combined_plot_groups) })
  output$centrality_group_plot    <- renderPlot({ req(rv$plot_centralidad_group); print(rv$plot_centralidad_group) })
  output$bridge_group_plot        <- renderPlot({ req(rv$bridge_plot_group);      print(rv$bridge_plot_group$plot) })
  output$combined_groups_plot     <- renderPlot({ req(rv$combined_groups_final);  print(rv$combined_groups_final) })
  
  output$download_groups <- downloadHandler(
    filename = function() paste0("grupos_",Sys.Date(),".jpg"),
    content  = function(file) {
      ggsave(file, rv$combined_groups_final, width=14, height=8, dpi=300)
    }
  )
}

# ---- 6. Lanzar la app ----
shinyApp(ui, server)