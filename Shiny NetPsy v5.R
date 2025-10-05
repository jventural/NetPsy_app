# app.R - NetPsy Professional con Diseño Mejorado
# PARTE 1: LIBRERÍAS Y FUNCIONES AUXILIARES

# ---- 1. Librerías ----
library(shiny)
library(shinydashboard)
library(shinyWidgets)    
library(shinyjs)         
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
library(tidyverse)       
library(forcats)         
library(tibble)          
library(parallel)        
library(cowplot)         
library(promises)        
library(future)          
library(NetworkComparisonTest)  
library(png)             
library(gridExtra)       
library(patchwork)       
library(scales)          

options(future.globals.maxSize = 891289600)  # ~850 MB
plan(multisession, workers = min(4, detectCores() - 1))

# ---- 2. Detectar núcleos físicos disponibles ----
max_cores <- detectCores(logical = FALSE)
if (is.na(max_cores) || max_cores < 1) max_cores <- 1

calc_auto_cores <- function(n_boot, max_cores) {
  cores <- floor(n_boot / 50)
  if (cores < 1) cores <- 1
  min(cores, max_cores)
}

auto_bootstrap_config <- function(n_vars, n_bootstraps = 50) {
  bootstraps <- as.numeric(n_bootstraps)
  
  return(list(
    case_boots = bootstraps,
    nonparam_boots = bootstraps,
    cores = min(detectCores(logical = FALSE) - 1, 4)
  ))
}

# ---- 3. Función centrality_plots2_fixed ----
centrality_plots2_fixed <- function(qgraph_obj,
                                    network,
                                    groups = NULL,
                                    measure0 = "ExpectedInfluence",
                                    measure1 = NULL,
                                    color_palette = c("#6C5CE7", "#00B894"),
                                    use_abbrev = TRUE) {
  library(qgraph); library(dplyr); library(networktools)
  library(tibble); library(ggplot2); library(tidyr)
  library(stringr); library(forcats); library(purrr)
  
  label0 <- measure0  
  label1 <- measure1  
  
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
  
  if (!is.null(measure1)) {
    
    if (is.null(groups)) {
      warning("No se pueden calcular medidas de puente sin definir grupos. Se omitirá measure1.")
      measure1 <- NULL
      label1 <- NULL
    } else {
      
      tryCatch({
        b_obj <- bridge(qgraph_obj,
                        communities    = groups,
                        useCommunities = "all",
                        normalize      = FALSE)
        
        if (is.null(b_obj)) {
          warning("La función bridge() devolvió NULL. Se omitirá measure1.")
          measure1 <- NULL
          label1 <- NULL
        } else if (!measure1 %in% names(b_obj)) {
          warning(paste("La medida", measure1, "no se encontró en los resultados de bridge. Medidas disponibles:", 
                        paste(names(b_obj), collapse = ", "), ". Se omitirá measure1."))
          measure1 <- NULL
          label1 <- NULL
        } else {
          
          bridge_values <- b_obj[[measure1]]
          
          if (is.null(bridge_values) || length(bridge_values) == 0) {
            warning(paste("La medida", measure1, "está vacía o es NULL. Se omitirá measure1."))
            measure1 <- NULL
            label1 <- NULL
          } else {
            
            qgraph_labels <- qgraph_obj$graphAttributes$Nodes$labels
            
            if (length(bridge_values) != length(qgraph_labels)) {
              warning(paste("Longitud de bridge_values (", length(bridge_values), 
                            ") no coincide con qgraph_labels (", length(qgraph_labels), 
                            "). Se omitirá measure1."))
              measure1 <- NULL
              label1 <- NULL
            } else {
              
              bridge_data <- tibble(
                qgraph_label = qgraph_labels,
                raw_bridge   = as.numeric(bridge_values)
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
                } else {
                  warning("No se pueden unir centralidad y puente por diferencia en número de filas. Se omitirá measure1.")
                  measure1 <- NULL
                  label1 <- NULL
                }
              }
            }
          }
        }
      }, error = function(e) {
        warning(paste("Error al calcular medidas de puente:", e$message, ". Se omitirá measure1."))
        measure1 <- NULL
        label1 <- NULL
      })
    }
  }
  
  if (!is.null(measure1) && !is.null(label1) && exists("cents2")) {
    y_var     <- if (use_abbrev) quo(Abrev) else quo(full_name)
    cents_long <- cents2 %>%
      pivot_longer(
        cols      = c(!!sym(measure0), !!sym(measure1)),
        names_to  = "Measure",
        values_to = "Value"
      )
    
    pal    <- setNames(color_palette, c(label0, label1))
    breaks <- c(label0, label1)
    labels <- c(label0, label1)
    
    Figura <- ggplot(cents_long,
                     aes(x = Value,
                         y = fct_reorder(!!y_var, Value),
                         color = Measure, group = Measure
                     )) +
      geom_point(size = 4, alpha = 0.8) +
      geom_line(size = 1.2, alpha = 0.7) +
      theme_minimal() +
      labs(x = "z-score", y = "Nodos", color = "Métrica") +
      scale_color_manual(
        values = pal,
        breaks = breaks,
        labels = labels
      ) +
      theme(
        axis.text.y  = element_text(size = 12, face = "bold"),
        axis.text.x  = element_text(size = 12),
        axis.title   = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_line(color = "grey95", size = 0.3),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    class(Figura) <- c("silent_gg", class(Figura))
    assign("print.silent_gg",
           function(x, ...) suppressWarnings(NextMethod()),
           envir = .GlobalEnv)
    
    return(list(
      table = cents2 %>% arrange(desc(!!sym(measure0))),
      plot  = Figura
    ))
  } else {
    cents_single <- cents_expect %>%
      select(full_name, Abrev, !!sym(measure0)) %>%
      arrange(desc(!!sym(measure0))) %>%
      rename(Value = !!sym(measure0))
    
    y_var <- if (use_abbrev) quo(Abrev) else quo(full_name)
    pal_single <- setNames(color_palette[1], label0)
    
    Figura <- ggplot(cents_single,
                     aes(x = Value,
                         y = fct_reorder(!!y_var, Value),
                         group = 1)) +
      geom_line(color = pal_single, size = 1.2, alpha = 0.7) +
      geom_point(color = pal_single, size = 4, alpha = 0.8) +
      theme_minimal() +
      labs(x = "z-score", y = "Nodos", color = "Métrica") +
      scale_color_manual(
        values = pal_single,
        breaks = label0,
        labels = label0
      ) +
      theme(
        axis.text.y  = element_text(size = 12, face = "bold"),
        axis.text.x  = element_text(size = 12),
        axis.title   = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_line(color = "grey95", size = 0.3),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    class(Figura) <- c("silent_gg", class(Figura))
    assign("print.silent_gg",
           function(x, ...) suppressWarnings(NextMethod()),
           envir = .GlobalEnv)
    
    return(list(table = cents_single, plot = Figura))
  }
}

# ---- 3.5 Funciones auxiliares ----
validate_groups <- function(group_names, group_values, n_vars) {
  nm   <- trimws(strsplit(group_names, ",")[[1]])
  vals_text <- trimws(strsplit(group_values, ",")[[1]])
  vals <- suppressWarnings(as.numeric(vals_text))
  
  errors <- c()
  
  if (any(is.na(vals))) {
    errors <- c(errors, "Los tamaños de grupos deben ser números válidos")
  }
  
  if (length(nm) != length(vals)) {
    errors <- c(errors, "El número de nombres debe coincidir con el número de tamaños")
  }
  
  if (sum(vals) != n_vars) {
    errors <- c(errors, paste("La suma de tamaños (", sum(vals), ") debe ser igual al número de variables (", n_vars, ")"))
  }
  
  if (any(vals <= 0)) {
    errors <- c(errors, "Todos los tamaños de grupos deben ser positivos")
  }
  
  return(list(valid = length(errors) == 0, errors = errors))
}

combine_network_centrality_custom <- function(qgraph_obj, 
                                              network_obj, 
                                              groups, 
                                              error_model, 
                                              centrality_plot,
                                              width_network = 0.65,
                                              width_centrality = 0.35,
                                              height_total = 600,
                                              dpi = 300) {
  
  library(cowplot)
  library(ggplot2)
  library(patchwork)
  
  tryCatch({
    temp_file <- tempfile(fileext = ".png")
    
    png(temp_file, width = 800, height = 600, res = 150)
    plot(qgraph_obj)
    dev.off()
    
    network_gg <- cowplot::ggdraw() + 
      cowplot::draw_image(temp_file) +
      ggtitle("A. Network Structure") +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#2D3748"),
            plot.margin = margin(5, 5, 5, 5))
    
    centrality_gg <- centrality_plot + 
      ggtitle("B. Node Centrality") +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#2D3748"),
            plot.margin = margin(5, 5, 5, 5),
            legend.position = "bottom")
    
    combined_plot <- network_gg + centrality_gg + 
      patchwork::plot_layout(
        ncol = 2, 
        widths = c(width_network, width_centrality)
      ) +
      patchwork::plot_annotation(
        title = "Network Analysis: Structure and Node Centrality",
        theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "#1A202C"))
      )
    
    unlink(temp_file)
    
    return(combined_plot)
    
  }, error = function(e) {
    tryCatch({
      library(gridExtra)
      library(grid)
      
      temp_file2 <- tempfile(fileext = ".png")
      png(temp_file2, width = 600, height = 400, res = 100)
      plot(qgraph_obj)
      dev.off()
      
      network_grob <- rasterGrob(readPNG(temp_file2))
      centrality_grob <- ggplotGrob(centrality_plot)
      
      combined_grid <- grid.arrange(
        network_grob, centrality_grob,
        ncol = 2,
        widths = c(width_network, width_centrality),
        top = "Network Analysis: Structure and Node Centrality"
      )
      
      unlink(temp_file2)
      return(combined_grid)
      
    }, error = function(e2) {
      warning("No se pudo combinar las figuras. Mostrando solo centralidad.")
      return(centrality_plot + 
               ggtitle("Node Centrality\n(Network plot available in 'Grafo' tab)") +
               theme(plot.title = element_text(hjust = 0.5, size = 14)))
    })
  })
}

combine_network_centrality_simple <- function(centrality_plot) {
  library(ggplot2)
  
  enhanced_plot <- centrality_plot + 
    ggtitle("Node Centrality Analysis") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "#1A202C"),
      plot.margin = margin(15, 15, 15, 15),
      legend.position = "bottom",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.3),
      panel.background = element_rect(fill = "white", color = "grey80"),
      plot.background = element_rect(fill = "#F7FAFC")
    ) +
    labs(
      subtitle = "Standardized centrality measures (z-scores)",
      caption = "Note: Network structure available in 'Grafo' tab"
    )
  
  return(enhanced_plot)
}

create_combined_figure <- function(qgraph_obj, network_obj, groups, error_model, cent_data_plot) {
  
  tryCatch({
    result <- combine_network_centrality_custom(
      qgraph_obj = qgraph_obj,
      network_obj = network_obj, 
      groups = groups,
      error_model = error_model,
      centrality_plot = cent_data_plot,
      width_network = 0.6,
      width_centrality = 0.4
    )
    return(result)
    
  }, error = function(e1) {
    tryCatch({
      result <- combine_network_centrality_simple(cent_data_plot)
      return(result)
      
    }, error = function(e2) {
      return(cent_data_plot + ggtitle("Node Centrality"))
    })
  })
}

# PARTE 2: CSS PERSONALIZADO Y UI COMPLETA

# ---- 4. CSS PERSONALIZADO PROFESIONAL ----
# ==== UI COMPLETO (incluye CSS y pestañas) ====

# ---- 4. CSS PERSONALIZADO PROFESIONAL ----
professional_css <- "
/* Variables de colores profesionales */
:root {
  --primary: #1e3a5f;      /* Azul marino profesional */
  --secondary: #4a7c7e;     /* Verde azulado */
  --accent: #e8b04b;        /* Dorado suave */
  --success: #5cb85c;
  --warning: #f0ad4e;
  --danger: #d9534f;
  --info: #5bc0de;
  --dark: #2c3e50;
  --light: #f8f9fa;
  --gray: #6c757d;
  --white: #ffffff;
  --shadow-sm: 0 1px 3px rgba(0,0,0,0.12);
  --shadow-md: 0 4px 6px rgba(0,0,0,0.16);
  --shadow-lg: 0 10px 20px rgba(0,0,0,0.19);
  --radius: 6px;
  --transition: all 0.3s ease;
}

/* Fuentes profesionales */
body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  font-size: 14px;
  line-height: 1.6;
  color: var(--dark);
  background: linear-gradient(135deg, #f5f7fa 0%, #e9ecef 100%);
}

/* Header mejorado */
.main-header {
  background: var(--primary) !important;
  box-shadow: var(--shadow-md) !important;
  border-bottom: 3px solid var(--accent) !important;
}

.main-header .logo {
  background: var(--primary) !important;
  font-weight: 600 !important;
  letter-spacing: 0.5px !important;
  text-transform: uppercase !important;
  font-size: 18px !important;
}

.main-header .navbar {
  background: var(--primary) !important;
}

/* Sidebar profesional */
.main-sidebar {
  background: var(--white) !important;
  box-shadow: var(--shadow-lg) !important;
}

.sidebar {
  padding-top: 10px !important;
}

.sidebar-menu > li {
  margin: 5px 10px !important;
}

.sidebar-menu > li > a {
  color: var(--dark) !important;
  border-radius: var(--radius) !important;
  padding: 12px 15px !important;
  transition: var(--transition) !important;
  border-left: 3px solid transparent !important;
  font-weight: 500 !important;
}

.sidebar-menu > li > a:hover {
  background: rgba(30, 58, 95, 0.1) !important;
  color: var(--primary) !important;
  border-left-color: var(--accent) !important;
}

.sidebar-menu > li.active > a {
  background: var(--primary) !important;
  color: var(--white) !important;
  border-left-color: var(--accent) !important;
  box-shadow: var(--shadow-sm) !important;
}

.sidebar-menu > li > a > .fa {
  width: 20px !important;
  text-align: center !important;
  margin-right: 10px !important;
}

/* Contenido principal */
.content-wrapper {
  background: var(--light) !important;
  padding: 20px !important;
  min-height: calc(100vh - 50px) !important;
}

/* Cajas mejoradas */
.box {
  border-radius: var(--radius) !important;
  box-shadow: var(--shadow-md) !important;
  border: none !important;
  margin-bottom: 20px !important;
  background: var(--white) !important;
}

.box-header {
  background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%) !important;
  color: var(--white) !important;
  padding: 12px 20px !important;
  border-radius: var(--radius) var(--radius) 0 0 !important;
}

.box-header .box-title {
  font-size: 16px !important;
  font-weight: 600 !important;
  letter-spacing: 0.3px !important;
}

.box-body {
  padding: 20px !important;
}

.box-footer {
  background: var(--light) !important;
  padding: 15px 20px !important;
  border-radius: 0 0 var(--radius) var(--radius) !important;
  border-top: 1px solid #e9ecef !important;
}

/* Tabs mejorados */
.nav-tabs-custom {
  box-shadow: var(--shadow-md) !important;
  border-radius: var(--radius) !important;
}

.nav-tabs-custom > .nav-tabs {
  background: var(--light) !important;
  border-bottom: 2px solid var(--primary) !important;
}

.nav-tabs-custom > .nav-tabs > li > a {
  color: var(--gray) !important;
  font-weight: 500 !important;
  padding: 10px 20px !important;
  transition: var(--transition) !important;
}

.nav-tabs-custom > .nav-tabs > li.active > a {
  color: var(--primary) !important;
  border-top: 3px solid var(--accent) !important;
  background: var(--white) !important;
}

.nav-tabs-custom > .tab-content {
  padding: 20px !important;
  background: var(--white) !important;
}

/* Botones profesionales */
.btn {
  border-radius: var(--radius) !important;
  padding: 8px 20px !important;
  font-weight: 500 !important;
  transition: var(--transition) !important;
  border: none !important;
  text-transform: uppercase !important;
  letter-spacing: 0.5px !important;
  font-size: 13px !important;
  box-shadow: var(--shadow-sm) !important;
}

.btn:hover {
  transform: translateY(-2px) !important;
  box-shadow: var(--shadow-md) !important;
}

.btn-primary { background: var(--primary) !important; color: var(--white) !important; }
.btn-success { background: var(--success) !important; color: var(--white) !important; }
.btn-info    { background: var(--info) !important;    color: var(--white) !important; }
.btn-warning { background: var(--warning) !important; color: var(--white) !important; }
.btn-danger  { background: var(--danger) !important;  color: var(--white) !important; }

.btn-lg { padding: 15px 30px !important; font-size: 16px !important; }

/* Inputs mejorados */
.form-control {
  border: 1px solid #dce0e4 !important;
  border-radius: var(--radius) !important;
  padding: 8px 12px !important;
  font-size: 14px !important;
  transition: var(--transition) !important;
}

.form-control:focus {
  border-color: var(--primary) !important;
  box-shadow: 0 0 0 0.2rem rgba(30, 58, 95, 0.25) !important;
  outline: none !important;
}

.form-group label {
  font-weight: 600 !important;
  color: var(--dark) !important;
  margin-bottom: 5px !important;
  font-size: 13px !important;
  text-transform: uppercase !important;
  letter-spacing: 0.3px !important;
}

/* Select múltiple mejorado */
select[multiple] { min-height: 150px !important; border-radius: var(--radius) !important; }
select[multiple] option { padding: 5px 10px !important; }
select[multiple] option:checked { background: var(--primary) !important; color: var(--white) !important; }

/* Progress bars */
.progress {
  height: 20px !important;
  border-radius: var(--radius) !important;
  background: var(--light) !important;
  box-shadow: inset 0 1px 2px rgba(0,0,0,0.1) !important;
}
.progress-bar {
  background: linear-gradient(90deg, var(--primary) 0%, var(--secondary) 100%) !important;
  line-height: 20px !important;
  font-weight: 600 !important;
  font-size: 12px !important;
}

/* Tablas profesionales */
.table { font-size: 14px !important; }
.table thead th {
  background: var(--primary) !important;
  color: var(--white) !important;
  font-weight: 600 !important;
  padding: 12px !important;
  border: none !important;
  text-transform: uppercase !important;
  font-size: 12px !important;
  letter-spacing: 0.5px !important;
}
.table tbody td { padding: 10px 12px !important; border-bottom: 1px solid #e9ecef !important; }
.table tbody tr:hover { background: var(--light) !important; }

/* Info cards */
.info-card {
  background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%) !important;
  color: var(--white) !important;
  padding: 20px !important;
  border-radius: var(--radius) !important;
  margin-bottom: 20px !important;
  box-shadow: var(--shadow-md) !important;
}
.info-card h4 { margin-top: 0 !important; font-weight: 600 !important; }
.info-card p { margin-bottom: 0 !important; opacity: 0.95 !important; }

/* Alert boxes */
.alert {
  border-radius: var(--radius) !important;
  border: none !important;
  box-shadow: var(--shadow-sm) !important;
  padding: 15px 20px !important;
}
.alert-info    { background: #e3f2fd !important; color: #0d47a1 !important; border-left: 4px solid #2196f3 !important; }
.alert-success { background: #e8f5e9 !important; color: #1b5e20 !important; border-left: 4px solid #4caf50 !important; }
.alert-warning { background: #fff3e0 !important; color: #e65100 !important; border-left: 4px solid #ff9800 !important; }
.alert-danger  { background: #ffebee !important; color: #b71c1c !important; border-left: 4px solid #f44336 !important; }

/* Stat cards */
.stat-card {
  background: var(--white) !important;
  padding: 20px !important;
  border-radius: var(--radius) !important;
  box-shadow: var(--shadow-sm) !important;
  text-align: center !important;
  transition: var(--transition) !important;
  border-top: 3px solid var(--accent) !important;
}
.stat-card:hover { transform: translateY(-3px) !important; box-shadow: var(--shadow-md) !important; }
.stat-number { font-size: 28px !important; font-weight: 700 !important; color: var(--primary) !important; margin: 10px 0 !important; }
.stat-label { font-size: 12px !important; text-transform: uppercase !important; color: var(--gray) !important; letter-spacing: 0.5px !important; font-weight: 600 !important; }

/* Responsive */
@media (max-width: 768px) {
  .box { margin: 10px 5px !important; }
  .box-body { padding: 15px !important; }
  .btn { width: 100% !important; margin-bottom: 10px !important; }
}
"

# ---- 5. UI PROFESIONAL COMPLETO ----
ui <- dashboardPage(
  dashboardHeader(
    title = tags$span(
      icon("brain"),
      "NETPSY PROFESSIONAL",
      style = "font-weight: 600; font-size: 18px;"
    )
  ),
  
  # --- Sidebar: reemplaza el item "Descriptivos" por dos items ---
  dashboardSidebar(
    sidebarMenu(
      menuItem("Configuración", tabName = "config", icon = icon("sliders-h"),
               badgeLabel = "INICIO", badgeColor = "green"),
      menuItem("Datos", tabName = "data", icon = icon("table")),
      menuItem("Pre-estimación", tabName = "pre", icon = icon("clipboard-check")),
      menuItem("Topología de la red", tabName = "topology", icon = icon("project-diagram")),
      menuItem("Estimación de la red", tabName = "report", icon = icon("diagram-project")),
      menuItem("Estabilidad y Precisión", tabName = "boot", icon = icon("sync-alt")),
      menuItem("Comparación Grupos", tabName = "groups", icon = icon("users"),
               badgeLabel = "FINAL", badgeColor = "red")
    )
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML(professional_css))),
    useShinyjs(),
    
    tabItems(
      # ================= CONFIGURACIÓN =================
      tabItem("config",
              h2("Configuración del Análisis", style = "color: #1e3a5f; margin-bottom: 20px;"),
              tabBox(
                width = 12,
                id = "config_subtabs",
                title = tags$span(icon("cogs"), "Configuración del Análisis"),
                
                # --- SUBTAB 1: Datos y Variables ---
                tabPanel("Datos y Variables",
                         icon = icon("database"),
                         fluidRow(
                           column(8,
                                  box(width = NULL, status = "primary", solidHeader = TRUE, 
                                      title = tags$span(icon("upload"), "Carga de Datos"),
                                      fileInput("file", 
                                                label = tags$div(
                                                  tags$h5("📁 Selecciona tu archivo", style = "margin-bottom: 10px;"),
                                                  tags$p("Formatos soportados: Excel (.xlsx) o CSV (.csv)", style = "color: #718096; font-size: 14px;")
                                                ), 
                                                accept = c(".xlsx",".csv"),
                                                buttonLabel = "Buscar archivo...",
                                                placeholder = "Ningún archivo seleccionado")
                                  ),
                                  box(width = NULL, status = "info", solidHeader = TRUE,
                                      title = tags$span(icon("list"), "Variables de Análisis"),
                                      selectInput("vars", 
                                                  label = tags$div(
                                                    tags$h5("Variables para el análisis", style = "margin-bottom: 10px;"),
                                                    tags$p("Selecciona las variables que deseas incluir en la red", style = "color: #718096; font-size: 14px;")
                                                  ),
                                                  choices = NULL, 
                                                  multiple = TRUE,
                                                  size = 8,
                                                  selectize = FALSE)
                                  )
                           ),
                           column(4,
                                  box(width = NULL, status = "success", solidHeader = TRUE,
                                      title = tags$span(icon("info-circle"), "Guía Rápida"),
                                      div(class = "stat-card",
                                          h5("📋 Pasos para comenzar:", style = "color: #2D3748; margin-bottom: 15px;"),
                                          tags$ol(
                                            tags$li("Carga tu archivo de datos"),
                                            tags$li("Selecciona las variables de interés"),
                                            tags$li("Configura los grupos en la siguiente pestaña"),
                                            tags$li("Ejecuta el análisis")
                                          ),
                                          br(),
                                          div(class = "alert alert-info",
                                              icon("lightbulb"), 
                                              "Asegúrate de que tu archivo contenga datos numéricos para el análisis de red."
                                          )
                                      )
                                  )
                           )
                         )
                ),
                
                # --- SUBTAB 2: Grupos y Configuración (con PARCHE de estimateNetwork) ---
                tabPanel("Grupos y Configuración",
                         icon = icon("layer-group"),
                         fluidRow(
                           column(8,
                                  # Configuración de Grupos
                                  box(width = NULL, status = "warning", solidHeader = TRUE,
                                      title = tags$span(icon("layer-group"), "Configuración de Grupos"),
                                      fluidRow(
                                        column(6,
                                               textInput("group_names", 
                                                         label = tags$div(
                                                           tags$h5("🏷️ Nombres de grupos", style = "margin-bottom: 5px;"),
                                                           tags$small("Separados por comas", style = "color: #718096;")
                                                         ),
                                                         value = "",
                                                         placeholder = "Dependencia,Miedo Soledad,Bienestar")
                                        ),
                                        column(6,
                                               textInput("group_values", 
                                                         label = tags$div(
                                                           tags$h5("📊 Tamaños de grupos", style = "margin-bottom: 5px;"),
                                                           tags$small("Números separados por comas", style = "color: #718096;")
                                                         ),
                                                         value = "",
                                                         placeholder = "3,1,1")
                                        )
                                      ),
                                      div(class = "alert alert-info",
                                          icon("info-circle"), 
                                          "La suma de tamaños debe coincidir con el número de variables seleccionadas."
                                      )
                                  ),
                                  
                                  # NUEVO: Parámetros de Estimación (bootnet::estimateNetwork)
                                  box(width = NULL, status = "primary", solidHeader = TRUE,
                                      title = tags$span(icon("project-diagram"), "Parámetros de Estimación de la Red"),
                                      fluidRow(
                                        column(4,
                                               selectInput(
                                                 "est_default",
                                                 label = tags$strong("Método (default)"),
                                                 choices = c("ggmModSelect", "EBICglasso", "pcor", "IsingFit", "IsingSampler", "mgm"),
                                                 selected = "ggmModSelect"
                                               )
                                        ),
                                        column(4,
                                               checkboxInput(
                                                 "est_stepwise",
                                                 label = tags$strong("Stepwise"),
                                                 value = TRUE
                                               )
                                        ),
                                        column(4,
                                               selectInput(
                                                 "est_cor",
                                                 label = tags$strong("Correlación (corMethod)"),
                                                 choices = c("cor", "cov", "cor_auto", "npn", "spearman"),
                                                 selected = "spearman"
                                               )
                                        )
                                      ),
                                      div(class = "alert alert-info",
                                          icon("lightbulb"),
                                          "Estos parámetros se aplican al presionar “EJECUTAR ANÁLISIS” y también en “Comparación Grupos”."
                                      )
                                  )
                           ),
                           
                           column(4,
                                  # Configuración Avanzada (centralidad/visual)
                                  box(width = NULL, status = "success", solidHeader = TRUE,
                                      title = tags$span(icon("cogs"), "Configuración Avanzada"),
                                      
                                      div(class = "stat-card",
                                          h5("Configuración de Centralidad", style = "color: #2D3748; margin-bottom: 15px;"),
                                          
                                          selectInput("measure0", 
                                                      label = tags$strong("Medida principal"),
                                                      choices = c("Expected Influence" = "ExpectedInfluence",
                                                                  "Strength" = "Strength"),
                                                      selected = "Strength"),
                                          
                                          selectInput("measure1", 
                                                      label = tags$strong("Medida puente (opcional)"),
                                                      choices = c("Ninguna" = "None",
                                                                  "Bridge Strength" = "Bridge Strength",
                                                                  "Bridge Expected Influence (1-step)" = "Bridge Expected Influence (1-step)"),
                                                      selected = "None"),
                                          
                                          div(class = "alert alert-info", style = "margin-top: 15px;",
                                              icon("lightbulb"), 
                                              "Las medidas de puente requieren al menos 2 grupos válidos."
                                          )
                                      ),
                                      
                                      div(class = "stat-card", style = "margin-top: 15px;",
                                          h5("Configuración Visual", style = "color: #2D3748; margin-bottom: 15px;"),
                                          
                                          selectInput("network_palette", 
                                                      label = tags$strong("Paleta de colores de la red"),
                                                      choices = c("GGplot2 (Moderno)" = "ggplot2",
                                                                  "Arcoíris" = "rainbow", 
                                                                  "Daltónicos" = "colorblind",
                                                                  "Pastel" = "pastel",
                                                                  "Grises" = "gray",
                                                                  "R Base" = "R"),
                                                      selected = "ggplot2"),
                                          
                                          div(class = "alert alert-info", style = "margin-top: 10px;",
                                              icon("palette"), 
                                              "La paleta 'Daltónicos' es más accesible para personas con daltonismo."
                                          )
                                      )
                                  )
                           )
                         ),
                         
                         # Botón de ejecución
                         fluidRow(
                           column(12,
                                  box(width = NULL, background = "navy",
                                      div(style = "text-align: center;",
                                          actionButton("run", 
                                                       label = tags$div(
                                                         icon("play", style = "font-size: 18px;"),
                                                         tags$strong("EJECUTAR ANÁLISIS", style = "font-size: 16px; margin-left: 10px;")
                                                       ),
                                                       class = "btn-success btn-lg",
                                                       style = "width: 300px; height: 60px;")
                                      )
                                  )
                           )
                         )
                )
              )
      ),
      
      # ================= DATOS =================
      tabItem("data",
              h2("Vista de Datos", style = "color: #1e3a5f; margin-bottom: 20px;"),
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE, 
                    title = tags$span(icon("database"), "Vista Previa de Datos"),
                    div(style = "overflow-x: auto;", tableOutput("df_vars"))
                )
              )
      ),
      
      # ================= DESCRIPTIVOS =================
      # --- Tab PRE-ESTIMACIÓN (antes era "desc") ---
      tabItem("pre",
              h2("Pre-estimación: inspección de datos", style = "color: #1e3a5f; margin-bottom: 20px;"),
              fluidRow(
                tabBox(width = 12, title = tags$span(icon("chart-line"), "Chequeos preliminares"),
                       tabPanel("Estadísticos",
                                icon = icon("calculator"),
                                div(style = "overflow-x: auto; margin-bottom: 10px;", tableOutput("desc_table")),
                                fluidRow(
                                  column(
                                    width = 12,
                                    div(class = "box-footer",
                                        div(style = "display: flex; gap: 12px; align-items: center; flex-wrap: wrap;",
                                            selectInput("desc_format", label = tags$strong("Formato de descarga"),
                                                        choices = c("CSV","XLSX"), selected = "CSV", width = "200px"),
                                            numericInput("desc_decimals", label = tags$strong("Decimales"),
                                                         value = 2, min = 0, max = 6, step = 1, width = "160px"),
                                            downloadButton("download_desc",
                                                           label = tags$div(icon("download"), "Descargar estadísticos"),
                                                           class = "btn-success")
                                        )
                                    )
                                  )
                                )
                       ),
                       tabPanel("Goldbricker", icon = icon("search"), verbatimTextOutput("goldbr"))
                )
              )
      ),
      
      # --- NUEVA Tab TOPOLOGÍA (mueve aquí Edges y Density) ---
      tabItem("topology",
              h2("Topología de la red", style = "color: #1e3a5f; margin-bottom: 20px;"),
              fluidRow(
                tabBox(width = 12, title = tags$span(icon("sitemap"), "Propiedades del grafo"),
                       tabPanel("Matriz de adyacencia", icon = icon("th"),
                                box(width = 12, status = "primary", solidHeader = TRUE, 
                                    title = tags$span(icon("th-large"), "Matriz de adyacencia (network$graph)"),
                                    div(style = "overflow-x: auto;", tableOutput("adjacency_matrix")),
                                    br(),
                                    div(style = "display: flex; gap: 12px; align-items: center; flex-wrap: wrap;",
                                        numericInput("adj_decimals", label = tags$strong("Decimales"),
                                                     value = 2, min = 0, max = 6, step = 1, width = "140px"),
                                        checkboxInput("adj_abs", label = tags$strong("Valores absolutos |w_ij|"), value = FALSE),
                                        # downloadButton("download_adj_csv",
                                        #                label = tags$div(icon("download"), "Descargar matriz (CSV)"),
                                        #                class = "btn-success"),
                                        downloadButton("download_adj_xlsx",
                                                       label = tags$div(icon("file-excel"), "Descargar matriz (XLSX)"),
                                                       class = "btn-success")
                                    )
                                )
                       ),
                       tabPanel("Edges", icon = icon("table"),
                                box(width = 12, status = "primary", solidHeader = TRUE,
                                    title = tags$span(icon("link"), "Resumen de aristas"),
                                    # --- Controles ---
                                    fluidRow(
                                      column(3,
                                             shinyWidgets::materialSwitch(
                                               inputId = "edges_abs",
                                               label   = tags$strong("Usar |w_ij| (abs_weights)"),
                                               value   = TRUE,        # igual a tu comportamiento actual
                                               status  = "primary",
                                               right   = TRUE
                                             )
                                      ),
                                      column(3,
                                             numericInput(
                                               "edge_decimals",
                                               label = tags$strong("Decimales"),
                                               value = 2, min = 0, max = 6, step = 1
                                             )
                                      ),
                                      column(6,
                                             div(style = "display:flex; gap:12px; justify-content:flex-end; align-items:end; height:100%;",
                                                 # downloadButton(
                                                 #   "download_edges_xlsx",
                                                 #   label = tags$div(icon("file-excel"), "Descargar (XLSX)"),
                                                 #   class = "btn-success"
                                                 # )
                                             )
                                      )
                                    ),
                                    br(),
                                    # --- Tabla ---
                                    div(style = "overflow-x: auto;", tableOutput("edge_summary"))
                                )
                       ),
                       tabPanel("Density", icon = icon("info-circle"),
                                verbatimTextOutput("density_report")
                       )
                )
              )
      ),
      
      # ================= ESTIMACIÓN DE LA RED =================
      tabItem("report",
              h2("Estimación de la Red", style = "color: #1e3a5f; margin-bottom: 20px;"),
              fluidRow(
                tabBox(width = 12, id = "network_tabs", 
                       title = tags$span(icon("project-diagram"), "Estimación de la Red"),
                       tabPanel("Grafo",
                                icon = icon("project-diagram"),
                                box(width = 12, status = "info", solidHeader = TRUE, 
                                    title = tags$span(icon("eye"), "Visualización de la Red"),
                                    plotOutput("network_plot", height = "600px")
                                )
                       ),
                       tabPanel("Centralidad",
                                icon = icon("bullseye"),
                                box(width = 12, status = "success", solidHeader = TRUE, 
                                    title = tags$span(icon("bullseye"), "Centralidad de Nodos"),
                                    div(style = "overflow-x: auto;", tableOutput("cent_table")),
                                    plotOutput("cent_plot", height = "400px")
                                )
                       ),
                       tabPanel("Figura Combinada",
                                icon = icon("palette"),
                                fluidRow(
                                  column(9,
                                         box(width = NULL, status = "primary", solidHeader = TRUE, 
                                             title = tags$span(icon("palette"), "Red + Centralidad"),
                                             plotOutput("combined_plot", height = "600px"),
                                             verbatimTextOutput("combined_plot_status")
                                         )
                                  ),
                                  column(3,
                                         box(width = NULL, status = "warning", solidHeader = TRUE,
                                             title = tags$span(icon("download"), "Configuración de Descarga"),
                                             div(class = "stat-card", div(class = "stat-number", "📐"), div(class = "stat-label", "Dimensiones")),
                                             numericInput("network_width", "Ancho (pulgadas):", value = 10, min = 5, max = 20, step = 0.5),
                                             numericInput("network_height","Alto (pulgadas):", value = 6, min = 3, max = 15, step = 0.5),
                                             numericInput("network_dpi",   "DPI:", value = 300, min = 150, max = 1200, step = 50),
                                             br(),
                                             downloadButton("downloadReport", label = tags$div(icon("download"), "Descargar Figura"),
                                                            class = "btn-success btn-lg", style = "width: 100%;"),
                                             br(), br(),
                                             div(class = "alert alert-info", icon("info"), "Formato: JPG de alta calidad"),
                                             br(),
                                             div(class = "stat-card",
                                                 h5("Estado:", style = "margin-bottom: 10px;"),
                                                 verbatimTextOutput("network_render_status")
                                             )
                                         )
                                  )
                                )
                       )
                )
              )
      ),
      
      # ================= BOOTSTRAP =================
      tabItem("boot",
              h2("Análisis de Estabilidad y Precisión", style = "color: #1e3a5f; margin-bottom: 20px;"),
              fluidRow(
                column(4,
                       box(width = NULL, status = "primary", solidHeader = TRUE, 
                           title = tags$span(icon("rocket"), "Configuración Rápida"),
                           selectInput("analysis_type", 
                                       label = tags$div(icon("clock"), tags$strong("Número de remuestreos:")),
                                       choices = list("Rápido (50 bootstraps)" = "50",
                                                      "Estándar (250 bootstraps)" = "250", 
                                                      "Completo (500 bootstraps)" = "500",
                                                      "Investigación (1000 bootstraps)" = "1000"),
                                       selected = "50"),
                           selectInput("bootstrap_focus", 
                                       label = tags$div(icon("crosshairs"), tags$strong("Enfoque del análisis:")),
                                       choices = list("Solo Estabilidad (Case Dropping)" = "case_only",
                                                      "Solo Precisión (Non-Parametric)" = "nonparam_only",
                                                      "Ambos (Recomendado)" = "both"),
                                       selected = "both"),
                           div(class = "stat-card", verbatimTextOutput("bootstrap_info")),
                           actionButton("run_smart_bootstrap", 
                                        label = tags$div(icon("rocket", style = "font-size: 16px;"), tags$strong("EJECUTAR ANÁLISIS")),
                                        class = "btn-primary btn-lg", style = "width: 100%; margin-top: 15px;"),
                           br(), br(),
                           tags$details(
                             tags$summary(tags$span(icon("cog"), "Configuración Avanzada", 
                                                    style = "cursor: pointer; color: #1e3a5f; font-weight: 600;")),
                             br(),
                             numericInput("custom_case_boot", "Case Dropping:", 50, 10, 1000, 10),
                             numericInput("custom_nonparam_boot", "Non-Parametric:", 50, 10, 1000, 10),
                             selectInput("custom_statistics", "Estadísticos:",
                                         choices = list("Strength" = "strength",
                                                        "Expected Influence" = "expectedInfluence",
                                                        "Bridge Strength" = "bridgeStrength"),
                                         selected = "strength", multiple = TRUE),
                             checkboxInput("boot_labels", "Mostrar etiquetas (labels) en los gráficos", TRUE),
                             checkboxInput("use_custom", "Usar configuración personalizada", FALSE)
                           )
                       )
                ),
                column(8,
                       box(width = NULL, status = "info", solidHeader = TRUE, 
                           title = tags$span(icon("chart-line"), "Estado del Análisis"),
                           fluidRow(
                             column(6,
                                    div(class = "stat-card",
                                        h4(icon("sync-alt"), "Case Dropping Bootstrap", style = "color: #1e3a5f;"),
                                        shinyWidgets::progressBar(id = "case_progress_new", value = 0, total = 100, 
                                                                  display_pct = TRUE, size = "md", status = "primary",
                                                                  striped = TRUE, title = "Progreso"),
                                        verbatimTextOutput("case_status_new")
                                    )
                             ),
                             column(6,
                                    div(class = "stat-card",
                                        h4(icon("chart-bar"), "Non-Parametric Bootstrap", style = "color: #f0ad4e;"),
                                        shinyWidgets::progressBar(id = "nonparam_progress_new", value = 0, total = 100,
                                                                  display_pct = TRUE, size = "md", status = "warning",
                                                                  striped = TRUE, title = "Progreso"),
                                        verbatimTextOutput("nonparam_status_new")
                                    )
                             )
                           ),
                           hr(style = "border-color: #e2e8f0; margin: 20px 0;"),
                           fluidRow(
                             column(4, div(class = "stat-card",
                                           div(class = "stat-number", icon("clock", style = "font-size: 30px; color: #1e3a5f;")),
                                           div(class = "stat-label", tags$strong("Tiempo Estimado")),
                                           verbatimTextOutput("time_estimate_text"))),
                             column(4, div(class = "stat-card",
                                           div(class = "stat-number", icon("microchip", style = "font-size: 30px; color: #5cb85c;")),
                                           div(class = "stat-label", tags$strong("Recursos")),
                                           verbatimTextOutput("cores_info_text"))),
                             column(4, div(class = "stat-card",
                                           div(class = "stat-number", icon("tasks", style = "font-size: 30px; color: #d9534f;")),
                                           div(class = "stat-label", tags$strong("Estado General")),
                                           verbatimTextOutput("overall_status_text")))
                           )
                       )
                )
              ),
              fluidRow(
                tabBox(width = 12, id = "bootstrap_results", 
                       title = tags$span(icon("chart-area"), "Resultados del Análisis"),
                       tabPanel("Análisis Detallado",
                                icon = icon("microscope"),
                                fluidRow(
                                  column(5,
                                         box(width = 12, status = "primary", solidHeader = TRUE,
                                             title = tags$span(icon("bullseye"), "Estabilidad (Case Dropping)"),
                                             verbatimTextOutput("case_detailed"),
                                             plotOutput("case_plot_detailed", height = "400px")
                                         )
                                  ),
                                  column(7,
                                         box(width = 12, status = "warning", solidHeader = TRUE, 
                                             title = tags$span(icon("ruler"), "Precisión (Non-Parametric)"),
                                             plotOutput("nonparam_plot_detailed", height = "500px")
                                         )
                                  )
                                ),
                                # --- RESUMEN CS (dentro de Análisis Detallado) ---
                                fluidRow(
                                  column(12,
                                         box(width = NULL, status = "success", solidHeader = TRUE,
                                             title = tags$span(icon("shield-alt"), "Índice de Estabilidad de Correlación (CS)"),
                                             # Tarjetas-resumen (semáforo)
                                             uiOutput("cs_cards"),
                                             br(),
                                             # Tabla con interpretación
                                             div(style = "overflow-x: auto;", tableOutput("cs_table")),
                                             br(),
                                             tags$small(
                                               HTML("Criterios: <b>Excelente</b> ≥ 0.75, <b>Adecuado</b> 0.50–0.74, <b>Bajo</b> &lt; 0.50.")
                                             ),
                                             br(),
                                             downloadButton("download_cs_csv",
                                                            label = tags$div(icon("download"), "Descargar CS (CSV)"),
                                                            class = "btn-success")
                                         )
                                  )
                                )
                                
                       ),
                       tabPanel("Figura Combinada",
                                icon = icon("palette"),
                                fluidRow(
                                  column(9,
                                         box(width = NULL, status = "success", solidHeader = TRUE,
                                             title = tags$span(icon("palette"), "Figura Final para Publicación"),
                                             plotOutput("combined_bootstrap_plot", height = "600px")
                                         )
                                  ),
                                  column(3,
                                         box(width = NULL, status = "info", solidHeader = TRUE,
                                             title = tags$span(icon("download"), "Configuración de Descarga"),
                                             div(class = "stat-card", div(class = "stat-number", "📐"), div(class = "stat-label", "Dimensiones de la Figura")),
                                             numericInput("bootstrap_width", "Ancho (pulgadas):",  value = 9, min = 3, max = 20, step = 0.5),
                                             numericInput("bootstrap_height","Alto (pulgadas):",   value = 5, min = 3, max = 15, step = 0.5),
                                             numericInput("bootstrap_dpi",   "DPI:",               value = 600, min = 150, max = 1200, step = 50),
                                             br(),
                                             downloadButton("download_bootstrap_combined", 
                                                            label = tags$div(icon("download"), "Descargar Figura"),
                                                            class = "btn-success btn-lg", style = "width: 100%;"),
                                             br(), br(),
                                             div(class = "alert alert-info", icon("info"), "Formato: JPG de alta calidad")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Case Dropping (tabla)",
                                icon = icon("table"),
                                fluidRow(
                                  box(width = 12, status = "primary", solidHeader = TRUE,
                                      title = tags$span(icon("table"), "Resumen de correlaciones boot vs. original"),
                                      # Controles
                                      fluidRow(
                                        column(3,
                                               selectInput(
                                                 "cd_type_filter",
                                                 label = tags$strong("Estadístico"),
                                                 choices = c("Todos", "strength", "expectedInfluence", "bridgeStrength"),
                                                 selected = "Todos"
                                               )
                                        ),
                                        column(3,
                                               numericInput(
                                                 "cd_decimals",
                                                 label = tags$strong("Decimales"),
                                                 value = 3, min = 0, max = 6, step = 1
                                               )
                                        ),
                                        column(6,
                                               div(style = "display: flex; gap: 12px; justify-content: flex-end; flex-wrap: wrap; margin-top: 28px;",
                                                   # downloadButton(
                                                   #   "download_cd_csv",
                                                   #   label = tags$div(icon("download"), "Descargar (CSV)"),
                                                   #   class = "btn-success"
                                                   # ),
                                                   downloadButton(
                                                     "download_cd_xlsx",
                                                     label = tags$div(icon("file-excel"), "Descargar (XLSX)"),
                                                     class = "btn-success"
                                                   )
                                               )
                                        )
                                      ),
                                      br(),
                                      div(style = "overflow-x: auto;", tableOutput("cd_table")),
                                      br(),
                                      tags$small(
                                        HTML("La tabla resume la correlación entre cada estadístico/nodo estimado con muestra completa y los valores re-estimados al reducir el tamaño muestral (Case Dropping). <br>
             Columnas: <i>type</i> (estadístico), <i>nPerson</i> (tamaño de submuestra), <i>mean_cor</i> (correlación media),
             <i>lo/hi</i> (IC 95% por percentiles), <i>n_boot</i> (repeticiones), <i>p_sampled</i> (proporción de muestra).")
                                      )
                                  )
                                )
                       ),
                       tabPanel("Bootstrap no paramétrico (aristas)",
                                icon = icon("link"),
                                fluidRow(
                                  box(width = 12, status = "primary", solidHeader = TRUE,
                                      title = tags$span(icon("link"), "Resumen de aristas (nonparametric bootstrap)"),
                                      # Controles
                                      fluidRow(
                                        column(3,
                                               checkboxInput(
                                                 "np_only_sig",
                                                 label = tags$strong("Solo aristas con IC95% que excluye 0"),
                                                 value = FALSE
                                               )
                                        ),
                                        column(3,
                                               selectInput(
                                                 "np_order",
                                                 label = tags$strong("Orden"),
                                                 choices = c(
                                                   "Por |media bootstrap| (desc)" = "mag",
                                                   "Por nombre de arista (A–B)"  = "name",
                                                   "IC: NO incluye 0 primero"    = "zerofirst"
                                                 ),
                                                 selected = "mag"
                                               )
                                        ),
                                        column(3,
                                               numericInput(
                                                 "np_decimals",
                                                 label = tags$strong("Decimales"),
                                                 value = 3, min = 0, max = 6, step = 1
                                               )
                                        ),
                                        column(3,
                                               div(style = "display: flex; gap: 12px; justify-content: flex-end; flex-wrap: wrap; margin-top: 28px;",
                                                   # downloadButton(
                                                   #   "download_np_csv",
                                                   #   label = tags$div(icon("download"), "Descargar (CSV)"),
                                                   #   class = "btn-success"
                                                   # ),
                                                   downloadButton(
                                                     "download_np_xlsx",
                                                     label = tags$div(icon("file-excel"), "Descargar (XLSX)"),
                                                     class = "btn-success"
                                                   )
                                               )
                                        )
                                      ),
                                      br(),
                                      div(style = "overflow-x: auto;", tableOutput("np_edges_table")),
                                      br(),
                                      tags$small(
                                        HTML("Fuente: <code>InterconectaR::summarise_nonparametric_edges()</code>. 
             Columnas: <i>id</i> (A--B), <i>sample</i> (peso en muestra completa), 
             <i>boot_mean</i> (media bootstrap), <i>ci_lo</i>/<i>ci_hi</i> (IC95%), 
             <i>n_boot</i> (repeticiones), <i>ci_covers_zero</i> (¿IC incluye 0?).")
                                      )
                                  )
                                )
                       )
                       
                       
                )
              )
      ),
      
      # ================= COMPARACIÓN DE GRUPOS =================
      tabItem("groups",
              h2("Comparación de Grupos", style = "color: #1e3a5f; margin-bottom: 20px;"),
              fluidRow(
                column(6,
                       box(width = NULL, status = "primary", solidHeader = TRUE, 
                           title = tags$span(icon("cogs"), "Configuración"),
                           selectInput("group_var", 
                                       label = tags$div(icon("tag"), tags$strong("Variable de agrupación:")),
                                       choices = NULL),
                           fluidRow(
                             column(6, textInput("group1_name", label = tags$strong("Grupo 1:"), value = "Grupo 1")),
                             column(6, textInput("group2_name", label = tags$strong("Grupo 2:"), value = "Grupo 2"))
                           ),
                           textInput("group_colors", 
                                     label = tags$div(icon("palette"), tags$strong("Colores (separados por coma):")),
                                     value = "",  # sin valor por defecto
                                     placeholder = "#1e3a5f, #CC1228"),
                           hr(),
                           div(class = "stat-card",
                               h5(icon("microscope"), "Network Comparison Test (NCT)", style = "color: #2D3748;"),
                               numericInput("nct_iterations", "Iteraciones NCT:", value = 100, min = 50, max = 1000, step = 50),
                               div(class = "alert alert-info",
                                   icon("clock"), "Más iteraciones = mayor precisión pero más tiempo"
                               )
                           ),
                           actionButton("run_groups",
                                        label = tags$div(icon("users"), tags$strong("ANALIZAR GRUPOS")),
                                        class = "btn-primary btn-lg", style = "width: 100%; margin-top: 15px;")
                       )
                ),
                column(6,
                       box(width = NULL, status = "info", solidHeader = TRUE, 
                           title = tags$span(icon("chart-line"), "Estado del Análisis"),
                           div(class = "stat-card",
                               h5("Progreso General", style = "margin-bottom: 15px;"),
                               shinyWidgets::progressBar(id="groups_progress", value=0, total=100, display_pct=TRUE),
                               verbatimTextOutput("groups_status")
                           ),
                           hr(),
                           div(class = "stat-card",
                               h5(icon("microscope"), "NCT Status", style = "color: #f0ad4e; margin-bottom: 15px;"),
                               shinyWidgets::progressBar(id="nct_progress", value=0, total=100, display_pct=TRUE, status="warning"),
                               verbatimTextOutput("nct_status")
                           )
                       )
                )
              ),
              fluidRow(
                tabBox(width=12, id="groups_tabs", 
                       title = tags$span(icon("users"), "Análisis por Grupos"),
                       tabPanel("Redes", icon = icon("project-diagram"), plotOutput("networks_by_group_plot", height="600px")),
                       tabPanel("Centralidad", icon = icon("bullseye"), plotOutput("centrality_group_plot", height="500px")),
                       tabPanel("Bridge", icon = icon("bridge"), plotOutput("bridge_group_plot", height="500px")),
                       tabPanel("NCT - Comparación Estadística",
                                icon = icon("microscope"),
                                fluidRow(
                                  column(6,
                                         box(width = NULL, status = "primary", solidHeader = TRUE,
                                             title = tags$span(icon("microscope"), "Network Comparison Test - Resultados"),
                                             verbatimTextOutput("nct_results"),
                                             br(),
                                             div(class = "alert alert-info",
                                                 h5("📋 Interpretación:", style = "margin-bottom: 10px;"),
                                                 tags$ul(
                                                   tags$li("Network Invariance: ¿Son las redes globalmente diferentes?"),
                                                   tags$li("Global Strength: ¿Difiere la conectividad total?"),
                                                   tags$li("p < 0.05: Diferencias significativas"),
                                                   tags$li("p ≥ 0.05: No hay diferencias significativas")
                                                 )
                                             )
                                         )
                                  ),
                                  column(6,
                                         box(width = NULL, status = "info", solidHeader = TRUE,
                                             title = tags$span(icon("chart-bar"), "Resumen Visual NCT"),
                                             plotOutput("nct_plot", height = "400px"),
                                             br(),
                                             downloadButton("download_nct", 
                                                            label = tags$div(icon("download"), "Descargar NCT"),
                                                            class = "btn-info", style = "width: 100%;")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Combinada",
                                icon = icon("palette"),
                                fluidRow(
                                  column(9,
                                         box(width = NULL, status = "success", solidHeader = TRUE, 
                                             title = tags$span(icon("palette"), "Comparación Grupos Completa"),
                                             plotOutput("combined_groups_plot", height = "700px")
                                         )
                                  ),
                                  column(3,
                                         box(width = NULL, status = "info", solidHeader = TRUE,
                                             title = tags$span(icon("download"), "Configuración de Descarga"),
                                             div(class = "stat-card", div(class = "stat-number", "📐"), div(class = "stat-label", "Dimensiones de la Figura")),
                                             numericInput("groups_width",  "Ancho total (pulgadas):", value = 14, min = 8, max = 25, step = 0.5),
                                             numericInput("groups_height", "Alto (pulgadas):",        value = 8,  min = 4, max = 20, step = 0.5),
                                             numericInput("groups_dpi",    "DPI:",                    value = 300,min = 150, max = 1200, step = 50),
                                             hr(),
                                             div(class = "stat-card",
                                                 h5("📊 Proporciones de Paneles", style = "margin-bottom: 10px;"),
                                                 numericInput("groups_width_a",  "Ancho Panel A (Redes):", value = 15,  min = 8, max = 25, step = 0.5),
                                                 numericInput("groups_width_bc", "Ancho Paneles B+C:",     value = 4.5, min = 2, max = 10, step = 0.1),
                                                 div(class = "alert alert-info",
                                                     "Panel A: Redes por grupo", br(),
                                                     "Paneles B+C: Centralidad y Bridge"
                                                 )
                                             ),
                                             br(),
                                             downloadButton("download_groups", 
                                                            label = tags$div(icon("download"), "Descargar Comparación"),
                                                            class = "btn-success btn-lg", style = "width: 100%;"),
                                             br(), br(),
                                             div(class = "alert alert-info", icon("info"), "Formato: JPG de alta calidad"),
                                             br(),
                                             div(class = "stat-card",
                                                 h5("Estado:", style = "margin-bottom: 10px;"),
                                                 verbatimTextOutput("groups_render_status")
                                             )
                                         )
                                  )
                                )
                       ),
                       # En el tabItem("groups"), dentro del tabBox, agrega esta pestaña después de "NCT - Comparación Estadística":
                       
                       tabPanel("Tamaño del Efecto",
                                icon = icon("ruler-combined"),
                                fluidRow(
                                  column(12,
                                         box(width = NULL, status = "primary", solidHeader = TRUE,
                                             title = tags$span(icon("calculator"), "Métricas de Tamaño del Efecto"),
                                             
                                             # Tarjetas con métricas principales
                                             fluidRow(
                                               column(4,
                                                      div(class = "stat-card",
                                                          div(class = "stat-number", textOutput("effect_size_correlation")),
                                                          div(class = "stat-label", "Correlación de Spearman"),
                                                          tags$small("Similitud entre matrices de adyacencia")
                                                      )
                                               ),
                                               column(4,
                                                      div(class = "stat-card",
                                                          div(class = "stat-number", textOutput("effect_size_difference")),
                                                          div(class = "stat-label", "Diferencia Absoluta Media"),
                                                          tags$small("Promedio de diferencias |w₁ - w₂|")
                                                      )
                                               ),
                                               column(4,
                                                      div(class = "stat-card",
                                                          div(class = "stat-number", textOutput("effect_size_max_diff")),
                                                          div(class = "stat-label", "Diferencia Máxima"),
                                                          tags$small("Mayor diferencia entre aristas")
                                                      )
                                               )
                                             ),
                                             
                                             br(),
                                             
                                             # Interpretación y detalles
                                             box(width = NULL, status = "info",
                                                 title = tags$span(icon("info-circle"), "Interpretación Detallada"),
                                                 verbatimTextOutput("effect_size_interpretation"),
                                                 br(),
                                                 div(class = "alert alert-info",
                                                     h5("📊 Guía de Interpretación:", style = "margin-bottom: 10px;"),
                                                     tags$ul(
                                                       tags$li("Correlación > 0.90: Redes muy similares"),
                                                       tags$li("Correlación 0.70-0.90: Similitud moderada-alta"),
                                                       tags$li("Correlación 0.50-0.70: Similitud moderada"),
                                                       tags$li("Correlación < 0.50: Redes sustancialmente diferentes")
                                                     )
                                                 )
                                             ),
                                             
                                             # Gráfico de comparación
                                             box(width = NULL, status = "success",
                                                 title = tags$span(icon("chart-bar"), "Visualización de Diferencias"),
                                                 plotOutput("effect_size_plot", height = "400px"),
                                                 br(),
                                                 downloadButton("download_effect_size",
                                                                label = tags$div(icon("download"), "Descargar Análisis"),
                                                                class = "btn-success", style = "width: 100%;")
                                             )
                                         )
                                  )
                                )
                       )
                )
              )
      )
    )
  )
)

# Cierre de dashboardPage
# PARTE 4: SERVER - INICIO Y BOOTSTRAP CORREGIDO

# ---- 6. Server ----
server <- function(input, output, session) {
  rv <- reactiveValues(df_full = NULL, network_obj = NULL, groups = NULL,
                       errorMod = NULL, qgraph_obj = NULL,
                       caseDroppingBoot = NULL, nonParametricBoot = NULL,
                       combined_stability = NULL, combined_plot_final = NULL,
                       case_running = FALSE, nonparam_running = FALSE,
                       case_start_time = NULL, nonparam_start_time = NULL,
                       networks_groups = NULL, plot_centralidad_group = NULL,
                       bridge_plot_group = NULL, combined_groups_final = NULL,
                       errores_groups = NULL, pie_values = NULL, groups_running = FALSE,
                       nct_result = NULL, nct_running = FALSE, nct_start_time = NULL,
                       effect_size_corr = NULL,
                       effect_size_mean_diff = NULL,
                       effect_size_max_diff = NULL,
                       effect_size_n_diff = NULL,
                       effect_size_r_squared = NULL,
                       effect_size_diff_matrix = NULL)
  
  rv_bootstrap <- reactiveValues(
    case_running = FALSE,
    nonparam_running = FALSE, 
    both_running = FALSE,
    case_start_time = NULL,
    nonparam_start_time = NULL,
    overall_start_time = NULL,
    case_result = NULL,
    nonparam_result = NULL,
    config = NULL,
    selected_statistics = NULL  # AGREGADO para guardar estadísticos seleccionados
  )
  
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
  
  df_vars <- eventReactive(input$run, {
    req(rv$df_full, input$vars)
    rv$df_full[, input$vars, drop = FALSE]
  })
  output$df_vars <- renderTable(df_vars())
  
  desc_data <- reactive({
    req(df_vars())
    psych::describe(df_vars()) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Variable") %>%
      dplyr::select(
        Variable,
        Mean      = mean,
        sd        = sd,
        Min       = min,
        Max       = max,
        g1      = skew,
        g2  = kurtosis
      )
  })
  
  # ===== Matriz de adyacencia (network$graph) =====
  adj_matrix_data <- reactive({
    req(rv$network_obj)
    M <- rv$network_obj$graph          # equivalente a network$graph
    req(M)
    if (isTRUE(input$adj_abs)) M <- abs(M)
    dec <- if (is.null(input$adj_decimals)) 2 else input$adj_decimals
    M <- round(M, dec)
    df <- as.data.frame(M)
    df <- tibble::rownames_to_column(df, var = "Node")
    df
  })
  
  output$adjacency_matrix <- renderTable({
    adj_matrix_data()
  })
  
  # ---- Descarga CSV ----
  output$download_adj_csv <- downloadHandler(
    filename = function() paste0("Matriz_Adyacencia_", Sys.Date(), ".csv"),
    content  = function(file) {
      utils::write.csv(adj_matrix_data(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
  
  # ---- Descarga XLSX (con fallback) ----
  output$download_adj_xlsx <- downloadHandler(
    filename = function() paste0("Matriz_Adyacencia_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      d <- adj_matrix_data()
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(d, path = file)
      } else if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(d, file)
      } else {
        showNotification("No se encontró 'writexl' ni 'openxlsx'. Se exportará en CSV.",
                         type = "warning", duration = 5)
        utils::write.csv(d, file, row.names = FALSE, fileEncoding = "UTF-8")
      }
    }
  )
  
  
  output$desc_table <- renderTable({
    d <- desc_data()
    dec <- if (is.null(input$desc_decimals)) 3 else input$desc_decimals
    d %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, dec)))
  })
  
  output$download_desc <- downloadHandler(
    filename = function() {
      base <- paste0("Descriptivos_", Sys.Date())
      if (is.null(input$desc_format) || input$desc_format == "CSV") {
        paste0(base, ".csv")
      } else {
        paste0(base, ".xlsx")
      }
    },
    content = function(file) {
      d <- desc_data()
      dec <- if (is.null(input$desc_decimals)) 2 else input$desc_decimals
      d_out <- d %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, dec)))
      
      fmt <- if (is.null(input$desc_format)) "CSV" else input$desc_format
      
      if (fmt == "XLSX") {
        # Usa writexl si está disponible, si no, fallback a CSV
        if (requireNamespace("writexl", quietly = TRUE)) {
          writexl::write_xlsx(d_out, path = file)
        } else {
          # Fallback: CSV y aviso
          showNotification("No se encontró 'writexl'. Se exportará en CSV.", type = "warning", duration = 5)
          utils::write.csv(d_out, file = file, row.names = FALSE, fileEncoding = "UTF-8")
        }
      } else {
        utils::write.csv(d_out, file = file, row.names = FALSE, fileEncoding = "UTF-8")
      }
    }
  )
  
  
  output$goldbr       <- renderPrint(
    networktools::goldbricker(df_vars(), p=0.05,
                              method="hittner2003",
                              threshold=0.25,
                              corMin=0.50,
                              progressbar=FALSE)
  )
  # Datos para la tabla de aristas (reutilizable)
  edge_summary_data <- reactive({
    req(rv$network_obj)
    abs_arg <- isTRUE(input$edges_abs)                     # <- toggle del usuario
    dec     <- if (is.null(input$edge_decimals)) 2 else input$edge_decimals
    InterconectaR::get_edge_weights_summary(
      rv$network_obj,
      abs_weights  = abs_arg,
      round_digits = dec
    )
  })
  
  # Render de la tabla
  output$edge_summary <- renderTable({
    edge_summary_data()
  })
  
  output$density_report <- renderPrint({
    req(rv$network_obj)
    InterconectaR::Density_report(rv$network_obj$graph)
  })
  
  observeEvent(input$run, {
    req(df_vars())
    
    tryCatch({
      validation <- validate_groups(input$group_names, input$group_values, length(input$vars))
      
      if (!validation$valid) {
        showNotification(paste("Error en configuración de grupos:", paste(validation$errors, collapse = "; ")), 
                         type = "error", duration = 10)
        return()
      }
      
      nm   <- trimws(strsplit(input$group_names, ",")[[1]])
      vals <- as.numeric(trimws(strsplit(input$group_values, ",")[[1]]))
      # Vector de comunidades para bootnet (uno por variable, en el orden de las columnas)
rv$communities_vec <- rep(seq_along(vals), times = vals)
names(rv$communities_vec) <- colnames(df_vars())

      rv$groups      <- InterconectaR::structure_groups(nm, vals)
      rv$network_obj <- estimateNetwork(df_vars(),
                                        default   = "ggmModSelect",
                                        stepwise  = TRUE,
                                        corMethod = "spearman")
      
      n <- ncol(df_vars())
      rv$errorMod   <- mgm_error_metrics(data = df_vars(),
                                         type  = rep("g", n),
                                         level = rep(1, n))$errorCon$R2
      
      rv$qgraph_obj <- qgraph(rv$network_obj$graph,
                              groups      = rv$groups,
                              curveAll    = 2,
                              vsize       = 12,
                              esize       = 12,
                              palette     = if(is.null(input$network_palette)) "ggplot2" else input$network_palette,
                              layout      = "spring",
                              edge.labels = TRUE,
                              pie         = rv$errorMod,
                              layoutScale = c(0.8, 0.8),
                              legend.cex  = 0.5,
                              labels      = abbreviate(names(df_vars()), 3))
      
      showNotification("Red estimada exitosamente!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error al estimar la red:", e$message), type = "error", duration = 10)
    })
  })
  
  output$network_plot <- renderPlot({ 
    req(rv$qgraph_obj)
    plot(rv$qgraph_obj) 
  })
  
  cent <- eventReactive(input$run, {
    req(rv$qgraph_obj, rv$network_obj, rv$groups)
    
    isolate({
      rv$combined_plot_cache <- NULL
    })
    
    tryCatch({
      m0 <- if(is.null(input$measure0)) "ExpectedInfluence" else input$measure0  
      m1 <- if (input$measure1 == "None") NULL else input$measure1
      
      if (!is.null(m1) && (is.null(rv$groups) || length(rv$groups) < 2)) {
        showNotification("Advertencia: Se necesitan al menos 2 grupos para calcular medidas de puente. Se omitirá la medida de puente.", 
                         type = "warning", duration = 5)
        m1 <- NULL
      }
      
      result <- centrality_plots2_fixed(
        qgraph_obj    = rv$qgraph_obj,
        network       = rv$network_obj,
        groups        = rv$groups,
        measure0      = m0,  
        measure1      = m1,
        color_palette = c("#1e3a5f","#CC1228"),
        use_abbrev    = TRUE
      )
      
      result$timestamp <- Sys.time()
      
      return(result)
      
    }, error = function(e) {
      showNotification(paste("Error en cálculo de centralidad:", e$message), 
                       type = "error", duration = 10)
      
      m0 <- if(is.null(input$measure0)) "ExpectedInfluence" else input$measure0
      basic_result <- centrality_plots2_fixed(
        qgraph_obj    = rv$qgraph_obj,
        network       = rv$network_obj,
        groups        = rv$groups,
        measure0      = m0,  
        measure1      = NULL,
        color_palette = c("#1e3a5f","#CC1228"),
        use_abbrev    = TRUE
      )
      basic_result$timestamp <- Sys.time()
      return(basic_result)
    })
  })
  
  output$cent_table <- renderTable({
    result <- cent()
    if (!is.null(result$table)) {
      result$table %>%
        mutate(across(where(is.numeric), ~ round(.x, 3)))
    } else {
      data.frame(Mensaje = "No se pudieron calcular las métricas de centralidad")
    }
  })
  
  output$cent_plot <- renderPlot({
    result <- cent()
    if (!is.null(result$plot)) {
      print(result$plot)
    } else {
      plot.new()
      text(0.5, 0.5, "No se pudo generar el gráfico de centralidad", 
           cex = 1.2, col = "red", adj = c(0.5, 0.5))
    }
  })
  
  output$combined_plot <- renderPlot({
    req(rv$network_obj, rv$groups, rv$errorMod, rv$qgraph_obj)
    cent_data <- cent()
    
    isolate({
      tryCatch({
        m0 <- if(is.null(input$measure0)) "ExpectedInfluence" else input$measure0  
        m1 <- if (input$measure1 == "None") NULL else input$measure1
        
        if (!is.null(m1) && (is.null(rv$groups) || length(rv$groups) < 2)) {
          m1 <- NULL
        }
        
        fresh_cent_data <- centrality_plots2_fixed(
          qgraph_obj    = rv$qgraph_obj,
          network       = rv$network_obj,
          groups        = rv$groups,
          measure0      = m0,  
          measure1      = m1,
          color_palette = c("#1e3a5f", "#CC1228"),
          use_abbrev    = TRUE
        )
        
        combined_figure <- create_combined_figure(
          qgraph_obj = rv$qgraph_obj,
          network_obj = rv$network_obj,
          groups = rv$groups, 
          error_model = rv$errorMod,
          cent_data_plot = fresh_cent_data$plot
        )
        
        print(combined_figure)
        
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error al generar figura combinada:\n", e$message),
             cex = 1.2, col = "red", adj = c(0.5, 0.5))
      })
    })
  }, height = 600)
  
  observeEvent(c(input$measure0, input$measure1, input$network_palette), {
    rv$combined_plot_cache <- NULL
    
    if (!is.null(rv$network_obj) && !is.null(rv$groups) && !is.null(rv$errorMod) && !is.null(rv$df_full)) {
      tryCatch({
        current_df <- rv$df_full[, input$vars, drop = FALSE]
        
        rv$qgraph_obj <- qgraph(rv$network_obj$graph,
                                groups      = rv$groups,
                                curveAll    = 2,
                                vsize       = 12,
                                esize       = 12,
                                palette     = if(is.null(input$network_palette)) "ggplot2" else input$network_palette,
                                layout      = "spring",
                                edge.labels = TRUE,
                                pie         = rv$errorMod,
                                layoutScale = c(0.8, 0.8),
                                legend.cex  = 0.5,
                                labels      = abbreviate(names(current_df), 3))
        
      }, error = function(e) {
        showNotification(paste("Error al cambiar paleta:", e$message), type = "warning", duration = 3)
      })
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$run, {
    rv$combined_plot_cache <- NULL
    invalidateLater(200, session)
  }, ignoreInit = TRUE)
  
  output$combined_plot_status <- renderText({
    cent_data <- cent()
    
    if (is.null(rv$network_obj)) {
      "Esperando red estimada..."
    } else if (is.null(rv$groups)) {
      "Esperando configuración de grupos..."
    } else if (is.null(rv$errorMod)) {
      "Calculando errores del modelo..."
    } else if (is.null(rv$qgraph_obj)) {
      "Generando visualización..."
    } else if (!is.null(cent_data$timestamp)) {
      paste("Figura actualizada:", format(cent_data$timestamp, "%H:%M:%S"))
    } else {
      "Figura combinada lista"
    }
  })
  
  output$downloadReport <- downloadHandler(
    filename = function() paste0("Figura_Red_Centralidad_", Sys.Date(), ".jpg"),
    content  = function(file) {
      req(rv$network_obj, rv$groups, rv$errorMod, rv$qgraph_obj)
      
      tryCatch({
        m0 <- if(is.null(input$measure0)) "ExpectedInfluence" else input$measure0  
        m1 <- if (input$measure1 == "None") NULL else input$measure1
        
        if (!is.null(m1) && (is.null(rv$groups) || length(rv$groups) < 2)) {
          m1 <- NULL
        }
        
        fresh_cent_data <- centrality_plots2_fixed(
          qgraph_obj    = rv$qgraph_obj,
          network       = rv$network_obj,
          groups        = rv$groups,
          measure0      = m0,  
          measure1      = m1,
          color_palette = c("#1e3a5f", "#CC1228"),
          use_abbrev    = TRUE
        )
        
        combined_figure <- create_combined_figure(
          qgraph_obj = rv$qgraph_obj,
          network_obj = rv$network_obj,
          groups = rv$groups, 
          error_model = rv$errorMod,
          cent_data_plot = fresh_cent_data$plot
        )
        
        width_val <- if(is.null(input$network_width)) 10 else input$network_width
        height_val <- if(is.null(input$network_height)) 6 else input$network_height
        dpi_val <- if(is.null(input$network_dpi)) 300 else input$network_dpi
        
        ggsave(file, combined_figure, width = width_val, height = height_val, dpi = dpi_val, device = "jpeg")
        
      }, error = function(e) {
        jpeg(file, width = 1000, height = 600, res = 100)
        plot.new()
        text(0.5, 0.5, paste("Error al generar descarga:", e$message), 
             cex = 1.2, col = "red")
        dev.off()
      })
    }
  )
  
  output$network_render_status <- renderText({
    req(rv$network_obj, rv$groups, rv$errorMod, rv$qgraph_obj)
    
    tryCatch({
      m0 <- if(is.null(input$measure0)) "ExpectedInfluence" else input$measure0  
      m1 <- if (input$measure1 == "None") NULL else input$measure1
      if (!is.null(m1) && (is.null(rv$groups) || length(rv$groups) < 2)) {
        m1 <- NULL
      }
      
      cent_data <- centrality_plots2_fixed(
        qgraph_obj    = rv$qgraph_obj,
        network       = rv$network_obj,
        groups        = rv$groups,
        measure0      = m0,  
        measure1      = m1,
        color_palette = c("#1e3a5f", "#CC1228"),
        use_abbrev    = TRUE
      )
      
      "✅ Figura generada correctamente"
      
    }, error = function(e) {
      paste("❌ Error:", e$message)
    })
  })
  
  outputOptions(output, "combined_plot", suspendWhenHidden = FALSE)
  
  # ---- CONFIGURACIÓN BOOTSTRAP ----
  # Reemplaza la función finish_bootstrap_cycle existente con esta versión mejorada:
  finish_bootstrap_cycle <- function() {
    # Remover notificación si existe
    removeNotification("long_process")
    
    # Verificar si ambos procesos han terminado
    if (!rv_bootstrap$case_running && !rv_bootstrap$nonparam_running) {
      rv_bootstrap$both_running <- FALSE
      shinyjs::enable("run_smart_bootstrap")
      updateActionButton(
        session, "run_smart_bootstrap",
        label = tags$div(icon("redo"), tags$strong("VOLVER A EJECUTAR ANÁLISIS")),
        icon = NULL
      )
    }
  }
  
  
  
  
  observe({
    req(input$vars)
    if (length(input$vars) > 0) {
      config <- auto_bootstrap_config(length(input$vars), input$analysis_type)
      rv_bootstrap$config <- config
      if (!is.null(input$use_custom) && !input$use_custom) {
        updateNumericInput(session, "custom_case_boot", value = config$case_boots)
        updateNumericInput(session, "custom_nonparam_boot", value = config$nonparam_boots)
      }
    }
  })
  
  output$bootstrap_has_results <- reactive({
    !is.null(rv_bootstrap$case_result) || !is.null(rv_bootstrap$nonparam_result)
  })
  outputOptions(output, "bootstrap_has_results", suspendWhenHidden = FALSE)
  
  observeEvent(input$clear_bootstrap, {
    rv_bootstrap$case_result <- NULL
    rv_bootstrap$nonparam_result <- NULL
    rv_bootstrap$selected_statistics <- NULL
    rv_bootstrap$stats_effective <- NULL
    rv_bootstrap$stats_available <- NULL
    
    updateProgressBar(session, "case_progress_new", value = 0)
    updateProgressBar(session, "nonparam_progress_new", value = 0)
    
    updateActionButton(session, "run_smart_bootstrap",
                       label = tags$div(icon("rocket", style = "font-size: 16px;"), 
                                        tags$strong("EJECUTAR ANÁLISIS")))
    
    showNotification("✅ Resultados limpiados. Listo para nuevo análisis.", 
                     type = "message", duration = 3)
  })
  
  observe({
    has_results <- !is.null(rv_bootstrap$case_result) || !is.null(rv_bootstrap$nonparam_result)
    if (has_results) shinyjs::enable("clear_bootstrap") else shinyjs::disable("clear_bootstrap")
  })
  
  observe({
    case_done     <- !is.null(rv_bootstrap$case_result)     && !rv_bootstrap$case_running
    nonparam_done <- !is.null(rv_bootstrap$nonparam_result) && !rv_bootstrap$nonparam_running
    is_running <- rv_bootstrap$both_running || rv_bootstrap$case_running || rv_bootstrap$nonparam_running
    
    if (is_running) {
      shinyjs::disable("run_smart_bootstrap")
      updateActionButton(session, "run_smart_bootstrap",
                         label = "EJECUTANDO...", icon = icon("spinner", class = "fa-spin"))
    } else {
      shinyjs::enable("run_smart_bootstrap")
      if (case_done && nonparam_done) {
        updateActionButton(session, "run_smart_bootstrap",
                           label = "VOLVER A EJECUTAR ANÁLISIS", icon = icon("redo"))
      } else {
        updateActionButton(session, "run_smart_bootstrap",
                           label = "EJECUTAR ANÁLISIS", icon = icon("rocket"))
      }
    }
  })
  
  observe({
    if (rv_bootstrap$both_running || rv_bootstrap$case_running || rv_bootstrap$nonparam_running) {
      invalidateLater(500, session)
    }
  })
  
  output$bootstrap_info <- renderText({
    req(input$vars)
    if (is.null(rv_bootstrap$config)) return("Seleccione variables primero")
    
    config <- rv_bootstrap$config
    type_text <- switch(input$analysis_type,
                        "50" = "Rápido (50 bootstraps)",
                        "250" = "Estándar (250 bootstraps)",
                        "500" = "Completo (500 bootstraps)",
                        "1000" = "Investigación (1000 bootstraps)")
    focus_text <- switch(input$bootstrap_focus,
                         "case_only" = "Solo estabilidad",
                         "nonparam_only" = "Solo precisión",
                         "both" = "Estabilidad + Precisión")
    use_custom_val <- isTRUE(input$use_custom)
    
    case_stat <- if(use_custom_val && !is.null(input$custom_statistics) && length(input$custom_statistics)>0) {
      paste(input$custom_statistics, collapse = ", ")
    } else {
      if(is.null(input$measure0)) "expectedInfluence" else if (input$measure0 == "ExpectedInfluence") "expectedInfluence" else "strength"
    }
    
    case_boots_display <- if(use_custom_val) input$custom_case_boot else config$case_boots
    nonparam_boots_display <- if(use_custom_val) input$custom_nonparam_boot else config$nonparam_boots
    
    paste0("Configuración: ", type_text, "\n",
           "Análisis: ", focus_text, "\n",
           "Case Dropping: ", case_boots_display, " bootstraps (", case_stat, ")\n",
           "Non-Parametric: ", nonparam_boots_display, " bootstraps (edge)\n",
           "Estadísticos seleccionados: ", case_stat, "\n",
           "Núcleos: ", config$cores)
  })
  
  output$time_estimate_text <- renderText({
    req(input$vars)
    if (is.null(rv_bootstrap$config)) return("Configurando...")
    
    config <- rv_bootstrap$config
    use_custom_val <- isTRUE(input$use_custom)
    
    case_boots <- if(use_custom_val) input$custom_case_boot else config$case_boots
    nonparam_boots <- if(use_custom_val) input$custom_nonparam_boot else config$nonparam_boots
    
    case_time <- case_boots * 0.8 / config$cores
    nonparam_time <- nonparam_boots * 1.2 / config$cores
    
    total_time <- switch(input$bootstrap_focus,
                         "case_only" = case_time,
                         "nonparam_only" = nonparam_time,
                         "both" = case_time + nonparam_time)
    
    if (total_time < 60) paste0(round(total_time), " segundos")
    else paste0(round(total_time/60, 1), " minutos")
  })
  
  output$cores_info_text <- renderText({
    if (is.null(rv_bootstrap$config)) return("Preparando...")
    paste0(rv_bootstrap$config$cores, " núcleos\n", detectCores(logical = FALSE), " disponibles")
  })
  
  output$overall_status_text <- renderText({
    if (rv_bootstrap$both_running) "Ejecutando..."
    else if (rv_bootstrap$case_running) "Case Dropping..."
    else if (rv_bootstrap$nonparam_running) "Non-Parametric..."
    else if (!is.null(rv_bootstrap$case_result) || !is.null(rv_bootstrap$nonparam_result)) "Completado"
    else "Listo para ejecutar"
  })
  
  # ---- RUN_SMART_BOOTSTRAP ----
  observeEvent(input$run_smart_bootstrap, {
    req(rv$network_obj, rv_bootstrap$config)
    if (isTRUE(rv_bootstrap$both_running)) return()
    
    config <- rv_bootstrap$config
    use_custom_val <- isTRUE(input$use_custom)
    
    case_boots_val     <- if (use_custom_val && !is.null(input$custom_case_boot))     input$custom_case_boot     else config$case_boots
    nonparam_boots_val <- if (use_custom_val && !is.null(input$custom_nonparam_boot)) input$custom_nonparam_boot else config$nonparam_boots
    
    # estadísticas solicitadas (desde Configuración Avanzada o por defecto)
    case_statistics_val <- if (use_custom_val && !is.null(input$custom_statistics) && length(input$custom_statistics) > 0) {
      as.character(input$custom_statistics)
    } else {
      c("expectedInfluence")
    }
    
    # --- Manejo de bridgeStrength y comunidades
    want_bridge <- "bridgeStrength" %in% case_statistics_val
    communities_vec_val <- isolate(rv$communities_vec)
    communities_arg <- NULL
    if (want_bridge) {
      if (!is.null(communities_vec_val) && !any(is.na(communities_vec_val))) {
        communities_arg <- communities_vec_val
      } else {
        showNotification("⚠️ 'bridgeStrength' omitido: no hay comunidades válidas.", type = "warning", duration = 5)
        case_statistics_val <- setdiff(case_statistics_val, "bridgeStrength")
      }
    }
    if (length(case_statistics_val) == 0) case_statistics_val <- "strength"  # fallback seguro
    
    rv_bootstrap$selected_statistics <- case_statistics_val
    
    bootstrap_focus_val <- if (is.null(input$bootstrap_focus)) "both" else input$bootstrap_focus
    network_obj_val      <- isolate(rv$network_obj)
    cores_val            <- config$cores
    
    rv_bootstrap$both_running <- TRUE
    rv_bootstrap$case_running <- FALSE
    rv_bootstrap$nonparam_running <- FALSE
    # AGREGAR ESTE BLOQUE AQUÍ:
    if (case_boots_val >= 1000 || nonparam_boots_val >= 1000) {
      showNotification(
        "⚠️ Con 1000+ bootstraps el análisis puede tomar varios minutos. Por favor espere...", 
        type = "warning", 
        duration = NULL,  # No desaparece automáticamente
        id = "long_process"
      )
    }
    updateActionButton(session, "run_smart_bootstrap", label = "EJECUTANDO...", icon = icon("spinner", class = "fa-spin"))
    shinyjs::disable("run_smart_bootstrap")
    
    # -------- CASE DROPPING ----------
    if (bootstrap_focus_val %in% c("case_only","both")) {
      rv_bootstrap$case_running <- TRUE
      rv_bootstrap$case_start_time <- Sys.time()  # ← AGREGAR ESTA LÍNEA
      updateProgressBar(session, "case_progress_new", value = 10)
      
      future({
        bootnet::bootnet(
          network_obj_val,
          nBoots      = case_boots_val,
          type        = "case",
          nCores      = cores_val,
          statistics  = case_statistics_val,
          communities = communities_arg
        )
      }) %...>% (function(res) {
        rv_bootstrap$case_result  <- res
        rv_bootstrap$case_running <- FALSE
        updateProgressBar(session, "case_progress_new", value = 100)
        
        # AGREGAR ESTA LÍNEA:
        removeNotification("long_process")
        
        
        # Guardar las estadísticas efectivas: intersección de solicitadas con disponibles (si existen)
        actual_stats <- tryCatch(unique(as.character(res$bootTable$statistic)), error = function(e) character(0))
        eff <- if (length(actual_stats)) {
          intersect(rv_bootstrap$selected_statistics %||% character(0), actual_stats)
        } else {
          rv_bootstrap$selected_statistics
        }
        
        if (length(eff) == 0) eff <- rv_bootstrap$selected_statistics
        rv_bootstrap$stats_effective <- unique(eff)
        
        showNotification(paste("✅ Case Dropping completado. Stats:",
                               paste(rv_bootstrap$stats_effective, collapse=", ")),
                         type="message", duration=5)
        
        if (bootstrap_focus_val == "case_only" && !rv_bootstrap$nonparam_running) {
          rv_bootstrap$both_running <- FALSE
          shinyjs::enable("run_smart_bootstrap")
          updateActionButton(session, "run_smart_bootstrap", label = "VOLVER A EJECUTAR ANÁLISIS", icon = icon("redo"))
        }
        finish_bootstrap_cycle()  # Llamar siempre para verificar si ambos terminaron
        NULL
      }) %...!% (function(err) {
        rv_bootstrap$case_running <- FALSE
        rv_bootstrap$case_result  <- NULL
        rv_bootstrap$stats_effective <- NULL
        updateProgressBar(session, "case_progress_new", value = 0)
        # AGREGAR:
        removeNotification("long_process")
        showNotification(paste("❌ Error en Case Dropping:", err$message), type = "error", duration = 10)
        finish_bootstrap_cycle()      # <--- NUEVO
        rv_bootstrap$both_running <- FALSE
        shinyjs::enable("run_smart_bootstrap")
        updateActionButton(session, "run_smart_bootstrap", label = "EJECUTAR ANÁLISIS", icon = icon("rocket"))

      })
    }
    
    # -------- NON-PARAMETRIC ----------
    if (bootstrap_focus_val %in% c("nonparam_only","both")) {
      rv_bootstrap$nonparam_running <- TRUE
      rv_bootstrap$nonparam_start_time <- Sys.time()  # ← AGREGAR ESTA LÍNEA
      updateProgressBar(session, "nonparam_progress_new", value = 10)
      
      future({
        bootnet::bootnet(
          network_obj_val,
          nBoots      = nonparam_boots_val,
          type        = "nonparametric",
          nCores      = cores_val,
          statistics  = "edge",
          communities = NULL
        )
      }) %...>% (function(res) {
        rv_bootstrap$nonparam_result  <- res
        rv_bootstrap$nonparam_running <- FALSE
        updateProgressBar(session, "nonparam_progress_new", value = 100)
        # AGREGAR ESTA LÍNEA:
        removeNotification("long_process")
        showNotification("✅ Non-Parametric completado!", type="message", duration=3)
        
        if (bootstrap_focus_val == "nonparam_only" && !rv_bootstrap$case_running) {
          rv_bootstrap$both_running <- FALSE
          shinyjs::enable("run_smart_bootstrap")
          updateActionButton(session, "run_smart_bootstrap", label = "VOLVER A EJECUTAR ANÁLISIS", icon = icon("redo"))
        } else if (bootstrap_focus_val == "both" && !rv_bootstrap$case_running) {
          rv_bootstrap$both_running <- FALSE
          shinyjs::enable("run_smart_bootstrap")
          updateActionButton(session, "run_smart_bootstrap", label = "VOLVER A EJECUTAR ANÁLISIS", icon = icon("redo"))
        }
        finish_bootstrap_cycle()  # Llamar siempre para verificar si ambos terminaron
        NULL
      }) %...!% (function(err) {
        rv_bootstrap$nonparam_running <- FALSE
        rv_bootstrap$nonparam_result  <- NULL
        updateProgressBar(session, "nonparam_progress_new", value = 0)
        # AGREGAR:
        removeNotification("long_process")
        showNotification(paste("❌ Error en Non-Parametric:", err$message), type = "error", duration = 10)
        finish_bootstrap_cycle()      # <--- NUEVO
        rv_bootstrap$both_running <- FALSE
        shinyjs::enable("run_smart_bootstrap")
        updateActionButton(session, "run_smart_bootstrap", label = "EJECUTAR ANÁLISIS", icon = "rocket")

      })
    }
  })
  
  # ---- OBSERVADORES DE PROGRESO ----
  observe({
    invalidateLater(1000, session)
    if (rv_bootstrap$case_running && !is.null(rv_bootstrap$case_start_time) && !is.null(rv_bootstrap$config)) {
      elapsed <- as.numeric(difftime(Sys.time(), rv_bootstrap$case_start_time, units = "secs"))
      case_boots <- if(isTRUE(input$use_custom) && !is.null(input$custom_case_boot)) input$custom_case_boot else rv_bootstrap$config$case_boots
      estimated_total <- case_boots * 0.8 / rv_bootstrap$config$cores
      progress <- min(95, 10 + (elapsed / estimated_total) * 85)
      updateProgressBar(session, "case_progress_new", value = progress)
    }
  })
  
  observe({
    invalidateLater(1000, session)
    if (rv_bootstrap$nonparam_running && !is.null(rv_bootstrap$nonparam_start_time) && !is.null(rv_bootstrap$config)) {
      elapsed <- as.numeric(difftime(Sys.time(), rv_bootstrap$nonparam_start_time, units = "secs"))
      nonparam_boots <- if(isTRUE(input$use_custom) && !is.null(input$custom_nonparam_boot)) input$custom_nonparam_boot else rv_bootstrap$config$nonparam_boots
      estimated_total <- nonparam_boots * 1.2 / rv_bootstrap$config$cores
      progress <- min(95, 10 + (elapsed / estimated_total) * 85)
      updateProgressBar(session, "nonparam_progress_new", value = progress)
    }
  })
  
  # Agrega este observe después del observeEvent(input$run_smart_bootstrap, {...})
  observe({
    # Verificar cada 2 segundos si los procesos realmente terminaron
    invalidateLater(2000, session)
    
    if (rv_bootstrap$both_running && 
        !rv_bootstrap$case_running && 
        !rv_bootstrap$nonparam_running &&
        (!is.null(rv_bootstrap$case_result) || !is.null(rv_bootstrap$nonparam_result))) {
      
      # Forzar actualización si detectamos que terminó pero no se actualizó
      finish_bootstrap_cycle()
    }
  })
  
  output$case_status_new <- renderText({
    if (!is.null(rv_bootstrap$case_result)) "Completado"
    else if (rv_bootstrap$case_running) {
      elapsed <- if(!is.null(rv_bootstrap$case_start_time)) round(as.numeric(difftime(Sys.time(), rv_bootstrap$case_start_time, units = "secs"))) else 0
      paste("Ejecutando...", elapsed, "seg")
    } else "Listo"
  })
  
  output$nonparam_status_new <- renderText({
    if (!is.null(rv_bootstrap$nonparam_result)) "Completado"
    else if (rv_bootstrap$nonparam_running) {
      elapsed <- if(!is.null(rv_bootstrap$nonparam_start_time)) round(as.numeric(difftime(Sys.time(), rv_bootstrap$nonparam_start_time, units = "secs"))) else 0
      paste("Ejecutando...", elapsed, "seg")
    } else "Listo"
  })
  
  observe({
    invalidateLater(2000, session)
    if (!is.null(rv_bootstrap$case_result) && !rv_bootstrap$case_running) {
      updateProgressBar(session, "case_progress_new", value = 100)
    }
    if (!is.null(rv_bootstrap$nonparam_result) && !rv_bootstrap$nonparam_running) {
      updateProgressBar(session, "nonparam_progress_new", value = 100)
    }
  })
  
  # ---- PLOTS Y OUTPUTS BOOTSTRAP ----
  output$case_detailed <- renderPrint({
    req(rv_bootstrap$case_result)
    cat("🎯 ESTABILIDAD DE CENTRALIDAD - DETALLE\n")
    cat("======================================\n\n")
    cat("📊 Estadísticos configurados: ",
        if (is.null(rv_bootstrap$selected_statistics)) "[no disponible]" else paste(rv_bootstrap$selected_statistics, collapse=", "),
        "\n\n", sep = "")
    
    actual_stats <- tryCatch(unique(as.character(rv_bootstrap$case_result$bootTable$statistic)),
                             error = function(e) character(0))
    cat("📈 Estadísticos en el resultado: ",
        if (length(actual_stats)==0) "[No se encontraron estadísticas]" else paste(actual_stats, collapse=", "),
        "\n\n", sep = "")
    
    cat("💡 Interpretación:\n")
    cat("• CS > 0.75: Excelente\n• CS > 0.50: Aceptable\n• CS < 0.50: Baja\n• CS < 0.25: Insuficiente\n")
    invisible(NULL)
  })
  
  output$case_plot_detailed <- renderPlot({
    req(rv_bootstrap$case_result)
    stats_to_show <- rv_bootstrap$stats_effective
    if (is.null(stats_to_show) || length(stats_to_show) == 0) stats_to_show <- "strength"
    stats_to_show <- unique(as.character(stats_to_show))
    
    p <- try(plot(rv_bootstrap$case_result, plot = "area", statistics = stats_to_show), silent = TRUE)
    if (inherits(p, "ggplot")) {
      p <- p + ggplot2::ggtitle(paste("Case Dropping Stability:", paste(stats_to_show, collapse = ", ")))
      print(p)
    } else {
      try(title(main = paste("Case Dropping Stability:", paste(stats_to_show, collapse = ", ")),
                cex.main = 1.2, font.main = 2), silent = TRUE)
    }
  }, height = 400, width = 500)
  
  output$nonparam_plot_detailed <- renderPlot({
    req(rv_bootstrap$nonparam_result)
    plot(
      rv_bootstrap$nonparam_result,
      statistics = "edge",
      labels = isTRUE(input$boot_labels)  # <- usa el toggle
    )
  })
  

  
  output$combined_bootstrap_plot <- renderPlot({
    req(rv_bootstrap$case_result, rv_bootstrap$nonparam_result)
    tryCatch({
      stats_to_use <- if(!is.null(rv_bootstrap$selected_statistics)) rv_bootstrap$selected_statistics else c("strength")
      if (requireNamespace("InterconectaR", quietly = TRUE)) {
        combined_plot <- InterconectaR::plot_centrality_stability(
          caseDroppingBoot = rv_bootstrap$case_result,
          nonParametricBoot = rv_bootstrap$nonparam_result,
          statistics = stats_to_use,
          labels = isTRUE(input$boot_labels)  # <- toggle
        )
        print(combined_plot)
      } else {
        stop("InterconectaR no disponible")
      }
    }, error = function(e) {
      tryCatch({
        primary_stat <- if(!is.null(rv_bootstrap$selected_statistics)) rv_bootstrap$selected_statistics[1] else "strength"
        par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
        plot(rv_bootstrap$case_result, statistics = primary_stat, labels = isTRUE(input$boot_labels), plot = "area")
        title(main = "A. Case Dropping", cex.main = 1.2)
        plot(rv_bootstrap$nonparam_result, labels = isTRUE(input$boot_labels), order = "sample", statistics = "edge")
        title(main = "B. Non-Parametric Bootstrap", cex.main = 1.2)
        par(mfrow = c(1, 1))
      }, error = function(e2) {
        plot.new(); text(0.5, 0.5, paste("Error:", e$message), cex = 1.2, col = "red")
      })
    })
  }, height = 500)
  
  output$download_bootstrap_combined <- downloadHandler(
    filename = function() paste0("Bootstrap_Analysis_", Sys.Date(), ".jpg"),
    content = function(file) {
      req(rv_bootstrap$case_result, rv_bootstrap$nonparam_result)
      tryCatch({
        statistics_val <- if(!is.null(rv_bootstrap$selected_statistics)) rv_bootstrap$selected_statistics else {
          measure_val <- if(is.null(input$measure0)) "ExpectedInfluence" else input$measure0
          if(measure_val == "ExpectedInfluence") "expectedInfluence" else "strength"
        }
        width_val <- if(is.null(input$bootstrap_width)) 9 else input$bootstrap_width
        height_val <- if(is.null(input$bootstrap_height)) 5 else input$bootstrap_height
        dpi_val <- if(is.null(input$bootstrap_dpi)) 600 else input$bootstrap_dpi
        p <- InterconectaR::plot_centrality_stability(
          caseDroppingBoot = rv_bootstrap$case_result,
          nonParametricBoot = rv_bootstrap$nonparam_result,
          statistics = statistics_val,
          labels = isTRUE(input$boot_labels)
        )
        ggsave(file, p, width = width_val, height = height_val, dpi = dpi_val, device = "jpeg")
      }, error = function(e) {
        jpeg(file, width = 900, height = 500, res = 100)
        plot.new()
        text(0.5, 0.5, paste("Error al generar figura combinada:", e$message), cex = 1.2, col = "red")
        dev.off()
      })
    }
  )
  
  # ====== CS: cálculo, interpretación y UI ======
  cs_table_data <- reactive({
    req(rv_bootstrap$case_result)
    tryCatch({
      InterconectaR::filter_correlation_stability(rv_bootstrap$case_result) %>%
        tibble::as_tibble() %>%
        dplyr::rename(Métrica = rowname, `CS (cor)` = Index) %>%
        dplyr::mutate(
          `CS (cor)` = round(`CS (cor)`, 2),
          Interpretación = dplyr::case_when(
            `CS (cor)` >= 0.75 ~ "Excelente (≥ .75)",
            `CS (cor)` >= 0.50 ~ "Adecuado (.50–.74)",
            TRUE               ~ "Bajo (< .50)"
          )
        ) %>%
        dplyr::arrange(desc(`CS (cor)`))
    }, error = function(e) {
      tibble::tibble(Métrica = "—", `CS (cor)` = NA_real_, Interpretación = "—")
    })
  })
  
  output$cs_table <- renderTable({
    cs_table_data()
  })
  
  # Tarjetas-resumen (usa tu clase .stat-card del CSS)
  output$cs_cards <- renderUI({
    d <- cs_table_data()
    if (nrow(d) == 0 || is.na(d$`CS (cor)`[1])) return(NULL)
    
    make_card <- function(m, v, txt) {
      # color de borde superior según semáforo
      col <- if (is.na(v)) "#6c757d" else if (v >= 0.75) "#5cb85c" else if (v >= 0.50) "#f0ad4e" else "#d9534f"
      div(class = "stat-card",
          style = paste0("border-top: 3px solid ", col, ";"),
          div(class = "stat-number", sprintf("%.2f", v)),
          div(class = "stat-label", strong(m)),
          tags$small(txt)
      )
    }
    
    # genera una tarjeta por fila
    cards <- lapply(seq_len(nrow(d)), function(i) {
      make_card(d$Métrica[i], d$`CS (cor)`[i], d$Interpretación[i])
    })
    
    fluidRow(
      # hasta 3 tarjetas por fila; ajusta si agregas más métricas
      column(4, if (length(cards) >= 1) cards[[1]]),
      column(4, if (length(cards) >= 2) cards[[2]]),
      column(4, if (length(cards) >= 3) cards[[3]])
    )
  })
  
  # Descarga CSV del CS
  output$download_cs_csv <- downloadHandler(
    filename = function() paste0("CS_correlation_", Sys.Date(), ".csv"),
    content  = function(file) {
      utils::write.csv(cs_table_data(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
  
  # ===== Case Dropping: tabla (plot_df) =====
  cd_table_data_raw <- reactive({
    req(rv_bootstrap$case_result)
    tryCatch({
      InterconectaR::summarise_case_drop_stability(rv_bootstrap$case_result) %>%
        dplyr::arrange(type, nPerson)  # opcional, para mantener orden
    }, error = function(e) {
      showNotification(paste("Error al resumir Case Dropping:", e$message),
                       type = "error", duration = 6)
      tibble::tibble()
    })
  })
  
  # Versión filtrada y redondeada para UI/descargas
  cd_table_data <- reactive({
    df <- cd_table_data_raw()
    # Filtro por tipo (estadístico)
    if (!is.null(input$cd_type_filter) && input$cd_type_filter != "Todos") {
      df <- df %>% dplyr::filter(type == input$cd_type_filter)
    }
    dec <- if (is.null(input$cd_decimals)) 3 else input$cd_decimals
    df %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, dec)))
  })
  
  output$cd_table <- renderTable({
    req(cd_table_data())
    cd_table_data()
  })
  
  # ---- Descarga XLSX (con fallback) ----
  output$download_cd_xlsx <- downloadHandler(
    filename = function() paste0("CaseDropping_table_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      d <- cd_table_data()
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(d, path = file)
      } else if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(d, file)
      } else {
        showNotification("No se encontró 'writexl' ni 'openxlsx'. Se exportará en CSV.",
                         type = "warning", duration = 5)
        utils::write.csv(d, file, row.names = FALSE, fileEncoding = "UTF-8")
      }
    }
  )
  
  # ===== Nonparametric bootstrap: tabla de aristas =====
  np_edges_raw <- reactive({
    req(rv_bootstrap$nonparam_result)
    tryCatch({
      InterconectaR::summarise_nonparametric_edges(rv_bootstrap$nonparam_result) %>%
        dplyr::mutate(abs_mean = abs(boot_mean))
    }, error = function(e) {
      showNotification(paste("Error al resumir nonparametric bootstrap:", e$message),
                       type = "error", duration = 6)
      tibble::tibble()
    })
  })
  
  np_edges_table_data <- reactive({
    df <- np_edges_raw()
    if (nrow(df) == 0) return(df)
    
    # Filtro: solo aristas con IC que NO incluye 0
    if (isTRUE(input$np_only_sig)) {
      df <- df %>% dplyr::filter(ci_covers_zero == FALSE)
    }
    
    # Orden
    ord <- if (is.null(input$np_order)) "mag" else input$np_order
    df <- switch(ord,
                 "name"     = df %>% dplyr::arrange(id),
                 "zerofirst"= df %>% dplyr::arrange(ci_covers_zero, dplyr::desc(abs_mean)),
                 # default "mag"
                 df %>% dplyr::arrange(dplyr::desc(abs_mean))
    )
    
    # Redondeo y renombre columnas para UI
    dec <- if (is.null(input$np_decimals)) 3 else input$np_decimals
    df %>%
      dplyr::select(id, sample, boot_mean, ci_lo, ci_hi, n_boot, ci_covers_zero) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, dec))) %>%
      dplyr::rename(
        Arista           = id,
        `Muestra`        = sample,
        `Media bootstrap`= boot_mean,
        `CI 2.5%`        = ci_lo,
        `CI 97.5%`       = ci_hi,
        `Reps`           = n_boot,
        `IC incluye 0`   = ci_covers_zero
      )
  })
  
  output$np_edges_table <- renderTable({
    np_edges_table_data()
  })
  
  # ---- Descarga XLSX (con fallback) ----
  output$download_np_xlsx <- downloadHandler(
    filename = function() paste0("Nonparametric_edges_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      d <- np_edges_table_data()
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(d, path = file)
      } else if (requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(d, file)
      } else {
        showNotification("No se encontró 'writexl' ni 'openxlsx'. Se exportará en CSV.",
                         type = "warning", duration = 5)
        utils::write.csv(d, file, row.names = FALSE, fileEncoding = "UTF-8")
      }
    }
  )
  
  # ---- COMPARACIÓN DE GRUPOS ----
  parse_colors <- function(x, fallback = c("#1e3a5f", "#CC1228")) {
    if (is.null(x) || !nzchar(x)) return(fallback)
    cols <- unlist(strsplit(x, "\\s*,\\s*"))
    cols <- trimws(cols)
    is_hex <- function(z) grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", z)
    cols <- cols[is_hex(cols)]
    if (length(cols) < 2) fallback else cols[1:2]
  }
  
  
  observeEvent(input$run_groups, {
    req(rv$df_full, input$group_var, input$vars)
    rv$groups_running <- TRUE
    updateProgressBar(session, "groups_progress", value = 0)
    withProgress(message = "Analizando por grupos...", value = 0, {
      tryCatch({
        df_for_groups <- rv$df_full %>%
          select(all_of(c(input$vars, input$group_var))) %>%
          na.omit() %>%
          rename(.group = !!input$group_var)
        
        incProgress(0.1, detail = "Estimando redes...")
        rv$networks_groups <- InterconectaR::estimate_networks_by_group(
          data            = df_for_groups,
          group_var       = ".group",
          columns         = input$vars,
          default         = "ggmModSelect",
          stepwise        = TRUE,
          corMethod       = "spearman",
          abbreviate_vars = TRUE,
          abbr_minlength  = 3
        )
        
        rv$networks_groups <- lapply(rv$networks_groups, function(net) {
          net$graph <- abs(net$graph)
          net
        })
        
        incProgress(0.3, detail = "Centralidad por grupo...")
        colors <- parse_colors(input$group_colors, fallback = c("#1e3a5f", "#CC1228"))
        group_names <- c(input$group1_name, input$group2_name)
        measure_for_groups <- if (is.null(input$measure0)) "ExpectedInfluence" else input$measure0
        
        rv$plot_centralidad_group <- InterconectaR::plot_centrality_by_group(
          rv$networks_groups,
          replacements  = group_names,
          measure_spec  = c(measure_for_groups),
          color_palette = colors
        )
        
        rv$bridge_plot_group <- InterconectaR::centrality_bridge_plot(
          networks_groups = rv$networks_groups,
          group_names     = group_names,
          measure         = "Bridge Expected Influence (1-step)",
          color_palette   = colors
        )
        
        incProgress(0.5, detail = "Errores por grupo...")
        rv$errores_groups <- tryCatch({
          InterconectaR::mgm_errors_groups(
            data    = df_for_groups,
            type    = rep("g", length(input$vars)),
            level   = rep(1, length(input$vars)),
            group   = .group,
            columns = input$vars
          )
        }, error = function(e) {
          warning("No se pudieron calcular mgm_errors_groups: ", e$message)
          NULL
        })
        rv$pie_values <- if (!is.null(rv$errores_groups)) {
          lapply(rv$errores_groups, `[[`, "R2")
        } else {
          rep(list(NULL), length(rv$networks_groups))
        }
        
        incProgress(0.7, detail = "Redes por grupo...")
        rv$combined_plot_groups <- InterconectaR::plot_networks_by_group(
          res               = 300,
          networks_by_group = rv$networks_groups,
          groups            = rv$groups,
          pie               = rv$pie_values,
          legend.cex        = 0.8
        )
        
        incProgress(0.8, detail = "Bridge por grupo...")
        rv$bridge_plot_group <- InterconectaR::centrality_bridge_plot(
          networks_groups = rv$networks_groups,
          group_names     = group_names,
          measure         = "Bridge Expected Influence (1-step)",
          color_palette   = colors
        )
        
        incProgress(0.9, detail = "Combinando grupos...")
        rv$combined_groups_final <- InterconectaR::combine_groupBy(
          red_group               = rv$combined_plot_groups,
          plot_centralidad_group  = rv$plot_centralidad_group,
          bridge_plot_group       = rv$bridge_plot_group$plot,
          width_a                 = if (is.null(input$groups_width_a)) 15 else input$groups_width_a,
          width_bc                = if (is.null(input$groups_width_bc)) 4.5 else input$groups_width_bc,
          show_plot               = FALSE
        )
        
        incProgress(1, detail = "¡Listo!")
        rv$groups_running <- FALSE
        updateProgressBar(session, "groups_progress", value = 100)
        showNotification("Análisis por grupos completado!", type = "message")
        
        # ========== CALCULAR TAMAÑO DEL EFECTO ==========
        if (length(rv$networks_groups) == 2) {
          tryCatch({
            # Obtener matrices de adyacencia
            group_names_effect <- names(rv$networks_groups)
            network1_adjacency <- getWmat(rv$networks_groups[[group_names_effect[1]]])
            network2_adjacency <- getWmat(rv$networks_groups[[group_names_effect[2]]])
            
            # Correlación de Spearman
            rv$effect_size_corr <- cor(
              network1_adjacency[lower.tri(network1_adjacency)], 
              network2_adjacency[lower.tri(network2_adjacency)], 
              method = "spearman"
            )
            
            # Diferencia absoluta media
            diff_matrix <- abs(network1_adjacency - network2_adjacency)
            rv$effect_size_mean_diff <- mean(diff_matrix[lower.tri(diff_matrix)])
            
            # Diferencia máxima
            rv$effect_size_max_diff <- max(diff_matrix[lower.tri(diff_matrix)])
            
            # Número de aristas que difieren significativamente (> 0.1 en valor absoluto)
            rv$effect_size_n_diff <- sum(diff_matrix[lower.tri(diff_matrix)] > 0.1)
            
            # Proporción de varianza compartida (R²)
            rv$effect_size_r_squared <- rv$effect_size_corr^2
            
            # Guardar matrices para análisis adicionales
            rv$effect_size_diff_matrix <- diff_matrix
            
          }, error = function(e) {
            showNotification(paste("Error calculando tamaño del efecto:", e$message), 
                             type = "warning", duration = 5)
          })
        }
        
        # Network Comparison Test (NCT)
        group_names_nct <- names(rv$networks_groups)
        if (length(group_names_nct) == 2) {
          rv$nct_running    <- TRUE
          rv$nct_start_time <- Sys.time()
          updateProgressBar(session, "nct_progress", value = 10)
          showNotification("🔬 Iniciando Network Comparison Test...", type = "message", duration = 3)
          
          network1   <- rv$networks_groups[[group_names_nct[1]]]
          network2   <- rv$networks_groups[[group_names_nct[2]]]
          iterations <- if (is.null(input$nct_iterations)) 100 else input$nct_iterations
          
          library(igraph)
          g1 <- graph_from_adjacency_matrix(network1$graph, weighted = TRUE)
          
          comm <- cluster_spinglass(
            g1,
            weights        = E(g1)$weight,
            implementation = "neg"
          )
          
          future({
            set.seed(234)
            NCT(
              network1, network2,
              binary.data      = FALSE,
              test.edges       = FALSE,
              edges            = "all",
              it               = iterations,
              weighted         = TRUE,
              test.centrality  = FALSE,
              communities      = comm$membership,
              useCommunities   = "all"
            )
          }) %...>% (function(result) {
            rv$nct_result  <- result
            rv$nct_running <- FALSE
            updateProgressBar(session, "nct_progress", value = 100)
            showNotification("✅ Network Comparison Test completado!", type = "message", duration = 5)
          }) %...!% (function(error) {
            rv$nct_running <- FALSE
            updateProgressBar(session, "nct_progress", value = 0)
            showNotification(paste("❌ Error en NCT:", error$message), type = "error", duration = 10)
          })
          
        } else {
          showNotification("⚠️ NCT requiere exactamente 2 grupos para la comparación", type = "warning", duration = 5)
        }
        
      }, error = function(e) {
        rv$groups_running <- FALSE
        updateProgressBar(session, "groups_progress", value = 0)
        showNotification(paste("Error en análisis por grupos:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  # Outputs de tamaño del efecto
  output$effect_size_correlation <- renderText({
    req(rv$effect_size_corr)
    sprintf("%.3f", rv$effect_size_corr)
  })
  
  output$effect_size_difference <- renderText({
    req(rv$effect_size_mean_diff)
    sprintf("%.3f", rv$effect_size_mean_diff)
  })
  
  output$effect_size_max_diff <- renderText({
    req(rv$effect_size_max_diff)
    sprintf("%.3f", rv$effect_size_max_diff)
  })
  
  output$effect_size_interpretation <- renderPrint({
    req(rv$effect_size_corr, rv$effect_size_mean_diff, rv$effect_size_r_squared)
    
    cat("📊 ANÁLISIS DE TAMAÑO DEL EFECTO\n")
    cat("=====================================\n\n")
    
    cat("1. CORRELACIÓN DE SPEARMAN\n")
    cat(sprintf("   • Valor: %.4f\n", rv$effect_size_corr))
    cat(sprintf("   • R² (varianza compartida): %.2f%%\n", rv$effect_size_r_squared * 100))
    
    # Interpretación de la correlación
    if (rv$effect_size_corr > 0.90) {
      cat("   • Interpretación: Redes MUY SIMILARES\n")
    } else if (rv$effect_size_corr > 0.70) {
      cat("   • Interpretación: Similitud MODERADA-ALTA\n")
    } else if (rv$effect_size_corr > 0.50) {
      cat("   • Interpretación: Similitud MODERADA\n")
    } else {
      cat("   • Interpretación: Redes SUSTANCIALMENTE DIFERENTES\n")
    }
    
    cat("\n2. DIFERENCIAS EN PESOS DE ARISTAS\n")
    cat(sprintf("   • Diferencia absoluta media: %.4f\n", rv$effect_size_mean_diff))
    cat(sprintf("   • Diferencia máxima: %.4f\n", rv$effect_size_max_diff))
    cat(sprintf("   • Aristas con diferencia > 0.1: %d\n", rv$effect_size_n_diff))
    
    cat("\n3. TAMAÑO DEL EFECTO (Cohen's criteria adaptado)\n")
    # Interpretación basada en diferencia media
    if (rv$effect_size_mean_diff < 0.05) {
      cat("   • Efecto: PEQUEÑO (diferencia < 0.05)\n")
    } else if (rv$effect_size_mean_diff < 0.15) {
      cat("   • Efecto: MEDIANO (0.05 ≤ diferencia < 0.15)\n")
    } else {
      cat("   • Efecto: GRANDE (diferencia ≥ 0.15)\n")
    }
    
    cat("\n4. RESUMEN EJECUTIVO\n")
    if (rv$effect_size_corr > 0.90 && rv$effect_size_mean_diff < 0.05) {
      cat("   ✅ Las redes son prácticamente idénticas\n")
    } else if (rv$effect_size_corr > 0.70 && rv$effect_size_mean_diff < 0.15) {
      cat("   ⚠️ Las redes son similares con diferencias menores\n")
    } else {
      cat("   🔴 Las redes presentan diferencias sustanciales\n")
    }
  })
  
  # Gráfico de diferencias
  output$effect_size_plot <- renderPlot({
    req(rv$effect_size_diff_matrix, rv$networks_groups)
    
    # Convertir matriz de diferencias a formato largo
    diff_df <- rv$effect_size_diff_matrix %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Node1") %>%
      tidyr::pivot_longer(-Node1, names_to = "Node2", values_to = "Difference") %>%
      filter(Node1 < Node2)  # Solo triángulo inferior
    
    # Crear histograma de diferencias
    p1 <- ggplot(diff_df, aes(x = Difference)) +
      geom_histogram(bins = 30, fill = "#1e3a5f", alpha = 0.7, color = "white") +
      geom_vline(xintercept = rv$effect_size_mean_diff, 
                 color = "#CC1228", linetype = "dashed", size = 1.2) +
      labs(
        title = "Distribución de Diferencias entre Aristas",
        x = "Diferencia Absoluta |w₁ - w₂|",
        y = "Frecuencia"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      annotate("text", 
               x = rv$effect_size_mean_diff + 0.02, 
               y = Inf, 
               label = sprintf("Media = %.3f", rv$effect_size_mean_diff),
               vjust = 2, hjust = 0, color = "#CC1228", fontface = "bold")
    
    print(p1)
  })
  
  # Descarga del análisis
  output$download_effect_size <- downloadHandler(
    filename = function() paste0("Effect_Size_Analysis_", Sys.Date(), ".csv"),
    content = function(file) {
      # Crear data frame con todos los resultados
      results_df <- data.frame(
        Metrica = c("Correlación de Spearman", "R² (varianza compartida)", 
                    "Diferencia absoluta media", "Diferencia máxima", 
                    "Número de aristas diferentes (>0.1)"),
        Valor = c(rv$effect_size_corr, rv$effect_size_r_squared,
                  rv$effect_size_mean_diff, rv$effect_size_max_diff,
                  rv$effect_size_n_diff)
      )
      write.csv(results_df, file, row.names = FALSE)
    }
  )
  
  # Outputs de grupos
  output$groups_status <- renderText({
    if (rv$groups_running) "Analizando..."
    else if (!is.null(rv$networks_groups)) "Completado"
    else "Presiona 'Analizar'"
  })
  
  output$nct_status <- renderText({
    if (rv$nct_running) {
      elapsed <- if(!is.null(rv$nct_start_time)) {
        round(as.numeric(difftime(Sys.time(), rv$nct_start_time, units = "secs")))
      } else 0
      paste("Ejecutando NCT...", elapsed, "seg")
    } else if (!is.null(rv$nct_result)) {
      "NCT Completado"
    } else {
      "Esperando grupos"
    }
  })
  
  observe({
    invalidateLater(1000, session)
    
    if (rv$nct_running && !is.null(rv$nct_start_time)) {
      elapsed <- as.numeric(difftime(Sys.time(), rv$nct_start_time, units = "secs"))
      iterations <- if(is.null(input$nct_iterations)) 100 else input$nct_iterations
      
      estimated_total <- iterations * 0.5
      progress <- min(95, 10 + (elapsed / estimated_total) * 85)
      updateProgressBar(session, "nct_progress", value = progress)
    }
  })
  
  output$nct_results <- renderPrint({
    req(rv$nct_result)
    
    cat("🔬 NETWORK COMPARISON TEST - RESULTADOS\n")
    cat("======================================\n\n")
    
    cat("📊 NETWORK INVARIANCE TEST\n")
    cat("Test statistic M:", rv$nct_result$nwinv.real, "\n")
    cat("p-value:", rv$nct_result$nwinv.pval, "\n")
    
    if (rv$nct_result$nwinv.pval < 0.05) {
      cat("🔴 RESULTADO: Las redes son SIGNIFICATIVAMENTE DIFERENTES (p < 0.05)\n\n")
    } else {
      cat("🟢 RESULTADO: Las redes NO son significativamente diferentes (p ≥ 0.05)\n\n")
    }
    
    cat("💪 GLOBAL STRENGTH INVARIANCE TEST\n")
    cat("Test statistic S:", rv$nct_result$glstrinv.real, "\n")
    cat("p-value:", rv$nct_result$glstrinv.pval, "\n")
    
    if (rv$nct_result$glstrinv.pval < 0.05) {
      cat("🔴 RESULTADO: La conectividad global es SIGNIFICATIVAMENTE DIFERENTE (p < 0.05)\n\n")
    } else {
      cat("🟢 RESULTADO: La conectividad global NO es significativamente diferente (p ≥ 0.05)\n\n")
    }
    
    cat("📋 RESUMEN EJECUTIVO:\n")
    if (rv$nct_result$nwinv.pval < 0.05 || rv$nct_result$glstrinv.pval < 0.05) {
      cat("• Existe evidencia estadística de diferencias entre los grupos\n")
      if (rv$nct_result$nwinv.pval < 0.05) {
        cat("• Las estructuras de red difieren significativamente\n")
      }
      if (rv$nct_result$glstrinv.pval < 0.05) {
        cat("• La conectividad global difiere significativamente\n")
      }
    } else {
      cat("• No hay evidencia estadística de diferencias entre los grupos\n")
      cat("• Las redes pueden considerarse equivalentes estadísticamente\n")
    }
  })
  
  output$nct_plot <- renderPlot({
    req(rv$nct_result)
    
    data_nct <- data.frame(
      Test = c("Network\nInvariance", "Global\nStrength"),
      P_Value = c(rv$nct_result$nwinv.pval, rv$nct_result$glstrinv.pval),
      Significant = c(rv$nct_result$nwinv.pval < 0.05, rv$nct_result$glstrinv.pval < 0.05)
    )
    
    ggplot(data_nct, aes(x = Test, y = P_Value, fill = Significant)) +
      geom_col(width = 0.6, alpha = 0.8) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d9534f", size = 1.2) +
      scale_fill_manual(values = c("TRUE" = "#d9534f", "FALSE" = "#5cb85c")) +
      labs(
        title = "Network Comparison Test - p-values",
        y = "p-value",
        x = "Test Type",
        fill = "Significativo\n(p < 0.05)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#1A202C"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "#F7FAFC")
      ) +
      annotate("text", x = 1.5, y = 0.07, label = "α = 0.05", color = "#d9534f", size = 5, fontface = "bold")
  })
  
  output$download_nct <- downloadHandler(
    filename = function() paste0("NCT_Results_", Sys.Date(), ".jpg"),
    content = function(file) {
      req(rv$nct_result)
      
      data_nct <- data.frame(
        Test = c("Network\nInvariance", "Global\nStrength"),
        P_Value = c(rv$nct_result$nwinv.pval, rv$nct_result$glstrinv.pval),
        Significant = c(rv$nct_result$nwinv.pval < 0.05, rv$nct_result$glstrinv.pval < 0.05)
      )
      
      p_nct <- ggplot(data_nct, aes(x = Test, y = P_Value, fill = Significant)) +
        geom_col(width = 0.6, alpha = 0.8) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d9534f", size = 1.2) +
        scale_fill_manual(values = c("TRUE" = "#d9534f", "FALSE" = "#5cb85c")) +
        labs(
          title = "Network Comparison Test - p-values",
          y = "p-value",
          x = "Test Type",
          fill = "Significativo\n(p < 0.05)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#1A202C"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "#F7FAFC")
        ) +
        annotate("text", x = 1.5, y = 0.07, label = "α = 0.05", color = "#d9534f", size = 5, fontface = "bold")
      
      ggsave(file, p_nct, width = 8, height = 6, dpi = 300, device = "jpeg")
    }
  )
  
  output$networks_by_group_plot   <- renderPlot({ req(rv$combined_plot_groups);   print(rv$combined_plot_groups) })
  output$centrality_group_plot    <- renderPlot({ req(rv$plot_centralidad_group); print(rv$plot_centralidad_group) })
  output$bridge_group_plot        <- renderPlot({ req(rv$bridge_plot_group);      print(rv$bridge_plot_group$plot) })
  output$combined_groups_plot     <- renderPlot({ req(rv$combined_groups_final);  print(rv$combined_groups_final) })
  
  output$groups_render_status <- renderText({
    req(rv$combined_groups_final)
    
    tryCatch({
      if (!is.null(rv$combined_groups_final)) {
        "✅ Figura generada correctamente"
      } else {
        "⏳ Generando figura combinada..."
      }
      
    }, error = function(e) {
      paste("❌ Error:", e$message)
    })
  })
  
  output$download_groups <- downloadHandler(
    filename = function() paste0("Comparacion_Grupos_", Sys.Date(), ".jpg"),
    content  = function(file) {
      req(rv$combined_groups_final)
      
      tryCatch({
        width_val <- if(is.null(input$groups_width)) 14 else input$groups_width
        height_val <- if(is.null(input$groups_height)) 8 else input$groups_height
        dpi_val <- if(is.null(input$groups_dpi)) 300 else input$groups_dpi
        
        ggsave(file, rv$combined_groups_final, 
               width = width_val, height = height_val, dpi = dpi_val, device = "jpeg")
        
      }, error = function(e) {
        jpeg(file, width = 1400, height = 800, res = 100)
        plot.new()
        text(0.5, 0.5, paste("Error al generar descarga:", e$message), 
             cex = 1.2, col = "red")
        dev.off()
      })
    }
  )
} # Cierre del server

# ---- 7. Lanzar la app ----
shinyApp(ui, server)