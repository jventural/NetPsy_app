# ---- 1. INSTALACIÓN LIBRERÍAS ----
if (!require("devtools")) {
  install.packages("devtools")
}
# devtools::install_github("jventural/ThesiStats")
# devtools::install_github("jventural/PsyMetricTools")
# devtools::install_github("jventural/InterconectaR")
# remove.packages("InterconectaR")
# ---- 2. CARGA DE LIBRERÍAS ----
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
library(tidyverse)
library(bootnet)
library(mgm)
# ---- 3. CARGA Y REVISAR LA DATA ----
# Cargar la base de datos
df <- read_excel("df_completo.xlsx") %>% slice(1:1000)

# openxlsx::write.xlsx(Table_descriptive, file = "Table_descriptive.xlsx", overwrite = FALSE)
df_variables <-  df %>% select(Abandono:Bienestar)
df_variables
# Estadísticos descriptivos
Table_descriptive <- df_variables %>%
  psych::describe()

Table_descriptive

# Análisis de redundancia
goldbricker(df_variables, p = 0.05, method = "hittner2003", 
            threshold = 0.25, corMin = 0.50, progressbar = TRUE)

# Organización de las comunidades
names <- c("Dependencia", "Miedo Soledad", "Bienestar")
values <- c(3,1,1)
groups <- structure_groups(names, values)
groups

# ---- 5. ESTIMACION DE LA RED ----
network <- estimateNetwork(df_variables, 
                           default = "ggmModSelect", 
                           stepwise = T,
                           corMethod = "spearman")
centralityPlot()
network$graph
# Especificar los type y levels para estimar los R2
type <- c(rep("g", 5))
level <- c(rep(1, 5))

#Estimar los R2
error_Model <- mgm_error_metrics(data = df_variables,
                                 type = type,
                                 level = level)

# #Extraer solo los valores de interes
error_Model <- error_Model$errorCon$R2
error_Model

# Dibujar la red estimada
g1 <- qgraph(network$graph,
             groups = groups, 
             curveAll = 2,
             vsize = 12,
             esize = 12,
             palette = "ggplot2",
             layout = "spring", # Fruchterman-Reingold
             edge.labels = T, 
             pie = error_Model,
             layoutScale =c(0.8,0.8),
             legend.cex = 0.50,
             labels = c("Aba", "Auto", "Nec", "Mie", "Bie"))

# Generar el grafico de centralidad de la red estimada
# Caso A: solo ExpectedInfluence
Centralitys_A <- centrality_plots2(
  qgraph_obj    = g1,
  network       = network,
  groups        = groups,                       # <- aquí
  measure0      = "Strength",
  measure1      = NULL,
  color_palette = c("#FF0000"),
  labels        = NULL,
  legend_labels = NULL,
  use_abbrev = T
)
Centralitys_A$plot
# Caso B: ExpectedInfluence + Bridge Expected Influence (1-step)
centrality_plots2_fixed <- function(qgraph_obj,
                                    network,
                                    groups = NULL,
                                    measure0 = "ExpectedInfluence",
                                    measure1 = NULL,
                                    color_palette = c("#FF0000", "#00A08A"),
                                    labels = NULL,
                                    legend_labels = NULL,
                                    use_abbrev = TRUE) {
  
  # Cargar librerías necesarias
  library(qgraph)
  library(dplyr)
  library(networktools)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(purrr)
  
  # Extraer tabla de centralidad para measure0
  ctbl <- centralityTable(network)
  
  cents_expect <- ctbl %>%
    filter(measure == !!measure0) %>%
    select(full_name = node, value) %>%
    arrange(full_name) %>%
    rename(!!measure0 := value)
  
  # Etiquetas personalizadas (si las hay)
  if (!is.null(labels)) {
    cents_expect <- cents_expect %>%
      mutate(full_name = if_else(full_name %in% names(labels),
                                 labels[full_name],
                                 full_name))
  }
  
  # Abreviar "full_name" siempre (primeras 3 letras de la primera palabra)
  cents_expect <- cents_expect %>%
    mutate(
      Abrev = full_name %>%
        str_split(" ") %>%
        map_chr(~ .x[[1]]) %>%
        substr(1, 3) %>%
        str_to_title()
    )
  
  # Si measure1 no es NULL: calcular puente
  if (!is.null(measure1)) {
    # Ejecutar bridge()
    b_obj <- bridge(
      qgraph_obj,
      communities    = groups,
      useCommunities = "all",
      normalize      = FALSE
    )
    
    # CORRECCIÓN: Usar los labels del qgraph para crear el mapeo
    qgraph_labels <- qgraph_obj$graphAttributes$Nodes$labels
    
    # Crear data.frame bridge con el orden correcto
    bridge_data <- tibble(
      qgraph_label = qgraph_labels,
      raw_bridge = as.numeric(b_obj[[measure1]])
    ) %>%
      mutate(!!measure1 := as.numeric(scale(raw_bridge))) %>%
      select(-raw_bridge)
    
    # Crear un mapeo entre labels del qgraph y nombres completos
    # Esto asume que el orden es el mismo en ambos objetos
    full_names <- rownames(network$graph)
    if (is.null(full_names)) {
      full_names <- colnames(network$graph)
    }
    
    # Crear tabla de mapeo
    name_mapping <- tibble(
      full_name = full_names,
      qgraph_label = qgraph_labels
    )
    
    # Añadir nombres completos a bridge_data
    bridge_data <- bridge_data %>%
      left_join(name_mapping, by = "qgraph_label")
    
    # Crear abreviaturas consistentes para bridge_data
    bridge_data <- bridge_data %>%
      mutate(
        Abrev = full_name %>%
          str_split(" ") %>%
          map_chr(~ .x[[1]]) %>%
          substr(1, 3) %>%
          str_to_title()
      )
    
    # Unir con cents_expect por full_name
    cents2 <- inner_join(
      cents_expect,
      bridge_data %>% select(full_name, !!sym(measure1)),
      by = "full_name"
    )
    
    # Verificar que la unión fue exitosa
    if (nrow(cents2) == 0) {
      cat("ERROR: No se pudieron unir los datos. Intentando unión por posición...\n")
      
      # Plan B: unir por posición
      if (nrow(cents_expect) == nrow(bridge_data)) {
        cents2 <- bind_cols(
          cents_expect,
          bridge_data %>% select(!!sym(measure1))
        )
      } else {
        stop("No se pueden unir los datos de centralidad y puente")
      }
    }
    
    # Preparar datos en formato "largo" para ggplot (dos métricas)
    y_var <- if (use_abbrev) quo(Abrev) else quo(full_name)
    
    cents_long <- cents2 %>%
      pivot_longer(
        cols      = c(!!sym(measure0), !!sym(measure1)),
        names_to  = "Centrality",
        values_to = "Value"
      ) %>%
      rename(Measure = Centrality)
    
    # Definir paleta y etiquetas de leyenda
    if (is.null(legend_labels)) {
      pal         <- setNames(color_palette, c(measure1, measure0))
      legend_lbls <- setNames(c(measure1, measure0), c(measure1, measure0))
    } else {
      pal         <- setNames(color_palette, c(measure1, measure0))
      legend_lbls <- legend_labels
    }
    
    # Construir ggplot con dos líneas
    Figura <-  ggplot(
      cents_long,
      aes(
        x     = Value,
        y     = fct_reorder(!!y_var, Value),
        color = Measure,
        group = Measure
      )
    ) +
      geom_point(size = 3) +
      geom_line(size = 0.5) +
      theme_minimal() +
      labs(
        x     = "z-score",
        y     = "Nodos",
        color = "Métrica"
      ) +
      scale_color_manual(
        values = pal,
        labels = legend_lbls
      ) +
      theme(
        axis.text.y     = element_text(size = 12),
        axis.text.x     = element_text(size = 12),
        legend.text     = element_text(size = 12),
        legend.title    = element_text(size = 12),
        axis.title.y    = element_text(size = 12),
        axis.title.x    = element_text(size = 12),
        legend.position = "bottom"
      )
    
    # Asignar clase personalizada para suprimir advertencias al imprimir
    class(Figura) <- c("silent_gg", class(Figura))
    assign("print.silent_gg",
           function(x, ...) suppressWarnings(NextMethod()),
           envir = .GlobalEnv)
    
    # Devolver tabla y gráfico (caso dual)
    return(
      list(
        table = cents2 %>% arrange(desc(!!sym(measure0))),
        plot  = Figura
      )
    )
  }
  
  # Caso measure1 = NULL: solo measure0
  cents_single <- cents_expect %>%
    select(full_name, Abrev, !!sym(measure0)) %>%
    arrange(desc(!!sym(measure0))) %>%
    rename(Value = !!sym(measure0))
  
  # Preparar paleta y leyenda para una sola métrica
  pal_single         <- setNames(color_palette[1], measure0)
  legend_lbls_single <- setNames(measure0, measure0)
  
  # Construir ggplot para solo measure0
  y_var_single <- if (use_abbrev) quo(Abrev) else quo(full_name)
  
  Figura <- ggplot(
    cents_single,
    aes(
      x     = Value,
      y     = fct_reorder(!!y_var_single, Value),
      group = 1
    )
  ) +
    geom_line(color = pal_single, size = 0.5) +
    geom_point(color = pal_single, size = 3) +
    theme_minimal() +
    labs(
      x     = "z-score",
      y     = "Nodos",
      color = "Métrica"
    ) +
    scale_color_manual(
      values = pal_single,
      labels = legend_lbls_single
    ) +
    theme(
      axis.text.y     = element_text(size = 12),
      axis.text.x     = element_text(size = 12),
      legend.text     = element_text(size = 12),
      legend.title    = element_text(size = 12),
      axis.title.y    = element_text(size = 12),
      axis.title.x    = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Asignar clase personalizada para suprimir advertencias al imprimir
  class(Figura) <- c("silent_gg", class(Figura))
  assign("print.silent_gg",
         function(x, ...) suppressWarnings(NextMethod()),
         envir = .GlobalEnv)
  
  # Devolver tabla y gráfico (caso single)
  return(
    list(
      table = cents_single,
      plot  = Figura
    )
  )
}

Centralitys_B <- centrality_plots2_fixed(
  qgraph_obj    = g1,
  network       = network,
  groups        = groups,                       # <- aquí
  measure0      = "ExpectedInfluence",
  measure1      = "Bridge Expected Influence (1-step)",
  color_palette = c("#FF0000", "#00A08A"),
  labels        = NULL,
  legend_labels = NULL,
  use_abbrev = T
)
Centralitys_B$plot

# ---- 6. DESCRIPCIÓN DE LA RED ----

# # Describir la matriz adyacente
describe <- InterconectaR::get_edge_weights_summary(network) #corregir
describe
# 
# # Density
InterconectaR::Density_report(network$graph)

# Grafico combinado la red y centralidad
plot_combinet <- combine_graphs_centrality2(
  Figura1_Derecha = Centralitys_A$plot, 
  network = network, 
  groups = groups, 
  error_Model = error_Model, 
  ncol = 2, 
  widths = c(0.50, 0.25), 
  dpi = 300,
  legend.cex = 0.35,
  abbreviate_labels = T
)
plot_combinet
# Ahora puede guardar el resultado con ggsave:
ggsave(
  filename = "Figura_1_Final.jpg",
  plot     = plot_combinet,
  width    = 15,
  height   = 6.5,
  units    = "in",
  dpi      = 600
)
# ---- 7. EVALUACION DE LA ESTABILIDAD Y PRECISIÓN DE LA RED ----
# 1. Bootstrap network estimation , "expectedInfluence", "bridgeStrength"
caseDroppingBoot <- bootnet(network, 
                            boot = 10, 
                            type = "case", #Esto cambio
                            nCores = 12,
                            statistics = c("strength"),
                            communities=groups)

result_case <- InterconectaR::filter_correlation_stability(caseDroppingBoot)
result_case

# Combining Precision
nonParametricBoot <- bootnet(network, 
                             boot = 10, 
                             type = "nonparametric", #Esto cambio
                             nCores = 12,
                             statistics = "all",
                             communities=groups)

# Combining Stability and Accuracy
combined_plot <- InterconectaR::plot_centrality_stability(caseDroppingBoot, 
                                           nonParametricBoot,
                                           statistics = c("strength"))
combined_plot
ggsave("Figura_2.jpg", combined_plot, width = 9, height = 5, dpi = 600)

# ---- 8. COMPARACIÓN DE REDES ----
networks_groups <- InterconectaR::estimate_networks_by_group(data = df,
                                               group_var = "Sex",
                                               columns = df %>% select(Abandono:Bienestar) %>% names(),
                                               default = "ggmModSelect", 
                                               stepwise = TRUE,        
                                               corMethod = "spearman",
                                               abbreviate_vars  = T,
                                               abbr_minlength = 3)

networks_groups$Mujer$data
networks_groups$Varon$data

# Gráfico de centralidad por grupo
plot_centralidad_group <- InterconectaR::plot_centrality_by_group(networks_groups,
                                                   replacements  = c("Mujer", "Varon"),
                                                   measure_spec  = c("ExpectedInfluence"),
                                                   color_palette = c("#FF5733", "#33FFCE"))
plot_centralidad_group

# Cálculo de errores de redes por grupo
errores <- InterconectaR::mgm_errors_groups(
  data = df,
  type = c(rep("g", 5)),
  level = c(rep(1, 5)),
  group = Sex,
  columns = df %>% select(Abandono:Bienestar) %>% names()
)

errores$Mujer$R2
errores$Varon$R2

# Valores para gráficos de "pie"
pie_values <- list(
  "Varon" = errores$Varon$R2,
  "Mujer" = errores$Mujer$R2
)

# Gráfico de redes por grupo
# Supongamos que networks_by_group ya está definido, por ejemplo:
combined_plot22 <- InterconectaR::plot_networks_by_group(
  res = 300,
  networks_by_group = networks_groups,
  groups            = groups,
  pie               = pie_values,
  legend.cex        = 0.8
)
combined_plot22

# Guardar con ggsave (ejemplo: 14" × 8", 600 dpi)
ggsave(
  filename = "combined_networks_by_group.jpg",
  plot     = combined_plot22,
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 600
)

# Gráfico de centralidad y puentes
bridge_plot_group <- InterconectaR::centrality_bridge_plot(
  networks_groups = networks_groups,  
  group_names = c("Mujer", "Varon"), 
  measure = "Bridge Expected Influence (1-step)",
  color_palette = c("#FF5733", "#33FFCE"))

bridge_plot_group$plot

# Combinación de gráficos
combinado22 <- InterconectaR::combine_groupBy(
  red_group = combined_plot22,
  plot_centralidad_group = plot_centralidad_group,
  bridge_plot_group = bridge_plot_group$plot,
  width_a  = 12,   # darle más espacio horizontal a A
  width_bc = 4.5,
  show_plot = TRUE
)
combinado22
# Para guardarlo:
ggsave(
  filename = "figura_combinada.jpg",
  plot     = combinado22,
  width    = 14,
  height   = 8,
  dpi      = 300
)



set.seed(234) 
library(NetworkComparisonTest)
res <-NCT(networks_groups$Mujer, networks_groups$Varon, binary.data=F,
          test.edges=TRUE, edges="all", it = 100, test.centrality = TRUE)

#NETWORK INVARIANCE TEST 
res$nwinv.real # Test statistic M:  
res$nwinv.pval #p-value 

#GLOBAL STRENGTH INVARIANCE TEST 
res$glstrinv.real #Test statistic S
res$glstrinv.pval #p-value


#Calculation of effect size - Comparing adjacent matrices
network1_adjacency <- getWmat(networks_groups$Mujer)
network2_adjacency <- getWmat(networks_groups$Varon)

cor(network1_adjacency[lower.tri(network1_adjacency)], network2_adjacency[lower.tri(network2_adjacency)], method = "spearman")



res.mujer <- myboots2(networks_groups$Mujer$data, nBoots = 100, seed = 2024, estimator = "ggmModSelect", corMethod = "spearman")
res.varon <- myboots2(networks_groups$Varon$data, nBoots = 100, seed = 2024, estimator = "ggmModSelect", corMethod = "spearman")

tabla_resultados_bootstrap <- generate_table_bootstrap(res.mujer, res.varon)

# Ver los resultados
print(tabla_resultados_bootstrap)

# Ejecutar el cálculo y obtener el promedio de q de Cohen
result_qcohen_summary <- myboots_qcohen_summary(networks_groups$Mujer$data, networks_groups$Varon$data,
                                                nBoots = 100, 
                                                estimator = "ggmModSelect", corMethod = "spearman", 
                                                seed = 2024)

# Visualizar los resultados
result_qcohen_summary %>% arrange(Absolute_q)

result <- calculate_effect_sizes(networks_groups$Mujer$data, networks_groups$Varon$data, 
                                 nBoots = 100, seed = 123)

head(result$q_results)
head(result$delta_z_results)
head(result$absolute_diff_results)
head(result$relative_diff_results)
result$global_results$Metric

