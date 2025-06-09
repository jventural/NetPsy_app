# =============================================================================
# ANÁLISIS DE REDES PSICOLÓGICAS: DEPENDENCIA, MIEDO A LA SOLEDAD Y BIENESTAR
# Versión: optimizado y estructurado
# Autor: [Tu Nombre]
# Fecha: [Fecha]
# =============================================================================

# -----------------------------------------------------------------------------
# 1. CONFIGURACIÓN: INSTALACIÓN Y CARGA DE PAQUETES
# -----------------------------------------------------------------------------

# 1.1 Definir paquetes CRAN obligatorios
paquetes_cran <- c(
  "readxl", "dplyr",    "tidyverse",  "psych",      "huge",        "corrplot",
  "matrixcalc", "networktools", "qgraph", "bootnet",    "mgm",         "openxlsx",
  "NetworkComparisonTest"
)

# 1.2 Definir paquetes de GitHub (si no están instalados)
paquetes_github <- c(
  "jventural/ThesiStats",
  "jventural/PsyMetricTools",
  "jventural/InterconectaR"
)

# 1.3 Instalar paquetes CRAN faltantes
instalar_si_falta <- function(pkgs) {
  faltantes <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(faltantes)) {
    install.packages(faltantes, dependencies = TRUE)
  }
}
instalar_si_falta(paquetes_cran)

# 1.4 Instalar paquete 'devtools' si falta, para poder instalar paquetes GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", dependencies = TRUE)
}

# 1.5 Instalar paquetes GitHub (solo si aún no existen)
for (repo in paquetes_github) {
  nombre_pkg <- sub(".*/", "", repo)
  if (!(nombre_pkg %in% installed.packages()[, "Package"])) {
    devtools::install_github(repo)
  }
}

# 1.6 Cargar todas las librerías
library(readxl)
library(dplyr)
library(tidyverse)
library(psych)
library(huge)
library(corrplot)
library(matrixcalc)
library(networktools)
library(qgraph)
library(bootnet)
library(mgm)
library(openxlsx)
library(NetworkComparisonTest)
# Paquetes de GitHub
library(ThesiStats)         # funciones personalizadas para análisis estadístico
library(PsyMetricTools)     # herramientas psicométricas
library(InterconectaR)      # funciones para describir y comparar redes

# -----------------------------------------------------------------------------
# 2. CARGA Y PREPARACIÓN DE LOS DATOS
# -----------------------------------------------------------------------------

# 2.1 Leer el archivo Excel (solo primeras 1000 filas para propósitos de ejemplo)
df <- read_excel("df_completo.xlsx") %>%
  slice(1:1000)

# 2.2 Seleccionar solo las columnas de interés (Abandono a Bienestar)
vars_interes <- names(df)[which(names(df) == "Abandono") : which(names(df) == "Bienestar")]
df_variables <- df %>% select(all_of(vars_interes))

# 2.3 Mostrar las primeras filas y estatísticos descriptivos
print(head(df_variables))
Table_descriptive <- df_variables %>% psych::describe()
print(Table_descriptive)

# -----------------------------------------------------------------------------
# 3. ANÁLISIS DE REDUNDANCIA Y ESTRUCTURA DE GRUPOS
# -----------------------------------------------------------------------------

# 3.1 Identificar ítems redundantes con la función goldbricker
goldbricker(
  df_variables,
  p            = 0.05,
  method       = "hittner2003",
  threshold    = 0.25,
  corMin       = 0.50,
  progressbar  = TRUE
)

# 3.2 Definir agrupaciones teóricas de ítems (p. ej. Dependencia, Miedo Soledad, Bienestar)
nombres_grupos <- c("Dependencia", "Miedo Soledad", "Bienestar")
valores_grupos  <- c(3, 1, 1)
groups <- structure_groups(nombres_grupos, valores_grupos)
print(groups)

# -----------------------------------------------------------------------------
# 4. ESTIMACIÓN DE LA RED GLOBAL
# -----------------------------------------------------------------------------

# 4.1 Estimar la red usando ggmModSelect con correlación de Spearman
network <- estimateNetwork(
  df_variables,
  default   = "ggmModSelect",
  stepwise  = TRUE,
  corMethod = "spearman"
)

# 4.2 Extraer la matriz de adyacencia y mostrarla
matriz_adyacencia <- network$graph
print(matriz_adyacencia)

# 4.3 Calcular R2 por nodo usando mgm_error_metrics
#      (aquí asumimos variables continuas: type = "g", level = 1)
type  <- rep("g", ncol(df_variables))
level <- rep(1, ncol(df_variables))

error_Model_full <- mgm_error_metrics(
  data   = df_variables,
  type   = type,
  level  = level
)
error_Model <- error_Model_full$errorCon$R2
print(error_Model)

# 4.4 Dibujar la red estimada con qgraph
g1 <- qgraph(
  matriz_adyacencia,
  groups          = groups,
  curveAll        = 2,
  vsize           = 12,
  esize           = 12,
  palette         = "ggplot2",
  layout          = "spring",       # Fruchterman-Reingold
  edge.labels     = TRUE,
  pie             = error_Model,
  layoutScale     = c(0.8, 0.8),
  legend.cex      = 0.5,
  labels          = c("Aba", "Auto", "Nec", "Mie", "Bie")
)

# -----------------------------------------------------------------------------
# 5. GRÁFICOS DE CENTRALIDAD
# -----------------------------------------------------------------------------

# 5.1 ExpectedInfluence + Bridge Expected Influence (1-step)
Centralitys_B <- centrality_plots(
  qgraph_obj    = g1,
  network       = network,
  groups        = groups,
  measure0      = "ExpectedInfluence",
  measure1      = "Bridge Expected Influence (1-step)",
  color_palette = c("#FF0000", "#00A08A"),
  use_abbrev    = TRUE
)
print(Centralitys_B$plot)

# -----------------------------------------------------------------------------
# 6. DESCRIPCIÓN DE LA RED GLOBAL
# -----------------------------------------------------------------------------

# 6.1 Resumen de pesos de aristas
resumen_aristas <- InterconectaR::get_edge_weights_summary(network)
print(resumen_aristas)

# 6.2 Informe de densidad
densidad <- InterconectaR::Density_report(matriz_adyacencia)
print(densidad)

# 6.3 Gráfico combinado: red + centralidad
plot_combinet <- combine_graphs_centrality(
  Figura1_Derecha     = Centralitys_B$plot,
  network             = network,
  groups              = groups,
  error_Model         = error_Model,
  ncol                = 2,
  widths              = c(0.50, 0.25),
  dpi                 = 300,
  legend.cex          = 0.35,
  abbreviate_labels   = TRUE
)
print(plot_combinet)

# 6.4 Guardar figura final
ggsave(
  filename = "Figura_1_Final.jpg",
  plot     = plot_combinet,
  width    = 15,
  height   = 6.5,
  units    = "in",
  dpi      = 600
)

# -----------------------------------------------------------------------------
# 7. ESTABILIDAD Y PRECISIÓN DE LA RED GLOBAL
# -----------------------------------------------------------------------------

# 7.1 Bootstrap con remoción de casos (case-dropping)
caseDroppingBoot <- bootnet(
  network,
  boot       = 10,
  type       = "case",
  nCores     = 12,
  statistics = c("strength", "expectedInfluence", "bridgeStrength"),
  communities = groups
)
estabilidad_centralidad <- InterconectaR::filter_correlation_stability(caseDroppingBoot)
print(estabilidad_centralidad)

# 7.2 Bootstrap no paramétrico (precisión)
nonParametricBoot <- bootnet(
  network,
  boot       = 10,
  type       = "nonparametric",
  nCores     = 12,
  statistics = "all",
  communities = groups
)

# 7.3 Gráfico combinado de estabilidad y precisión
combined_plot <- plot_centrality_stability(
  caseDroppingBoot,
  nonParametricBoot,
  statistics = c("ExpectedInfluence", "strength", "bridgeStrength")
)
print(combined_plot)

# 7.4 Guardar figura
ggsave(
  filename = "Figura_2.jpg",
  plot     = combined_plot,
  width    = 9,
  height   = 5,
  dpi      = 600
)

# -----------------------------------------------------------------------------
# 8. COMPARACIÓN DE REDES POR GRUPO (p.ej. SEXO)
# -----------------------------------------------------------------------------

# 8.1 Estimar redes separadas por grupo (variable: Sex)
networks_groups <- estimate_networks_by_group(
  data            = df,
  group_var       = "Sex",
  columns         = vars_interes,
  default         = "ggmModSelect",
  stepwise        = TRUE,
  corMethod       = "spearman",
  abbreviate_vars = TRUE,
  abbr_minlength  = 3
)

# 8.2 Verificar estructura de datos en cada grupo
print(networks_groups$Mujer$data)
print(networks_groups$Varon$data)

# 8.3 Gráfico de centralidad por grupo (sólo ExpectedInfluence)
plot_centralidad_group <- plot_centrality_by_group(
  networks_groups = networks_groups,
  replacements    = c("Mujer", "Varon"),
  measure_spec    = c("ExpectedInfluence"),
  color_palette   = c("#FF5733", "#33FFCE")
)
print(plot_centralidad_group)

# 8.4 Cálculo de errores (R2) por grupo con mgm_errors_groups
errores <- mgm_errors_groups(
  data    = df,
  type    = rep("g", length(vars_interes)),
  level   = rep(1, length(vars_interes)),
  group   = Sex,
  columns = vars_interes
)

R2_Mujer <- errores$Mujer$R2
R2_Varon <- errores$Varon$R2
print(R2_Mujer)
print(R2_Varon)

# 8.5 Preparar valores "pie" para cada grupo
pie_values <- list(
  "Varon" = R2_Varon,
  "Mujer" = R2_Mujer
)

# 8.6 Gráfico de redes por grupo (usando plot_networks_by_group)
combined_plot22 <- plot_networks_by_group(
  res               = 300,
  networks_by_group = networks_groups,
  groups            = groups,
  pie               = pie_values,
  legend.cex        = 0.8
)
print(combined_plot22)

# 8.7 Guardar figura de redes por grupo
ggsave(
  filename = "combined_networks_by_group.jpg",
  plot     = combined_plot22,
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 600
)

# 8.8 Gráfico de centralidad de puentes por grupo
bridge_plot_group <- centrality_bridge_plot(
  networks_groups = networks_groups,
  group_names     = c("Mujer", "Varon"),
  measure         = "Bridge Expected Influence (1-step)",
  color_palette   = c("#FF5733", "#33FFCE")
)
print(bridge_plot_group$plot)

# 8.9 Combinar gráficos: redes, centralidad y puentes
combinado22 <- combine_groupBy(
  red_group               = combined_plot22,
  plot_centralidad_group  = plot_centralidad_group,
  bridge_plot_group       = NULL,           # Ajustar si se desea incluir bridge_plot_group
  width_a                 = 12,             # Espacio horizontal para el primer panel
  width_bc                = 4.5,
  show_plot               = TRUE
)
print(combinado22)

# 8.10 Guardar figura combinada
ggsave(
  filename = "figura_combinada.jpg",
  plot     = combinado22,
  width    = 14,
  height   = 8,
  dpi      = 300
)

# -----------------------------------------------------------------------------
# 9. PRUEBA DE COMPARACIÓN DE REDES (Network Comparison Test)
# -----------------------------------------------------------------------------

# 9.1 Configurar semilla para reproducibilidad
set.seed(234)

# 9.2 Ejecutar NCT entre grupo Mujer y Varon
res_nct <- NCT(
  networks_groups$Mujer,
  networks_groups$Varon,
  binary.data      = FALSE,
  test.edges       = TRUE,
  edges            = "all",
  it               = 100,
  test.centrality  = TRUE
)

# 9.3 Imprimir resultados de invariancia de red y de fuerza global
cat("TEST DE INVARIANCIA DE RED:\n")
cat("  Estadístico M =", res_nct$nwinv.real, "\n")
cat("  p-valor       =", res_nct$nwinv.pval, "\n\n")

cat("TEST DE FUERZA GLOBAL:\n")
cat("  Estadístico S =", res_nct$glstrinv.real, "\n")
cat("  p-valor       =", res_nct$glstrinv.pval, "\n\n")

# 9.4 Tamaño del efecto: correlación de matrices de adyacencia
network1_adjacency <- getWmat(networks_groups$Mujer)
network2_adjacency <- getWmat(networks_groups$Varon)

corr_adjacent <- cor(
  network1_adjacency[lower.tri(network1_adjacency)],
  network2_adjacency[lower.tri(network2_adjacency)],
  method = "spearman"
)
cat("Correlación entre matrices de adyacencia (efecto):", corr_adjacent, "\n")
