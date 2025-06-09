# ------------------------------------------------------------------------------
# 1. Instalar (si hace falta) y cargar paquetes
# ------------------------------------------------------------------------------
required_pkgs <- c(
  "ThesiStats",   # funciones psicométricas y utilidades
  "tidyverse",    # dplyr, tidyr, ggplot2, etc.
  "readxl",       # para leer .xlsx
  "openxlsx",     # para escribir .xlsx
  "pwr"           # cálculo de tamaño de muestra
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

if (!require("devtools")) {
  install.packages("devtools")
}
# devtools::install_github("jventural/ThesiStats")

# ------------------------------------------------------------------------------
# 2. Leer datos
# ------------------------------------------------------------------------------
input_file <- "Dependencia_new.xlsx"
df_raw <- read_excel(input_file)

# ------------------------------------------------------------------------------
# 3. Asignar y revisar etiquetas de columnas
# ------------------------------------------------------------------------------
labels(df_raw)

# ------------------------------------------------------------------------------
# 8. Lectura de texto y generación de nuevas variables
# ------------------------------------------------------------------------------
texto <- readLines("Texto.txt", encoding = "UTF-8")
generate_and_apply(df_raw, texto, new_name = "df_completo")

df_completo <- df_completo %>% select(Age,Sex,Abandono:Bienestar)

openxlsx::write.xlsx(df_completo, file = "df_completo.xlsx", overwrite = T)

