# librerias
library(sna)
library(qgraph)
library(MASS)
# Configurar semilla para reproducibilidad
set.seed(42)

# Definir matriz de correlación deseada (11 variables)
cor_matrix <- matrix(0.6, nrow = 11, ncol = 11) # Correlación base de 0.6
diag(cor_matrix) <- 1 # Correlación perfecta en la diagonal

# Generar datos con estructura correlacional
n <- 200 # Número de observaciones

# Aplicar descomposición de Cholesky para generar datos correlacionados
data_correlated <- mvrnorm(n = n, mu = rep(0, 11), Sigma = cor_matrix)

# Convertir a data frame y asignar nombres de columnas
datarem_11 <- as.data.frame(data_correlated)
colnames(datarem_11) <- c("agi", "con", "dep", "ene", "gui", "hyp", "ins", "int", "ret", "sui", "wap")

# Generar otro conjunto con las mismas correlaciones
datanonrem_11 <- as.data.frame(mvrnorm(n = n, mu = rep(0, 11), Sigma = cor_matrix))
colnames(datanonrem_11) <- c("agi", "con", "dep", "ene", "gui", "hyp", "ins", "int", "ret", "sui", "wap")

# Mostrar las primeras filas de los datos para verificar
head(datarem_11)
head(datanonrem_11)

myboots2 <- function(data, nBoots, estimator = "ggmModSelect", corMethod = "spearman", 
                     stepwise = NULL, seed = NULL) {
  library(bootnet)
  
  # Establecer la semilla si se proporciona
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  eigv <- matrix(NA, nBoots, ncol(data))
  colnames(eigv) <- colnames(data)
  
  for (b in seq_len(nBoots)) {
    # Crear muestra bootstrap
    bootData <- data[sample(seq_len(nrow(data)), nrow(data), replace = TRUE), ]
    
    # Preparar argumentos para estimateNetwork
    est_args <- list(data = bootData, default = estimator, corMethod = corMethod)
    if (!is.null(stepwise)) {
      est_args$stepwise <- stepwise
    }
    
    # Crear la red utilizando el estimador seleccionado
    network <- do.call(estimateNetwork, est_args)
    
    # Extraer matriz de pesos
    res <- getWmat(network)
    
    if (sum(res != 0) != 0) {
      # Calcular centralidad de eigenvector
      evc <- evcent(res, maxiter = 1e6)
      evc <- evc[order(colnames(res), decreasing = FALSE)]
      eigv[b, ] <- evc
    }
    print(paste("Bootstrap Iteration:", b))
  }
  
  return(list(eigv = eigv))
}

# Realizar bootstrap para ambos grupos
resbootr <- myboots2(datarem_11, nBoots = 100, seed = 2024, estimator = "EBICglasso", corMethod = "spearman")
resbootp <- myboots2(datanonrem_11, nBoots = 100, seed = 2024, estimator = "EBICglasso", corMethod = "spearman")

# Función para calcular tamaños del efecto y generar tabla
generate_table_bootstrap <- function(boots1, boots2) {
  results <- data.frame()
  
  for (var in colnames(boots1$eigv)) {
    # Extraer eigenvector centrality de cada grupo para la variable
    E1 <- boots1$eigv[, var]
    E2 <- boots2$eigv[, var]
    
    # Calcular medias y desviaciones estándar
    M1 <- mean(E1, na.rm = TRUE)
    DE1 <- sd(E1, na.rm = TRUE)
    M2 <- mean(E2, na.rm = TRUE)
    DE2 <- sd(E2, na.rm = TRUE)
    
    # Prueba t de Welch
    t_test <- t.test(E1, E2)
    t_val <- t_test$statistic
    gl <- t_test$parameter
    p_val <- t_test$p.value
    
    # Calcular d de Cohen
    pooled_sd <- sqrt(((length(E1) - 1) * DE1^2 + (length(E2) - 1) * DE2^2) / 
                        (length(E1) + length(E2) - 2))
    d_cohen <- (M1 - M2) / pooled_sd
    
    # Interpretación de d de Cohen
    interpret_d <- ifelse(abs(d_cohen) < 0.2, "Trivial",
                          ifelse(abs(d_cohen) < 0.5, "Small",
                                 ifelse(abs(d_cohen) < 0.8, "Medium", "Large")))
    
    # Agregar fila a la tabla
    results <- rbind(results, data.frame(
      Variable = var,
      M_Grupo1 = round(M1, 2),
      DE_Grupo1 = round(DE1, 2),
      M_Grupo2 = round(M2, 2),
      DE_Grupo2 = round(DE2, 2),
      t = round(t_val, 2),
      gl = round(gl, 2),
      p = round(p_val, 3),
      d = round(d_cohen, 2),
      Interpretation = interpret_d
    ))
  }
  
  return(results)
}

tabla_resultados_bootstrap <- generate_table_bootstrap(resbootr, resbootp)

# Ver los resultados
print(tabla_resultados_bootstrap)

datarem_11
datanonrem_11
###############
calculate_qcohen <- function(res1, res2, n1, n2) {
  q_matrix <- matrix(NA, nrow = ncol(res1), ncol = ncol(res1),
                     dimnames = list(colnames(res1), colnames(res1)))
  
  for (i in seq_len(ncol(res1))) {
    for (j in seq_len(ncol(res1))) {
      # Correlaciones entre pares de variables en ambos grupos
      r1 <- res1[i, j]
      r2 <- res2[i, j]
      
      # Evitar problemas numéricos si r = 1 o r = -1
      r1 <- ifelse(abs(r1) == 1, sign(r1) * 0.9999, r1)
      r2 <- ifelse(abs(r2) == 1, sign(r2) * 0.9999, r2)
      
      # Transformación Z de Fisher
      Z1 <- if (!is.na(r1)) 0.5 * log((1 + r1) / (1 - r1)) else NA
      Z2 <- if (!is.na(r2)) 0.5 * log((1 + r2) / (1 - r2)) else NA
      
      # Calcular q de Cohen
      q_matrix[i, j] <- if (!is.na(Z1) && !is.na(Z2)) Z1 - Z2 else NA
    }
  }
  
  return(q_matrix)
}

myboots_qcohen_summary <- function(data1, data2, nBoots, estimator = "ggmModSelect", corMethod = "spearman", 
                                   stepwise = NULL, seed = NULL) {
  library(bootnet)
  var_names <- colnames(data1)
  
  # Establecer la semilla si se proporciona
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Inicializar array para guardar valores de q en cada iteración
  q_values <- array(NA, dim = c(nBoots, length(var_names), length(var_names)),
                    dimnames = list(NULL, var_names, var_names))
  
  for (b in seq_len(nBoots)) {
    # Crear muestras bootstrap para ambos grupos
    bootData1 <- data1[sample(seq_len(nrow(data1)), nrow(data1), replace = TRUE), ]
    bootData2 <- data2[sample(seq_len(nrow(data2)), nrow(data2), replace = TRUE), ]
    
    # Preparar argumentos para estimateNetwork
    est_args <- list(data = bootData1, default = estimator, corMethod = corMethod)
    if (!is.null(stepwise)) {
      est_args$stepwise <- stepwise
    }
    
    # Crear redes utilizando bootnet
    network1 <- do.call(estimateNetwork, est_args)
    
    est_args$data <- bootData2
    network2 <- do.call(estimateNetwork, est_args)
    
    # Extraer matrices de correlaciones de ambas redes
    res1 <- getWmat(network1)
    res2 <- getWmat(network2)
    
    if (sum(res1 != 0) != 0 && sum(res2 != 0) != 0) {
      for (i in seq_len(ncol(res1))) {
        for (j in seq_len(ncol(res1))) {
          # Correlaciones entre pares de variables en ambos grupos
          r1 <- res1[i, j]
          r2 <- res2[i, j]
          
          # Evitar problemas numéricos si r = 1 o r = -1
          r1 <- ifelse(abs(r1) == 1, sign(r1) * 0.9999, r1)
          r2 <- ifelse(abs(r2) == 1, sign(r2) * 0.9999, r2)
          
          # Transformación Z de Fisher
          Z1 <- if (!is.na(r1)) 0.5 * log((1 + r1) / (1 - r1)) else NA
          Z2 <- if (!is.na(r2)) 0.5 * log((1 + r2) / (1 - r2)) else NA
          
          # Calcular q de Cohen
          q_values[b, i, j] <- if (!is.na(Z1) && !is.na(Z2)) Z1 - Z2 else NA
        }
      }
    }
    print(paste("Bootstrap Iteration:", b))
  }
  
  # Calcular el promedio de q de Cohen para cada par de variables
  avg_q <- apply(q_values, c(2, 3), mean, na.rm = TRUE)
  
  # Convertir el resultado en un data.frame
  result <- as.data.frame(as.table(avg_q))
  colnames(result) <- c("Variable1", "Variable2", "Mean_q")
  
  # Eliminar combinaciones repetidas (ejemplo: agi -- agi)
  result <- result[result$Variable1 != result$Variable2, ]
  
  # Separar el signo y el valor absoluto de Mean_q
  result$Sign <- ifelse(result$Mean_q >= 0, "+", "-")
  result$Absolute_q <- abs(result$Mean_q)
  
  # Agregar interpretación basada en los valores de Absolute_q
  result$Interpretation <- cut(
    result$Absolute_q,
    breaks = c(-Inf, 0.10, 0.30, 0.50, Inf),
    labels = c("Trivial", "Small", "Medium", "Large"),
    right = FALSE
  )
  
  # Reordenar columnas para mayor claridad
  result <- result[, c("Variable1", "Variable2", "Sign", "Absolute_q", "Interpretation")]
  
  return(result)
}


# Ejecutar el cálculo y obtener el promedio de q de Cohen
result_qcohen_summary <- myboots_qcohen_summary(datarem_11, datanonrem_11, nBoots = 100, 
                                                estimator = "EBICglasso", corMethod = "cor_auto", 
                                                seed = 2024)

# Visualizar los resultados
print(result_qcohen_summary)


calculate_effect_sizes <- function(data1, data2, nBoots, estimator = "ggmModSelect", 
                                   corMethod = "spearman", stepwise = NULL, seed = NULL) {
  library(bootnet)
  library(dplyr)
  
  var_names <- colnames(data1)
  
  # Establecer la semilla si se proporciona
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Inicializar estructuras para guardar los resultados de cada medida
  q_values <- array(NA, dim = c(nBoots, length(var_names), length(var_names)),
                    dimnames = list(NULL, var_names, var_names))
  frobenius_distances <- numeric(nBoots)
  tucker_congruences <- numeric(nBoots)
  delta_z_matrices <- array(NA, dim = c(nBoots, length(var_names), length(var_names)),
                            dimnames = list(NULL, var_names, var_names))
  absolute_diffs <- array(NA, dim = c(nBoots, length(var_names), length(var_names)),
                          dimnames = list(NULL, var_names, var_names))
  relative_diffs <- array(NA, dim = c(nBoots, length(var_names), length(var_names)),
                          dimnames = list(NULL, var_names, var_names))
  procrustes_similarities <- numeric(nBoots)
  pearson_correlations <- numeric(nBoots)
  spearman_correlations <- numeric(nBoots)
  
  for (b in seq_len(nBoots)) {
    # Crear muestras bootstrap para ambos grupos
    bootData1 <- data1[sample(seq_len(nrow(data1)), nrow(data1), replace = TRUE), ]
    bootData2 <- data2[sample(seq_len(nrow(data2)), nrow(data2), replace = TRUE), ]
    
    # Preparar argumentos para estimateNetwork
    est_args <- list(data = bootData1, default = estimator, corMethod = corMethod)
    if (!is.null(stepwise)) {
      est_args$stepwise <- stepwise
    }
    
    # Crear redes utilizando bootnet
    network1 <- do.call(estimateNetwork, est_args)
    est_args$data <- bootData2
    network2 <- do.call(estimateNetwork, est_args)
    
    # Extraer matrices de correlaciones de ambas redes
    res1 <- getWmat(network1)
    res2 <- getWmat(network2)
    network1_vector <- as.vector(res1)
    network2_vector <- as.vector(res2)
    
    if (sum(res1 != 0) != 0 && sum(res2 != 0) != 0) {
      # Calcular q de Cohen, ΔZ, y diferencias para cada par de variables
      for (i in seq_len(ncol(res1))) {
        for (j in seq_len(ncol(res1))) {
          r1 <- res1[i, j]
          r2 <- res2[i, j]
          
          # Evitar problemas numéricos si r = 1 o r = -1
          r1 <- ifelse(abs(r1) == 1, sign(r1) * 0.9999, r1)
          r2 <- ifelse(abs(r2) == 1, sign(r2) * 0.9999, r2)
          
          # Transformación Z de Fisher
          Z1 <- if (!is.na(r1)) 0.5 * log((1 + r1) / (1 - r1)) else NA
          Z2 <- if (!is.na(r2)) 0.5 * log((1 + r2) / (1 - r2)) else NA
          
          # Calcular q de Cohen
          q_values[b, i, j] <- if (!is.na(Z1) && !is.na(Z2)) Z1 - Z2 else NA
          
          # Calcular ΔZ de Fisher
          delta_z_matrices[b, i, j] <- if (!is.na(Z1) && !is.na(Z2)) abs(Z1 - Z2) else NA
          
          # Diferencia absoluta
          absolute_diffs[b, i, j] <- abs(r1 - r2)
          
          # Diferencia relativa
          relative_diffs[b, i, j] <- if (!is.na(r1) && !is.na(r2)) abs(r1 - r2) / mean(c(r1, r2)) else NA
        }
      }
      
      # Calcular Distancia de Frobenius normalizada
      num_elements <- ncol(res1) * nrow(res1) # Número total de elementos en la matriz
      frobenius_distances[b] <- sqrt(sum((res1 - res2)^2, na.rm = TRUE)) / num_elements
      
      # Calcular Coeficiente de Congruencia Tucker
      tucker_congruences[b] <- sum(res1 * res2) / (sqrt(sum(res1^2)) * sqrt(sum(res2^2)))
      
      # Calcular similaridad de Procrustes
      procrustes_similarities[b] <- sum(network1_vector * network2_vector) /
        (sqrt(sum(network1_vector^2)) * sqrt(sum(network2_vector^2)))
      
      # Calcular correlaciones globales
      pearson_correlations[b] <- cor(network1_vector, network2_vector, method = "pearson", use = "complete.obs")
      spearman_correlations[b] <- cor(network1_vector, network2_vector, method = "spearman", use = "complete.obs")
    }
    print(paste("Bootstrap Iteration:", b))
  }
  
  # Resumir resultados
  frobenius_mean <- mean(frobenius_distances, na.rm = TRUE)
  tucker_mean <- mean(tucker_congruences, na.rm = TRUE)
  procrustes_mean <- mean(procrustes_similarities, na.rm = TRUE)
  pearson_mean <- mean(pearson_correlations, na.rm = TRUE)
  spearman_mean <- mean(spearman_correlations, na.rm = TRUE)
  
  # Agregar interpretaciones
  frobenius_interpretation <- if (frobenius_mean < 0.1) "Small" else if (frobenius_mean < 0.5) "Medium" else "Large"
  tucker_interpretation <- if (tucker_mean > 0.90) "Highly Similar" else if (tucker_mean > 0.60) "Moderately Similar" else "Low Similarity"
  procrustes_interpretation <- if (procrustes_mean > 0.90) "Highly Similar" else if (procrustes_mean > 0.60) "Moderately Similar" else "Low Similarity"
  pearson_interpretation <- if (pearson_mean < 0.10) "Trivial" else if (pearson_mean < 0.30) "Small" else if (pearson_mean < 0.50) "Medium" else "Large"
  spearman_interpretation <- if (spearman_mean < 0.10) "Trivial" else if (spearman_mean < 0.30) "Small" else if (spearman_mean < 0.50) "Medium" else "Large"
  
  # Resultados globales
  global_results <- data.frame(
    Metric = c("Frobenius", "Tucker", "Procrustes", "Pearson", "Spearman"),
    Value = c(frobenius_mean, tucker_mean, procrustes_mean, pearson_mean, spearman_mean),
    Interpretation = c(frobenius_interpretation, tucker_interpretation, procrustes_interpretation, pearson_interpretation, spearman_interpretation)
  )
  
  # Crear data.frames para q, ΔZ, diferencias absolutas y relativas con interpretaciones
  q_results <- as.data.frame(as.table(apply(q_values, c(2, 3), mean, na.rm = TRUE)))
  colnames(q_results) <- c("Variable1", "Variable2", "Mean_q")
  q_results <- q_results[q_results$Variable1 != q_results$Variable2, ]
  
  delta_z_results <- as.data.frame(as.table(apply(delta_z_matrices, c(2, 3), mean, na.rm = TRUE)))
  colnames(delta_z_results) <- c("Variable1", "Variable2", "Mean_Delta_Z")
  delta_z_results <- delta_z_results[delta_z_results$Variable1 != delta_z_results$Variable2, ]
  
  absolute_diff_results <- as.data.frame(as.table(apply(absolute_diffs, c(2, 3), mean, na.rm = TRUE)))
  colnames(absolute_diff_results) <- c("Variable1", "Variable2", "Mean_Absolute_Diff")
  absolute_diff_results$Interpretation <- cut(
    absolute_diff_results$Mean_Absolute_Diff,
    breaks = c(-Inf, 0.10, 0.30, 0.50, Inf),
    labels = c("Trivial", "Small", "Medium", "Large"),
    right = FALSE
  )
  
  relative_diff_results <- as.data.frame(as.table(apply(relative_diffs, c(2, 3), mean, na.rm = TRUE)))
  colnames(relative_diff_results) <- c("Variable1", "Variable2", "Mean_Relative_Diff")
  
  # Retornar resultados
  results <- list(
    global_results = global_results,
    q_results = q_results,
    delta_z_results = delta_z_results,
    absolute_diff_results = absolute_diff_results,
    relative_diff_results = relative_diff_results
  )
  
  return(results)
}



result <- calculate_effect_sizes(datarem_11, datanonrem_11, nBoots = 10, seed = 123)

result$q_results
result$delta_z_results
result$absolute_diff_results
result$relative_diff_results
result$global_results$Metric
