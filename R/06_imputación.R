# ============================================================
# 06_imputacion.R (Versión final — todo con PMM)
#   - Dengue      -> mice + method = "pmm"
#   - Chikungunya -> mice + method = "pmm"
#   - Numéricas y categóricas (binarias y nominales) con pmm
#   - Respeta NA estructurales (no_aplica) con máscara
# Librerías: tidyverse, here, janitor, naniar, mice, scales
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(janitor)
  library(naniar)
  library(mice)
  library(scales)
})

# --- Cargar datos preparados ---
dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

dir.create(here("data", "imputed"),
           recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0. FUNCIONES AUXILIARES
# ============================================================

clasificar_celdas <- function(df) {
  df_analisis <- df %>% select(-starts_with("na_estructural_"))
  vars_analisis <- names(df_analisis)
  
  estados <- df_analisis %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "faltante", "presente")))
  
  for (v in vars_analisis) {
    col_est <- paste0("na_estructural_", v)
    if (col_est %in% names(df)) {
      idx <- which(!is.na(df[[col_est]]) & df[[col_est]])
      if (length(idx) > 0) estados[[v]][idx] <- "no_aplica"
    }
  }
  estados
}

# ============================================================
# 1. SELECCIÓN DE VARIABLES A IMPUTAR
# ============================================================

vars_excluir <- c(
  "consecutive", "cod_eve", "nombre_evento", "evento",
  "fec_not", "fecha_epi", "ini_sin", "fec_con",
  "nombre_upgd", "va_sispro", "fuente", "confirmados",
  "nom_est_f_caso", "estado_final_de_caso",
  "cod_pais_o", "cod_dpto_o", "cod_mun_o",
  "cod_pais_r", "cod_dpto_r", "cod_mun_r",
  "cod_dpto_n", "cod_mun_n",
  "ano", "anio", "semana"
)

seleccionar_vars_imputar <- function(df, vars_excluir) {
  estados <- clasificar_celdas(df)
  vars_candidatas <- setdiff(names(estados), vars_excluir)
  
  faltantes_reales <- sapply(vars_candidatas, function(v) {
    sum(estados[[v]] == "faltante", na.rm = TRUE)
  })
  
  tibble(
    variable        = vars_candidatas,
    n_faltante_real = as.integer(faltantes_reales),
    imputar         = faltantes_reales > 0
  ) %>%
    arrange(desc(n_faltante_real))
}

cat("\n===== VARIABLES A IMPUTAR — DENGUE =====\n")
selec_dengue <- seleccionar_vars_imputar(dengue, vars_excluir)
print(as.data.frame(filter(selec_dengue, imputar)), row.names = FALSE)

cat("\n===== VARIABLES A IMPUTAR — CHIKUNGUNYA =====\n")
selec_chik <- seleccionar_vars_imputar(chikungunya, vars_excluir)
print(as.data.frame(filter(selec_chik, imputar)), row.names = FALSE)

# ============================================================
# 2. PREPARAR DATA FRAME PARA IMPUTACIÓN
# ============================================================

preparar_para_imputar <- function(df, vars_imputar) {
  
  estados <- clasificar_celdas(df)
  
  df_imp <- df %>%
    select(all_of(vars_imputar)) %>%
    mutate(across(where(is.character), as.factor))
  
  mask_no_aplica <- list()
  for (v in vars_imputar) {
    idx_no_aplica <- which(estados[[v]] == "no_aplica")
    mask_no_aplica[[v]] <- idx_no_aplica
    if (length(idx_no_aplica) > 0) {
      df_imp[[v]][idx_no_aplica] <- NA
    }
  }
  
  list(df_imp = df_imp, mask = mask_no_aplica)
}

# ============================================================
# 3. ASIGNAR MÉTODOS — TODO CON PMM
# ============================================================

asignar_metodos <- function(df_imp) {
  
  metodos <- make.method(df_imp)
  
  es_binaria    <- sapply(df_imp, function(x) is.factor(x) && nlevels(x) == 2)
  es_categorica <- sapply(df_imp, function(x) is.factor(x) && nlevels(x) > 2)
  es_numerica   <- sapply(df_imp, is.numeric)
  es_entero     <- sapply(df_imp, is.integer)
  es_fecha      <- sapply(df_imp, function(x) inherits(x, "Date") || inherits(x, "POSIXct"))
  
  # TODO con pmm (numéricas y categóricas)
  metodos[es_numerica]   <- "pmm"
  metodos[es_entero]     <- "pmm"
  metodos[es_binaria]    <- "pmm"
  metodos[es_categorica] <- "pmm"
  
  # Fechas no se imputan (mice no maneja bien Date con pmm)
  if (any(es_fecha)) metodos[es_fecha] <- ""
  
  metodos
}

# ============================================================
# 4. FUNCIÓN MAESTRA DE IMPUTACIÓN (todo pmm)
# ============================================================

imputar_evento <- function(df, nombre_evento, vars_imputar,
                           m = 5, maxit = 5, seed = 123) {
  
  if (length(vars_imputar) == 0) {
    message("[skip] ", nombre_evento, ": nada que imputar")
    return(invisible(NULL))
  }
  
  cat("\n")
  cat("=========================================================\n")
  cat("IMPUTANDO:", nombre_evento, " (todo con pmm)\n")
  cat("=========================================================\n")
  cat("Variables a imputar:  ", length(vars_imputar), "\n")
  cat("Método:               pmm\n")
  cat("m (imputaciones):     ", m, "\n")
  cat("maxit (iteraciones):  ", maxit, "\n")
  cat("=========================================================\n")
  
  # --- Preparar datos ---
  prep <- preparar_para_imputar(df, vars_imputar)
  df_imp <- prep$df_imp
  mask_no_aplica <- prep$mask
  
  # --- Asignar métodos (todo pmm) ---
  metodos <- asignar_metodos(df_imp)
  
  cat("\nMétodos asignados:\n")
  print(table(metodos))
  
  # --- Ejecutar mice con tiempo medido ---
  t0 <- Sys.time()
  imp <- mice(
    df_imp,
    m         = m,
    maxit     = maxit,
    method    = metodos,
    seed      = seed,
    printFlag = TRUE
  )
  t1 <- Sys.time()
  
  cat("\nTiempo de imputación:",
      round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")
  
  # --- Extraer la primera imputación completada ---
  df_imputado <- complete(imp, 1)
  
  # --- Restaurar no_aplica (volver a NA en esas posiciones) ---
  for (v in vars_imputar) {
    idx <- mask_no_aplica[[v]]
    if (length(idx) > 0) {
      df_imputado[[v]][idx] <- NA
    }
  }
  
  # --- Comparación NA antes/después ---
  comparacion <- tibble(
    variable     = vars_imputar,
    n_na_antes   = sapply(df_imp,      function(x) sum(is.na(x))),
    n_na_despues = sapply(df_imputado, function(x) sum(is.na(x)))
  ) %>%
    mutate(
      n_imputado   = n_na_antes - n_na_despues,
      pct_imputado = round(n_imputado / pmax(n_na_antes, 1) * 100, 1)
    )
  
  cat("\n--- Comparación NA antes/después ---\n")
  print(as.data.frame(comparacion), row.names = FALSE)
  
  list(
    imp            = imp,
    df_imputado    = df_imputado,
    comparacion    = comparacion,
    mask_no_aplica = mask_no_aplica,
    metodos        = metodos,
    tiempo_min     = round(as.numeric(difftime(t1, t0, units = "mins")), 2)
  )
}

# ============================================================
# 5. EJECUTAR IMPUTACIÓN POR EVENTO
# ============================================================

# --- DENGUE: pmm, m=3, maxit=2 (rápido y suficiente) ---
vars_imp_dengue <- selec_dengue %>% filter(imputar) %>% pull(variable)

cat("\n>>> Dengue: ", length(vars_imp_dengue),
    "variables a imputar con pmm\n", sep = "")

res_dengue <- imputar_evento(
  df            = dengue,
  nombre_evento = "Dengue",
  vars_imputar  = vars_imp_dengue,
  m             = 3,
  maxit         = 2,
  seed          = 123
)

# --- CHIKUNGUNYA: pmm, m=5, maxit=5 (base chica, no hay prisa) ---
vars_imp_chik <- selec_chik %>% filter(imputar) %>% pull(variable)

cat("\n>>> Chikungunya: ", length(vars_imp_chik),
    "variables a imputar con pmm\n", sep = "")

res_chik <- imputar_evento(
  df            = chikungunya,
  nombre_evento = "Chikungunya",
  vars_imputar  = vars_imp_chik,
  m             = 5,
  maxit         = 5,
  seed          = 123
)

# ============================================================
# 6. RECONSTRUIR DATA FRAMES COMPLETOS
# ============================================================

reconstruir_df <- function(df_original, df_imputado, vars_imputar) {
  
  df_final <- df_original
  
  for (v in vars_imputar) {
    df_final[[v]] <- df_imputado[[v]]
  }
  
  estados_orig <- clasificar_celdas(df_original)
  
  for (v in vars_imputar) {
    col_flag <- paste0("imputado_", v)
    df_final[[col_flag]] <- as.integer(estados_orig[[v]] == "faltante")
  }
  
  df_final
}

dengue_final      <- reconstruir_df(dengue,      res_dengue$df_imputado,
                                    vars_imp_dengue)
chikungunya_final <- reconstruir_df(chikungunya, res_chik$df_imputado,
                                    vars_imp_chik)

# ============================================================
# 7. GUARDAR RESULTADOS
# ============================================================

saveRDS(dengue_final,
        here("data", "imputed", "dengue_imputado.rds"))
saveRDS(chikungunya_final,
        here("data", "imputed", "chikungunya_imputado.rds"))

saveRDS(res_dengue$imp,
        here("data", "imputed", "dengue_mice_mids.rds"))
saveRDS(res_chik$imp,
        here("data", "imputed", "chikungunya_mice_mids.rds"))

write.csv(res_dengue$comparacion,
          here("data", "imputed", "comparacion_na_dengue.csv"),
          row.names = FALSE)
write.csv(res_chik$comparacion,
          here("data", "imputed", "comparacion_na_chikungunya.csv"),
          row.names = FALSE)

# ============================================================
# 8. DIAGNÓSTICO FINAL
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("8. DIAGNÓSTICO FINAL DE IMPUTACIÓN\n")
cat("=========================================================\n")

diagnostico_imputacion <- function(df_orig, df_final, nombre, res) {
  cat("\n---", nombre, "---\n")
  
  estados_orig <- clasificar_celdas(df_orig)
  na_orig  <- sum(sapply(estados_orig, function(x) sum(x == "faltante")))
  
  cols_data <- df_final %>%
    select(-starts_with("na_estructural_"), -starts_with("imputado_"))
  na_final <- sum(is.na(cols_data))
  
  cat("Celdas faltante real (original): ", comma(na_orig), "\n")
  cat("Celdas NA en df final:            ", comma(na_final), "\n")
  cat("Método usado:                     pmm\n")
  cat("Tiempo total:                     ", res$tiempo_min, "min\n")
  
  cols_con_na <- names(cols_data)[sapply(cols_data, function(x) any(is.na(x)))]
  
  if (length(cols_con_na) > 0) {
    cat("Variables excluidas con NA:      ", length(cols_con_na), "\n")
    cat("  (identificadores, fechas de proceso, administrativas)\n")
    print(head(cols_con_na, 20))
  } else {
    cat("Todas las variables candidatas fueron imputadas.\n")
  }
}

diagnostico_imputacion(dengue,      dengue_final,      "DENGUE",      res_dengue)
diagnostico_imputacion(chikungunya, chikungunya_final, "CHIKUNGUNYA", res_chik)

cat("\n")
cat("=========================================================\n")
cat("✅ IMPUTACIÓN COMPLETADA (todo con pmm)\n")
cat("=========================================================\n")
cat("Archivos en data/imputed/:\n")
cat("  - dengue_imputado.rds\n")
cat("  - chikungunya_imputado.rds\n")
cat("  - dengue_mice_mids.rds\n")
cat("  - chikungunya_mice_mids.rds\n")
cat("  - comparacion_na_dengue.csv\n")
cat("  - comparacion_na_chikungunya.csv\n")