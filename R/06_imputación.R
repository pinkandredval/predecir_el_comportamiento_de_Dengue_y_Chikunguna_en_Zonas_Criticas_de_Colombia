# ============================================================
# 06_imputacion.R
# Imputación coherente con el diagnóstico del 05_faltantes.R
#   - Dengue      -> mice + method = "rf"  (MAR, base grande)
#   - Chikungunya -> mice + method = "pmm" (MCAR, base chica)
#   - Respeta NA estructurales (no_aplica) con máscara
# Librerías: tidyverse, here, janitor, naniar, mice, ranger, scales
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(janitor)
  library(naniar)
  library(mice)
  library(ranger)     # requerido por mice cuando method = "rf"
  library(scales)
})

# --- Cargar datos preparados ---
dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

dir.create(here("data", "imputed"),
           recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0. FUNCIONES AUXILIARES (reutilizadas del 05)
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
  # Identificadores
  "consecutive", "cod_eve", "nombre_evento", "evento",
  # Fechas de proceso
  "fec_not", "fecha_epi", "ini_sin", "fec_con",
  # Administrativas / institucionales
  "nombre_upgd", "va_sispro", "fuente", "confirmados",
  "nom_est_f_caso", "estado_final_de_caso",
  # Geográficas rígidas (códigos)
  "cod_pais_o", "cod_dpto_o", "cod_mun_o",
  "cod_pais_r", "cod_dpto_r", "cod_mun_r",
  "cod_dpto_n", "cod_mun_n",
  # Temporales
  "ano", "anio", "semana",
  
  # >>> NUEVO: NA estructurales de subpoblación <<<
  # Militares (solo aplican a personal militar)
  "fm_fuerza", "fm_unidad", "fm_grado",
  
  # Grupos poblacionales específicos
  "gp_gestan", "sem_ges",           # solo gestantes
  "gp_discapa", "gp_desplaz",       # solo discapacitados/desplazados
  "gp_migrant", "gp_carcela",       # solo migrantes/carcelarios
  "gp_indigen", "gp_pobicfb",       # solo indígenas/afro
  "gp_mad_com", "gp_desmovi",       # solo madres comunitarias/desmovilizados
  "gp_psiquia", "gp_vic_vio",       # solo psiquiátricos/víctimas violencia
  "gp_otros", "gru_pob", "nom_grupo",  # grupo poblacional general
  
  # Otras específicas
  "fec_hos", "fec_def",             # solo hospitalizados/fallecidos
  "fec_aju", "ajuste",              # solo casos ajustados
  "fecha_nto",                      # fecha de nacimiento (puede ser específica)
  "con_fin",                        # condición final (solo fin de caso)
  "ocupacion"                       # si es muy específica por grupo
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
#    - Convertir caracteres a factor
#    - Convertir no_aplica -> NA temporal (con máscara)
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
# 3. SELECCIÓN AUTOMÁTICA DE MÉTODOS POR TIPO DE VARIABLE
# ============================================================

asignar_metodos <- function(df_imp, metodo_principal = c("rf", "pmm")) {
  
  metodo_principal <- match.arg(metodo_principal)
  
  metodos <- make.method(df_imp)
  
  es_binaria    <- sapply(df_imp, function(x) is.factor(x) && nlevels(x) == 2)
  es_categorica <- sapply(df_imp, function(x) is.factor(x) && nlevels(x) > 2)
  es_numerica   <- sapply(df_imp, is.numeric)
  es_entero     <- sapply(df_imp, is.integer)
  es_fecha      <- sapply(df_imp, function(x) inherits(x, "Date") || inherits(x, "POSIXct"))
  
  if (metodo_principal == "rf") {
    # Random Forest para todo lo que pueda (numéricas y categóricas)
    metodos[es_numerica]   <- "rf"
    metodos[es_entero]     <- "rf"
    metodos[es_categorica] <- "rf"
    metodos[es_binaria]    <- "rf"
  } else {
    # pmm para numéricas, logreg/polyreg para categóricas
    metodos[es_numerica]   <- "pmm"
    metodos[es_entero]     <- "pmm"
    metodos[es_binaria]    <- "logreg"
    metodos[es_categorica] <- "polyreg"
  }
  
  # Las fechas o columnas problemáticas -> dejar sin imputar
  if (any(es_fecha)) metodos[es_fecha] <- ""
  
  metodos
}

# ============================================================
# 4. FUNCIÓN MAESTRA DE IMPUTACIÓN
# ============================================================

imputar_evento <- function(df, nombre_evento, vars_imputar,
                           metodo_principal = "pmm",
                           m = 5, maxit = 5, seed = 123,
                           ntree = 10) {
  
  if (length(vars_imputar) == 0) {
    message("[skip] ", nombre_evento, ": nada que imputar")
    return(invisible(NULL))
  }
  
  cat("\n")
  cat("=========================================================\n")
  cat("IMPUTANDO:", nombre_evento, "\n")
  cat("=========================================================\n")
  cat("Variables a imputar:  ", length(vars_imputar), "\n")
  cat("Método principal:     ", metodo_principal, "\n")
  cat("m (imputaciones):     ", m, "\n")
  cat("maxit (iteraciones):  ", maxit, "\n")
  if (metodo_principal == "rf") {
    cat("ntree (árboles RF):   ", ntree, "\n")
  }
  cat("=========================================================\n")
  
  # --- Preparar datos ---
  prep <- preparar_para_imputar(df, vars_imputar)
  df_imp <- prep$df_imp
  mask_no_aplica <- prep$mask
  
  # --- Asignar métodos ---
  metodos <- asignar_metodos(df_imp, metodo_principal)
  
  cat("\nMétodos asignados:\n")
  print(table(metodos))
  
  # --- Parámetros extra para mice ---
  args_mice <- list(
    data      = df_imp,
    m         = m,
    maxit     = maxit,
    method    = metodos,
    seed      = seed,
    printFlag = TRUE
  )
  
  if (metodo_principal == "rf") {
    # Parámetros RF: ntree, etc.
    args_mice$ntree <- ntree
    # rf necesita al menos 2 variables predictoras
    if (ncol(df_imp) < 2) {
      warning("Se necesitan >=2 variables para RF. Cambiando a pmm.")
      metodos[metodos == "rf"] <- "pmm"
      metodos[metodos == ""]   <- ""  # nada
      args_mice$method <- metodos
      args_mice$ntree  <- NULL
    }
  }
  
  # --- Ejecutar mice con tiempo medido ---
  t0 <- Sys.time()
  imp <- do.call(mice, args_mice)
  t1 <- Sys.time()
  
  cat("\nTiempo de ejecución:",
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
  
  # --- Diagnóstico de convergencia ---
  cat("\n--- Convergencia (primeras 5 variables) ---\n")
  vars_diag <- head(vars_imputar, 5)
  tryCatch({
    print(imp$loggedEvents)
  }, error = function(e) NULL)
  
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

# --- DENGUE: rf, MAR, base grande ---
vars_imp_dengue <- selec_dengue %>% filter(imputar) %>% pull(variable)

cat("\n>>> Dengue: ", length(vars_imp_dengue),
    "variables a imputar con rf\n", sep = "")

res_dengue <- imputar_evento(
  df               = dengue,
  nombre_evento    = "Dengue",
  vars_imputar     = vars_imp_dengue,
  metodo_principal = "rf",
  m                = 3,
  maxit            = 2,
  ntree            = 5,
  seed             = 123
)

# --- CHIKUNGUNYA: pmm, MCAR, base chica ---
vars_imp_chik <- selec_chik %>% filter(imputar) %>% pull(variable)

cat("\n>>> Chikungunya: ", length(vars_imp_chik),
    "variables a imputar con pmm\n", sep = "")

res_chik <- imputar_evento(
  df               = chikungunya,
  nombre_evento    = "Chikungunya",
  vars_imputar     = vars_imp_chik,
  metodo_principal = "pmm",
  m                = 5,
  maxit            = 5,
  seed             = 123
)

# ============================================================
# 6. RECONSTRUIR DATA FRAMES COMPLETOS
# ============================================================

reconstruir_df <- function(df_original, df_imputado, vars_imputar) {
  
  df_final <- df_original
  
  # Reemplazar SOLO las variables imputadas
  for (v in vars_imputar) {
    df_final[[v]] <- df_imputado[[v]]
  }
  
  # Banderas de imputación (1 si fue faltante real, 0 si no)
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

# Data frames finales (con variables imputadas + banderas)
saveRDS(dengue_final,
        here("data", "imputed", "dengue_imputado.rds"))
saveRDS(chikungunya_final,
        here("data", "imputed", "chikungunya_imputado.rds"))

# Objetos mids (por si quieres pooling después)
saveRDS(res_dengue$imp,
        here("data", "imputed", "dengue_mice_mids.rds"))
saveRDS(res_chik$imp,
        here("data", "imputed", "chikungunya_mice_mids.rds"))

# Tablas de comparación NA antes/después
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
  cat("Método usado:                     ",
      ifelse(nombre == "DENGUE", "rf", "pmm"), "\n")
  cat("Tiempo total:                     ", res$tiempo_min, "min\n")
  
  cols_con_na <- names(cols_data)[sapply(cols_data, function(x) any(is.na(x)))]
  
  if (length(cols_con_na) > 0) {
    cat("Variables excluidas con NA:      ", length(cols_con_na), "\n")
    cat("  (son las que no entraron a imputar: identificadores, fechas de proceso, etc.)\n")
    print(head(cols_con_na, 20))
  } else {
    cat("Todas las variables candidatas fueron imputadas.\n")
  }
}

diagnostico_imputacion(dengue,      dengue_final,      "DENGUE",      res_dengue)
diagnostico_imputacion(chikungunya, chikungunya_final, "CHIKUNGUNYA", res_chik)

cat("\n")
cat("=========================================================\n")
cat("✅ IMPUTACIÓN COMPLETADA\n")
cat("=========================================================\n")
cat("Archivos en data/imputed/:\n")
cat("  - dengue_imputado.rds\n")
cat("  - chikungunya_imputado.rds\n")
cat("  - dengue_mice_mids.rds        (objeto mids para pooling)\n")
cat("  - chikungunya_mice_mids.rds\n")
cat("  - comparacion_na_dengue.csv\n")
cat("  - comparacion_na_chikungunya.csv\n")