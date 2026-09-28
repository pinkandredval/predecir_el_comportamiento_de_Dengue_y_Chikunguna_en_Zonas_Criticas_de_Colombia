# ============================================================
# 05_faltantes.R
# Análisis de faltantes (estilo taller)
# ============================================================

library(tidyverse)
library(here)
library(janitor)
library(naniar)
library(mice)
library(scales)

# --- Cargar datos preparados ---
dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

# ============================================================
# 1. DESCRIPCIÓN DE FALTANTES
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("1. DESCRIPCIÓN DE FALTANTES\n")
cat("=========================================================\n")

# --- naniar: miss_var_summary() ---
cat("\n--- DENGUE: resumen de faltantes por variable ---\n")
miss_dengue <- miss_var_summary(dengue)
print(as.data.frame(miss_dengue), row.names = FALSE)

cat("\n--- CHIKUNGUNYA: resumen de faltantes por variable ---\n")
miss_chik <- miss_var_summary(chikungunya)
print(as.data.frame(miss_chik), row.names = FALSE)

# --- Resumen de casos y variables ---
cat("\n--- Resumen de casos (Dengue) ---\n")
print(miss_case_summary(dengue))

cat("\n--- Resumen de casos (Chikungunya) ---\n")
print(miss_case_summary(chikungunya))

# ============================================================
# 2. VARIABLES 100% VACÍAS Y CONSTANTES
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("2. VARIABLES 100% VACÍAS O CONSTANTES\n")
cat("=========================================================\n")

# --- Función para detectar problemas ---
detectar_problematicas <- function(df, nombre) {
  cat("\n---", nombre, "---\n")
  
  # 100% vacías
  vacias <- df %>%
    summarise(across(everything(), ~ all(is.na(.)))) %>%
    pivot_longer(everything()) %>%
    filter(value) %>%
    pull(name)
  
  cat("100% vacías (", length(vacias), "):\n", sep = "")
  print(vacias)
  
  # Constantes
  constantes <- df %>%
    summarise(across(everything(), ~ n_distinct(., na.rm = TRUE) == 1)) %>%
    pivot_longer(everything()) %>%
    filter(value) %>%
    pull(name)
  
  cat("\nConstantes (", length(constantes), "):\n", sep = "")
  print(constantes)
  
  # Retornar
  list(vacias = vacias, constantes = constantes)
}

prob_dengue <- detectar_problematicas(dengue, "DENGUE")
prob_chik   <- detectar_problematicas(chikungunya, "CHIKUNGUNYA")

# ============================================================
# 3. LIMPIEZA PREVIA (eliminar filas totalmente vacías)
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("3. LIMPIEZA PREVIA\n")
cat("=========================================================\n")

limpiar <- function(df, evento) {
  # Columnas que NO son identificadores (variables de análisis)
  vars_analisis <- setdiff(names(df),
                           c("consecutive", "cod_eve", "nombre_evento",
                             "evento", "va_sispro", "ano"))
  
  # 1. Eliminar filas totalmente vacías (en variables de análisis)
  n_antes <- nrow(df)
  df <- df %>%
    filter(rowSums(!is.na(select(., all_of(vars_analisis)))) > 0)
  n_despues <- nrow(df)
  
  cat("\n---", evento, "---\n")
  cat("Filas antes:", n_antes, "\n")
  cat("Filas después:", n_despues, "\n")
  cat("Filas eliminadas (totalmente vacías):", n_antes - n_despues, "\n")
  
  # 2. Eliminar columnas 100% vacías
  cols_vacias <- df %>%
    summarise(across(everything(), ~ all(is.na(.)))) %>%
    pivot_longer(everything()) %>%
    filter(value) %>%
    pull(name)
  
  if (length(cols_vacias) > 0) {
    cat("Columnas eliminadas (100% vacías):", length(cols_vacias), "\n")
    print(cols_vacias)
    df <- df %>% select(-all_of(cols_vacias))
  } else {
    cat("Sin columnas 100% vacías.\n")
  }
  
  df
}

dengue_limpio <- limpiar(dengue, "DENGUE")
chik_limpio   <- limpiar(chikungunya, "CHIKUNGUNYA")

# ============================================================
# 4. PATRÓN DE FALTANTES
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("4. PATRÓN DE FALTANTES\n")
cat("=========================================================\n")


# --- A) Tabla de patrones con TODAS las variables (texto) ---

cat("\n--- Patrones de faltantes (Dengue) ---\n")
patrones_dengue <- dengue_limpio %>%
  mutate(patron = apply(is.na(.), 1, paste, collapse = "")) %>%
  count(patron, sort = TRUE) %>%
  head(15)

# Convertir el patrón a lista de variables faltantes (más legible)
patrones_dengue_legible <- patrones_dengue %>%
  mutate(
    n_faltantes = str_count(patron, "T"),
    variables_faltantes = sapply(patron, function(p) {
      cols <- names(dengue_limpio)
      faltan <- cols[which(strsplit(p, "")[[1]] == "T")]
      if (length(faltan) == 0) return("(ninguna)")
      paste(faltan, collapse = ", ")
    })
  ) %>%
  select(-patron)

print(as.data.frame(patrones_dengue_legible), row.names = FALSE)

cat("\n--- Patrones de faltantes (Chikungunya) ---\n")
patrones_chik <- chik_limpio %>%
  mutate(patron = apply(is.na(.), 1, paste, collapse = "")) %>%
  count(patron, sort = TRUE) %>%
  head(15)

patrones_chik_legible <- patrones_chik %>%
  mutate(
    n_faltantes = str_count(patron, "T"),
    variables_faltantes = sapply(patron, function(p) {
      cols <- names(chik_limpio)
      faltan <- cols[which(strsplit(p, "")[[1]] == "T")]
      if (length(faltan) == 0) return("(ninguna)")
      paste(faltan, collapse = ", ")
    })
  ) %>%
  select(-patron)

print(as.data.frame(patrones_chik_legible), row.names = FALSE)

# --- B) gg_miss_upset (visualización moderna) ---

cat("\n--- Combinaciones de faltantes más frecuentes (Dengue) ---\n")
gg_miss_upset(dengue_limpio, nsets = 15)

cat("\n--- Combinaciones de faltantes más frecuentes (Chikungunya) ---\n")
gg_miss_upset(chik_limpio, nsets = 15)

# --- C) Mapa de faltantes (cuadrícula azul/rojo) ---

cat("\n--- Mapa de faltantes (Dengue) ---\n")
mapa_faltantes(dengue_limpio, "Dengue")

cat("\n--- Mapa de faltantes (Chikungunya) ---\n")
mapa_faltantes(chik_limpio, "Chikungunya")

# ============================================================
# 5. MAPA DE FALTANTES (cuadrícula azul/rojo)
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("5. MAPA DE FALTANTES\n")
cat("=========================================================\n")

mapa_faltantes <- function(df, nombre_evento, max_filas = 5000) {
  if (nrow(df) > max_filas) {
    df_muestra <- df %>% slice_sample(n = max_filas)
    subtitulo <- paste0("(muestra de ", comma(max_filas), " filas)")
  } else {
    df_muestra <- df
    subtitulo <- paste0("(", comma(nrow(df)), " filas)")
  }
  
  df_long <- df_muestra %>%
    mutate(fila = row_number()) %>%
    mutate(across(-fila, is.na)) %>%
    pivot_longer(-fila, names_to = "variable", values_to = "faltante")
  
  orden_vars <- df_long %>%
    group_by(variable) %>%
    summarise(pct_falt = mean(faltante) * 100, .groups = "drop") %>%
    arrange(desc(pct_falt)) %>%
    pull(variable)
  
  df_long <- df_long %>%
    mutate(variable = factor(variable, levels = orden_vars))
  
  ggplot(df_long, aes(x = variable, y = fila, fill = faltante)) +
    geom_raster() +
    scale_fill_manual(
      values = c("FALSE" = "#2166AC", "TRUE" = "#B2182B"),
      labels = c("FALSE" = "Presente", "TRUE" = "Faltante"),
      name   = ""
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title    = paste0("Mapa de faltantes - ", nombre_evento),
      subtitle = subtitulo,
      x        = NULL,
      y        = "Registros (muestra)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_line(color = "white", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top"
    )
}

mapa_faltantes(dengue_limpio, "Dengue")
mapa_faltantes(chik_limpio,   "Chikungunya")

# ============================================================
# 6. FALTANTES AGRUPADOS POR TIPO DE VARIABLE
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("6. FALTANTES POR GRUPO DE VARIABLES\n")
cat("=========================================================\n")

vars_socio    <- c("edad", "sexo", "tip_ss", "per_etn", "area",
                   "estrato", "ocupacion", "cod_ase")
vars_clinicas <- c("tip_cas", "pac_hos", "con_fin", "ajuste")
vars_fechas   <- c("fec_not", "ini_sin", "fec_hos", "fec_def",
                   "fecha_nto", "fec_con")

resumen_grupo <- function(df, vars, nombre_grupo, evento) {
  vars_presentes <- intersect(vars, names(df))
  if (length(vars_presentes) == 0) return(NULL)
  
  df %>%
    select(all_of(vars_presentes)) %>%
    summarise(across(everything(), ~ round(mean(is.na(.)) * 100, 2))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct") %>%
    mutate(grupo = nombre_grupo, evento = evento)
}

falt_grupos <- bind_rows(
  resumen_grupo(dengue_limpio, vars_socio,    "Sociodemográficas", "Dengue"),
  resumen_grupo(dengue_limpio, vars_clinicas, "Clínicas",          "Dengue"),
  resumen_grupo(dengue_limpio, vars_fechas,   "Fechas",            "Dengue"),
  resumen_grupo(chik_limpio,   vars_socio,    "Sociodemográficas", "Chikungunya"),
  resumen_grupo(chik_limpio,   vars_clinicas, "Clínicas",          "Chikungunya"),
  resumen_grupo(chik_limpio,   vars_fechas,   "Fechas",            "Chikungunya")
)

cat("\n--- Faltantes por grupo (Dengue) ---\n")
print(as.data.frame(falt_grupos %>% filter(evento == "Dengue")))

cat("\n--- Faltantes por grupo (Chikungunya) ---\n")
print(as.data.frame(falt_grupos %>% filter(evento == "Chikungunya")))

cat("\n--- Promedio de faltantes por grupo ---\n")
print(
  falt_grupos %>%
    group_by(evento, grupo) %>%
    summarise(promedio_pct = round(mean(pct), 2), .groups = "drop") %>%
    as.data.frame()
)

# ============================================================
# 7. COMPARACIÓN ENTRE EVENTOS
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("7. COMPARACIÓN DE FALTANTES ENTRE EVENTOS\n")
cat("=========================================================\n")

columnas_comunes <- intersect(names(dengue_limpio), names(chik_limpio))

comparacion <- bind_rows(
  miss_var_summary(dengue_limpio %>% select(all_of(columnas_comunes))) %>%
    mutate(evento = "Dengue"),
  miss_var_summary(chik_limpio %>% select(all_of(columnas_comunes))) %>%
    mutate(evento = "Chikungunya")
) %>%
  select(variable, evento, pct_miss) %>%
  pivot_wider(names_from = evento, values_from = pct_miss, values_fill = 0) %>%
  mutate(diferencia = abs(Dengue - Chikungunya)) %>%
  arrange(desc(diferencia))

cat("\n--- Top 20 variables con más diferencia ---\n")
print(head(as.data.frame(comparacion), 20), row.names = FALSE)

# ============================================================
# 8. DIAGNÓSTICO FINAL
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("8. DIAGNÓSTICO FINAL\n")
cat("=========================================================\n")

diagnostico <- function(df_original, df_limpio, nombre) {
  n_orig   <- nrow(df_original)
  n_limpio <- nrow(df_limpio)
  
  miss_total <- sum(is.na(df_limpio))
  celdas_totales <- nrow(df_limpio) * ncol(df_limpio)
  
  cat("\n---", nombre, "---\n")
  cat("Filas originales:", comma(n_orig), "\n")
  cat("Filas después de limpieza:", comma(n_limpio), "\n")
  cat("Filas eliminadas:", comma(n_orig - n_limpio), "\n")
  cat("Columnas:", ncol(df_limpio), "\n")
  cat("Celdas totales:", comma(celdas_totales), "\n")
  cat("Celdas faltantes:", comma(miss_total), "\n")
  cat("Porcentaje faltante global:", round(miss_total / celdas_totales * 100, 2), "%\n")
}

diagnostico(dengue, dengue_limpio, "DENGUE")
diagnostico(chikungunya, chik_limpio, "CHIKUNGUNYA")

cat("\n✅ Análisis de faltantes completado.\n")
cat("Los objetos 'dengue_limpio' y 'chik_limpio' están listos para imputación.\n")
