# ============================================================
# 02_preparacion.R
# Limpieza básica de nombres y tipos + marcado de NA estructurales
# ============================================================

library(tidyverse)
library(here)
library(janitor)
library(lubridate)

# --- Cargar crudos ---
dengue      <- readRDS(here("data", "processed", "dengue_crudo.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_crudo.rds"))

# --- Función para parsear fechas del SIVIGILA ---
parsear_fecha <- function(x) {
  x[x == "" | x == "NA"] <- NA
  solo_fecha <- substr(x, 1, 10)
  as.Date(solo_fecha, format = "%d/%m/%Y")
}

# --- Función de preparación ---
preparar <- function(df, evento) {
  df %>%
    clean_names() %>%
    mutate(
      evento    = evento,
      anio      = as.integer(anio),
      edad      = as.numeric(edad),
      semana    = as.integer(semana),
      
      fec_not   = parsear_fecha(fec_not),
      ini_sin   = parsear_fecha(ini_sin),
      fec_hos   = parsear_fecha(fec_hos),
      fec_def   = parsear_fecha(fec_def),
      fecha_nto = parsear_fecha(fecha_nto),
      fec_con   = parsear_fecha(fec_con),
      
      fecha_epi = make_date(anio) + weeks(semana - 1)
    ) %>%
    select(-any_of(c("gru_pob", "cbmete")))
}

# --- Aplicar preparación básica ---
dengue      <- preparar(dengue, "Dengue")
chikungunya <- preparar(chikungunya, "Chikungunya")

# ============================================================
# MARCADO DE FALTANTES ESTRUCTURALES
# Basado en el diccionario oficial del INS
# ============================================================
# Un NA es estructural cuando la variable no aplica para ese registro
# según las condiciones de la ficha SIVIGILA.
#
# Reglas:
#   - sem_ges:   solo aplica si gp_gestan = 1
#   - fec_hos:   solo aplica si pac_hos = 1
#   - fec_def:   solo aplica si con_fin = 2
#   - cbmte:     solo aplica si con_fin = 2
#   - fec_aju:   solo aplica si ajuste != 0
#   - fm_fuerza: aplica a militares (las 3 NA = no militar = estructural)
#   - fm_unidad: aplica a militares (las 3 NA = no militar = estructural)
#   - fm_grado:  aplica a militares (las 3 NA = no militar = estructural)
#   - nom_grupo: solo aplica si hay grupo poblacional

marcar_estructurales <- function(df) {
  
  # Detectar qué columnas existen en la base
  tiene <- function(v) v %in% names(df)
  
  tiene_fm_todas <- tiene("fm_fuerza") & tiene("fm_unidad") & tiene("fm_grado")
  
  df %>%
    mutate(
      # --- sem_ges: solo si es gestante ---
      na_estructural_sem_ges = if (tiene("sem_ges")) {
        is.na(sem_ges) & gp_gestan == "2"
      } else FALSE,
      
      # --- fec_hos: solo si fue hospitalizado ---
      na_estructural_fec_hos = is.na(fec_hos) & pac_hos == "2",
      
      # --- fec_def: solo si murió ---
      na_estructural_fec_def = if (tiene("fec_def")) {
        is.na(fec_def) & con_fin == "1"
      } else FALSE,
      
      # --- cbmte: solo si murió ---
      na_estructural_cbmte = if (tiene("cbmte")) {
        is.na(cbmte) & con_fin == "1"
      } else FALSE,
      
      # --- fec_aju: solo si hubo ajuste ---
      na_estructural_fec_aju = if (tiene("fec_aju")) {
        is.na(fec_aju) & ajuste == "0"
      } else FALSE,
      
      # --- fm_*: aplica solo a militares ---
      # Es estructural SOLO si las TRES están NA (no es militar).
      # Si alguna tiene dato, es militar y los NA son faltantes reales.
      na_estructural_fm_fuerza = if (tiene_fm_todas) {
        is.na(fm_fuerza) & is.na(fm_unidad) & is.na(fm_grado)
      } else FALSE,
      
      na_estructural_fm_unidad = if (tiene_fm_todas) {
        is.na(fm_unidad) & is.na(fm_fuerza) & is.na(fm_grado)
      } else FALSE,
      
      na_estructural_fm_grado = if (tiene_fm_todas) {
        is.na(fm_grado) & is.na(fm_fuerza) & is.na(fm_unidad)
      } else FALSE,
      
      # --- nom_grupo: solo si hay grupo poblacional ---
      na_estructural_nom_grupo = if (tiene("nom_grupo")) {
        is.na(nom_grupo)
      } else FALSE
    )
}

dengue      <- marcar_estructurales(dengue)
chikungunya <- marcar_estructurales(chikungunya)

# --- Verificar ---
cat("\n--- Columnas na_estructural en dengue ---\n")
print(names(dengue)[grepl("^na_estructural", names(dengue))])
cat("\nResumen de NA estructurales (dengue):\n")
dengue %>%
  summarise(across(starts_with("na_estructural"), sum)) %>%
  print()

cat("\n--- Columnas na_estructural en chikungunya ---\n")
print(names(chikungunya)[grepl("^na_estructural", names(chikungunya))])
cat("\nResumen de NA estructurales (chikungunya):\n")
chikungunya %>%
  summarise(across(starts_with("na_estructural"), sum)) %>%
  print()

# --- Verificar que los militares se conservan ---
cat("\n--- Verificación de militares en DENGUE ---\n")
cat("fm_fuerza no-NA:", sum(!is.na(dengue$fm_fuerza)), "\n")
cat("fm_unidad no-NA:", sum(!is.na(dengue$fm_unidad)), "\n")
cat("fm_grado  no-NA:", sum(!is.na(dengue$fm_grado)), "\n")
cat("Faltantes reales fm_unidad (militares sin unidad):",
    sum(dengue$na_estructural_fm_unidad == FALSE & is.na(dengue$fm_unidad)), "\n")
cat("Faltantes reales fm_grado (militares sin grado):",
    sum(dengue$na_estructural_fm_grado == FALSE & is.na(dengue$fm_grado)), "\n")

cat("\n--- Verificación de militares en CHIKUNGUNYA ---\n")
cat("fm_fuerza no-NA:", sum(!is.na(chikungunya$fm_fuerza)), "\n")
cat("fm_unidad no-NA:", sum(!is.na(chikungunya$fm_unidad)), "\n")
cat("fm_grado  no-NA:", sum(!is.na(chikungunya$fm_grado)), "\n")
cat("Faltantes reales fm_unidad (militares sin unidad):",
    sum(chikungunya$na_estructural_fm_unidad == FALSE & is.na(chikungunya$fm_unidad)), "\n")
cat("Faltantes reales fm_grado (militares sin grado):",
    sum(chikungunya$na_estructural_fm_grado == FALSE & is.na(chikungunya$fm_grado)), "\n")

# --- Guardar ---
saveRDS(dengue,      here("data", "processed", "dengue_prep.rds"))
saveRDS(chikungunya, here("data", "processed", "chikungunya_prep.rds"))

cat("\n✅ Bases procesadas y guardadas con marcado de NA estructurales.\n")