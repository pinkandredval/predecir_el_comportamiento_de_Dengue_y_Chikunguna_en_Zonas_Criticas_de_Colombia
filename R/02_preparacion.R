# ============================================================
# 02_preparacion.R
# Limpieza básica de nombres y tipos
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
      # --- Identificador del evento ---
      evento    = evento,
      
      # --- Numéricas ---
      anio      = as.integer(anio),
      edad      = as.numeric(edad),
      semana    = as.integer(semana),
      
      # --- Fechas ---
      fec_not   = parsear_fecha(fec_not),
      ini_sin   = parsear_fecha(ini_sin),
      fec_hos   = parsear_fecha(fec_hos),
      fec_def   = parsear_fecha(fec_def),
      fecha_nto = parsear_fecha(fecha_nto),
      fec_con   = parsear_fecha(fec_con),
      
      # --- Semana epidemiológica (lunes de la semana ISO) ---
      fecha_epi = make_date(anio) + weeks(semana - 1)
    ) %>%
    select(-any_of(c("gru_pob", "cbmete")))
}

# --- Aplicar ---
dengue      <- preparar(dengue, "Dengue")
chikungunya <- preparar(chikungunya, "Chikungunya")

# --- Guardar ---
saveRDS(dengue,      here("data", "processed", "dengue_prep.rds"))
saveRDS(chikungunya, here("data", "processed", "chikungunya_prep.rds"))

cat("\n✅ Bases procesadas y guardadas.\n")