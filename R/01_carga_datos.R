# ============================================================
# 01_carga_datos.R
# Carga y unión de archivos crudos de SIVIGILA
# ============================================================

library(readxl)
library(tidyverse)
library(here)

# --- Listar archivos ---
archivos_dengue <- list.files(
  here("data", "raw", "dengue"),
  pattern = "\\.xlsx$",
  full.names = TRUE
)

archivos_chik <- list.files(
  here("data", "raw", "chikungunya"),
  pattern = "\\.xlsx$",
  full.names = TRUE
)

# --- Leer y unir dengue ---
dengue <- archivos_dengue |>
  set_names(str_extract(basename(archivos_dengue), "\\d{4}")) |>
  map(read_xlsx) |>
  bind_rows(.id = "anio")

# --- Leer y unir chikungunya ---
chikungunya <- archivos_chik |>
  set_names(str_extract(basename(archivos_chik), "\\d{4}")) |>
  map(read_xlsx) |>
  bind_rows(.id = "anio")

# --- Verificación ---
glimpse(dengue)
glimpse(chikungunya)
table(dengue$anio)
table(chikungunya$anio)

# --- Guardar como RDS para no repetir la lectura ---
dir.create(here("data", "processed"), showWarnings = FALSE)
saveRDS(dengue, here("data", "processed", "dengue_crudo.rds"))
saveRDS(chikungunya, here("data", "processed", "chikungunya_crudo.rds"))
