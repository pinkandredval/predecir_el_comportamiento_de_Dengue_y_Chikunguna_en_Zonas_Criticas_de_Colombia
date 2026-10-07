# ============================================================
# 03_auditoria.R
# Auditoría después de la preparación
# ============================================================

library(tidyverse)
library(here)
library(janitor)

dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

auditar <- function(df, nombre) {
  cat("\n==================================================\n")
  cat("AUDITORÍA:", nombre, "\n")
  cat("==================================================\n")
  
  # --- Excluir columnas auxiliares na_estructural_* ---
  df_analisis <- df %>% select(-starts_with("na_estructural"))
  
  cat("\nDimensiones:", nrow(df), "filas x", ncol(df), "columnas\n")
  cat("Columnas de análisis (sin auxiliares):", ncol(df_analisis), "\n")
  
  # --- Tipos ---
  cat("\n--- TIPOS DE COLUMNAS ---\n")
  df_analisis %>%
    summarise(across(everything(), ~ class(.)[1])) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "tipo") %>%
    count(tipo) %>%
    print()
  
  # --- 100% vacías ---
  cat("\n--- COLUMNAS 100% VACÍAS ---\n")
  vacias <- df_analisis %>%
    summarise(across(everything(), ~ all(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "vacia") %>%
    filter(vacia) %>% pull(variable)
  if (length(vacias) == 0) cat("Ninguna. ✅\n") else print(vacias)
  
  # --- Constantes ---
  cat("\n--- COLUMNAS CONSTANTES ---\n")
  constantes <- df_analisis %>%
    summarise(across(everything(), ~ n_distinct(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_unicos") %>%
    filter(n_unicos == 1)
  if (nrow(constantes) == 0) cat("Ninguna. ✅\n") else print(constantes)
  
  # --- Cobertura temporal ---
  cat("\n--- COBERTURA POR AÑO ---\n")
  print(table(df$anio, useNA = "ifany"))
  
  cat("\n--- RANGO fecha_epi ---\n")
  print(range(df$fecha_epi, na.rm = TRUE))
  
  # --- Fuera de rango ---
  cat("\n--- VALORES FUERA DE RANGO ---\n")
  cat("Edad mín/máx:", min(df$edad, na.rm = TRUE), "/",
      max(df$edad, na.rm = TRUE), "\n")
  cat("Edad < 0:", sum(df$edad < 0, na.rm = TRUE), "\n")
  cat("Edad > 110:", sum(df$edad > 110, na.rm = TRUE), "\n")
  cat("Semana < 1:", sum(df$semana < 1, na.rm = TRUE), "\n")
  cat("Semana > 53:", sum(df$semana > 53, na.rm = TRUE), "\n")
  
  # --- Categóricas clave ---
  cat("\n--- CATEGÓRICAS CLAVE ---\n")
  vars_clave <- c("evento", "sexo", "tip_cas", "pac_hos", "con_fin",
                  "ajuste", "tip_ss", "area", "per_etn")
  for (v in vars_clave) {
    if (v %in% names(df)) {
      cat("\n", v, ":\n")
      print(table(df[[v]], useNA = "ifany"))
    }
  }
}

auditar(dengue, "DENGUE")
auditar(chikungunya, "CHIKUNGUNYA")