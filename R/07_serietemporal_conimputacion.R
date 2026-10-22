# ============================================================
# 07_serie_temporal.R
# Construye series temporales semanales y mensuales
# desde datos imputados (nacional + por departamento + por grupo etario)
#
# Corrección Tabla 2 (#1):
#   - 33 códigos DIVIPOLA (32 departamentos + Bogotá D.C.) x 314 semanas.
#   - expand_grid() + left_join() + replace_na(0) para hacer explícitos
#     los ceros. No se usa unique() para construir la rejilla.
#   - El código 01 (EXTERIOR, no DIVIPOLA) se conserva como serie aparte
#     para que la serie NACIONAL siga cuadrando con la suma de sus partes.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(scales)
})

dengue_imp      <- readRDS(here("data", "imputed", "dengue_imputado.rds"))
chikungunya_imp <- readRDS(here("data", "imputed", "chikungunya_imputado.rds"))

dir.create(here("data", "series"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("graficas", "series"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0. Rejilla temporal: 33 códigos DIVIPOLA x 314 semanas
# ============================================================

# 32 departamentos + Bogotá D.C. (11). El código 01 no es DIVIPOLA:
# en esta base corresponde a casos importados (EXTERIOR).
codigos_divipola <- c(
  "05", "08", "11", "13", "15", "17", "18", "19", "20", "23",
  "25", "27", "41", "44", "47", "50", "52", "54", "63", "66",
  "68", "70", "73", "76", "81", "85", "86", "88", "91", "94",
  "95", "97", "99"
)
stopifnot(length(codigos_divipola) == 33L, !anyDuplicated(codigos_divipola))

# Semanas ISO/MMWR presentes en el calendario del proyecto:
# 2020 y 2025 tienen 53 semanas; 2021-2024 tienen 52. Total = 314.
semanas_por_anio <- tibble(
  anio      = 2020:2025,
  n_semanas = c(53L, 52L, 52L, 52L, 52L, 53L)
)

calendario_semanal <- expand_grid(anio = 2020:2025, semana = 1:53) %>%
  left_join(semanas_por_anio, by = "anio") %>%
  filter(semana <= n_semanas) %>%
  mutate(
    fecha_sem = as.Date(paste0(anio, "-01-01")) + (semana - 1) * 7
  ) %>%
  select(anio, semana, fecha_sem) %>%
  arrange(anio, semana)

stopifnot(nrow(calendario_semanal) == 314L)

cat("\n===== REJILLA TEMPORAL =====\n")
cat("Semanas del calendario:", nrow(calendario_semanal), "\n")
cat("Códigos DIVIPOLA:", length(codigos_divipola), "\n")
cat("Filas esperadas por código:", nrow(calendario_semanal), "\n")

# ============================================================
# 1. FUNCIÓN: construir serie semanal
# ============================================================

construir_serie <- function(df, nombre_evento) {

  col_dpto <- intersect(
    c("cod_dpto_o", "cod_dpto_r", "cod_dpto_n",
      "departamento_ocurrencia", "departamento_residencia",
      "departamento_notificacion"),
    names(df)
  )[1]

  if (is.na(col_dpto)) {
    warning("No se encontró columna de departamento en ", nombre_evento)
    return(NULL)
  }

  df_serie <- df %>%
    mutate(
      anio_num   = as.integer(anio),
      semana_num = as.integer(semana),
      dpto       = as.character(.data[[col_dpto]]),
      grupo_edad = ifelse(!is.na(edad) & edad < 18, "Joven (0-17)",
                          ifelse(!is.na(edad), "Adulto (18+)", NA_character_))
    ) %>%
    filter(!is.na(anio_num), !is.na(semana_num),
           anio_num >= 2020, anio_num <= 2025) %>%
    mutate(
      fecha_sem = as.Date(paste0(anio_num, "-01-01")) +
        (semana_num - 1) * 7
    )

  # --- Serie por departamento: 33 códigos DIVIPOLA x 314 semanas ---
  conteos_dpto <- df_serie %>%
    filter(dpto %in% codigos_divipola) %>%
    count(dpto, anio = anio_num, semana = semana_num, name = "casos")

  serie_dpto <- expand_grid(dpto = codigos_divipola, calendario_semanal) %>%
    left_join(conteos_dpto, by = c("dpto", "anio", "semana")) %>%
    mutate(casos = tidyr::replace_na(casos, 0L)) %>%
    select(dpto, anio, semana, fecha_sem, casos)

  # --- Serie EXTERIOR (código 01): se conserva fuera del panel DIVIPOLA ---
  conteos_exterior <- df_serie %>%
    filter(!dpto %in% codigos_divipola) %>%
    count(anio = anio_num, semana = semana_num, name = "casos")

  serie_exterior <- calendario_semanal %>%
    left_join(conteos_exterior, by = c("anio", "semana")) %>%
    mutate(
      dpto  = "EXTERIOR",
      casos = tidyr::replace_na(casos, 0L)
    ) %>%
    select(dpto, anio, semana, fecha_sem, casos)

  # --- Serie nacional (incluye todos los códigos, también EXTERIOR) ---
  conteos_nacional <- df_serie %>%
    count(anio = anio_num, semana = semana_num, name = "casos")

  serie_nacional <- calendario_semanal %>%
    left_join(conteos_nacional, by = c("anio", "semana")) %>%
    mutate(
      dpto  = "NACIONAL",
      casos = tidyr::replace_na(casos, 0L)
    ) %>%
    select(dpto, anio, semana, fecha_sem, casos)

  # --- Serie por grupo etario (nacional), también con ceros explícitos ---
  conteos_edad <- df_serie %>%
    filter(!is.na(grupo_edad)) %>%
    count(grupo_edad, anio = anio_num, semana = semana_num, name = "casos")

  serie_edad <- expand_grid(
    grupo_edad = c("Joven (0-17)", "Adulto (18+)"),
    calendario_semanal
  ) %>%
    left_join(conteos_edad, by = c("grupo_edad", "anio", "semana")) %>%
    mutate(
      dpto  = paste0("NACIONAL_", grupo_edad),
      casos = tidyr::replace_na(casos, 0L)
    ) %>%
    select(dpto, anio, semana, fecha_sem, casos)

  bind_rows(serie_dpto, serie_exterior, serie_nacional, serie_edad) %>%
    mutate(
      evento = nombre_evento,
      casos  = as.integer(casos)
    )
}

# ============================================================
# 2. Construir series
# ============================================================

serie_dengue <- construir_serie(dengue_imp, "Dengue")
serie_chik   <- construir_serie(chikungunya_imp, "Chikungunya")

# Validación explícita de la rejilla y de los meses
validar_serie <- function(serie, nombre) {
  div <- serie %>% filter(dpto %in% codigos_divipola)

  stopifnot(n_distinct(div$dpto) == 33L)
  stopifnot(nrow(div) == 33L * 314L)
  stopifnot(all(count(div, dpto)$n == 314L))
  stopifnot(nrow(filter(serie, dpto == "NACIONAL")) == 314L)
  stopifnot(nrow(filter(serie, dpto == "EXTERIOR")) == 314L)
  stopifnot(!anyNA(serie$casos), all(serie$casos >= 0))

  mensual <- serie %>%
    filter(dpto == "NACIONAL") %>%
    mutate(mes = floor_date(fecha_sem, "month")) %>%
    group_by(mes) %>%
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") %>%
    arrange(mes)

  stopifnot(nrow(mensual) == 72L)

  cat("\n===== VALIDACIÓN:", nombre, "=====\n")
  cat("Filas departamentales (33 x 314):", nrow(div), "\n")
  cat("Semanas por departamento:", paste(range(count(div, dpto)$n), collapse = "-"), "\n")
  cat("Filas NACIONAL:", nrow(filter(serie, dpto == "NACIONAL")), "\n")
  cat("Meses nacionales:", nrow(mensual),
      "| rango:", format(min(mensual$mes)), "→", format(max(mensual$mes)), "\n")
  cat("Total de casos:", sum(serie$casos[serie$dpto == "NACIONAL"]), "\n")

  invisible(mensual)
}

mensual_dengue <- validar_serie(serie_dengue, "Dengue")
mensual_chik   <- validar_serie(serie_chik,   "Chikungunya")

# ============================================================
# 3. Guardar
# ============================================================

saveRDS(serie_dengue, here("data", "series", "serie_dengue.rds"))
saveRDS(serie_chik,   here("data", "series", "serie_chik.rds"))

write.csv(serie_dengue, here("data", "series", "serie_dengue.csv"), row.names = FALSE)
write.csv(serie_chik,   here("data", "series", "serie_chik.csv"),   row.names = FALSE)

# ============================================================
# 4. Gráfico exploratorio nacional
# ============================================================

grafico_serie <- function(serie, nombre, color) {
  serie %>%
    filter(dpto == "NACIONAL") %>%
    ggplot(aes(x = fecha_sem, y = casos)) +
    geom_line(color = color, linewidth = 0.6) +
    geom_point(size = 0.5, alpha = 0.4, color = color) +
    labs(
      title = paste0("Serie semanal nacional — ", nombre),
      subtitle = paste0("Total: ",
                        comma(sum(serie$casos[serie$dpto == "NACIONAL"])),
                        " casos"),
      x = "Fecha", y = "Casos"
    ) +
    theme_minimal(base_size = 11)
}

p_d <- grafico_serie(serie_dengue, "Dengue",      "#08519C")
p_c <- grafico_serie(serie_chik,   "Chikungunya", "#B2182B")

print(p_d); print(p_c)

ggsave(here("graficas", "series", "serie_nacional_dengue.png"),
       p_d, width = 10, height = 4, dpi = 150, bg = "white")
ggsave(here("graficas", "series", "serie_nacional_chikungunya.png"),
       p_c, width = 10, height = 4, dpi = 150, bg = "white")

cat("\n✅ Series guardadas en data/series/\n")

