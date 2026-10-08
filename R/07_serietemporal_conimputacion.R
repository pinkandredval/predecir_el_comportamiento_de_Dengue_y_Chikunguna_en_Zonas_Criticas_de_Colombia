# ============================================================
# 07_serie_temporal.R
# Construye series temporales semanales y mensuales
# desde datos imputados (nacional + por departamento + por grupo etario)
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
           anio_num >= 2013, anio_num <= 2025) %>%
    mutate(
      fecha_sem = as.Date(paste0(anio_num, "-01-01")) +
        (semana_num - 1) * 7
    )
  
  # Serie por departamento
  serie_dpto <- df_serie %>%
    count(dpto, anio = anio_num, semana = semana_num, fecha_sem,
          name = "casos") %>%
    arrange(dpto, anio, semana)
  
  # Serie nacional
  serie_nacional <- df_serie %>%
    count(anio = anio_num, semana = semana_num, fecha_sem,
          name = "casos") %>%
    mutate(dpto = "NACIONAL") %>%
    select(dpto, everything())
  
  # Serie por grupo etario (nacional)
  serie_edad <- df_serie %>%
    filter(!is.na(grupo_edad)) %>%
    count(grupo_edad, anio = anio_num, semana = semana_num, fecha_sem,
          name = "casos") %>%
    mutate(dpto = paste0("NACIONAL_", grupo_edad)) %>%
    select(dpto, everything())
  
  bind_rows(serie_dpto, serie_nacional, serie_edad) %>%
    mutate(evento = nombre_evento)
}

# ============================================================
# 2. Construir series
# ============================================================

serie_dengue <- construir_serie(dengue_imp, "Dengue")
serie_chik   <- construir_serie(chikungunya_imp, "Chikungunya")

cat("\n===== SERIE DENGUE =====\n")
cat("Filas:", nrow(serie_dengue), "\n")
cat("Series únicas (dpto):", n_distinct(serie_dengue$dpto), "\n")
cat("Rango:", format(min(serie_dengue$fecha_sem)),
    "→", format(max(serie_dengue$fecha_sem)), "\n")

cat("\n===== SERIE CHIKUNGUNYA =====\n")
cat("Filas:", nrow(serie_chik), "\n")
cat("Series únicas (dpto):", n_distinct(serie_chik$dpto), "\n")

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
