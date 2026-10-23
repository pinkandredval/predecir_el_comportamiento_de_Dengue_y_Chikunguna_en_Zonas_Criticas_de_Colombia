# ============================================================
# 12_indice_atencion.R
# Índice de atención en Google Trends (regiones de la EMA)
# Dengue y Chikungunya | 2020-01-01 a 2025-12-31
# ============================================================

library(gtrendsR)
library(dplyr)
library(readr)
library(here)

# --- 1. Regiones de la EMA (códigos geo de Google) ---
# Departamentos principales que componen las 12 regiones de la EMA (DANE).
regiones_ema <- c(
  "CO-DC", "CO-ANT", "CO-VAC", "CO-ATL", "CO-BOL", "CO-SAN",
  "CO-NSA", "CO-CAL", "CO-QUI", "CO-RIS", "CO-TOL", "CO-SAP"
)

# --- 2. Términos de búsqueda ---
# Google Trends admite máximo 5 términos por consulta y normaliza la escala
# según la consulta; por eso se usa una sola consulta de 5 términos por
# enfermedad (así los hits de los 5 términos comparten la misma escala).
terminos <- list(
  dengue = c(
    "brote dengue", "epidemia dengue", "dengue",
    "casos dengue", "síntomas dengue"
  ),
  chikungunya = c(
    "chikungunya", "chicungunya", "chikunguya",
    "chicunguya", "síntomas chikungunya"
  )
)

# --- 3. Parámetros de la descarga ---
periodo <- "2020-01-01 2025-12-31"

# --- 4. Carpeta de salida ---
dir_trends <- here("data", "trends")
if (!dir.exists(dir_trends)) {
  dir.create(dir_trends, recursive = TRUE)
}

# --- 5. Función auxiliar: una consulta + limpieza ---
# Devuelve un data.frame con date y hits (suma de los 5 términos) o NULL
# si la consulta falla o no trae datos. Nunca lanza un error: así el bucle
# continúa con la siguiente consulta.
descargar_trends <- function(terminos_enfermedad, region, periodo) {
  consulta <- tryCatch(
    gtrends(
      keyword = terminos_enfermedad,
      geo = region,
      time = periodo
    ),
    error = function(e) {
      # Se devuelve el error para poder informar la causa exacta.
      e
    }
  )

  # Pausa obligatoria después de cada consulta (evita el error 429 de Google).
  Sys.sleep(5)

  # El mensaje de advertencia se emite aquí, una sola vez por consulta fallida.
  if (inherits(consulta, "error")) {
    warning(
      "Error en Google Trends para ", region, ": ",
      conditionMessage(consulta),
      call. = FALSE
    )
    return(NULL)
  }

  if (is.null(consulta)) {
    warning("Google Trends no devolvió datos para ", region, ".", call. = FALSE)
    return(NULL)
  }

  serie <- consulta$interest_over_time

  if (is.null(serie) || nrow(serie) == 0) {
    warning("Google Trends no devolvió datos para ", region, ".", call. = FALSE)
    return(NULL)
  }

  if (!all(c("date", "hits") %in% names(serie))) {
    warning(
      "La respuesta para ", region, " no tiene las columnas date/hits.",
      call. = FALSE
    )
    return(NULL)
  }

  serie |>
    mutate(
      # gtrendsR puede devolver "<1" o texto en hits: se fuerza a numérico
      # y esos valores sin volumen se llevan a 0.
      hits = suppressWarnings(as.numeric(hits)),
      hits = if_else(is.na(hits), 0, hits),
      date = as.Date(date)
    ) |>
    group_by(date) |>
    summarise(hits = sum(hits), .groups = "drop") |>
    arrange(date)
}

# --- 6. Descarga por región y enfermedad ---
for (region in regiones_ema) {
  for (enfermedad in names(terminos)) {
    message("Descargando ", enfermedad, " para ", region, "...")

    datos <- descargar_trends(
      terminos_enfermedad = terminos[[enfermedad]],
      region = region,
      periodo = periodo
    )

    # Si la consulta falló o no trajo datos (ya se avisó con warning() dentro
    # de descargar_trends()), se salta esta consulta sin romper el bucle.
    if (is.null(datos)) {
      next
    }

    # archivo <- file.path(dir_trends, paste0(enfermedad, "_", region, ".csv"))
    write_csv(datos, archivo)
  }
}

message("✅ Descarga de Google Trends finalizada")
