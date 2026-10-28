# ============================================================
# 12_indice_atencion.R
# Índice de atención (Google Trends) y miedo excedente por
# regiones de la EMA - Dengue y Chikungunya, 2020-2025
# ============================================================

library(gtrendsR)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(here)

# --- 1. Departamentos (códigos ISO de Google Trends) ---
departamentos <- c(
  "CO-DC", "CO-ANT", "CO-VAC", "CO-ATL", "CO-BOL", "CO-SAN",
  "CO-NSA", "CO-CAL", "CO-QUI", "CO-RIS", "CO-TOL", "CO-SAP",
  "CO-BOY", "CO-CUN", "CO-HUI", "CO-CES", "CO-LAG", "CO-MAG",
  "CO-CAU", "CO-NAR", "CO-CHO", "CO-MET", "CO-CAS", "CO-ARA",
  "CO-VID", "CO-COR", "CO-SUC", "CO-AMA", "CO-CAQ", "CO-GUA",
  "CO-PUT", "CO-GUV", "CO-VAU"
)

# --- 2. Términos de búsqueda (máx. 5 por consulta) ---
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

periodo <- "2020-01-01 2025-12-31"

# --- 3. Diccionario departamento -> región EMA ---
# CO-BOL aparece dos veces: pertenece a "Cartagena" y también a
# "Golfo Morrosquillo y Sabana", según el mapeo de la EMA.
mapeo_ema <- tribble(
  ~codigo, ~region_ema,
  "CO-DC",  "Bogotá D.C.",
  "CO-ANT", "Antioquia",
  "CO-SAP", "San Andrés y Providencia",
  "CO-BOL", "Cartagena",
  "CO-BOY", "Región Central",
  "CO-CUN", "Región Central",
  "CO-HUI", "Región Central",
  "CO-TOL", "Región Central",
  "CO-ATL", "Región Costa Caribe",
  "CO-CES", "Región Costa Caribe",
  "CO-LAG", "Región Costa Caribe",
  "CO-MAG", "Región Costa Caribe",
  "CO-CAL", "Región Eje Cafetero",
  "CO-QUI", "Región Eje Cafetero",
  "CO-RIS", "Región Eje Cafetero",
  "CO-VAC", "Región Pacífico",
  "CO-CAU", "Región Pacífico",
  "CO-NAR", "Región Pacífico",
  "CO-CHO", "Región Pacífico",
  "CO-SAN", "Región Santanderes",
  "CO-NSA", "Región Santanderes",
  "CO-MET", "Región Llanos Orinoquía",
  "CO-CAS", "Región Llanos Orinoquía",
  "CO-ARA", "Región Llanos Orinoquía",
  "CO-VID", "Región Llanos Orinoquía",
  "CO-BOL", "Golfo Morrosquillo y Sabana",
  "CO-COR", "Golfo Morrosquillo y Sabana",
  "CO-SUC", "Golfo Morrosquillo y Sabana",
  "CO-AMA", "Región Amazonía",
  "CO-CAQ", "Región Amazonía",
  "CO-GUA", "Región Amazonía",
  "CO-PUT", "Región Amazonía",
  "CO-GUV", "Región Amazonía",
  "CO-VAU", "Región Amazonía"
)

# Validación del diccionario: todos los departamentos deben estar mapeados.
faltan_mapeo <- setdiff(departamentos, mapeo_ema$codigo)
if (length(faltan_mapeo) > 0) {
  warning(
    "Departamentos sin región asignada: ",
    paste(faltan_mapeo, collapse = ", "),
    call. = FALSE
  )
}

# --- 4. Carpeta de salida ---
dir_trends <- here("data", "D_12")
if (!dir.exists(dir_trends)) {
  dir.create(dir_trends, recursive = TRUE)
}

# --- 5. Descarga y limpieza de una consulta ---
# Devuelve un data.frame con date y hits (suma de los términos) o NULL
# si la consulta falla o no trae datos. Nunca lanza un error: así el
# bucle continúa con el siguiente departamento.
descargar_trends <- function(terminos_enfermedad, region, periodo) {
  consulta <- tryCatch(
    gtrends(
      keyword = terminos_enfermedad,
      geo = region,
      time = periodo
    ),
    error = function(e) {
      # Se devuelve el error para informar la causa exacta.
      e
    }
  )

  # Pausa obligatoria después de cada consulta (evita el error 429 de Google).
  Sys.sleep(5)

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

# --- 6. Descarga por departamento y enfermedad ---
exito_dengue <- setNames(rep(FALSE, length(departamentos)), departamentos)
exito_chik <- setNames(rep(FALSE, length(departamentos)), departamentos)

for (codigo in departamentos) {
  for (enfermedad in names(terminos)) {
    message("Descargando ", enfermedad, " para ", codigo, "...")

    datos <- descargar_trends(
      terminos_enfermedad = terminos[[enfermedad]],
      region = codigo,
      periodo = periodo
    )

    # Si la consulta falló o no trajo datos (ya se avisó con warning()
    # dentro de descargar_trends()), se salta sin romper el bucle.
    if (is.null(datos)) {
      next
    }

    archivo <- file.path(
      dir_trends,
      paste0("trends_", enfermedad, "_", codigo, ".csv")
    )
    write_csv(datos, archivo)

    if (enfermedad == "dengue") {
      exito_dengue[[codigo]] <- TRUE
    } else {
      exito_chik[[codigo]] <- TRUE
    }
  }
}

# --- 7. Resumen de descargas ---
fallidos_dengue <- departamentos[!exito_dengue]
fallidos_chik <- departamentos[!exito_chik]

cat("\n=== Resumen de descargas ===\n")
cat(sprintf(
  "Dengue:      %d/%d departamentos\n",
  sum(exito_dengue), length(departamentos)
))
cat(sprintf(
  "Chikungunya: %d/%d departamentos\n",
  sum(exito_chik), length(departamentos)
))
cat(sprintf(
  "Ambas:       %d/%d departamentos\n",
  sum(exito_dengue & exito_chik), length(departamentos)
))
if (length(fallidos_dengue) > 0) {
  cat("Fallaron en dengue: ", paste(fallidos_dengue, collapse = ", "), "\n", sep = "")
}
if (length(fallidos_chik) > 0) {
  cat("Fallaron en chikungunya: ", paste(fallidos_chik, collapse = ", "), "\n", sep = "")
}

# --- 8. Agregación por región EMA (media, no suma) ---
# Google Trends normaliza cada departamento a 0-100, así que no se pueden
# sumar los hits: la región se construye con la media aritmética.
construir_atencion <- function(archivos, mapeo) {
  if (length(archivos) == 0) {
    return(tibble(
      region_ema = character(),
      date = as.Date(character()),
      atencion = numeric()
    ))
  }

  series <- lapply(archivos, function(archivo) {
    codigo <- str_extract(basename(archivo), "CO-[A-Z]{2,3}")
    datos <- tryCatch(
      read_csv(archivo, show_col_types = FALSE),
      error = function(e) NULL
    )

    if (is.null(datos) || !all(c("date", "hits") %in% names(datos))) {
      warning("Archivo ilegible, se omite: ", basename(archivo), call. = FALSE)
      return(NULL)
    }

    datos |>
      transmute(
        codigo = codigo,
        date = floor_date(as.Date(date), "month"),
        hits = suppressWarnings(as.numeric(hits))
      ) |>
      filter(!is.na(date)) |>
      mutate(hits = if_else(is.na(hits), 0, hits))
  })

  series <- bind_rows(series)

  if (nrow(series) == 0 || !all(c("codigo", "date", "hits") %in% names(series))) {
    return(tibble(
      region_ema = character(),
      date = as.Date(character()),
      atencion = numeric()
    ))
  }

  series |>
    inner_join(mapeo, by = "codigo", relationship = "many-to-many") |>
    group_by(region_ema, date) |>
    summarise(atencion = mean(hits), .groups = "drop") |>
    arrange(region_ema, date)
}

archivos_dengue <- list.files(
  dir_trends,
  pattern = "^trends_dengue_CO-[A-Z]{2,3}\\.csv$",
  full.names = TRUE
)
archivos_chik <- list.files(
  dir_trends,
  pattern = "^trends_chikungunya_CO-[A-Z]{2,3}\\.csv$",
  full.names = TRUE
)

atencion_dengue <- construir_atencion(archivos_dengue, mapeo_ema)
atencion_chik <- construir_atencion(archivos_chik, mapeo_ema)

if (nrow(atencion_dengue) > 0) {
  ruta_atencion <- here("data", "D_12", "atencion_dengue_EMA.csv")
  write_csv(atencion_dengue, ruta_atencion)
  # Se relee el archivo guardado para el cálculo del miedo excedente.
  atencion_dengue <- read_csv(ruta_atencion, show_col_types = FALSE)
}

if (nrow(atencion_chik) > 0) {
  ruta_atencion_chik <- here("data", "D_12", "atencion_chikungunya_EMA.csv")
  write_csv(atencion_chik, ruta_atencion_chik)
  # Se relee el archivo guardado para el cálculo del miedo excedente.
  atencion_chik <- read_csv(ruta_atencion_chik, show_col_types = FALSE)
}

# --- 9. Incidencia mensual por región EMA ---
# Las series vienen a nivel municipal-semanal: se suman los casos de los
# municipios de cada región y se agregan por mes.
leer_incidencia <- function(ruta) {
  read_csv(ruta, show_col_types = FALSE) |>
    mutate(mes = floor_date(fecha_sem + 3, "month")) |>  # jueves de la semana (ISO)
    group_by(region_ema, mes) |>
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop")
}

ruta_incidencia_dengue <- here("data", "D_10", "serie_dengue_tasas.csv")
ruta_incidencia_chik <- here("data", "D_10", "serie_chik_tasas.csv")

# --- 10. Miedo excedente: regresión auxiliar por región ---
# A_t = b0 + b1*Casos_{t-1} + b2*Casos_{t-2} + b3*Casos_{t-3}
#       + efectos fijos de mes + e_t,  con M_t = e_t.
calcular_miedo <- function(atencion, incidencia) {
  datos <- atencion |>
    left_join(incidencia, by = c("region_ema", "date" = "mes")) |>
    arrange(region_ema, date) |>
    group_by(region_ema) |>
    mutate(
      casos_lag1 = lag(casos, 1),
      casos_lag2 = lag(casos, 2),
      casos_lag3 = lag(casos, 3),
      mes_num = month(date)
    ) |>
    ungroup()

  regiones <- sort(unique(datos$region_ema))
  salidas <- vector("list", length(regiones))
  r2 <- rep(NA_real_, length(regiones))

  for (i in seq_along(regiones)) {
    region <- regiones[i]
    d <- datos[datos$region_ema == region, ]
    d$M <- NA_real_

    completo <- stats::complete.cases(
      d[, c("atencion", "casos_lag1", "casos_lag2", "casos_lag3")]
    )

    # Se exigen suficientes observaciones para estimar 3 rezagos + 11 meses.
    if (sum(completo) >= 15) {
      modelo <- tryCatch(
        stats::lm(
          atencion ~ casos_lag1 + casos_lag2 + casos_lag3 + factor(mes_num),
          data = d[completo, ]
        ),
        error = function(e) NULL
      )

      if (!is.null(modelo)) {
        d$M[completo] <- stats::resid(modelo)
        r2[i] <- summary(modelo)$r.squared
        if (!is.finite(r2[i])) {
          r2[i] <- NA_real_
        }
      }
    }

    salidas[[i]] <- d
  }

  list(
    datos = bind_rows(salidas) |> arrange(region_ema, date),
    r2 = tibble(region_ema = regiones, r2 = r2)
  )
}

# Calcula y guarda el miedo excedente de una enfermedad.
# Devuelve la tabla de R-cuadrado (vacía si no hay datos).
procesar_miedo <- function(atencion, ruta_incidencia, ruta_salida, etiqueta) {
  if (nrow(atencion) == 0) {
    warning(
      "No hay datos de atención de ", etiqueta,
      "; no se calcula su miedo excedente.",
      call. = FALSE
    )
    return(tibble(region_ema = character(), r2 = numeric()))
  }

  if (!file.exists(ruta_incidencia)) {
    warning(
      "No se encontró ", ruta_incidencia,
      "; no se calcula el miedo excedente de ", etiqueta, ".",
      call. = FALSE
    )
    return(tibble(region_ema = character(), r2 = numeric()))
  }

  resultado <- calcular_miedo(atencion, leer_incidencia(ruta_incidencia))
  write_csv(resultado$datos, ruta_salida)

  resultado$r2
}

r2_dengue <- procesar_miedo(
  atencion = atencion_dengue,
  ruta_incidencia = ruta_incidencia_dengue,
  ruta_salida = here("data", "D_12", "miedo_excedente_EMA.csv"),
  etiqueta = "dengue"
)

r2_chik <- procesar_miedo(
  atencion = atencion_chik,
  ruta_incidencia = ruta_incidencia_chik,
  ruta_salida = here("data", "D_12", "miedo_excedente_chik_EMA.csv"),
  etiqueta = "chikungunya"
)

# --- 11. R-cuadrado de las regresiones auxiliares ---
imprimir_r2 <- function(r2_tabla, etiqueta) {
  cat(
    "\n=== R-cuadrado de la regresión auxiliar: ", etiqueta, " ===\n",
    sep = ""
  )

  if (nrow(r2_tabla) == 0) {
    cat("No disponible: no hay datos de atención e incidencia para estimar.\n")
    return(invisible(NULL))
  }

  print(as.data.frame(r2_tabla), row.names = FALSE)

  media_r2 <- mean(r2_tabla$r2, na.rm = TRUE)
  if (is.finite(media_r2)) {
    cat(sprintf(
      "R-cuadrado promedio: %.3f (%.1f%% de la atención explicada por casos rezagados y estacionalidad)\n",
      media_r2, 100 * media_r2
    ))
  } else {
    cat("R-cuadrado promedio: no disponible\n")
  }
}

imprimir_r2(r2_dengue, "dengue")
imprimir_r2(r2_chik, "chikungunya")

message("✅ Proceso finalizado.")
