# ============================================================
# 13_panel_ema.R
# Panel región-mes (12 regiones EMA × 48 meses, 2022-2025)
# Dengue y chikungunya: incidencia, miedo excedente y clima
# ============================================================
#
# OBJETIVO (plan, §4.4 «Vía principal: panel con efectos fijos»)
#
#   Y_rt = a_r + g_t + s_{r,m(t)}
#          + sum_{k=1..3} b_k  I_{r,t-k}
#          + sum_{k=1..3} v_k  M_{r,t-k}
#          + delta' C_{r,t-1} + e_rt
#
#   Y  ocupación hotelera (EMA, porcentaje)
#   I  incidencia (casos por 100.000 habitantes)
#   M  miedo excedente (residuo de la regresión auxiliar, script 12)
#   C  clima local (lluvia y temperatura), rezagado un mes
#   g_t efectos fijos de año-mes; s_{r,m(t)} estacionalidad por región
#
# Unidad de análisis: 12 regiones EMA × 48 meses (2022-01 a 2025-12)
# = 576 filas. 2020-2021 queda fuera del panel por el COVID.
#
# INSUMOS
#   data/D_10/serie_dengue_tasas.csv        (municipio-semana)
#   data/D_10/serie_chik_tasas.csv          (municipio-semana)
#   data/D_12/miedo_excedente_EMA.csv       (M dengue, script 12)
#   data/D_12/miedo_excedente_chik_EMA.csv  (M chikungunya, script 12)
#   data/raw/EMA/anex-EMA-ago2026.xlsx        (DANE, ocupación)
#   NASA POWER (lluvia y temperatura de las capitales): vía nasapower
#   si está instalado; si no, la misma API con httr2 (no es obligatorio
#   instalar nasapower para que el script funcione).
#
# SALIDAS
#   data/D_13/panel_ema_den_chik.rds
#   data/D_13/panel_ema_den_chik.csv
#   data/D_13/clima/clima_<capital>.csv       (caché de la descarga)
#
# NOTA SOBRE LOS REZAGOS
# Los rezagos se calculan sobre la serie mensual 2020-2025 y el panel se
# recorta después a 2022-2025. Así enero-marzo de 2022 usan diciembre de
# 2021 como insumo (pre-muestra) y ninguna fila 2020-2021 entra al panel.
# ============================================================

library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(here)

# --- 0. Rutas y constantes ----------------------------------------

ruta_ema <- here("data", "raw", "EMA", "anex-EMA-ago2026.xlsx")
ruta_dengue <- here("data", "D_10", "serie_dengue_tasas.csv")
ruta_chik <- here("data", "D_10", "serie_chik_tasas.csv")
dir_trends <- here("data", "D_12")
dir_panel <- here("data", "D_13")
dir_clima <- here("data", "D_13", "clima")

dir.create(dir_panel, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_clima, showWarnings = FALSE, recursive = TRUE)

# Las 12 regiones de la EMA, con los nombres canónicos que ya usan los
# scripts 10 y 12. Todo lo que entre al panel se homologa a esta lista.
regiones_ema <- c(
  "Bogotá D.C.",
  "Cartagena",
  "Región Central",
  "Región Costa Caribe",
  "Región Eje Cafetero",
  "Antioquia",
  "San Andrés y Providencia",
  "Región Pacífico",
  "Región Santanderes",
  "Región Llanos Orinoquía",
  "Golfo Morrosquillo y Sabana",
  "Región Amazonía"
)

# Periodo del panel principal y periodo de cálculo de rezagos.
desde_panel <- as.Date("2022-01-01")
hasta_panel <- as.Date("2025-12-01")
desde_series <- as.Date("2020-01-01")
hasta_series <- as.Date("2025-12-01")

# Clima: única descarga del script (12 series de NASA POWER). Si la caché
# de data/D_13/clima ya existe no se repite la petición; ponga TRUE aquí
# para forzar una descarga nueva.
forzar_descarga_clima <- FALSE

# --- 1. Utilidades -------------------------------------------------

# Texto comparable: sin tildes, en minúsculas y sin espacios sobrantes.
normalizar <- function(x) {
  x |>
    as.character() |>
    iconv(to = "ASCII//TRANSLIT") |>
    str_to_lower() |>
    str_squish()
}

# Número tolerante a separadores de miles y coma decimal.
a_numero <- function(x) {
  x <- str_remove_all(as.character(x), "\\s")
  x <- if_else(
    str_detect(x, ",") & str_detect(x, "\\."),
    str_replace_all(str_replace_all(x, "\\.", ""), ",", "."),
    str_replace_all(x, ",", ".")
  )
  suppressWarnings(as.numeric(x))
}

# Mes en español (texto del DANE) a número 1-12.
meses_es <- c(
  enero = 1, febrero = 2, marzo = 3, abril = 4, mayo = 5, junio = 6,
  julio = 7, agosto = 8, septiembre = 9, setiembre = 9, octubre = 10,
  noviembre = 11, diciembre = 12
)

parsear_mes <- function(x) {
  unname(meses_es[normalizar(x)])
}

# Nombre de archivo seguro para la caché del clima.
slugificar <- function(x) {
  x |>
    iconv(to = "ASCII//TRANSLIT") |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_+|_+$", "")
}

# Homologa el nombre de región de la EMA (o de cualquier fuente) al
# nombre canónico. Devuelve NA para columnas que no son región
# («Total nacional», «Año», «Meses», columnas vacías).
homologar_region_ema <- function(x) {
  # clean_names() deja los nombres con guion bajo («san_andres»); se
  # pasan a espacios para que las expresiones regulares sean legibles.
  x_limpio <- normalizar(x) |>
    str_replace_all("_", " ") |>
    str_squish()

  case_when(
    is.na(x_limpio) | x_limpio == "" ~ NA_character_,
    str_detect(x_limpio, "^total") ~ NA_character_,
    str_detect(x_limpio, "bogota") ~ "Bogotá D.C.",
    str_detect(x_limpio, "cartagena") ~ "Cartagena",
    str_detect(x_limpio, "san andres") ~ "San Andrés y Providencia",
    str_detect(x_limpio, "antioquia") ~ "Antioquia",
    str_detect(x_limpio, "golfo|morrosquillo|sabana") ~
      "Golfo Morrosquillo y Sabana",
    str_detect(x_limpio, "central") ~ "Región Central",
    str_detect(x_limpio, "caribe") ~ "Región Costa Caribe",
    str_detect(x_limpio, "eje") ~ "Región Eje Cafetero",
    str_detect(x_limpio, "pacifico") ~ "Región Pacífico",
    str_detect(x_limpio, "santander") ~ "Región Santanderes",
    str_detect(x_limpio, "llanos") ~ "Región Llanos Orinoquía",
    str_detect(x_limpio, "amazon") ~ "Región Amazonía",
    TRUE ~ NA_character_
  )
}

# Recorta una serie mensual al periodo de trabajo.
recortar_periodo <- function(datos, desde = desde_series, hasta = hasta_series) {
  datos |>
    filter(!is.na(mes), mes >= desde, mes <= hasta)
}

# --- 2. Ocupación hotelera (Y) desde el anexo de la EMA ------------
#
# El anexo del DANE trae una hoja por cuadro y, dentro de cada hoja,
# varios títulos antes del encabezado real. Se listan las hojas con
# excel_sheets(), se inspeccionan las primeras filas con n_max = 15 y
# se ajusta skip a la fila donde aparecen «Año» y «Meses».

detectar_encabezado_ema <- function(ruta, hoja, n_max = 15) {
  vista <- tryCatch(
    suppressWarnings(
      suppressMessages(
        readxl::read_excel(ruta, sheet = hoja, col_names = FALSE, n_max = n_max)
      )
    ),
    error = function(e) NULL
  )

  if (is.null(vista) || nrow(vista) == 0) {
    return(list(vista = NULL, fila = NA_integer_))
  }

  fila <- NA_integer_
  for (i in seq_len(nrow(vista))) {
    celdas <- vista[i, ] |> as.character() |> str_squish()
    celdas <- celdas[!is.na(celdas) & celdas != ""]
    tiene_anio <- any(str_detect(celdas, regex("^a(ñ|n)o$", ignore_case = TRUE)))
    tiene_mes <- any(str_detect(celdas, regex("^mes(es)?$", ignore_case = TRUE)))
    if (tiene_anio && tiene_mes && length(celdas) >= 10) {
      fila <- i
      break
    }
  }

  list(vista = vista, fila = fila)
}

# Busca la hoja de ocupación mensual por región y devuelve los datos ya
# leídos con el skip correcto. Se apoya en el nombre de la hoja, pero
# valida el contenido: debe tener encabezado y al menos 10 regiones.
leer_ocupacion_ema <- function(ruta) {
  if (!file.exists(ruta)) {
    stop("No existe el anexo de la EMA: ", ruta, call. = FALSE)
  }

  hojas <- readxl::excel_sheets(ruta)
  message("Hojas del anexo EMA (", length(hojas), "):")
  message("  ", paste(hojas, collapse = " | "))

  # Prioridad: la hoja de ocupación mensual por región («4.2 Porc Mens
  # Ocupación.reg»), luego el resto de hojas de ocupación.
  candidatas <- hojas[str_detect(normalizar(hojas), "ocupaci")]
  if (length(candidatas) == 0) {
    candidatas <- hojas
  }
  preferidas <- candidatas[
    str_detect(normalizar(candidatas), "mens\\s+ocupacion")
  ]
  orden <- unique(c(preferidas, candidatas))

  for (hoja in orden) {
    info <- detectar_encabezado_ema(ruta, hoja, n_max = 15)
    if (is.na(info$fila)) next

    datos <- suppressWarnings(
      suppressMessages(
        readxl::read_excel(ruta, sheet = hoja, skip = info$fila - 1)
      )
    )
    datos <- janitor::clean_names(datos)

    columnas_region <- names(datos)[
      !is.na(homologar_region_ema(names(datos)))
    ]

    if (length(unique(homologar_region_ema(columnas_region))) >= 10) {
      message(
        "Hoja seleccionada: '", hoja, "' (encabezado en la fila ",
        info$fila, "; skip = ", info$fila - 1, ")."
      )
      return(list(datos = datos, hoja = hoja, fila = info$fila))
    }
  }

  stop(
    "Ninguna hoja de '", basename(ruta), "' tiene el cuadro de ocupación ",
    "mensual por región. Hojas disponibles: ",
    paste(hojas, collapse = ", "), ".",
    call. = FALSE
  )
}

construir_y_ema <- function(ruta) {
  leido <- leer_ocupacion_ema(ruta)
  datos <- leido$datos

  # No todos los años vienen como número limpio: el DANE marca los
  # provisionales como «2023p», «2024p», etc.
  col_anio <- names(datos)[str_detect(names(datos), "^(ano|anio)$")]
  col_meses <- names(datos)[str_detect(names(datos), "^mes(es)?$")]

  if (length(col_anio) == 0 || length(col_meses) == 0) {
    stop(
      "No se encontraron las columnas de año y mes en la hoja '",
      leido$hoja, "'. Columnas: ", paste(names(datos), collapse = ", "), ".",
      call. = FALSE
    )
  }
  col_anio <- col_anio[1]
  col_meses <- col_meses[1]

  # Las columnas de región llegan con tipos mezclados (unas numéricas y
  # otras como texto por las notas del DANE); se homogenizan antes del
  # pivot para que pivot_longer() no falle.
  cols_valor <- setdiff(names(datos), c(col_anio, col_meses))
  datos <- datos |>
    mutate(across(all_of(cols_valor), a_numero)) |>
    # El DANE deja «Año» solo en la fila de enero: se arrastra hacia abajo.
    fill(all_of(col_anio), .direction = "down")

  datos |>
    mutate(
      anio = as.integer(str_extract(as.character(.data[[col_anio]]), "[0-9]{4}")),
      mes_num = parsear_mes(.data[[col_meses]])
    ) |>
    filter(!is.na(anio), !is.na(mes_num)) |>
    pivot_longer(
      cols = all_of(cols_valor),
      names_to = "region_raw",
      values_to = "Y_raw"
    ) |>
    mutate(
      region_ema = homologar_region_ema(region_raw),
      Y = Y_raw,
      mes = make_date(anio, mes_num, 1)
    ) |>
    filter(!is.na(region_ema)) |>
    dplyr::select(region_ema, mes, Y) |>
    distinct(region_ema, mes, .keep_all = TRUE) |>
    arrange(region_ema, mes)
}

message("\nLeyendo ocupación hotelera (EMA)...")
y_ema <- construir_y_ema(ruta_ema)

faltan_y_panel <- sum(is.na(
  y_ema$Y[y_ema$mes >= desde_panel & y_ema$mes <= hasta_panel]
))

if (faltan_y_panel > 0) {
  warning(
    "La ocupación hotelera tiene ", faltan_y_panel,
    " valores faltantes dentro de 2022-2025.",
    call. = FALSE
  )
}

if (anyNA(y_ema$Y)) {
  message(
    "  (", sum(is.na(y_ema$Y)),
    " NA en la serie completa; no afectan al panel si están fuera de 2022-2025)"
  )
}

faltan_y <- setdiff(regiones_ema, unique(y_ema$region_ema))
if (length(faltan_y) > 0) {
  stop(
    "Faltan regiones EMA en la ocupación hotelera: ",
    paste(faltan_y, collapse = ", "), ".",
    call. = FALSE
  )
}

message(
  "Ocupación: ", nrow(y_ema), " filas, ", n_distinct(y_ema$region_ema),
  " regiones, ", format(min(y_ema$mes)), " a ", format(max(y_ema$mes)), "."
)

# --- 3. Incidencia (I) desde las tasas municipales -----------------
#
# Las series vienen por municipio-semana. Se suman los casos de los
# municipios de cada región y se agregan por mes; la tasa mensual usa
# la población regional del año (constante dentro del año).

leer_incidencia <- function(ruta, nombre) {
  if (!file.exists(ruta)) {
    stop(
      "No existe la serie de tasas de ", nombre, ": ", ruta, ". ",
      "Ejecute primero R/10_tasas_incidencia.R.",
      call. = FALSE
    )
  }

  datos <- read_csv(ruta, show_col_types = FALSE)

  requeridas <- c("fecha_sem", "casos", "region_ema", "poblacion_region")
  faltantes <- setdiff(requeridas, names(datos))
  if (length(faltantes) > 0) {
    stop(
      "Faltan columnas en la serie de ", nombre, ": ",
      paste(faltantes, collapse = ", "), ".",
      call. = FALSE
    )
  }

  sin_region <- setdiff(unique(datos$region_ema), regiones_ema)
  if (length(sin_region) > 0) {
    warning(
      "Regiones fuera de la EMA en la serie de ", nombre, ": ",
      paste(sin_region, collapse = ", "), ". Se omiten.",
      call. = FALSE
    )
  }

  datos |>
    filter(region_ema %in% regiones_ema) |>
    mutate(mes = floor_date(as.Date(fecha_sem) + 3, "month")) |>  # jueves de la semana (ISO)
    group_by(region_ema, mes) |>
    summarise(
      casos = sum(casos, na.rm = TRUE),
      poblacion = dplyr::first(poblacion_region),
      .groups = "drop"
    ) |>
    mutate(I = casos / poblacion * 100000) |>
    arrange(region_ema, mes)
}

message("\nLeyendo incidencia de dengue y chikungunya...")
i_dengue <- leer_incidencia(ruta_dengue, "dengue")
i_chik <- leer_incidencia(ruta_chik, "chikungunya")

names(i_dengue)[names(i_dengue) == "casos"] <- "casos_dengue"
names(i_dengue)[names(i_dengue) == "poblacion"] <- "poblacion_dengue"
names(i_dengue)[names(i_dengue) == "I"] <- "I_dengue"

names(i_chik)[names(i_chik) == "casos"] <- "casos_chik"
names(i_chik)[names(i_chik) == "poblacion"] <- "poblacion_chik"
names(i_chik)[names(i_chik) == "I"] <- "I_chik"

message(
  "Dengue: ", nrow(i_dengue), " regiones-mes, ",
  formatC(
    sum(i_dengue$casos_dengue),
    format = "d", big.mark = ".", decimal.mark = ","
  ),
  " casos. ",
  "Chikungunya: ", nrow(i_chik), " regiones-mes, ",
  formatC(
    sum(i_chik$casos_chik),
    format = "d", big.mark = ".", decimal.mark = ","
  ),
  " casos."
)

# --- 4. Miedo excedente (M) desde el script 12 ---------------------
#
# El script 12 deja un archivo por enfermedad:
#   miedo_excedente_EMA.csv       (dengue)
#   miedo_excedente_chik_EMA.csv  (chikungunya)
# Este script no estima M: lo único que hace es leerlo y homologarlo.
# Se aceptan también un único archivo con columna `evento` (se separa con
# pivot_wider) o nombres de archivo que indiquen la enfermedad.

asignar_evento <- function(archivo) {
  nombre <- normalizar(basename(archivo))
  case_when(
    str_detect(nombre, "chik|chicung") ~ "chik",
    str_detect(nombre, "dengue") ~ "dengue",
    TRUE ~ NA_character_
  )
}

# Regresión auxiliar del script 12 (referencia, no se recalcula aquí):
#   A_t = b0 + b1*casos_{t-1} + b2*casos_{t-2} + b3*casos_{t-3}
#         + efectos fijos de mes + e_t,  con M_t = e_t.
# Si M falta, el panel se construye con M = NA y se avisa: el miedo
# excedente es responsabilidad del script 12, no de este.

cargar_miedo <- function(dir_trends) {
  archivos <- list.files(
    dir_trends,
    pattern = "^miedo_excedente.*\\.csv$",
    full.names = TRUE
  )

  # Si conviven una versión oficial (script 12) y una reconstruida, la
  # oficial se lee primero para que sea la que entra al panel.
  archivos <- archivos[
    order(str_detect(normalizar(basename(archivos)), "reconstruido"))
  ]

  if (length(archivos) == 0) {
    warning(
      "No hay archivos de miedo excedente en ", dir_trends, ". ",
      "M_dengue y M_chik quedarán como NA.",
      call. = FALSE
    )
    return(list(dengue = NULL, chik = NULL))
  }

  message(
    "Archivos de miedo excedente: ",
    paste(basename(archivos), collapse = ", ")
  )

  piezas <- lapply(archivos, function(archivo) {
    d <- read_csv(archivo, show_col_types = FALSE)
    if (!"M" %in% names(d)) {
      warning("Archivo sin columna M, se omite: ", basename(archivo),
              call. = FALSE)
      return(NULL)
    }
    if ("date" %in% names(d) && !"mes" %in% names(d)) {
      d <- rename(d, mes = date)
    }
    d$archivo <- basename(archivo)
    d
  })
  piezas <- bind_rows(piezas)

  if (nrow(piezas) == 0) {
    return(list(dengue = NULL, chik = NULL))
  }

  piezas <- piezas |>
    mutate(
      region_ema = as.character(region_ema),
      mes = floor_date(as.Date(mes), "month"),
      M = as.numeric(M)
    ) |>
    filter(region_ema %in% regiones_ema, !is.na(mes))

  # (a) Columna `evento` o `enfermedad`: se separa con pivot_wider().
  col_evento <- names(piezas)[
    str_detect(normalizar(names(piezas)), "^(evento|enfermedad)$")
  ]

  salida <- list(dengue = NULL, chik = NULL)

  if (length(col_evento) > 0) {
    evento <- normalizar(piezas[[col_evento[1]]])
    ancho <- piezas |>
      mutate(evento = if_else(str_detect(evento, "chik|chicung"),
                              "chik", "dengue")) |>
      dplyr::select(region_ema, mes, evento, M) |>
      pivot_wider(
        names_from = evento,
        values_from = M,
        values_fn = dplyr::first
      )

    if ("dengue" %in% names(ancho)) {
      salida$dengue <- ancho |>
        dplyr::select(region_ema, mes, M_dengue = dengue) |>
        filter(!is.na(M_dengue))
    }
    if ("chik" %in% names(ancho)) {
      salida$chik <- ancho |>
        dplyr::select(region_ema, mes, M_chik = chik) |>
        filter(!is.na(M_chik))
    }
    return(salida)
  }

  # (b) Sin columna de evento: el archivo se asigna por enfermedad.
  # Primero por el nombre del archivo; si no dice nada, comparando la
  # columna `casos` con el total mensual de cada serie.
  total_dengue <- sum(i_dengue$casos_dengue, na.rm = TRUE)
  total_chik <- sum(i_chik$casos_chik, na.rm = TRUE)

  por_archivo <- split(piezas, piezas$archivo)

  for (nombre in names(por_archivo)) {
    d <- por_archivo[[nombre]]
    evento <- asignar_evento(nombre)

    if (is.na(evento) && "casos" %in% names(d)) {
      total_archivo <- sum(d$casos, na.rm = TRUE)
      evento <- if (
        abs(total_archivo - total_chik) < abs(total_archivo - total_dengue)
      ) "chik" else "dengue"
      message(
        "  '", nombre, "' no indica la enfermedad; por la columna casos se ",
        "asigna a ", evento, "."
      )
    }

    if (is.na(evento) || is.null(salida[[evento]])) {
      if (is.na(evento)) {
        evento <- "dengue"
        message("  '", nombre, "' se asigna a dengue por defecto.")
      }
      salida[[evento]] <- d |>
        transmute(region_ema, mes, M) |>
        filter(!is.na(M))
      message("  M_", evento, " desde '", nombre, "'.")
    }
  }

  salida
}

message("\nLeyendo miedo excedente...")
miedo <- cargar_miedo(dir_trends)

if (is.null(miedo$chik)) {
  warning(
    "No se encontró ", here("data", "D_12", "miedo_excedente_chik_EMA.csv"),
    ". Ejecute R/12_indice_atencion.R antes de este script; ",
    "M_chik quedará como NA.",
    call. = FALSE
  )
}

m_dengue <- if (is.null(miedo$dengue)) {
  NULL
} else {
  rename(miedo$dengue, M_dengue = M)
}

m_chik <- if (is.null(miedo$chik)) {
  NULL
} else {
  rename(miedo$chik, M_chik = M)
}

# --- 5. Clima (C) de las capitales regionales ----------------------
#
# No todas las regiones EMA coinciden con una ciudad: se usa la capital
# más representativa de cada agrupación (la misma que reporta el DANE en
# la nota del cuadro 4.2). La descarga va por nasapower si el paquete
# está instalado (no viene en el renv.lock de este proyecto); si no, se
# consulta la misma API mensual de NASA POWER con httr2, que sí está.
# Cada respuesta se cachea en data/D_13/clima para no repetir la descarga.

capitales_ema <- tribble(
  ~region_ema,                     ~capital,        ~lat,    ~lon,
  "Bogotá D.C.",                   "Bogotá",         4.7110, -74.0721,
  "Antioquia",                     "Medellín",       6.2442, -75.5812,
  "Cartagena",                     "Cartagena",     10.3910, -75.4794,
  "San Andrés y Providencia",      "San Andrés",    12.5310, -81.7020,
  "Región Central",                "Ibagué",         4.4389, -75.2322,
  "Región Costa Caribe",           "Barranquilla",  10.9685, -74.7813,
  "Región Eje Cafetero",           "Pereira",        4.8133, -75.6961,
  "Región Pacífico",               "Cali",           3.4516, -76.5320,
  "Región Santanderes",            "Bucaramanga",    7.1193, -73.1227,
  "Región Llanos Orinoquía",       "Villavicencio",  4.1420, -73.6266,
  "Golfo Morrosquillo y Sabana",   "Montería",       8.7480, -75.8814,
  "Región Amazonía",               "Florencia",      1.6144, -75.6062
)

# Petición directa a la API mensual de NASA POWER (misma fuente que
# nasapower). Devuelve yyyymm, PRECTOTCORR (mm/día) y T2M (°C).
consultar_power_api <- function(lat, lon, anio_inicio, anio_fin) {
  url <- "https://power.larc.nasa.gov/api/temporal/monthly/point"
  parametros <- list(
    parameters = "PRECTOTCORR,T2M",
    community = "AG",
    latitude = lat,
    longitude = lon,
    start = anio_inicio,
    end = anio_fin,
    format = "JSON"
  )

  if (requireNamespace("httr2", quietly = TRUE)) {
    peticion <- httr2::request(url)
    peticion <- do.call(httr2::req_url_query, c(list(peticion), parametros))
    peticion <- httr2::req_timeout(peticion, 120)
    peticion <- httr2::req_retry(peticion, max_tries = 3)
    respuesta <- httr2::req_perform(peticion)
    cuerpo <- httr2::resp_body_json(respuesta, simplifyVector = FALSE)
  } else {
    consulta <- paste0(
      url, "?",
      paste(names(parametros), parametros, sep = "=", collapse = "&")
    )
    cuerpo <- jsonlite::fromJSON(consulta, simplifyVector = FALSE)
  }

  parametro <- cuerpo$properties$parameter
  lluvia <- unlist(parametro$PRECTOTCORR, use.names = TRUE)
  temp <- unlist(parametro$T2M, use.names = TRUE)

  tibble(
    yyyymm = names(lluvia),
    PRECTOTCORR = as.numeric(lluvia),
    T2M = as.numeric(temp)
  )
}

# La API mensual de POWER devuelve, además de los 12 meses, una fila
# anual por año con el mes 13 (p. ej. «202113»). No es un mes: se descarta
# antes de convertir a fecha (si no, ym() avisa «5 failed to parse»).
solo_meses_power <- function(datos) {
  datos |>
    mutate(yyyymm = as.character(yyyymm)) |>
    filter(str_detect(yyyymm, "^[0-9]{4}(0[1-9]|1[0-2])$"))
}

descargar_clima_region <- function(lat, lon, anio_inicio, anio_fin,
                                   ruta_cache, forzar = FALSE) {
  if (file.exists(ruta_cache) && !forzar) {
    cache <- tryCatch(
      read_csv(ruta_cache, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (!is.null(cache) &&
        all(c("yyyymm", "PRECTOTCORR", "T2M") %in% names(cache))) {
      return(solo_meses_power(cache))
    }
  }

  clima <- tryCatch(
    {
      if (requireNamespace("nasapower", quietly = TRUE)) {
        # nasapower pide c(lon, lat): longitud primero.
        nasapower::get_power(
          community = "AG",
          lonlat = c(lon, lat),
          dates = c(
            paste0(anio_inicio, "-01-01"),
            paste0(anio_fin, "-12-31")
          ),
          temporal_api = "monthly",
          pars = c("PRECTOTCORR", "T2M")
        ) |>
          transmute(
            yyyymm = as.character(YYYYMM),
            PRECTOTCORR = as.numeric(PRECTOTCORR),
            T2M = as.numeric(T2M)
          )
      } else {
        consultar_power_api(lat, lon, anio_inicio, anio_fin)
      }
    },
    error = function(e) e
  )

  if (inherits(clima, "error")) {
    warning(
      "No se pudo descargar el clima (", lat, ", ", lon, "): ",
      conditionMessage(clima),
      call. = FALSE
    )
    return(NULL)
  }

  clima <- solo_meses_power(clima)

  write_csv(clima, ruta_cache)
  clima
}

construir_clima <- function(capitales, anio_inicio, anio_fin, forzar = FALSE) {
  if (!requireNamespace("nasapower", quietly = TRUE)) {
    message(
      "nasapower no está instalado; se consulta la API mensual de NASA ",
      "POWER directamente (mismos datos)."
    )
  }

  lista <- vector("list", length = nrow(capitales))

  for (i in seq_len(nrow(capitales))) {
    region <- capitales$region_ema[i]
    capital <- capitales$capital[i]

    ruta_cache <- file.path(
      dir_clima,
      paste0("clima_", slugificar(capital), ".csv")
    )

    origen <- if (file.exists(ruta_cache) && !forzar) {
      " [caché]"
    } else {
      " [descarga NASA POWER]"
    }
    message("  Clima ", region, " (", capital, ")...", origen)

    datos <- descargar_clima_region(
      lat = capitales$lat[i],
      lon = capitales$lon[i],
      anio_inicio = anio_inicio,
      anio_fin = anio_fin,
      ruta_cache = ruta_cache,
      forzar = forzar
    )

    if (is.null(datos)) next

    lista[[i]] <- datos |>
      mutate(
        region_ema = region,
        mes = ym(yyyymm),
        # PRECTOTCORR viene en mm/día (promedio mensual): se pasa a mm/mes.
        lluvia_mm = PRECTOTCORR * days_in_month(mes),
        temp_c = T2M
      ) |>
      filter(!is.na(mes), mes >= desde_series, mes <= hasta_series) |>
      filter(lluvia_mm > -100, temp_c > -100) |>
      dplyr::select(region_ema, mes, lluvia_mm, temp_c)
  }

  clima <- bind_rows(lista)

  if (nrow(clima) == 0) {
    return(NULL)
  }

  sin_clima <- setdiff(regiones_ema, unique(clima$region_ema))
  if (length(sin_clima) > 0) {
    warning(
      "Sin clima para: ", paste(sin_clima, collapse = ", "), ".",
      call. = FALSE
    )
  }

  clima |>
    distinct(region_ema, mes, .keep_all = TRUE) |>
    arrange(region_ema, mes)
}

message("\nClima de las capitales (NASA POWER; se usa la caché si existe)...")
# Se pide desde 2021 para que el rezago de enero de 2022 exista.
clima <- construir_clima(
  capitales_ema, 2021, 2025,
  forzar = forzar_descarga_clima
)

if (is.null(clima)) {
  warning(
    "No se pudo descargar el clima y no hay caché en ", dir_clima, ". ",
    "El panel se guarda con las columnas de clima en NA: vuelva a ejecutar ",
    "el script con conexión a internet para completarlas.",
    call. = FALSE
  )
}

# --- 6. Unión del panel --------------------------------------------
#
# Se parte de la rejilla completa región × mes (2020-2025) para que los
# rezagos no salten meses, y se hace full_join() de las seis fuentes.

rejilla <- expand_grid(
  region_ema = regiones_ema,
  mes = seq(desde_series, hasta_series, by = "1 month")
)

message("\nUniendo las bases...")

panel_ext <- rejilla |>
  full_join(recortar_periodo(y_ema), by = c("region_ema", "mes")) |>
  full_join(recortar_periodo(i_dengue), by = c("region_ema", "mes")) |>
  full_join(recortar_periodo(i_chik), by = c("region_ema", "mes"))

if (!is.null(m_dengue)) {
  panel_ext <- full_join(
    panel_ext, recortar_periodo(m_dengue), by = c("region_ema", "mes")
  )
} else {
  panel_ext <- mutate(panel_ext, M_dengue = NA_real_)
}

if (!is.null(m_chik)) {
  panel_ext <- full_join(
    panel_ext, recortar_periodo(m_chik), by = c("region_ema", "mes")
  )
} else {
  panel_ext <- mutate(panel_ext, M_chik = NA_real_)
}

if (!is.null(clima)) {
  panel_ext <- full_join(panel_ext, clima, by = c("region_ema", "mes"))
} else {
  panel_ext <- mutate(panel_ext, lluvia_mm = NA_real_, temp_c = NA_real_)
}

panel_ext <- arrange(panel_ext, region_ema, mes)

# Columnas que deben existir aunque la fuente no haya traído nada.
for (col in c("Y", "casos_dengue", "I_dengue", "casos_chik", "I_chik",
              "M_dengue", "M_chik", "lluvia_mm", "temp_c")) {
  if (!col %in% names(panel_ext)) {
    panel_ext[[col]] <- NA_real_
  }
}

# --- 7. Rezagos y estacionalidad -----------------------------------

panel_ext <- panel_ext |>
  group_by(region_ema) |>
  mutate(
    casos_dengue_lag1 = lag(casos_dengue, 1),
    casos_dengue_lag2 = lag(casos_dengue, 2),
    casos_dengue_lag3 = lag(casos_dengue, 3),
    I_dengue_lag1 = lag(I_dengue, 1),
    I_dengue_lag2 = lag(I_dengue, 2),
    I_dengue_lag3 = lag(I_dengue, 3),
    M_dengue_lag1 = lag(M_dengue, 1),
    M_dengue_lag2 = lag(M_dengue, 2),
    M_dengue_lag3 = lag(M_dengue, 3),
    casos_chik_lag1 = lag(casos_chik, 1),
    casos_chik_lag2 = lag(casos_chik, 2),
    casos_chik_lag3 = lag(casos_chik, 3),
    I_chik_lag1 = lag(I_chik, 1),
    I_chik_lag2 = lag(I_chik, 2),
    I_chik_lag3 = lag(I_chik, 3),
    M_chik_lag1 = lag(M_chik, 1),
    M_chik_lag2 = lag(M_chik, 2),
    M_chik_lag3 = lag(M_chik, 3),
    # El plan usa C_{r,t-1}: el clima del mes anterior.
    lluvia_lag1 = lag(lluvia_mm, 1),
    temp_lag1 = lag(temp_c, 1)
  ) |>
  ungroup()

# --- 8. Filtro temporal estricto (2022-2025) -----------------------

panel <- panel_ext |>
  filter(mes >= desde_panel, mes <= hasta_panel) |>
  mutate(
    anio = year(mes),
    mes_num = month(mes),
    anio_mes = format(mes, "%Y-%m")
  ) |>
  dplyr::select(
    region_ema, mes, anio, mes_num, anio_mes, Y,
    casos_dengue, casos_dengue_lag1, casos_dengue_lag2, casos_dengue_lag3,
    I_dengue, I_dengue_lag1, I_dengue_lag2, I_dengue_lag3,
    M_dengue, M_dengue_lag1, M_dengue_lag2, M_dengue_lag3,
    casos_chik, casos_chik_lag1, casos_chik_lag2, casos_chik_lag3,
    I_chik, I_chik_lag1, I_chik_lag2, I_chik_lag3,
    M_chik, M_chik_lag1, M_chik_lag2, M_chik_lag3,
    lluvia_mm, temp_c, lluvia_lag1, temp_lag1
  ) |>
  arrange(region_ema, mes)

# --- 9. Validación --------------------------------------------------

if (nrow(panel) != 576) {
  stop(
    "El panel debe tener 576 filas (12 regiones × 48 meses) y tiene ",
    nrow(panel), ". Revisa la rejilla y el filtro temporal.",
    call. = FALSE
  )
}
stopifnot(nrow(panel) == 576)

if (anyDuplicated(panel[, c("region_ema", "mes")])) {
  stop("Hay regiones-mes duplicadas en el panel.", call. = FALSE)
}

if (n_distinct(panel$region_ema) != 12) {
  stop(
    "El panel debe tener 12 regiones EMA y tiene ",
    n_distinct(panel$region_ema), ".",
    call. = FALSE
  )
}

if (n_distinct(panel$mes) != 48) {
  stop(
    "El panel debe tener 48 meses (2022-01 a 2025-12) y tiene ",
    n_distinct(panel$mes), ".",
    call. = FALSE
  )
}

if (anyNA(panel$Y)) {
  warning(
    "La ocupación hotelera (Y) tiene ", sum(is.na(panel$Y)),
    " valores faltantes en el panel.",
    call. = FALSE
  )
}

# --- 10. Guardado ---------------------------------------------------

ruta_rds <- here("data", "D_13", "panel_ema_den_chik.rds")
ruta_csv <- here("data", "D_13", "panel_ema_den_chik.csv")

saveRDS(panel, ruta_rds)
write_csv(panel, ruta_csv)

message("\nPanel guardado en:")
message("  ", ruta_rds)
message("  ", ruta_csv)

# --- 11. Resumen en consola ----------------------------------------

cat("\n============================================================\n")
cat(" RESUMEN DEL PANEL REGIÓN-MES EMA (2022-2025)\n")
cat("============================================================\n")
cat(sprintf("Filas:    %d\n", nrow(panel)))
cat(sprintf("Columnas: %d\n", ncol(panel)))
cat(sprintf(
  "Regiones: %d | Meses: %d | Periodo: %s a %s\n",
  n_distinct(panel$region_ema), n_distinct(panel$mes),
  format(min(panel$mes), "%Y-%m"), format(max(panel$mes), "%Y-%m")
))

cat(sprintf(
  "\nValores faltantes de Y (ocupación hotelera): %d de %d (%.1f%%)\n",
  sum(is.na(panel$Y)), nrow(panel), 100 * mean(is.na(panel$Y))
))

cat(sprintf(
  "Y: mínimo %.2f | mediana %.2f | media %.2f | máximo %.2f\n",
  min(panel$Y, na.rm = TRUE),
  stats::median(panel$Y, na.rm = TRUE),
  mean(panel$Y, na.rm = TRUE),
  max(panel$Y, na.rm = TRUE)
))

cat("\nOcupación media por región (%):\n")
resumen_y <- panel |>
  group_by(region_ema) |>
  summarise(
    n = n(),
    faltantes_Y = sum(is.na(Y)),
    Y_media = mean(Y, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(Y_media))

for (i in seq_len(nrow(resumen_y))) {
  cat(sprintf(
    "  %-32s %3d obs  %2d NA  %6.2f\n",
    resumen_y$region_ema[i], resumen_y$n[i], resumen_y$faltantes_Y[i],
    resumen_y$Y_media[i]
  ))
}

cat("\nValores faltantes por columna (solo columnas con NA):\n")
faltantes <- colSums(is.na(panel))
faltantes <- faltantes[faltantes > 0]

if (length(faltantes) == 0) {
  cat("  Ninguno: el panel está completo.\n")
} else {
  for (nombre in names(faltantes)) {
    cat(sprintf("  %-24s %4d\n", nombre, faltantes[[nombre]]))
  }
}

cat("\nCobertura de las variables del modelo:\n")
cobertura <- function(x) sprintf("%d/%d", sum(!is.na(x)), length(x))
cat(sprintf("  Y (ocupación)     %s\n", cobertura(panel$Y)))
cat(sprintf("  I_dengue          %s\n", cobertura(panel$I_dengue)))
cat(sprintf("  M_dengue          %s\n", cobertura(panel$M_dengue)))
cat(sprintf("  I_chik            %s\n", cobertura(panel$I_chik)))
cat(sprintf("  M_chik            %s\n", cobertura(panel$M_chik)))
cat(sprintf("  lluvia_mm         %s\n", cobertura(panel$lluvia_mm)))
cat(sprintf("  temp_c            %s\n", cobertura(panel$temp_c)))

if (anyNA(panel$lluvia_mm) || anyNA(panel$temp_c)) {
  cat(
    "\nAviso: el clima quedó incompleto. Vuelva a ejecutar este script con\n",
    "conexión para descargar NASA POWER; la caché evita repetir la descarga.\n",
    sep = ""
  )
}

cat("\nPrimeras filas del panel:\n")
print(
  as.data.frame(
    panel[1:4, c("region_ema", "mes", "Y", "I_dengue", "M_dengue",
                 "I_chik", "M_chik", "lluvia_lag1", "temp_lag1")]
  ),
  row.names = FALSE
)

message("\n✅ Script 13 finalizado.")

