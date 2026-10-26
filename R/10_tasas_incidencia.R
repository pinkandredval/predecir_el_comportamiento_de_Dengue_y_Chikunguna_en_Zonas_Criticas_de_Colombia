# ============================================================
# 10_tasas_incidencia.R
# Tasas acumuladas por 100.000 habitantes y 12 regiones EMA
# Dengue y Chikungunya | 2020-2025
# ============================================================
#
# CASOS
# Series municipales generadas por 07_serietemporal_conimputacion.R:
#   data/series/serie_dengue_municipal.rds
#   data/series/serie_chik_municipal.rds
# columnas: cod_mun_o (DIVIPOLA de 5 dígitos), anio, semana, fecha_sem, casos.
#
# POBLACIÓN
# Proyecciones municipales DANE 2018-2042:
#   data/raw/proyecciones_dane_2018-2042/PPED-AreaMun-2018-2042_VP.xlsx
# Se filtra ÁREA GEOGRÁFICA == "Total" y se suman los municipios de cada
# región EMA. Cartagena (13001) se separa del resto de Bolívar; Bogotá
# (11001), Antioquia (05) y San Andrés (88) se conservan como regiones.
#
# REGLA DE BOGOTÁ
# Bogotá se mantiene con los casos que traiga la serie municipal. En
# dengue la serie tiene 0 casos por decisión de diseño del análisis
# (placebo de miedo sin enfermedad); no se imputa ni se corrige.
# ============================================================

library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(here)

# --- 0. Rutas -----------------------------------------------------

ruta_series <- here("data", "series")
ruta_graficas <- here("graficas")

dir.create(ruta_series, showWarnings = FALSE, recursive = TRUE)
dir.create(ruta_graficas, showWarnings = FALSE, recursive = TRUE)

ruta_dengue_mun <- here("data", "series", "serie_dengue_municipal.rds")
ruta_chik_mun <- here("data", "series", "serie_chik_municipal.rds")

ruta_pob_mun <- here(
  "data", "raw", "proyecciones_dane_2018-2042",
  "PPED-AreaMun-2018-2042_VP.xlsx"
)

# --- 1. Regiones EMA ---------------------------------------------

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

# Asigna el municipio (DIVIPOLA de 5 dígitos) a una de las 12 regiones.
# La regla especial es Bolívar: 13001 es Cartagena y el resto de los
# municipios 13xxx van a Golfo Morrosquillo y Sabana.
asignar_region_ema <- function(cod_mun_o) {
  cod_mun_o <- str_pad(as.character(cod_mun_o), 5, pad = "0")
  dpto <- substr(cod_mun_o, 1, 2)

  case_when(
    cod_mun_o == "11001" ~ "Bogotá D.C.",
    cod_mun_o == "13001" ~ "Cartagena",
    dpto == "05" ~ "Antioquia",
    dpto == "88" ~ "San Andrés y Providencia",
    dpto == "13" ~ "Golfo Morrosquillo y Sabana",
    dpto %in% c("08", "20", "44", "47") ~ "Región Costa Caribe",
    dpto %in% c("15", "25", "41", "73") ~ "Región Central",
    dpto %in% c("17", "63", "66") ~ "Región Eje Cafetero",
    dpto %in% c("76", "19", "52", "27") ~ "Región Pacífico",
    dpto %in% c("68", "54") ~ "Región Santanderes",
    dpto %in% c("50", "85", "81", "99") ~ "Región Llanos Orinoquía",
    dpto %in% c("23", "70") ~ "Golfo Morrosquillo y Sabana",
    dpto %in% c("91", "18", "94", "86", "95", "97") ~ "Región Amazonía",
    dpto == "11" ~ "Bogotá D.C.",
    TRUE ~ "Otras regiones"
  )
}

# --- 2. Lectura de las series municipales -------------------------

leer_serie_municipal <- function(ruta, nombre) {
  if (!file.exists(ruta)) {
    warning(
      "No existe la serie municipal de ", nombre, ": ", ruta, ". ",
      "Ejecute primero R/07_serietemporal_conimputacion.R.",
      call. = FALSE
    )
    return(NULL)
  }

  serie <- readRDS(ruta)

  if (!"cod_mun_o" %in% names(serie)) {
    warning(
      "La serie municipal de ", nombre, " no tiene la columna 'cod_mun_o'. ",
      "Vuelva a ejecutar R/07_serietemporal_conimputacion.R.",
      call. = FALSE
    )
    return(NULL)
  }

  requeridas <- c("cod_mun_o", "anio", "semana", "fecha_sem", "casos")
  faltantes <- setdiff(requeridas, names(serie))

  if (length(faltantes) > 0) {
    warning(
      "Faltan columnas en la serie municipal de ", nombre, ": ",
      paste(faltantes, collapse = ", "), ".",
      call. = FALSE
    )
    return(NULL)
  }

  serie |>
    mutate(
      cod_mun_o = str_pad(str_trim(as.character(cod_mun_o)), 5, pad = "0"),
      anio = as.integer(anio),
      semana = as.integer(semana),
      fecha_sem = as.Date(fecha_sem),
      casos = as.integer(casos)
    ) |>
    filter(
      !is.na(anio),
      anio >= 2020,
      anio <= 2025,
      !is.na(cod_mun_o),
      grepl("^[0-9]{5}$", cod_mun_o),
      substr(cod_mun_o, 1, 2) != "01"
    )
}

# --- 3. Lectura de la población municipal -------------------------

# Busca la hoja municipal con el encabezado real en la fila 8 (skip = 7).
leer_hoja_poblacion_municipal <- function(ruta) {
  if (!file.exists(ruta)) {
    warning(
      "No existe el archivo municipal del DANE: ", ruta,
      call. = FALSE
    )
    return(NULL)
  }

  hojas <- readxl::excel_sheets(ruta)

  for (hoja in hojas) {
    encabezado <- tryCatch(
      suppressWarnings(
        suppressMessages(
          readxl::read_excel(ruta, sheet = hoja, skip = 7, n_max = 1)
        )
      ),
      error = function(e) NULL
    )

    if (is.null(encabezado)) next
    nombres <- names(janitor::clean_names(encabezado))

    if (all(c("mpio", "ano", "area_geografica", "total") %in% nombres)) {
      return(hoja)
    }
  }

  warning(
    "El archivo '", basename(ruta), "' no tiene una hoja con el formato ",
    "municipal esperado (mpio, año, área geográfica, total). ",
    "Hojas disponibles: ", paste(hojas, collapse = ", "), ".",
    call. = FALSE
  )
  NULL
}

cargar_poblacion_municipal <- function(ruta) {
  hoja <- leer_hoja_poblacion_municipal(ruta)
  if (is.null(hoja)) return(NULL)

  pob <- suppressWarnings(
    suppressMessages(
      readxl::read_excel(ruta, sheet = hoja, skip = 7)
    )
  ) |>
    janitor::clean_names()

  requeridas <- c("mpio", "ano", "area_geografica", "total")
  faltantes <- setdiff(requeridas, names(pob))

  if (length(faltantes) > 0) {
    warning(
      "Faltan columnas en la población municipal: ",
      paste(faltantes, collapse = ", "), ".",
      call. = FALSE
    )
    return(NULL)
  }

  pob <- pob |>
    mutate(
      mpio = str_pad(str_trim(as.character(mpio)), 5, pad = "0"),
      anio = as.integer(ano),
      total = as.numeric(total),
      area_geografica = str_trim(as.character(area_geografica))
    ) |>
    filter(
      !is.na(anio),
      anio >= 2020,
      anio <= 2025,
      !is.na(area_geografica),
      area_geografica == "Total",
      !is.na(mpio),
      grepl("^[0-9]{5}$", mpio),
      substr(mpio, 1, 2) != "01",
      !is.na(total),
      total > 0
    ) |>
    mutate(region_ema = asignar_region_ema(mpio))

  sin_region_pob <- pob |>
    filter(region_ema == "Otras regiones") |>
    distinct(mpio)

  if (nrow(sin_region_pob) > 0) {
    warning(
      "Hay municipios sin región EMA en la población: ",
      paste(sin_region_pob$mpio, collapse = ", "), ".",
      call. = FALSE
    )
  }

  pob <- pob |>
    filter(region_ema != "Otras regiones")

  if (nrow(pob) == 0) {
    warning(
      "La población municipal quedó vacía después de filtrar área == 'Total' ",
      "y 2020-2025. Revise el archivo del DANE.",
      call. = FALSE
    )
    return(NULL)
  }

  duplicados <- pob |>
    count(mpio, anio) |>
    filter(n > 1)

  if (nrow(duplicados) > 0) {
    warning(
      "Hay municipios repetidos por año en la población; se conserva la ",
      "primera fila: ",
      paste(sprintf("%s-%s", duplicados$mpio, duplicados$anio), collapse = ", "),
      call. = FALSE
    )
    pob <- pob |>
      distinct(mpio, anio, .keep_all = TRUE)
  }

  faltan_region <- setdiff(regiones_ema, unique(pob$region_ema))
  if (length(faltan_region) > 0) {
    warning(
      "Faltan regiones EMA en la población municipal: ",
      paste(faltan_region, collapse = ", "), ".",
      call. = FALSE
    )
  }

  region_anual <- pob |>
    group_by(region_ema, anio) |>
    summarise(total = sum(total, na.rm = TRUE), .groups = "drop")

  region_promedio <- region_anual |>
    group_by(region_ema) |>
    summarise(
      poblacion_promedio = mean(total, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    municipal = pob,
    region_anual = region_anual,
    region_promedio = region_promedio,
    hoja = hoja
  )
}

# --- 4. Cargar insumos --------------------------------------------

message("Leyendo series municipales...")

dengue_mun <- leer_serie_municipal(ruta_dengue_mun, "dengue")
chik_mun <- leer_serie_municipal(ruta_chik_mun, "chikungunya")

if (is.null(dengue_mun) || is.null(chik_mun)) {
  stop(
    "No se pueden calcular las tasas sin las series municipales. ",
    "Ejecute primero R/07_serietemporal_conimputacion.R.",
    call. = FALSE
  )
}

message("Leyendo población municipal del DANE...")

poblacion_lista <- cargar_poblacion_municipal(ruta_pob_mun)

if (is.null(poblacion_lista)) {
  stop(
    "No se pueden calcular las tasas sin la población municipal. ",
    "Verifique el archivo: ", ruta_pob_mun,
    call. = FALSE
  )
}

poblacion_region_anual <- poblacion_lista$region_anual
poblacion_region_promedio <- poblacion_lista$region_promedio

message(
  "Población municipal cargada desde '", basename(ruta_pob_mun),
  "' (hoja '", poblacion_lista$hoja, "'): ",
  n_distinct(poblacion_lista$municipal$mpio), " municipios, ",
  n_distinct(poblacion_region_anual$region_ema), " regiones, ",
  n_distinct(poblacion_region_anual$anio), " años."
)

# --- 5. Tasas por municipio y región ------------------------------

procesar_serie_municipal <- function(serie, nombre) {
  serie_reg <- serie |>
    mutate(region_ema = asignar_region_ema(cod_mun_o))

  sin_region <- serie_reg |>
    filter(region_ema == "Otras regiones") |>
    distinct(cod_mun_o)

  if (nrow(sin_region) > 0) {
    warning(
      "Hay municipios sin región EMA en la serie de ", nombre, ": ",
      paste(sin_region$cod_mun_o, collapse = ", "), ".",
      call. = FALSE
    )
  }

  serie_reg <- serie_reg |>
    filter(region_ema != "Otras regiones")

  sin_poblacion <- serie_reg |>
    anti_join(
      poblacion_region_anual,
      by = c("region_ema", "anio")
    ) |>
    distinct(region_ema, anio)

  if (nrow(sin_poblacion) > 0) {
    warning(
      "Hay regiones/años sin denominador de población en ", nombre, ": ",
      paste(
        sprintf("%s-%s", sin_poblacion$region_ema, sin_poblacion$anio),
        collapse = ", "
      ),
      ". Esas tasas quedarán como NA.",
      call. = FALSE
    )
  }

  serie_reg |>
    left_join(
      poblacion_region_anual,
      by = c("region_ema", "anio")
    ) |>
    rename(poblacion_region = total) |>
    mutate(
      tasa = casos / poblacion_region * 100000
    ) |>
    select(
      cod_mun_o, anio, semana, fecha_sem, casos,
      region_ema, poblacion_region, tasa
    ) |>
    arrange(region_ema, cod_mun_o, anio, semana)
}

message("Calculando tasas de dengue...")
dengue_tasas <- procesar_serie_municipal(dengue_mun, "dengue")

message("Calculando tasas de chikungunya...")
chik_tasas <- procesar_serie_municipal(chik_mun, "chikungunya")

# --- 6. Ranking de las 12 regiones EMA ----------------------------

casos_region <- dengue_tasas |>
  group_by(region_ema) |>
  summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop")

# Se parte de las 12 regiones para que las que no tengan casos queden en 0.
ranking_dengue <- tibble(region_ema = regiones_ema) |>
  left_join(casos_region, by = "region_ema") |>
  mutate(casos = tidyr::replace_na(casos, 0L)) |>
  left_join(poblacion_region_promedio, by = "region_ema") |>
  mutate(
    tasa_acumulada = casos / poblacion_promedio * 100000
  ) |>
  arrange(desc(tasa_acumulada)) |>
  mutate(
    posicion = row_number(),
    periodo = "2020-2025"
  ) |>
  select(
    posicion, periodo, region_ema,
    casos, poblacion_promedio, tasa_acumulada
  )

if (anyNA(ranking_dengue$tasa_acumulada)) {
  warning(
    "Hay regiones sin población promedio; revise la cobertura municipal.",
    call. = FALSE
  )
}

ranking_csv <- ranking_dengue |>
  mutate(
    poblacion_promedio = round(poblacion_promedio),
    tasa_acumulada = round(tasa_acumulada, 2)
  )

# --- 7. Guardar CSV -----------------------------------------------

write_csv(
  dengue_tasas |> mutate(tasa = round(tasa, 4)),
  here("data", "series", "serie_dengue_tasas.csv")
)

write_csv(
  chik_tasas |> mutate(tasa = round(tasa, 4)),
  here("data", "series", "serie_chik_tasas.csv")
)

write_csv(
  ranking_csv,
  here("data", "series", "ranking_tasas_dengue.csv")
)

message(
  "CSV guardados: serie_dengue_tasas.csv (", nrow(dengue_tasas), " filas), ",
  "serie_chik_tasas.csv (", nrow(chik_tasas), " filas), ",
  "ranking_tasas_dengue.csv (", nrow(ranking_csv), " filas)."
)

# --- 8. Top 12 en consola -----------------------------------------

ranking_consola <- ranking_csv |>
  select(posicion, region_ema, casos, poblacion_promedio, tasa_acumulada)

cat("\n--- Top 12 regiones EMA por tasa acumulada (2020-2025) ---\n")
cat(sprintf(
  "%2s  %-34s  %9s  %12s  %11s\n",
  "#", "REGIÓN EMA", "CASOS", "POBL. PROM.", "TASA"
))
cat(sprintf(
  "%2d  %-34s  %9s  %12s  %11s\n",
  ranking_consola$posicion,
  str_trunc(ranking_consola$region_ema, 34),
  formatC(
    ranking_consola$casos,
    format = "d", big.mark = ".", decimal.mark = ","
  ),
  formatC(
    ranking_consola$poblacion_promedio,
    format = "d", big.mark = ".", decimal.mark = ","
  ),
  formatC(
    ranking_consola$tasa_acumulada,
    format = "f", digits = 2, big.mark = ".", decimal.mark = ","
  )
))
cat("\n")

# --- 9. Gráfico de barras de las 12 regiones ----------------------

datos_plot <- ranking_dengue |>
  mutate(
    region_ema = factor(
      region_ema,
      levels = region_ema[order(tasa_acumulada)]
    )
  )

if (nrow(datos_plot) != 12) {
  warning(
    "El gráfico debería tener 12 regiones EMA, pero tiene ",
    nrow(datos_plot), ".",
    call. = FALSE
  )
}

p_ranking <- ggplot(
  datos_plot,
  aes(x = region_ema, y = tasa_acumulada, fill = tasa_acumulada)
) +
  geom_col(width = 0.72, color = NA, show.legend = FALSE) +
  geom_text(
    aes(
      label = formatC(
        tasa_acumulada,
        format = "f",
        digits = 0,
        big.mark = ".",
        decimal.mark = ","
      )
    ),
    hjust = -0.12,
    size = 3.2,
    color = "grey20"
  ) +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) str_wrap(x, 22)
  ) +
  scale_fill_gradient(low = "#C6DBEF", high = "#08306B", guide = "none") +
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ",",
      accuracy = 1
    ),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Tasa de incidencia acumulada de dengue por regiones EMA",
    subtitle = "Casos por 100.000 habitantes, 2020-2025",
    x = NULL,
    y = "Tasa por 100.000 habitantes",
    caption = paste(
      "Fuente: SIVIGILA/INS (series municipales) y proyecciones municipales DANE 2018-2042.",
      "\nCartagena se separa del resto de Bolívar. La tasa acumulada usa la población anual promedio."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "grey20"),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey35"),
    plot.caption = element_text(hjust = 0, color = "grey45", size = 8),
    plot.margin = margin(10, 20, 10, 10)
  )

ruta_grafica <- here("graficas", "10_ranking_regiones.png")

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = ruta_grafica,
    plot = p_ranking,
    width = 10,
    height = 6.5,
    dpi = 300,
    bg = "white",
    device = ragg::agg_png
  )
} else {
  ggsave(
    filename = ruta_grafica,
    plot = p_ranking,
    width = 10,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

message("Gráfico guardado en: ", ruta_grafica)
message("✅ Script 10 finalizado.")
