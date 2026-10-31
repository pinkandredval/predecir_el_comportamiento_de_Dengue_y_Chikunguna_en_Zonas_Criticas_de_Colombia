# ============================================================
# 22_riesgo_credito.R
# Teoría del riesgo: pérdida esperada de la cartera de crédito
# y causalidad de Granger dengue -> morosidad (PD)
# ============================================================
#
# OBJETIVO (módulo de Teoría del Riesgo)
#
#   Traducir el canal epidemiológico a riesgo de crédito:
#
#     EL = PD * LGD * EAD
#
#   EAD  exposición al incumplimiento: saldo total de la cartera del
#        departamento en el trimestre (Superfinanciera, Socrata u2wk-tfe3).
#   PD   probabilidad de incumplimiento: saldo en categorías de riesgo
#        C, D y E sobre el EAD (categorías B, C, D y E de la SFC; la C
#        es «riesgo apreciable», la D «riesgo significativo» y la E
#        «riesgo de incobrabilidad»).
#   LGD  severidad: 45% estándar de Basilea.
#
#   Luego se cruza la EL con el panel epidémico de la EMA (script 13)
#   para comparar trimestres epidémicos (M_dengue > percentil 75) con
#   trimestres normales, y se estima un VAR nacional entre casos de
#   dengue y PD para contrastar la causalidad de Granger.
#
# INSUMOS
#   https://www.datos.gov.co/resource/u2wk-tfe3.json
#     «Saldo de las captaciones y colocaciones por municipios»,
#     Superintendencia Financiera de Colombia.
#   data/D_13/panel_ema_den_chik.rds
#     576 filas: 12 regiones EMA x 48 meses (2022-2025).
#
# SALIDAS
#   data/D_22/perdida_esperada.csv
#   data/D_22/var_granger_credito.csv
#   graficas/G_22/el_estres.png
#
# NOTAS METODOLÓGICAS
#   - La API entrega filas por entidad-municipio y también filas
#     «TOTAL <departamento>» y «TOTAL NACIONAL». Estas últimas se
#     excluyen para no duplicar saldos; los totales departamentales se
#     reconstruyen sumando los municipios.
#   - El código de departamento de Socrata NO es el DIVIPOLA: es un
#     código secuencial alfabético (1 Antioquia, 2 Atlántico,
#     3 Bogotá D.C., 4 Bolívar, ...). El diccionario de abajo traduce
#     ese código a las 12 regiones EMA, con la misma regla del script
#     10: Cartagena se separa del resto de Bolívar y va a su propia
#     región; el resto de Bolívar va a Golfo Morrosquillo y Sabana.
#   - El umbral epidémico es un único percentil 75 calculado sobre
#     todas las observaciones región-trimestre de M_dengue. Cada celda
#     región-trimestre se clasifica como epidémica o normal.
#   - El VAR nacional tiene solo 16 trimestres. Con p = 2 quedan 14
#     observaciones efectivas; la prueba de Granger tiene baja potencia
#     y debe leerse como evidencia descriptiva, no estructural.
#   - Si la API falla, el script genera datos simulados con rnorm(),
#     lanza un warning() y continúa. El resto del pipeline no cambia.
#
# REGLAS DEL PROYECTO
#   - Todas las rutas usan here::here().
#   - Todas las funciones no base se llaman con paquete::funcion().
#   - No se ejecuta git add ni git commit.
# ============================================================

# --- 0. Entorno, paquetes y rutas ---------------------------------

# En Windows, el backend Schannel de curl puede fallar en sesiones sin
# credenciales del almacén de certificados. OpenSSL evita ese error y
# usa el mismo almacén de CA. Solo se fija si el usuario no lo definió.
if (.Platform$OS.type == "windows" && Sys.getenv("CURL_SSL_BACKEND") == "") {
  Sys.setenv(CURL_SSL_BACKEND = "openssl")
}

paquetes <- c(
  "httr2", "jsonlite", "janitor", "dplyr", "tidyr", "lubridate",
  "stringr", "tibble", "readr", "ggplot2", "vars", "here"
)

for (paquete in paquetes) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop(
      "Falta el paquete '", paquete, "'. Instálelo con renv::install(\"",
      paquete, "\").",
      call. = FALSE
    )
  }
}

dir.create(here::here("data", "D_22"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("graficas", "G_22"), recursive = TRUE, showWarnings = FALSE)

ruta_panel <- here::here("data", "D_13", "panel_ema_den_chik.rds")
ruta_perdida <- here::here("data", "D_22", "perdida_esperada.csv")
ruta_var <- here::here("data", "D_22", "var_granger_credito.csv")
ruta_grafica <- here::here("graficas", "G_22", "el_estres.png")

url_cartera <- "https://www.datos.gov.co/resource/u2wk-tfe3.json"

# La API devuelve ~94.000 filas para 2022-2025; se pide de a 50.000
# (el límite que exige el plan) y se pagina con $offset.
lgd <- 0.45
anio_desde <- 2022L
anio_hasta <- 2025L
semilla <- 2025L

# --- 1. Utilidades ------------------------------------------------

# Texto comparable: mayúsculas, sin tildes y sin espacios sobrantes.
normalizar <- function(x) {
  x |>
    as.character() |>
    iconv(to = "ASCII//TRANSLIT") |>
    stringr::str_to_upper() |>
    stringr::str_squish()
}

# Fecha tolerante a POSIXct, Date, ISO con hora y epoch numérico.
parsear_fecha <- function(x) {
  if (inherits(x, "Date")) {
    return(as.Date(x))
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  x_chr <- as.character(x)
  salida <- suppressWarnings(as.Date(substr(x_chr, 1, 10)))

  if (all(is.na(salida))) {
    salida <- as.Date(suppressWarnings(as.numeric(x_chr)), origin = "1970-01-01")
  }

  salida
}

# Devuelve el primer nombre de columna que cumple alguno de los patrones.
# Se usa grepl() para no depender de que los nombres sean exactos.
detectar_columna <- function(nombres, patrones, etiqueta, obligatoria = TRUE) {
  for (patron in patrones) {
    coincide <- grepl(patron, nombres)
    if (any(coincide)) {
      return(nombres[which(coincide)[1L]])
    }
  }

  if (obligatoria) {
    stop(
      "No se encontró la columna de ", etiqueta, ". Patrones probados: ",
      paste(patrones, collapse = ", "), ". Columnas disponibles: ",
      paste(nombres, collapse = ", "), ".",
      call. = FALSE
    )
  }

  NA_character_
}

# Porcentaje con coma decimal y punto de miles, para consola y gráfico.
formatear_pct <- function(x, decimales = 2L) {
  paste0(
    format(
      round(x, decimales),
      nsmall = decimales,
      decimal.mark = ",",
      big.mark = "."
    ),
    "%"
  )
}

# Miles de millones de pesos con coma decimal.
formatear_miles_millones <- function(x) {
  format(
    round(x / 1e9, 1),
    nsmall = 1,
    decimal.mark = ",",
    big.mark = ".",
    scientific = FALSE,
    trim = TRUE
  )
}

# --- 2. Diccionarios de departamentos y regiones EMA --------------

# Código secuencial de Socrata -> nombre del departamento.
# El código 3 es Bogotá D.C. y no tiene fila «TOTAL BOGOTA» porque es
# a la vez departamento y municipio.
deptos_socrata <- tibble::tribble(
  ~codigo, ~departamento,
  1L,  "Antioquia",
  2L,  "Atlántico",
  3L,  "Bogotá D.C.",
  4L,  "Bolívar",
  5L,  "Boyacá",
  6L,  "Caldas",
  7L,  "Caquetá",
  8L,  "Cauca",
  9L,  "Cesar",
  10L, "Córdoba",
  11L, "Cundinamarca",
  12L, "Chocó",
  13L, "Huila",
  14L, "La Guajira",
  15L, "Magdalena",
  16L, "Meta",
  17L, "Nariño",
  18L, "Norte de Santander",
  19L, "Quindío",
  20L, "Risaralda",
  21L, "Santander",
  22L, "Sucre",
  23L, "Tolima",
  24L, "Valle del Cauca",
  25L, "Arauca",
  26L, "Casanare",
  27L, "Putumayo",
  28L, "San Andrés y Providencia",
  29L, "Amazonas",
  30L, "Guainía",
  31L, "Guaviare",
  32L, "Vaupés",
  33L, "Vichada"
)

# Las 12 regiones canónicas de la EMA (las mismas del script 13).
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

# Diccionario del script 10 traducido al código de Socrata.
# Bolívar va a Golfo Morrosquillo y Sabana; Cartagena se separa después
# con la columna de municipio.
region_por_codigo <- c(
  "1"  = "Antioquia",
  "2"  = "Región Costa Caribe",
  "3"  = "Bogotá D.C.",
  "4"  = "Golfo Morrosquillo y Sabana",
  "5"  = "Región Central",
  "6"  = "Región Eje Cafetero",
  "7"  = "Región Amazonía",
  "8"  = "Región Pacífico",
  "9"  = "Región Costa Caribe",
  "10" = "Golfo Morrosquillo y Sabana",
  "11" = "Región Central",
  "12" = "Región Pacífico",
  "13" = "Región Central",
  "14" = "Región Costa Caribe",
  "15" = "Región Costa Caribe",
  "16" = "Región Llanos Orinoquía",
  "17" = "Región Pacífico",
  "18" = "Región Santanderes",
  "19" = "Región Eje Cafetero",
  "20" = "Región Eje Cafetero",
  "21" = "Región Santanderes",
  "22" = "Golfo Morrosquillo y Sabana",
  "23" = "Región Central",
  "24" = "Región Pacífico",
  "25" = "Región Llanos Orinoquía",
  "26" = "Región Llanos Orinoquía",
  "27" = "Región Amazonía",
  "28" = "San Andrés y Providencia",
  "29" = "Región Amazonía",
  "30" = "Región Amazonía",
  "31" = "Región Amazonía",
  "32" = "Región Amazonía",
  "33" = "Región Llanos Orinoquía"
)

# --- 3. Descarga de la cartera (con respaldo simulado) -------------

# Descarga paginada. El $where recorta el periodo en el servidor para no
# bajar los 323.000 registros históricos completos; $order = :id da un
# orden estable entre páginas.
descargar_cartera <- function() {
  peticion <- httr2::request(url_cartera)
  peticion <- httr2::req_url_query(
    peticion,
    "$where" = paste0(
      "fecha_de_corte >= '", anio_desde, "-01-01T00:00:00' AND ",
      "fecha_de_corte <= '", anio_hasta, "-12-31T23:59:59'"
    ),
    "$order" = ":id"
  )
  peticion <- httr2::req_timeout(peticion, 180)
  peticion <- httr2::req_retry(peticion, max_tries = 3)

  paginas <- list()
  desplazamiento <- 0L

  for (pagina in seq_len(10L)) {
    consulta <- httr2::req_url_query(
      peticion,
      "$limit" = 50000L,
      "$offset" = desplazamiento
    )
    respuesta <- httr2::req_perform(consulta)
    cuerpo <- httr2::resp_body_json(respuesta, simplifyVector = TRUE)

    if (is.null(cuerpo) || length(cuerpo) == 0L) {
      break
    }
    if (!is.data.frame(cuerpo)) {
      cuerpo <- as.data.frame(cuerpo, stringsAsFactors = FALSE)
    }

    paginas[[length(paginas) + 1L]] <- cuerpo
    message(sprintf(
      "  Descarga Socrata: página %d con %s filas.",
      pagina,
      format(nrow(cuerpo), big.mark = ".", decimal.mark = ",")
    ))

    if (nrow(cuerpo) < 50000L) {
      break
    }
    desplazamiento <- desplazamiento + nrow(cuerpo)
  }

  if (length(paginas) == 0L) {
    stop("La API no devolvió registros de cartera.", call. = FALSE)
  }

  dplyr::bind_rows(paginas)
}

# Respaldo reproducible si la API falla. Imita el esquema municipal de
# Socrata para que prepare_cartera() no tenga una rama distinta.
simular_cartera <- function() {
  set.seed(semilla)

  trimestres <- seq(
    as.Date("2022-01-01"), as.Date("2025-10-01"), by = "quarter"
  )

  base <- tidyr::expand_grid(
    codigo = deptos_socrata$codigo,
    trimestre = trimestres
  ) |>
    dplyr::left_join(deptos_socrata, by = "codigo") |>
    dplyr::mutate(
      nombre_del_municipio = dplyr::if_else(
        codigo == 4L,
        "RESTO BOLIVAR",
        normalizar(departamento)
      )
    )

  # Cartagena se duplica como municipio aparte para conservar la regla
  # de las 12 regiones EMA.
  cartagena <- base |>
    dplyr::filter(codigo == 4L) |>
    dplyr::mutate(nombre_del_municipio = "CARTAGENA")

  dplyr::bind_rows(base, cartagena) |>
    dplyr::mutate(
      c_digo_del_departamento = codigo,
      c_digo_del_municipio = dplyr::row_number(),
      fecha_de_corte = lubridate::ceiling_date(trimestre, "quarter") - 1,
      cartera_de_cr_ditos = pmax(
        1e8,
        stats::rnorm(dplyr::n(), mean = 2e11, sd = 6e10)
      ),
      categor_a_c_riesgo_apreciable = cartera_de_cr_ditos * pmax(
        0, stats::rnorm(dplyr::n(), mean = 0.035, sd = 0.012)
      ),
      categor_a_d_riesgo = cartera_de_cr_ditos * pmax(
        0, stats::rnorm(dplyr::n(), mean = 0.018, sd = 0.008)
      ),
      categor_a_e_riesgo_de = cartera_de_cr_ditos * pmax(
        0, stats::rnorm(dplyr::n(), mean = 0.008, sd = 0.004)
      )
    ) |>
    dplyr::select(-codigo, -trimestre, -departamento)
}

# --- 4. Preparación de la cartera ----------------------------------

# Limpia, detecta columnas, filtra 2022-2025, excluye totales, mapea a
# regiones EMA y agrega por región-trimestre.
preparar_cartera <- function(datos) {
  datos <- janitor::clean_names(datos)
  nombres <- names(datos)

  col_fecha <- detectar_columna(
    nombres, c("fecha"), "la fecha de corte"
  )
  col_depto <- detectar_columna(
    nombres,
    c(
      "(codigo|c_digo).*departamento",
      "departamento.*(codigo|c_digo)",
      "departamento"
    ),
    "el código del departamento"
  )
  col_municipio <- detectar_columna(
    nombres,
    c("nombre.*municipio", "municipio.*nombre"),
    "el nombre del municipio",
    obligatoria = FALSE
  )
  col_saldo <- detectar_columna(
    nombres,
    c("^cartera.*cr.dito", "^cartera"),
    "el saldo total de la cartera"
  )

  cols_cde <- nombres[
    grepl("^categor.*_[cde]_riesgo", nombres)
  ]

  if (length(cols_cde) == 0L) {
    stop(
      "No se encontraron columnas de categorías de riesgo C, D y E. ",
      "Columnas disponibles: ", paste(nombres, collapse = ", "), ".",
      call. = FALSE
    )
  }

  message("  Columnas detectadas:")
  message("    fecha       -> ", col_fecha)
  message("    departamento-> ", col_depto)
  message(
    "    municipio   -> ",
    if (is.na(col_municipio)) "(no disponible)" else col_municipio
  )
  message("    saldo       -> ", col_saldo)
  message("    categorías C/D/E -> ", length(cols_cde), " columnas")

  datos$fecha <- parsear_fecha(datos[[col_fecha]])
  datos$depto_cod <- suppressWarnings(as.integer(datos[[col_depto]]))
  datos$saldo <- suppressWarnings(as.numeric(datos[[col_saldo]]))

  if (is.na(col_municipio)) {
    datos$municipio <- NA_character_
    warning(
      "La base no trae nombre de municipio; no se podrá separar Cartagena ",
      "del resto de Bolívar. Cartagena quedará dentro de Golfo ",
      "Morrosquillo y Sabana.",
      call. = FALSE
    )
  } else {
    datos$municipio <- as.character(datos[[col_municipio]])
  }

  # Suma de saldos en categorías C, D y E de todas las modalidades.
  matriz_cde <- as.data.frame(
    lapply(datos[cols_cde], function(x) suppressWarnings(as.numeric(x)))
  )
  datos$cde <- rowSums(matriz_cde, na.rm = TRUE)

  # Se excluyen «TOTAL <departamento>» y «TOTAL NACIONAL» para no
  # duplicar la cartera; el total se reconstruye sumando municipios.
  datos <- datos |>
    dplyr::filter(
      is.na(municipio) |
        !grepl("^TOTAL", normalizar(municipio)),
      !is.na(fecha),
      lubridate::year(fecha) >= anio_desde,
      lubridate::year(fecha) <= anio_hasta
    ) |>
    dplyr::select(depto_cod, municipio, fecha, saldo, cde)

  if (nrow(datos) == 0L) {
    stop("La cartera quedó vacía después del filtro 2022-2025.", call. = FALSE)
  }

  datos <- datos |>
    dplyr::left_join(deptos_socrata, by = c("depto_cod" = "codigo")) |>
    dplyr::mutate(
      region_ema = unname(region_por_codigo[as.character(depto_cod)]),
      region_ema = dplyr::if_else(
        !is.na(municipio) &
          grepl("^CARTAGENA", normalizar(municipio)),
        "Cartagena",
        region_ema
      )
    )

  sin_region <- datos |>
    dplyr::filter(is.na(region_ema)) |>
    dplyr::distinct(depto_cod, departamento)

  if (nrow(sin_region) > 0L) {
    warning(
      "Departamentos sin equivalencia EMA (se omiten): ",
      paste(
        sprintf("%s (%s)", sin_region$departamento, sin_region$depto_cod),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  datos <- datos |>
    dplyr::filter(!is.na(region_ema), !is.na(saldo)) |>
    dplyr::mutate(trimestre = lubridate::floor_date(fecha, "quarter"))

  cartera_region <- datos |>
    dplyr::group_by(region_ema, trimestre) |>
    dplyr::summarise(
      EAD = sum(saldo, na.rm = TRUE),
      saldo_cde = sum(cde, na.rm = TRUE),
      n_registros = dplyr::n(),
      n_municipios = dplyr::n_distinct(municipio),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      PD_cruda = dplyr::if_else(EAD > 0, saldo_cde / EAD, NA_real_),
      PD = pmin(pmax(PD_cruda, 0), 1),
      LGD = lgd,
      EL = PD * LGD * EAD
    )

  fuera_rango <- sum(
    cartera_region$PD_cruda < 0 | cartera_region$PD_cruda > 1,
    na.rm = TRUE
  )
  if (fuera_rango > 0L) {
    warning(
      "Había ", fuera_rango, " observaciones región-trimestre con PD ",
      "fuera de [0, 1]; se acotaron a ese intervalo.",
      call. = FALSE
    )
  }

  list(
    region_trimestre = cartera_region,
    n_registros = nrow(datos),
    n_departamentos = dplyr::n_distinct(datos$depto_cod),
    n_municipios = dplyr::n_distinct(datos$municipio),
    columnas = list(
      fecha = col_fecha,
      departamento = col_depto,
      municipio = col_municipio,
      saldo = col_saldo,
      cde = cols_cde
    )
  )
}

# --- 5. Descarga y preparación con respaldo ------------------------

message("\nDescargando cartera de la Superfinanciera (Socrata)...")

cartera_cruda <- tryCatch(
  descargar_cartera(),
  error = function(e) {
    warning(
      "No se pudo descargar la cartera desde Socrata: ",
      conditionMessage(e), ". Se generan datos simulados con rnorm() y ",
      "el script continúa.",
      call. = FALSE
    )
    NULL
  }
)

if (is.null(cartera_cruda) || nrow(cartera_cruda) == 0L) {
  if (!is.null(cartera_cruda)) {
    warning(
      "La API respondió sin filas; se generan datos simulados con rnorm().",
      call. = FALSE
    )
  }
  cartera_cruda <- simular_cartera()
  origen_cartera <- "simulada"
} else {
  origen_cartera <- "Socrata u2wk-tfe3"
}

message("\nPreparando la cartera...")

prep_cartera <- tryCatch(
  preparar_cartera(cartera_cruda),
  error = function(e) {
    warning(
      "No se pudo preparar la cartera descargada: ", conditionMessage(e),
      ". Se usa el respaldo simulado con rnorm().",
      call. = FALSE
    )
    preparar_cartera(simular_cartera())
  }
)

cartera_region <- prep_cartera$region_trimestre

faltan_regiones <- setdiff(regiones_ema, unique(cartera_region$region_ema))
if (length(faltan_regiones) > 0L) {
  warning(
    "Faltan regiones EMA en la cartera: ",
    paste(faltan_regiones, collapse = ", "), ".",
    call. = FALSE
  )
}

message(sprintf(
  "Cartera preparada: %s registros municipales, %d departamentos, ",
  format(prep_cartera$n_registros, big.mark = ".", decimal.mark = ","),
  prep_cartera$n_departamentos
))
message(sprintf(
  "  %d municipios, %d celdas región-trimestre, origen: %s.",
  prep_cartera$n_municipios, nrow(cartera_region), origen_cartera
))

# --- 6. Panel epidémico trimestral ---------------------------------

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel epidémico: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

panel <- readRDS(ruta_panel)

requeridas_panel <- c("region_ema", "mes", "casos_dengue", "M_dengue")
faltantes_panel <- setdiff(requeridas_panel, names(panel))

if (length(faltantes_panel) > 0L) {
  stop(
    "Faltan columnas en el panel epidémico: ",
    paste(faltantes_panel, collapse = ", "), ".",
    call. = FALSE
  )
}

panel_trimestral <- panel |>
  dplyr::mutate(
    trimestre = lubridate::floor_date(as.Date(mes), "quarter")
  ) |>
  dplyr::group_by(region_ema, trimestre) |>
  dplyr::summarise(
    casos_dengue = sum(casos_dengue, na.rm = TRUE),
    M_dengue = if (all(is.na(M_dengue))) {
      NA_real_
    } else {
      mean(M_dengue, na.rm = TRUE)
    },
    n_meses = dplyr::n(),
    .groups = "drop"
  )

if (nrow(panel_trimestral) != 192L) {
  warning(
    "El panel trimestral tiene ", nrow(panel_trimestral),
    " filas; se esperaban 192 (12 regiones x 16 trimestres).",
    call. = FALSE
  )
}

# --- 7. Cruce cartera x epidemiología ------------------------------

cartera_epi <- cartera_region |>
  dplyr::inner_join(panel_trimestral, by = c("region_ema", "trimestre")) |>
  dplyr::arrange(region_ema, trimestre)

if (nrow(cartera_epi) == 0L) {
  stop(
    "El cruce cartera-epidemiología quedó vacío. Revise los nombres de ",
    "las regiones EMA y el periodo del panel.",
    call. = FALSE
  )
}

if (nrow(cartera_epi) < 192L) {
  warning(
    "El cruce tiene ", nrow(cartera_epi), " celdas región-trimestre; ",
    "se esperaban 192. Faltan combinaciones de región o trimestre.",
    call. = FALSE
  )
}

# Umbral epidémico: percentil 75 de M_dengue en todas las celdas
# región-trimestre. Un trimestre epidémico es el que supera ese umbral.
umbral_epidemico <- stats::quantile(
  cartera_epi$M_dengue,
  probs = 0.75,
  na.rm = TRUE,
  names = FALSE
)

if (is.na(umbral_epidemico)) {
  stop(
    "No hay valores válidos de M_dengue para calcular el percentil 75.",
    call. = FALSE
  )
}

cartera_epi <- cartera_epi |>
  dplyr::mutate(
    umbral_epidemico = umbral_epidemico,
    epidemico = !is.na(M_dengue) & M_dengue > umbral_epidemico,
    escenario = dplyr::if_else(epidemico, "Epidémico", "Normal"),
    anio = lubridate::year(trimestre),
    trimestre_num = lubridate::quarter(trimestre),
    trimestre_id = sprintf("%d-T%d", anio, trimestre_num),
    tasa_EL = dplyr::if_else(EAD > 0, EL / EAD, NA_real_)
  )

# --- 8. Pérdida esperada: epidémico vs. normal ---------------------

resumen_el <- cartera_epi |>
  dplyr::group_by(epidemico) |>
  dplyr::summarise(
    n_region_trimestre = dplyr::n(),
    EL_promedio = mean(EL, na.rm = TRUE),
    EL_total = sum(EL, na.rm = TRUE),
    EAD_promedio = mean(EAD, na.rm = TRUE),
    PD_promedio = mean(PD, na.rm = TRUE),
    tasa_EL_promedio = mean(tasa_EL, na.rm = TRUE),
    casos_promedio = mean(casos_dengue, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    escenario = dplyr::if_else(epidemico, "Epidémico", "Normal")
  )

el_normal <- resumen_el$EL_promedio[!resumen_el$epidemico]
el_epi <- resumen_el$EL_promedio[resumen_el$epidemico]
tasa_normal <- resumen_el$tasa_EL_promedio[!resumen_el$epidemico]
tasa_epi <- resumen_el$tasa_EL_promedio[resumen_el$epidemico]

aumento_pct <- if (
  length(el_normal) == 1L && length(el_epi) == 1L &&
    !is.na(el_normal) && !is.na(el_epi) && el_normal > 0
) {
  (el_epi - el_normal) / el_normal * 100
} else {
  NA_real_
}

aumento_tasa_pct <- if (
  length(tasa_normal) == 1L && length(tasa_epi) == 1L &&
    !is.na(tasa_normal) && !is.na(tasa_epi) && tasa_normal > 0
) {
  (tasa_epi - tasa_normal) / tasa_normal * 100
} else {
  NA_real_
}

# --- 9. Gráfico de barras: EL promedio normal vs. epidémico --------

datos_grafico <- resumen_el |>
  dplyr::mutate(
    escenario = factor(escenario, levels = c("Normal", "Epidémico")),
    EL_miles_millones = EL_promedio / 1e9
  )

subtitulo_grafico <- if (is.na(aumento_pct)) {
  paste0(
    "EL = PD × LGD × EAD | LGD = 45%. ",
    "No fue posible comparar ambos escenarios."
  )
} else {
  sprintf(
    paste0(
      "EL = PD × LGD × EAD | LGD = 45%% | ",
      "Aumento en trimestres epidémicos: %s"
    ),
    formatear_pct(aumento_pct)
  )
}

p_el <- ggplot2::ggplot(
  datos_grafico,
  ggplot2::aes(x = escenario, y = EL_miles_millones, fill = escenario)
) +
  ggplot2::geom_col(width = 0.58, colour = NA) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = format(
        round(EL_miles_millones, 1),
        nsmall = 1,
        decimal.mark = ",",
        big.mark = ".",
        scientific = FALSE
      )
    ),
    vjust = -0.45,
    size = 3.6,
    colour = "grey20"
  ) +
  ggplot2::scale_fill_manual(
    values = c("Normal" = "#4C72B0", "Epidémico" = "#C44E52"),
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::scale_y_continuous(
    labels = function(x) {
      format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
    },
    expand = ggplot2::expansion(mult = c(0, 0.16))
  ) +
  ggplot2::labs(
    title = "Pérdida esperada de la cartera: trimestres normales vs. epidémicos",
    subtitle = subtitulo_grafico,
    x = NULL,
    y = "EL promedio (miles de millones de COP)",
    caption = paste0(
      "EAD: saldo total de cartera; PD: saldo C+D+E / EAD; LGD = 45%. ",
      "Umbral epidémico: P75 de M_dengue = ",
      format(round(umbral_epidemico, 2), nsmall = 2, decimal.mark = ","),
      ".\n",
      "Fuente: Superfinanciera (Socrata u2wk-tfe3) y panel EMA (D_13)."
    )
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_grafica,
  plot = p_el,
  width = 8,
  height = 5,
  dpi = 200,
  bg = "white"
)

# --- 10. VAR nacional y Causalidad de Granger ----------------------

nacional_q <- cartera_epi |>
  dplyr::group_by(trimestre) |>
  dplyr::summarise(
    casos_dengue = sum(casos_dengue, na.rm = TRUE),
    EAD = sum(EAD, na.rm = TRUE),
    saldo_cde = sum(saldo_cde, na.rm = TRUE),
    M_dengue = mean(M_dengue, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(PD = dplyr::if_else(EAD > 0, saldo_cde / EAD, NA_real_)) |>
  dplyr::arrange(trimestre)

periodo_nacional <- paste(
  sprintf(
    "%d-T%d",
    lubridate::year(min(nacional_q$trimestre)),
    lubridate::quarter(min(nacional_q$trimestre))
  ),
  "a",
  sprintf(
    "%d-T%d",
    lubridate::year(max(nacional_q$trimestre)),
    lubridate::quarter(max(nacional_q$trimestre))
  )
)

# Extrae el htest de vars::causality() a una fila ordenada.
extraer_granger <- function(objeto, causa, efecto, prueba) {
  gl <- unname(objeto$parameter)

  tibble::tibble(
    prueba = prueba,
    causa = causa,
    efecto = efecto,
    estadistico = unname(objeto$statistic)[1L],
    gl1 = gl[1L],
    gl2 = if (length(gl) >= 2L) gl[2L] else NA_real_,
    p_valor = unname(objeto$p.value)[1L]
  )
}

# Prueba F de exclusión de los rezagos de `causa` en la ecuación
# `ecuacion`. Es la versión ecuación por ecuación de Granger y responde
# directamente por la PD, mientras vars::causality() contrasta todas las
# ecuaciones del sistema a la vez.
probar_granger_ecuacion <- function(modelo, ecuacion, causa) {
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    warning(
      "lmtest no está disponible; se omite la prueba de Granger por ",
      "ecuación.",
      call. = FALSE
    )
    return(NULL)
  }

  ajuste <- modelo$varresult[[ecuacion]]
  if (is.null(ajuste)) {
    return(NULL)
  }

  marco <- stats::model.frame(ajuste)
  terminos <- names(marco)[-1L]
  rezagos_causa <- paste0(causa, ".l", seq_len(modelo$p))

  if (!all(rezagos_causa %in% terminos)) {
    warning(
      "La ecuación ", ecuacion, " no contiene los rezagos de ", causa, ".",
      call. = FALSE
    )
    return(NULL)
  }

  conservar <- setdiff(terminos, rezagos_causa)
  formula_restringida <- stats::reformulate(
    conservar,
    response = names(marco)[1L],
    intercept = FALSE
  )
  restringido <- stats::lm(formula_restringida, data = marco)
  tabla <- as.data.frame(lmtest::waldtest(ajuste, restringido, test = "F"))

  tibble::tibble(
    prueba = "Granger por ecuación (Wald F)",
    causa = causa,
    efecto = ecuacion,
    estadistico = as.numeric(tabla[2L, "F"]),
    gl1 = abs(as.numeric(tabla[2L, "Df"])),
    gl2 = as.numeric(tabla[1L, "Res.Df"]),
    p_valor = as.numeric(tabla[2L, "Pr(>F)"])
  )
}

# Estima el VAR con lag.max = 2 y devuelve los contrastes de Granger.
ajustar_var_granger <- function(datos_nacional) {
  serie_var <- as.matrix(datos_nacional[, c("casos_dengue", "PD")])

  if (anyNA(serie_var)) {
    stop("La serie nacional tiene valores faltantes.", call. = FALSE)
  }

  seleccion <- vars::VARselect(
    serie_var,
    lag.max = 2,
    type = "const"
  )

  p_aic <- as.integer(seleccion$selection[["AIC(n)"]])
  if (is.na(p_aic)) {
    p_aic <- 2L
  }

  # p = 0 no permite contrastar Granger; el tope es 2 por grados de
  # libertad (16 trimestres).
  p_var <- max(1L, min(2L, p_aic))

  modelo <- vars::VAR(serie_var, p = p_var, type = "const")

  granger_dengue_pd <- vars::causality(
    modelo, cause = "casos_dengue"
  )$Granger

  granger_pd_dengue <- vars::causality(
    modelo, cause = "PD"
  )$Granger

  tabla <- dplyr::bind_rows(
    extraer_granger(
      granger_dengue_pd,
      causa = "casos_dengue",
      efecto = "PD",
      prueba = "vars::causality (sistema)"
    ),
    probar_granger_ecuacion(modelo, "PD", "casos_dengue"),
    extraer_granger(
      granger_pd_dengue,
      causa = "PD",
      efecto = "casos_dengue",
      prueba = "vars::causality (sistema)"
    )
  )

  list(
    p = p_var,
    modelo = modelo,
    seleccion = seleccion,
    tabla = tabla,
    granger_dengue_pd = granger_dengue_pd
  )
}

ajuste_var <- tryCatch(
  ajustar_var_granger(nacional_q),
  error = function(e) {
    warning(
      "No se pudo estimar el VAR/Granger: ", conditionMessage(e), ". ",
      "Los resultados se guardan con NA y el script continúa.",
      call. = FALSE
    )
    NULL
  }
)

if (is.null(ajuste_var)) {
  p_var <- NA_integer_
  p_granger <- NA_real_
  p_ecuacion <- NA_real_
  efecto_acumulado <- NA_real_
  raiz_maxima <- NA_real_

  var_granger <- tibble::tibble(
    modelo = "VAR(casos_dengue, PD)",
    variables = "casos_dengue, PD",
    rezagos = NA_integer_,
    criterio_rezago = "AIC (lag.max = 2)",
    observaciones = NA_integer_,
    n_trimestres = nrow(nacional_q),
    periodo = periodo_nacional,
    prueba = "vars::causality (sistema)",
    causa = "casos_dengue",
    efecto = "PD",
    estadistico = NA_real_,
    gl1 = NA_real_,
    gl2 = NA_real_,
    p_valor = NA_real_,
    raiz_maxima = NA_real_,
    significativo_5pct = NA
  )
} else {
  p_var <- ajuste_var$p
  p_granger <- ajuste_var$tabla$p_valor[
    ajuste_var$tabla$prueba == "vars::causality (sistema)" &
      ajuste_var$tabla$causa == "casos_dengue" &
      ajuste_var$tabla$efecto == "PD"
  ][1]

  fila_ecuacion <- ajuste_var$tabla$p_valor[
    ajuste_var$tabla$prueba == "Granger por ecuación (Wald F)" &
      ajuste_var$tabla$causa == "casos_dengue" &
      ajuste_var$tabla$efecto == "PD"
  ]
  p_ecuacion <- if (length(fila_ecuacion) > 0L) fila_ecuacion[1] else NA_real_

  # Signo acumulado de los rezagos del dengue en la ecuación de PD.
  coefs_pd <- stats::coef(ajuste_var$modelo$varresult$PD)
  nombres_dengue <- grep("^casos_dengue\\.l", names(coefs_pd), value = TRUE)
  efecto_acumulado <- if (length(nombres_dengue) > 0L) {
    sum(coefs_pd[nombres_dengue])
  } else {
    NA_real_
  }

  raiz_maxima <- tryCatch(
    max(vars::roots(ajuste_var$modelo, modulus = TRUE)),
    error = function(e) NA_real_
  )

  var_granger <- ajuste_var$tabla |>
    dplyr::mutate(
      modelo = "VAR(casos_dengue, PD)",
      variables = "casos_dengue, PD",
      rezagos = p_var,
      criterio_rezago = "AIC (lag.max = 2)",
      observaciones = ajuste_var$modelo$obs,
      n_trimestres = nrow(nacional_q),
      periodo = periodo_nacional,
      raiz_maxima = raiz_maxima,
      significativo_5pct = !is.na(p_valor) & p_valor < 0.05
    ) |>
    dplyr::relocate(
      modelo, variables, rezagos, criterio_rezago, observaciones,
      n_trimestres, periodo, .before = 1
    )
}

# --- 11. Guardado de resultados ------------------------------------

perdida_salida <- cartera_epi |>
  dplyr::select(
    region_ema, trimestre, trimestre_id, anio, trimestre_num,
    EAD, saldo_cde, PD, LGD, EL, tasa_EL,
    casos_dengue, M_dengue, umbral_epidemico, epidemico, escenario
  ) |>
  dplyr::arrange(region_ema, trimestre)

readr::write_csv(perdida_salida, ruta_perdida)
readr::write_csv(var_granger, ruta_var)

rutas_salida <- c(ruta_perdida, ruta_var, ruta_grafica)
for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 12. Resumen en consola ----------------------------------------

cat("\n============================================================\n")
cat(" RIESGO DE CRÉDITO, DENGUE Y PÉRDIDA ESPERADA\n")
cat("============================================================\n")
cat(sprintf(
  "Cartera: %s | %d regiones | %d trimestres | origen: %s\n",
  periodo_nacional,
  dplyr::n_distinct(cartera_epi$region_ema),
  dplyr::n_distinct(cartera_epi$trimestre),
  origen_cartera
))
cat(sprintf(
  "Celdas región-trimestre: %d | umbral epidémico (P75 M_dengue): %.3f\n",
  nrow(cartera_epi), umbral_epidemico
))

cat("\n--- PÉRDIDA ESPERADA (EL = PD * LGD * EAD) ---\n")
cat(sprintf(
  "EL promedio en trimestres normales:   %s miles de millones COP\n",
  formatear_miles_millones(el_normal)
))
cat(sprintf(
  "EL promedio en trimestres epidémicos: %s miles de millones COP\n",
  formatear_miles_millones(el_epi)
))
cat(sprintf(
  "Tasa de EL (PD * LGD) normal: %s | epidémica: %s\n",
  formatear_pct(100 * tasa_normal),
  formatear_pct(100 * tasa_epi)
))

cat("\n>>> AUMENTO PORCENTUAL DE LA EL EN EL ESCENARIO EPIDÉMICO: ")
if (is.na(aumento_pct)) {
  cat("no calculable\n")
} else {
  cat(formatear_pct(aumento_pct), "\n", sep = "")
}

if (!is.na(aumento_tasa_pct)) {
  cat(sprintf(
    "    (La tasa de EL, que descuenta el crecimiento del EAD, cambia %s.)\n",
    formatear_pct(aumento_tasa_pct)
  ))
}

cat("\n--- VAR NACIONAL Y CAUSALIDAD DE GRANGER ---\n")
if (is.null(ajuste_var)) {
  cat("VAR no estimado: no hay p-valor de Granger disponible.\n")
} else {
  cat(sprintf(
    "VAR(casos_dengue, PD): p = %d | %d observaciones efectivas | %s\n",
    p_var, ajuste_var$modelo$obs, periodo_nacional
  ))
  cat(sprintf(
    "Raíz máxima del polinomio característico: %.4f%s\n",
    raiz_maxima,
    if (!is.na(raiz_maxima) && raiz_maxima < 1) {
      " (estable)"
    } else {
      " (revisar estabilidad)"
    }
  ))

  cat("\n>>> P-VALOR DE CAUSALIDAD DE GRANGER (casos_dengue -> PD): ")
  if (is.na(p_granger)) {
    cat("no disponible\n")
  } else {
    cat(sprintf("%.4f\n", p_granger))
  }

  if (!is.na(p_granger)) {
    if (p_granger < 0.05) {
      cat(
        "    Se rechaza H0 al 5%: el dengue causa en sentido de Granger\n",
        "    la PD del sistema financiero.\n"
      )
    } else {
      cat(
        "    No se rechaza H0 al 5%: no hay evidencia suficiente de que\n",
        "    el dengue cause en sentido de Granger la PD.\n"
      )
    }
  }

  if (!is.na(p_ecuacion)) {
    cat(sprintf(
      "    Prueba por ecuación (PD): p = %.4f | efecto acumulado de los\n",
      p_ecuacion
    ))
    cat(sprintf(
      "    rezagos de dengue sobre PD: %.4e\n",
      efecto_acumulado
    ))
  }

  cat(
    "\nNota: con 16 trimestres y p = 2 la prueba tiene baja potencia;\n",
    "el resultado es evidencia descriptiva, no una relación estructural.\n",
    sep = ""
  )
}

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

message("\n✅ Script 22 finalizado.")
