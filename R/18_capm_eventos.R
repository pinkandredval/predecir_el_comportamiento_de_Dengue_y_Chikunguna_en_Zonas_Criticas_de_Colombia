# ============================================================
# 18_capm_eventos.R
# Teoría del Riesgo — CAPM regional y estudio de eventos del miedo
# ============================================================
#
# OBJETIVO
#   Tratar el rendimiento de la ocupación hotelera de cada región EMA
#   como un activo y el rendimiento de las demás regiones como el
#   "mercado" (leave-one-out). A partir del CAPM regional se obtienen
#   el alfa, el beta y los rendimientos anormales (AR) de las 12
#   regiones. Con los AR se
#   construye un estudio de eventos alrededor de los picos de miedo
#   (M_dengue nacional por encima de su percentil 90).
#
# MODELO
#   Rendimiento logarítmico mensual de la ocupación:
#
#       r_{r,t} = log(Y_{r,t}) - log(Y_{r,t-1})
#               = diff(log(Y_{r,t}))
#
#   Rendimiento nacional de referencia (media simple de las 12 regiones):
#
#       r_{n,t} = (1/12) sum_r r_{r,t}
#
#   Mercado CAPM de cada región (leave-one-out):
#
#       r_{m,r,t} = (1/11) sum_{s != r} r_{s,t}
#
#   CAPM regional, una regresión por región:
#
#       r_{r,t} = alfa_r + beta_r r_{m,r,t} + eps_{r,t}
#
#   Rendimiento anormal (AR), residuo del CAPM:
#
#       AR_{r,t} = r_{r,t} - (alfa_r + beta_r r_{m,r,t}) = eps_{r,t}
#
#   Evento: mes e en el que el miedo nacional supera su percentil 90.
#   Ventana de 7 meses: k = -3, -2, ..., 0, ..., +3, con k = 0 en e.
#
#   CAR de la región r en el evento e:
#
#       CAR_{r,e} = sum_{k=-3}^{3} AR_{r,e+k}
#
#   CAR del evento (promedio de las 12 regiones):
#
#       CAR_e = (1/12) sum_r CAR_{r,e}
#
#   Trayectoria del CAR promedio por mes relativo k:
#
#       CAR(k) = (1/E) sum_e sum_{j=-3}^{k} AR_{e+j}
#
#   donde AR_{e+j} es el AR promedio de las 12 regiones en el mes
#   relativo j. En k = +3, CAR(+3) coincide con el promedio de los
#   CAR_e de todos los eventos.
#
# INSUMOS
#   data/D_13/panel_ema_den_chik.rds
#       Panel región-mes EMA (12 regiones × 48 meses, 2022-2025).
#       La columna económica puede llamarse Y u ocupacion; el script
#       la detecta y la homologa a Y.
#
# SALIDAS
#   data/D_18/capm_resultados.csv
#       region_ema, Alfa, Beta, R2, p_valor (mercado leave-one-out)
#   data/D_18/capm_resultados_media_12.csv
#       La misma tabla con el mercado de media simple de las 12 regiones.
#       Es un diagnóstico: por la identidad algebraica, su AR promedio
#       entre regiones es cero en todos los meses.
#   data/D_18/event_study_resultados.csv
#       mes_relativo, ar_promedio, car_promedio, ee_car, n_eventos
#   graficas/G_18/event_study_car.png
#
# NOTAS METODOLÓGICAS
#   1. El filtro 2022-2025 se aplica antes de calcular los rendimientos,
#      tal como pide el guion. Por eso enero de 2022 no tiene rendimiento
#      (no hay diciembre de 2021 en el panel) y cada región aporta 47
#      rendimientos, de febrero de 2022 a diciembre de 2025.
#   2. ELECCIÓN DEL MERCADO (leave-one-out). Si el mercado fuera la
#      media simple de las 12 regiones, cada región se regresaría sobre
#      un promedio que incluye su propio rendimiento. Las ecuaciones
#      normales implican entonces sum_r alfa_r = 0 y sum_r beta_r = 12;
#      por tanto, el AR promedio entre regiones es exactamente 0 en
#      cada mes y el CAR promedio de cualquier evento es 0, sin importar
#      la ventana ni la muestra de estimación. Para que el estudio de
#      eventos tenga contenido, el mercado de cada región excluye su
#      propio rendimiento: media de las otras 11. La media simple se
#      conserva como referencia en r_nacional y su CAPM se guarda en
#      capm_resultados_media_12.csv.
#   3. M_dengue es un índice regional (residuo del script 12). El miedo
#      nacional se obtiene promediando las 12 regiones de cada mes, no
#      sumándolas: sumar cambiaría la escala y la interpretación.
#   4. El percentil 90 se calcula sobre los 48 meses nacionales de
#      2022-2025. "Superó" se implementa con > (estrictamente mayor).
#   5. Los eventos pueden ser meses contiguos y, por tanto, sus ventanas
#      de 7 meses pueden solaparse. El guion define cada mes por encima
#      del p90 como un evento, así que se conservan todos y se advierte
#      en consola; no se fuerza que las ventanas sean disjuntas.
#   6. Los AR son residuos del CAPM estimado con toda la muestra
#      2022-2025 (event study "in-sample"). No se reestima el CAPM
#      excluyendo la ventana de cada evento.
#   7. El p-valor de Beta contrasta H0: beta_r = 0 con el estadístico t
#      habitual de la regresión. El alfa se reporta con su error estándar
#      en la tabla interna, pero el CSV final conserva las cinco columnas
#      pedidas: Región, Alfa, Beta, R² y p-valor.
#
# USO DE ::
#   Todas las funciones no base se llaman con paquete::funcion() para
#   evitar el enmascaramiento entre dplyr, tidyr, stats y ggplot2.
#   No se adjuntan paquetes con library().
# ============================================================

# --- 0. Paquetes, rutas y carpetas ---------------------------------

paquetes <- c(
  "dplyr", "tidyr", "ggplot2", "lubridate", "readr", "tibble", "here"
)
faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]
if (length(faltantes) > 0L) {
  stop(
    "Faltan paquetes: ", paste(faltantes, collapse = ", "), ". ",
    "Instálelos con renv::install(c(",
    paste0("\"", faltantes, "\"", collapse = ", "), ")).",
    call. = FALSE
  )
}

dir.create(here::here("data", "D_18"), recursive = TRUE, showWarnings = FALSE)
dir.create(
  here::here("graficas", "G_18"),
  recursive = TRUE, showWarnings = FALSE
)

ruta_panel        <- here::here("data", "D_13", "panel_ema_den_chik.rds")
ruta_capm         <- here::here("data", "D_18", "capm_resultados.csv")
ruta_capm_media12 <- here::here(
  "data", "D_18", "capm_resultados_media_12.csv"
)
ruta_eventos      <- here::here("data", "D_18", "event_study_resultados.csv")
ruta_png          <- here::here("graficas", "G_18", "event_study_car.png")

desde_anio      <- 2022L
hasta_anio      <- 2025L
ventana_evento  <- -3:3
prob_percentil  <- 0.90

# --- 1. Utilidades -------------------------------------------------

# Fecha desplazada k meses. Las fechas del panel son siempre el día 1,
# así que seq() evita los problemas de fin de mes de la aritmética de
# meses. k = 0 devuelve la misma fecha.
sumar_meses <- function(fecha, k) {
  if (k == 0L) {
    return(fecha)
  }
  seq(fecha, by = paste(k, "months"), length.out = 2L)[2L]
}

# p-valor legible para consola.
formatear_p <- function(p) {
  if (length(p) == 0L || is.na(p)) {
    return("NA")
  }
  if (p < 0.0001) {
    return("<0.0001")
  }
  sprintf("%.4f", p)
}

# Detecta la columna de economía: Y (prioritaria) u ocupacion/ocupación,
# sin distinguir mayúsculas ni tildes.
detectar_columna_economia <- function(datos) {
  nombres <- names(datos)

  if ("Y" %in% nombres) {
    return("Y")
  }

  normalizados <- tolower(iconv(nombres, to = "ASCII//TRANSLIT"))
  candidatas <- which(
    normalizados == "y" | startsWith(normalizados, "ocupacion")
  )

  if (length(candidatas) == 0L) {
    return(NA_character_)
  }
  nombres[candidatas[1L]]
}

# --- 2. Lectura y validación del panel -----------------------------

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

panel <- readRDS(ruta_panel)

if (!inherits(panel, "data.frame")) {
  stop("El objeto leído no es un data.frame: ", ruta_panel, ".", call. = FALSE)
}

requeridas <- c("region_ema", "mes", "M_dengue")
faltan_columnas <- setdiff(requeridas, names(panel))
if (length(faltan_columnas) > 0L) {
  stop(
    "Faltan columnas en el panel: ",
    paste(faltan_columnas, collapse = ", "), ". Columnas disponibles: ",
    paste(names(panel), collapse = ", "), ".",
    call. = FALSE
  )
}

col_economia <- detectar_columna_economia(panel)
if (is.na(col_economia)) {
  stop(
    "No se encontró la columna de economía. Se buscó 'Y' y ",
    "'ocupacion'/'ocupación'. Columnas disponibles: ",
    paste(names(panel), collapse = ", "), ".",
    call. = FALSE
  )
}

panel <- panel |>
  dplyr::mutate(
    region_ema = as.character(.data[["region_ema"]]),
    mes = as.Date(.data[["mes"]]),
    Y = as.numeric(.data[[col_economia]]),
    M_dengue = as.numeric(.data[["M_dengue"]])
  ) |>
  dplyr::filter(
    mes >= as.Date(sprintf("%d-01-01", desde_anio)),
    mes <= as.Date(sprintf("%d-12-01", hasta_anio))
  ) |>
  dplyr::select(region_ema, mes, Y, M_dengue) |>
  dplyr::arrange(region_ema, mes)

if (nrow(panel) == 0L) {
  stop("El panel quedó vacío después del filtro 2022-2025.", call. = FALSE)
}
if (anyNA(panel$mes)) {
  stop("Hay fechas faltantes (NA) en la columna mes.", call. = FALSE)
}
if (anyNA(panel$Y)) {
  stop(
    "La ocupación (Y) tiene ", sum(is.na(panel$Y)),
    " valores faltantes en 2022-2025; los rendimientos no se pueden ",
    "calcular sin una serie completa.",
    call. = FALSE
  )
}
if (anyNA(panel$M_dengue)) {
  stop(
    "M_dengue tiene ", sum(is.na(panel$M_dengue)),
    " valores faltantes en 2022-2025.",
    call. = FALSE
  )
}
if (any(panel$Y <= 0)) {
  stop(
    "Y debe ser estrictamente positiva para calcular log(Y); hay ",
    sum(panel$Y <= 0), " observaciones no positivas.",
    call. = FALSE
  )
}

regiones <- sort(unique(panel$region_ema))
if (length(regiones) != 12L) {
  stop(
    "El panel debe tener 12 regiones EMA y tiene ", length(regiones), ": ",
    paste(regiones, collapse = ", "), ".",
    call. = FALSE
  )
}

if (anyDuplicated(panel[, c("region_ema", "mes")])) {
  stop("Hay pares región-mes duplicados en el panel.", call. = FALSE)
}

meses <- sort(unique(panel$mes))
meses_esperados <- seq(
  as.Date(sprintf("%d-01-01", desde_anio)),
  as.Date(sprintf("%d-12-01", hasta_anio)),
  by = "1 month"
)

if (length(meses) != length(meses_esperados) ||
    any(meses != meses_esperados)) {
  stop(
    "El panel no cubre la rejilla mensual completa 2022-2025. ",
    "Meses encontrados: ", length(meses), "; esperados: ",
    length(meses_esperados), ".",
    call. = FALSE
  )
}

conteo_panel <- panel |>
  dplyr::count(region_ema, name = "n_meses")

if (any(conteo_panel$n_meses != length(meses_esperados))) {
  stop(
    "No todas las regiones tienen los 48 meses: ",
    paste(
      conteo_panel$region_ema[conteo_panel$n_meses != length(meses_esperados)],
      collapse = ", "
    ), ".",
    call. = FALSE
  )
}

# --- 3. Rendimientos regionales y rendimiento nacional -------------

# diff(log(Y)) dentro de cada región. La primera observación de cada
# región genera NA porque no tiene mes anterior; se elimina con
# filter(!is.na(r_region)). No debe quedar ningún otro NA.
rendimientos <- panel |>
  dplyr::arrange(region_ema, mes) |>
  dplyr::group_by(region_ema) |>
  dplyr::mutate(
    r_region = c(NA_real_, diff(log(Y)))
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(r_region))

if (anyNA(rendimientos$r_region) ||
    any(!is.finite(rendimientos$r_region))) {
  stop(
    "Quedaron rendimientos no finitos después de eliminar el primer ",
    "mes de cada región.",
    call. = FALSE
  )
}

n_ret_esperados <- length(meses_esperados) - 1L
if (nrow(rendimientos) != 12L * n_ret_esperados) {
  stop(
    "Se esperaban ", 12L * n_ret_esperados, " rendimientos región-mes y ",
    "hay ", nrow(rendimientos), ".",
    call. = FALSE
  )
}

conteo_ret <- rendimientos |>
  dplyr::count(region_ema, name = "n_ret")

if (any(conteo_ret$n_ret != n_ret_esperados)) {
  stop(
    "No todas las regiones tienen ", n_ret_esperados,
    " rendimientos mensuales.",
    call. = FALSE
  )
}

primer_mes <- rendimientos |>
  dplyr::group_by(region_ema) |>
  dplyr::summarise(mes_min = min(mes), .groups = "drop")

if (any(primer_mes$mes_min != meses_esperados[2L])) {
  stop(
    "El primer rendimiento no es el de ", format(meses_esperados[2L], "%Y-%m"),
    " en todas las regiones; revise los huecos de la serie.",
    call. = FALSE
  )
}

# Rendimiento nacional de referencia: promedio simple de las 12 regiones.
# Se conserva para la tabla de diagnóstico y para documentar la identidad
# algebraica que obliga a usar el mercado leave-one-out en el CAPM.
mercado <- rendimientos |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    r_nacional = mean(r_region),
    suma_r = sum(r_region),
    n_regiones = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(mes)

if (any(mercado$n_regiones != 12L)) {
  stop(
    "El rendimiento nacional no se calculó con 12 regiones en todos los ",
    "meses; revise la rejilla del panel.",
    call. = FALSE
  )
}
if (anyNA(mercado$r_nacional) || any(!is.finite(mercado$r_nacional))) {
  stop("El rendimiento nacional tiene valores no finitos.", call. = FALSE)
}

datos_capm <- rendimientos |>
  dplyr::left_join(
    mercado,
    by = "mes",
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(
    # Mercado de la región r: media de las otras 11 regiones. Es el
    # leave-one-out que evita que la región aparezca en su propio
    # regresor y que el AR promedio sea cero por construcción.
    r_mercado = (suma_r - r_region) / (n_regiones - 1L)
  ) |>
  dplyr::select(region_ema, mes, r_region, r_nacional, r_mercado) |>
  dplyr::arrange(region_ema, mes)

if (nrow(datos_capm) != nrow(rendimientos)) {
  stop(
    "La unión con el rendimiento nacional cambió el número de filas: ",
    nrow(datos_capm), " frente a ", nrow(rendimientos), ".",
    call. = FALSE
  )
}
if (anyNA(datos_capm)) {
  stop(
    "Hay valores faltantes en datos_capm después de unir el mercado.",
    call. = FALSE
  )
}
if (any(!is.finite(datos_capm$r_mercado))) {
  stop("El mercado leave-one-out tiene valores no finitos.", call. = FALSE)
}

# Comprobación: r_mercado = (12 * r_nacional - r_region) / 11.
r_mercado_implicito <- (
  12 * datos_capm$r_nacional - datos_capm$r_region
) / 11

if (any(abs(datos_capm$r_mercado - r_mercado_implicito) > 1e-10)) {
  stop(
    "El mercado leave-one-out no coincide con la media de las otras 11 ",
    "regiones.",
    call. = FALSE
  )
}

# --- 4. CAPM regional: alfa, beta, R², p-valor y AR ----------------
#
# El CAPM principal usa r_mercado (media de las otras 11 regiones). La
# regresión con r_nacional (media simple de las 12) se calcula en 4.1
# solo como diagnóstico de la identidad algebraica.

modelos <- vector("list", length(regiones))
names(modelos) <- regiones
tabla_capm <- vector("list", length(regiones))
residuos_lista <- vector("list", length(regiones))

for (i in seq_along(regiones)) {
  region <- regiones[i]
  datos_region <- dplyr::filter(datos_capm, region_ema == region)

  modelo <- stats::lm(r_region ~ r_mercado, data = datos_region)
  resumen <- summary(modelo)
  coefs <- resumen$coefficients

  if (nrow(coefs) < 2L) {
    stop(
      "La regresión de ", region, " no tiene coeficiente para r_mercado. ",
      "Revise si el mercado tiene variación.",
      call. = FALSE
    )
  }
  if (length(stats::residuals(modelo)) != nrow(datos_region)) {
    stop(
      "Los residuos de ", region, " no coinciden con el número de ",
      "observaciones (¿hay NA ocultos?).",
      call. = FALSE
    )
  }

  alfa <- unname(coefs[1L, 1L])
  beta <- unname(coefs[2L, 1L])
  r2 <- resumen$r.squared
  p_beta <- unname(coefs[2L, 4L])

  if (!is.finite(alfa) || !is.finite(beta) || !is.finite(r2) ||
      !is.finite(p_beta)) {
    stop(
      "La regresión de ", region, " produjo alfa, beta, R2 o p-valor ",
      "no finitos.",
      call. = FALSE
    )
  }

  tabla_capm[[i]] <- tibble::tibble(
    region_ema = region,
    Alfa = alfa,
    Beta = beta,
    R2 = r2,
    p_valor = p_beta,
    n_obs = nrow(datos_region),
    ee_alfa = unname(coefs[1L, 2L]),
    ee_beta = unname(coefs[2L, 2L])
  )

  residuos_lista[[i]] <- tibble::tibble(
    region_ema = region,
    mes = datos_region$mes,
    ar = as.numeric(stats::residuals(modelo))
  )

  modelos[[i]] <- modelo
}

tabla_capm <- dplyr::bind_rows(tabla_capm)
residuos <- dplyr::bind_rows(residuos_lista)

if (nrow(tabla_capm) != 12L) {
  stop("La tabla del CAPM debe tener 12 filas.", call. = FALSE)
}
if (anyNA(tabla_capm[, c("Alfa", "Beta", "R2", "p_valor")])) {
  stop("La tabla del CAPM tiene valores faltantes.", call. = FALSE)
}
if (any(!is.finite(
  as.matrix(tabla_capm[, c("Alfa", "Beta", "R2", "p_valor")])
))) {
  stop("La tabla del CAPM tiene valores no finitos.", call. = FALSE)
}
if (nrow(residuos) != nrow(datos_capm)) {
  stop(
    "Se esperaban ", nrow(datos_capm), " residuos y hay ",
    nrow(residuos), ".",
    call. = FALSE
  )
}
if (anyNA(residuos$ar) || any(!is.finite(residuos$ar))) {
  stop("Los residuos del CAPM tienen valores no finitos.", call. = FALSE)
}

# El residuo de una regresión con intercepto suma cero en cada región.
medias_ar <- residuos |>
  dplyr::group_by(region_ema) |>
  dplyr::summarise(media_ar = mean(ar), .groups = "drop")

if (any(abs(medias_ar$media_ar) > 1e-8)) {
  warning(
    "La media de los residuos no es cero en alguna región; revise el ",
    "ajuste del CAPM.",
    call. = FALSE
  )
}

# Tabla final del CAPM: exactamente las cinco columnas pedidas.
capm_salida <- dplyr::select(tabla_capm, region_ema, Alfa, Beta, R2, p_valor)
readr::write_csv(capm_salida, ruta_capm)

# --- 4.1 Diagnóstico: CAPM con la media simple de las 12 regiones ---
#
# Esta es la regresión literal del enunciado. Se guarda para comparar,
# pero no se usa en el estudio de eventos: por la identidad algebraica
# su AR promedio entre regiones es cero en todos los meses.
#
#   sum_r alfa_r = 0,  sum_r beta_r = 12  =>  sum_r AR_{r,t} = 0
#
# porque r_nacional_t es exactamente el promedio de los r_region_t y
# cada región se regresa sobre ese promedio incluyendo su propio dato.

tabla_media12 <- vector("list", length(regiones))
residuos_media12_lista <- vector("list", length(regiones))

for (i in seq_along(regiones)) {
  region <- regiones[i]
  datos_region <- dplyr::filter(datos_capm, region_ema == region)
  modelo <- stats::lm(r_region ~ r_nacional, data = datos_region)
  resumen <- summary(modelo)
  coefs <- resumen$coefficients

  tabla_media12[[i]] <- tibble::tibble(
    region_ema = region,
    Alfa = unname(coefs[1L, 1L]),
    Beta = unname(coefs[2L, 1L]),
    R2 = resumen$r.squared,
    p_valor = unname(coefs[2L, 4L])
  )

  residuos_media12_lista[[i]] <- tibble::tibble(
    region_ema = region,
    mes = datos_region$mes,
    ar = as.numeric(stats::residuals(modelo))
  )
}

tabla_media12 <- dplyr::bind_rows(tabla_media12)
residuos_media12 <- dplyr::bind_rows(residuos_media12_lista)

if (nrow(tabla_media12) != 12L ||
    anyNA(tabla_media12[, c("Alfa", "Beta", "R2", "p_valor")]) ||
    any(!is.finite(
      as.matrix(tabla_media12[, c("Alfa", "Beta", "R2", "p_valor")])
    ))) {
  stop(
    "La tabla del CAPM con media simple de 12 regiones no es válida.",
    call. = FALSE
  )
}
if (nrow(residuos_media12) != nrow(datos_capm)) {
  stop(
    "Los residuos del CAPM con media simple no cubren todas las ",
    "observaciones.",
    call. = FALSE
  )
}

suma_alfa_media12 <- sum(tabla_media12$Alfa)
suma_beta_media12 <- sum(tabla_media12$Beta)
media_ar_mensual_media12 <- residuos_media12 |>
  dplyr::group_by(mes) |>
  dplyr::summarise(media_ar = mean(ar), .groups = "drop")
max_abs_media_ar_media12 <- max(abs(media_ar_mensual_media12$media_ar))

if (abs(suma_alfa_media12) > 1e-8 ||
    abs(suma_beta_media12 - 12) > 1e-8 ||
    max_abs_media_ar_media12 > 1e-8) {
  warning(
    "La identidad del mercado de media simple no se cumple ",
    "numéricamente: sum(alfa) = ", suma_alfa_media12,
    ", sum(beta) = ", suma_beta_media12,
    ", máximo |AR promedio mensual| = ", max_abs_media_ar_media12, ".",
    call. = FALSE
  )
}

capm_media12_salida <- dplyr::select(
  tabla_media12, region_ema, Alfa, Beta, R2, p_valor
)
readr::write_csv(capm_media12_salida, ruta_capm_media12)

# --- 5. Miedo nacional y definición de eventos ---------------------

miedo_nacional <- panel |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    M_dengue = mean(M_dengue),
    n_regiones = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(mes)

if (nrow(miedo_nacional) != length(meses_esperados)) {
  stop(
    "El miedo nacional no cubre los 48 meses: tiene ",
    nrow(miedo_nacional), ".",
    call. = FALSE
  )
}
if (any(miedo_nacional$n_regiones != 12L)) {
  stop(
    "El miedo nacional no se promedió con 12 regiones en todos los meses.",
    call. = FALSE
  )
}
if (anyNA(miedo_nacional$M_dengue)) {
  stop("El miedo nacional tiene valores faltantes.", call. = FALSE)
}

umbral_p90 <- stats::quantile(
  miedo_nacional$M_dengue,
  probs = prob_percentil,
  na.rm = TRUE,
  names = FALSE
)

if (!is.finite(umbral_p90)) {
  stop("No se pudo calcular el percentil 90 del miedo nacional.",
       call. = FALSE)
}

eventos <- miedo_nacional$mes[miedo_nacional$M_dengue > umbral_p90]
eventos <- sort(unique(eventos))

if (length(eventos) == 0L) {
  stop(
    "Ningún mes de 2022-2025 supera el percentil 90 del miedo nacional (",
    sprintf("%.4f", umbral_p90), ").",
    call. = FALSE
  )
}

# Rejilla completa evento × mes relativo × región. Incluir la región
# desde el inicio permite unir por región-mes (many-to-one) y evita
# cualquier ambigüedad de many-to-many. rowwise() permite usar
# sumar_meses(), que trabaja con una fecha y un escalar por fila.
ventanas <- tidyr::expand_grid(
  evento = eventos,
  mes_relativo = ventana_evento,
  region_ema = regiones
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    mes_objetivo = sumar_meses(evento, mes_relativo)
  ) |>
  dplyr::ungroup()

if (!inherits(ventanas$mes_objetivo, "Date")) {
  ventanas$mes_objetivo <- as.Date(
    ventanas$mes_objetivo,
    origin = "1970-01-01"
  )
}

# --- 6. AR en las ventanas de evento -------------------------------
#
# Se usa left_join() y se descartan explícitamente los meses sin AR.
# residuos debe tener una sola fila por región-mes; se valida antes del
# cruce para que no haya productos cartesianos. Las fechas se homologan
# a Date para que el join no falle por horas o zonas horarias.

residuos_join <- residuos |>
  dplyr::mutate(mes = as.Date(mes)) |>
  dplyr::select(region_ema, mes, ar)

if (anyDuplicated(residuos_join[, c("region_ema", "mes")])) {
  duplicados <- residuos_join |>
    dplyr::count(region_ema, mes, name = "n") |>
    dplyr::filter(n > 1L)

  stop(
    "Hay regiones-mes duplicadas en los residuos; el cruce produciría ",
    "filas repetidas. Ejemplos: ",
    paste(
      sprintf(
        "%s/%s (%d filas)",
        duplicados$region_ema,
        format(duplicados$mes, "%Y-%m"),
        duplicados$n
      ),
      collapse = "; "
    ), ".",
    call. = FALSE
  )
}

ventanas_join <- ventanas |>
  dplyr::mutate(mes_objetivo = as.Date(mes_objetivo))

ar_ventanas <- ventanas_join |>
  dplyr::left_join(
    residuos_join,
    by = c("region_ema" = "region_ema", "mes_objetivo" = "mes")
  ) |>
  dplyr::filter(!is.na(ar))

if (nrow(ar_ventanas) == 0L) {
  stop(
    "Ninguna ventana de evento coincide con los meses con rendimiento. ",
    "Revise las fechas de los eventos y del panel.",
    call. = FALSE
  )
}

if (anyDuplicated(ar_ventanas[, c("evento", "region_ema", "mes_relativo")])) {
  stop(
    "El cruce generó más de una observación por evento, región y mes ",
    "relativo. Revise que residuos tenga una sola fila por región-mes.",
    call. = FALSE
  )
}

# Meses del panel con rendimiento: sirven para separar los eventos que
# caen en el borde de la muestra (ventana incompleta) de los huecos
# internos de datos.
meses_panel <- sort(unique(residuos_join$mes))

ventanas_temporales <- ventanas_join |>
  dplyr::distinct(evento, mes_relativo, mes_objetivo) |>
  dplyr::mutate(
    mes_disponible = mes_objetivo %in% meses_panel
  ) |>
  dplyr::group_by(evento) |>
  dplyr::summarise(
    n_meses_disponibles = sum(mes_disponible),
    .groups = "drop"
  )

eventos_validos <- ventanas_temporales$evento[
  ventanas_temporales$n_meses_disponibles == length(ventana_evento)
]

if (length(eventos_validos) == 0L) {
  stop(
    "Ningún evento tiene los 7 meses de la ventana dentro del panel ",
    "2022-2025. Los eventos están demasiado cerca del borde de la muestra.",
    call. = FALSE
  )
}

eventos_borde <- setdiff(eventos, eventos_validos)
if (length(eventos_borde) > 0L) {
  message(
    "Se descartan ", length(eventos_borde),
    " evento(s) con ventana incompleta por el borde de la muestra: ",
    paste(format(eventos_borde, "%Y-%m"), collapse = ", "), "."
  )
}

ar_ventanas <- dplyr::filter(ar_ventanas, evento %in% eventos_validos)

# Verificación estricta: cada evento × región debe tener exactamente 7
# filas (una por mes relativo). Si falta alguna, se informa el evento,
# la región y los meses relativos faltantes.
conteo_region <- ar_ventanas |>
  dplyr::count(evento, region_ema, name = "n_meses")

grid_region <- tidyr::expand_grid(
  evento = eventos_validos,
  region_ema = regiones
)

diagnostico_region <- grid_region |>
  dplyr::left_join(conteo_region, by = c("evento", "region_ema")) |>
  dplyr::mutate(n_meses = dplyr::coalesce(n_meses, 0L)) |>
  dplyr::filter(n_meses != length(ventana_evento))

if (nrow(diagnostico_region) > 0L) {
  detalle <- vapply(
    seq_len(nrow(diagnostico_region)),
    function(i) {
      evento_i <- diagnostico_region$evento[i]
      region_i <- diagnostico_region$region_ema[i]

      presentes <- ar_ventanas$mes_relativo[
        ar_ventanas$evento == evento_i &
          ar_ventanas$region_ema == region_i
      ]
      faltantes <- setdiff(ventana_evento, presentes)

      sprintf(
        "%s / %s: %d de %d meses (faltan %s)",
        format(evento_i, "%Y-%m"),
        region_i,
        diagnostico_region$n_meses[i],
        length(ventana_evento),
        paste(faltantes, collapse = ", ")
      )
    },
    character(1)
  )

  stop(
    "Ventanas incompletas por región. Cada evento debe tener 7 meses ",
    "por región:\n  ", paste(detalle, collapse = "\n  "),
    call. = FALSE
  )
}

if (any(conteo_region$n_meses != length(ventana_evento))) {
  stop(
    "Alguna región no tiene los 7 meses de la ventana en algún evento.",
    call. = FALSE
  )
}

# Cada evento × mes relativo debe tener las 12 regiones.
conteo_mes <- ar_ventanas |>
  dplyr::count(evento, mes_relativo, name = "n_regiones")

if (any(conteo_mes$n_regiones != 12L)) {
  stop(
    "Algún evento no tiene las 12 regiones en todos los meses relativos.",
    call. = FALSE
  )
}

# --- 7. CAR por evento y trayectoria promedio ----------------------

# AR promedio de las 12 regiones en cada evento y mes relativo.
aar_evento <- ar_ventanas |>
  dplyr::group_by(evento, mes_relativo) |>
  dplyr::summarise(
    ar_promedio = mean(ar),
    n_regiones = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(evento, mes_relativo)

# CAR de cada región en cada evento: suma de sus 7 AR.
car_region_evento <- ar_ventanas |>
  dplyr::group_by(evento, region_ema) |>
  dplyr::summarise(
    n_meses = dplyr::n(),
    car_region = sum(ar),
    .groups = "drop"
  )

if (any(car_region_evento$n_meses != length(ventana_evento))) {
  malas <- dplyr::filter(
    car_region_evento,
    n_meses != length(ventana_evento)
  )

  stop(
    "CAR regional incompleto: ",
    paste(
      sprintf(
        "%s / %s (%d de %d meses)",
        format(malas$evento, "%Y-%m"),
        malas$region_ema,
        malas$n_meses,
        length(ventana_evento)
      ),
      collapse = "; "
    ), ".",
    call. = FALSE
  )
}

# CAR del evento: promedio simple de los CAR de las 12 regiones.
car_por_evento <- car_region_evento |>
  dplyr::group_by(evento) |>
  dplyr::summarise(
    n_regiones = dplyr::n(),
    car_evento = mean(car_region),
    .groups = "drop"
  ) |>
  dplyr::arrange(evento)

if (any(car_por_evento$n_regiones != 12L)) {
  stop(
    "Algún evento no tiene el CAR de las 12 regiones.",
    call. = FALSE
  )
}

car_promedio_evento <- mean(car_por_evento$car_evento)

# Trayectoria: CAR acumulado desde -3 hasta cada mes relativo.
car_evento <- aar_evento |>
  dplyr::group_by(evento) |>
  dplyr::mutate(
    car_acumulado = cumsum(ar_promedio)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(evento, mes_relativo)

conteo_car <- car_evento |>
  dplyr::count(evento, name = "n_meses")

if (any(conteo_car$n_meses != length(ventana_evento))) {
  stop(
    "La trayectoria del CAR no tiene 7 meses en todos los eventos.",
    call. = FALSE
  )
}

car_promedio <- car_evento |>
  dplyr::group_by(mes_relativo) |>
  dplyr::summarise(
    ar_promedio = mean(ar_promedio),
    car_promedio = mean(car_acumulado),
    ee_car = if (dplyr::n() > 1L) {
      stats::sd(car_acumulado) / sqrt(dplyr::n())
    } else {
      NA_real_
    },
    n_eventos = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(mes_relativo)

if (nrow(car_promedio) != length(ventana_evento)) {
  stop(
    "La tabla del CAR debe tener 7 filas (meses -3 a +3) y tiene ",
    nrow(car_promedio), ".",
    call. = FALSE
  )
}
if (!all(car_promedio$mes_relativo == ventana_evento)) {
  stop("Los meses relativos del CAR no son -3 a +3.", call. = FALSE)
}
if (any(car_promedio$n_eventos != length(eventos_validos))) {
  stop(
    "El número de eventos no es constante en la trayectoria del CAR.",
    call. = FALSE
  )
}

# Coherencia: el CAR acumulado promedio debe ser la suma acumulada del
# AR promedio, y en +3 debe coincidir con el promedio de los CAR_e.
car_implicito <- cumsum(car_promedio$ar_promedio)
if (any(abs(car_implicito - car_promedio$car_promedio) > 1e-10)) {
  stop(
    "El CAR acumulado promedio no coincide con la suma acumulada del ",
    "AR promedio.",
    call. = FALSE
  )
}

car_final_acumulado <- car_evento |>
  dplyr::filter(mes_relativo == max(ventana_evento)) |>
  dplyr::select(evento, car_acumulado)

comparacion_car <- dplyr::inner_join(
  car_por_evento,
  car_final_acumulado,
  by = "evento"
)

if (any(abs(comparacion_car$car_evento -
            comparacion_car$car_acumulado) > 1e-10)) {
  stop(
    "El CAR del evento no coincide con el CAR acumulado en +3; revise ",
    "el cálculo de la ventana.",
    call. = FALSE
  )
}

car_en_mas_3 <- car_promedio$car_promedio[
  car_promedio$mes_relativo == max(ventana_evento)
]

if (abs(car_en_mas_3 - car_promedio_evento) > 1e-10) {
  stop(
    "El CAR promedio en +3 no coincide con el promedio de los CAR por ",
    "evento.",
    call. = FALSE
  )
}

# --- 8. Guardado de la tabla del CAR -------------------------------

eventos_salida <- dplyr::select(
  car_promedio,
  mes_relativo, ar_promedio, car_promedio, ee_car, n_eventos
)

readr::write_csv(eventos_salida, ruta_eventos)

# --- 9. Gráfico del CAR promedio -----------------------------------

subtitulo <- sprintf(
  paste0(
    "Ventana de 7 meses (-3 a +3) | %d eventos | mercado leave-one-out | ",
    "M_dengue nacional > p90 (%.3f)"
  ),
  length(eventos_validos), umbral_p90
)

caption <- sprintf(
  paste0(
    "CAR: suma de los AR (residuos del CAPM con mercado leave-one-out) ",
    "desde -3 hasta cada mes relativo,\n",
    "promediada entre las 12 regiones y los %d eventos. En +3 equivale ",
    "al CAR promedio de la ventana completa.\n",
    "Eventos: meses con miedo por encima del percentil 90; ventanas de ",
    "eventos contiguos pueden solaparse."
  ),
  length(eventos_validos)
)

p_car <- ggplot2::ggplot(
  car_promedio,
  ggplot2::aes(x = mes_relativo, y = car_promedio)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey55",
    linewidth = 0.3,
    linetype = "dotted"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = "#B2182B",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(colour = "#2166AC", linewidth = 0.9) +
  ggplot2::geom_point(colour = "#2166AC", size = 2.4) +
  ggplot2::scale_x_continuous(
    breaks = ventana_evento,
    limits = c(min(ventana_evento) - 0.4, max(ventana_evento) + 0.4)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.3f", x),
    expand = ggplot2::expansion(mult = c(0.10, 0.14))
  ) +
  ggplot2::labs(
    title = "Estudio de eventos: rendimiento anormal acumulado del turismo",
    subtitle = subtitulo,
    x = "Mes relativo al pico de miedo",
    y = "CAR promedio (rendimiento anormal acumulado)",
    caption = caption
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 9),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7.5, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_png,
  plot = p_car,
  width = 9.5,
  height = 5.4,
  dpi = 200,
  bg = "white"
)

# --- 10. Validación de salidas -------------------------------------

rutas_salida <- c(ruta_capm, ruta_capm_media12, ruta_eventos, ruta_png)
for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 11. Resumen en consola ----------------------------------------

cat("\n============================================================\n")
cat(" 18 — CAPM REGIONAL Y ESTUDIO DE EVENTOS DEL MIEDO\n")
cat("============================================================\n")
cat(sprintf(
  "Panel: %d regiones-mes | %d regiones | %s a %s\n",
  nrow(panel), length(regiones),
  format(min(panel$mes), "%Y-%m"), format(max(panel$mes), "%Y-%m")
))
cat(sprintf(
  "Columna económica detectada: %s (homologada a Y)\n",
  col_economia
))
cat(sprintf(
  "Rendimientos: %d regiones-mes | %d meses por región\n",
  nrow(datos_capm), as.integer(nrow(datos_capm) / 12L)
))
cat(
  "Mercado CAPM: leave-one-out (media de las otras 11 regiones) | ",
  "M_dengue nacional: promedio simple de las 12 regiones\n",
  sep = ""
)
cat(
  "Diagnóstico con la media simple de 12: ",
  sprintf("sum(alfa) = %.2e", suma_alfa_media12),
  sprintf(" | sum(beta) = %.6f", suma_beta_media12),
  sprintf(
    " | máx |AR medio mensual| = %.2e\n",
    max_abs_media_ar_media12
  ),
  sep = ""
)

cat(
  "\n--- CAPM REGIONAL (leave-one-out): ",
  "r_region = Alfa + Beta * r_mercado + AR ---\n",
  sep = ""
)
cat(sprintf(
  "%-34s %10s %10s %8s %10s %5s\n",
  "Región", "Alfa", "Beta", "R2", "p(Beta)", "n"
))

for (i in seq_len(nrow(tabla_capm))) {
  cat(sprintf(
    "%-34s %10.4f %10.4f %8.4f %10s %5d\n",
    tabla_capm$region_ema[i],
    tabla_capm$Alfa[i],
    tabla_capm$Beta[i],
    tabla_capm$R2[i],
    formatear_p(tabla_capm$p_valor[i]),
    tabla_capm$n_obs[i]
  ))
}

cat("\n--- BETAS DE LAS 12 REGIONES (de mayor a menor) ---\n")
orden_beta <- order(tabla_capm$Beta, decreasing = TRUE)

for (i in orden_beta) {
  cat(sprintf(
    "  %-34s %8.4f\n",
    tabla_capm$region_ema[i],
    tabla_capm$Beta[i]
  ))
}

cat("\n--- ESTUDIO DE EVENTOS ---\n")
cat("AR: residuos del CAPM con mercado leave-one-out.\n")
cat(sprintf(
  "Umbral p90 de M_dengue nacional: %.4f\n",
  umbral_p90
))
cat(sprintf(
  "Meses con M_dengue > p90 (%d): %s\n",
  length(eventos),
  paste(format(eventos, "%Y-%m"), collapse = ", ")
))
cat(sprintf(
  "Ventana: %d a %d (%d meses) | eventos con ventana completa: %d\n",
  min(ventana_evento), max(ventana_evento), length(ventana_evento),
  length(eventos_validos)
))

cat("\nCAR por evento (suma de los AR leave-one-out de -3 a +3):\n")
for (i in seq_len(nrow(car_por_evento))) {
  cat(sprintf(
    "  %s  %10.4f\n",
    format(car_por_evento$evento[i], "%Y-%m"),
    car_por_evento$car_evento[i]
  ))
}

cat("\nCAR promedio por mes relativo:\n")
for (i in seq_len(nrow(car_promedio))) {
  cat(sprintf(
    "  %+d  AR = %8.4f  CAR = %8.4f  ee = %8s  (n = %d)\n",
    car_promedio$mes_relativo[i],
    car_promedio$ar_promedio[i],
    car_promedio$car_promedio[i],
    if (is.na(car_promedio$ee_car[i])) {
      "NA"
    } else {
      sprintf("%.4f", car_promedio$ee_car[i])
    },
    car_promedio$n_eventos[i]
  ))
}

cat(sprintf(
  "\nCAR acumulado promedio en la ventana del evento (-3 a +3): %.4f\n",
  car_promedio_evento
))
cat("  (equivale al CAR promedio en el mes relativo +3)\n")

if (car_promedio_evento < 0) {
  cat(
    "Lectura: en promedio, el turismo acumula un rendimiento anormal ",
    "negativo alrededor de los picos de miedo. Es una asociación ",
    "descriptiva, no un efecto causal.\n",
    sep = ""
  )
} else {
  cat(
    "Lectura: en promedio, el turismo acumula un rendimiento anormal ",
    "positivo alrededor de los picos de miedo. Es una asociación ",
    "descriptiva, no un efecto causal.\n",
    sep = ""
  )
}

indice_mes_eventos <- lubridate::year(eventos) * 12L +
  lubridate::month(eventos)

if (any(diff(indice_mes_eventos) == 1L)) {
  cat(
    "\nNota: hay eventos en meses contiguos; sus ventanas de 7 meses se ",
    "solapan. Cada mes que supera el p90 recibe el mismo peso en el ",
    "promedio.\n",
    sep = ""
  )
}

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

message("\n✅ Script 18 finalizado.")
