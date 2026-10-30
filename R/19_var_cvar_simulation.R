# ============================================================
# 19_var_cvar_simulation.R
# Teoría del Riesgo — VaR y CVaR de la pérdida por dengue
# ============================================================
#
# OBJETIVO
#   Cuantificar la pérdida semanal asociada al dengue en Colombia
#   (2022-2025) y estimar su Value at Risk (VaR) y su Conditional
#   Value at Risk (CVaR) al 95% y 99% con tres métodos: histórico,
#   Montecarlo y bootstrap por bloques. El VaR histórico se somete
#   a backtesting fuera de muestra con ventana móvil de 52 semanas.
#
# DEFINICIÓN DE LA PÉRDIDA
#
#   Loss_t = casos_dengue_t * 1.000.000 + M_dengue_t * 10.000
#
#   - 1.000.000 (COP por caso) es un COSTO FICTICIO DE SIMULACIÓN
#     inspirado en la literatura de carga económica del dengue
#     (costos directos e indirectos por caso). No se estima en este
#     script: es un supuesto para llevar casos y miedo a una misma
#     escala monetaria.
#   - 10.000 (COP por punto de M_dengue) es un FACTOR DE ESCALA PARA
#     EL MIEDO. M_dengue es un índice estandarizado (residuo de la
#     regresión auxiliar del script 12) que puede ser negativo; por
#     eso el término del miedo suma o resta una fracción pequeña
#     frente al término de casos, que domina la pérdida.
#
#   La pérdida está en pesos colombianos (COP). En la gráfica y en
#   la consola se expresa en miles de millones de COP (1e9).
#
# MÉTODOS
#
#   1. Histórico (incondicional): percentil empírico de la serie de
#      pérdidas al 95% y 99%. El CVaR es la media de las pérdidas
#      que caen en la cola derecha (>= percentil).
#
#   2. Montecarlo (condicional, horizonte 1 semana): se ajusta una
#      t-Student a los rendimientos logarítmicos de la pérdida,
#
#          r_t = log(Loss_t) - log(Loss_{t-1}),
#
#      y se simulan 5.000 choques semanales. Cada choque se aplica
#      a la última pérdida observada:
#
#          Loss_sim = Loss_T * exp(r_sim).
#
#      Los percentiles 95% y 99% de Loss_sim son el VaR; el CVaR es
#      la media de la cola. Si el ajuste t falla, se usa una normal
#      y se avisa en consola.
#
#   3. Bootstrap por bloques: boot::tsboot con bloques fijos de 8
#      semanas (l = 8, sim = "fixed"), que conservan la
#      autocorrelación y la estacionalidad local. En cada réplica se
#      recalculan VaR y CVaR; el estimador bootstrap es la mediana
#      de las 5.000 réplicas y se reporta también su intervalo del
#      95% (percentiles 2,5 y 97,5 de la distribución bootstrap).
#
#   Los tres métodos usan el percentil empírico tipo 7 (el default
#   de stats::quantile), para que sean comparables.
#
# BACKTESTING (VENTANA MÓVIL)
#
#   Para que Kupiec y Christoffersen sean válidos se hace backtesting
#   fuera de muestra. Con una ventana inicial de 52 semanas se calcula
#   el VaR al 95% de la semana 53 como el percentil 95% de las 52
#   pérdidas previas; luego la ventana se mueve una semana y se repite
#   hasta el final. Se obtienen 156 pronósticos (208 - 52).
#
#   rugarch::VaRTest(alpha = 0.05, actual = -loss_obs,
#                    VaR = -VaR_pronosticado) define la violación como
#   actual < VaR, es decir, -loss_obs < -VaR_pronosticado, que equivale
#   a loss_obs > VaR_pronosticado. Se niegan ambas series para llevar
#   la pérdida a la convención de retornos de rugarch; el signo se
#   cancela en la comparación.
#
#   De la prueba se extraen uc.LRp (Kupiec, cobertura incondicional) y
#   cc.LRp (Christoffersen, cobertura e independencia). El backtesting
#   se aplica únicamente al VaR histórico rodante; por eso las filas
#   de Montecarlo y bootstrap dejan esos p-valores en NA.
#
# INSUMOS
#   data/D_13/panel_ema_den_chik.rds   (panel región-mes EMA)
#   data/series/serie_dengue.rds       (serie nacional semanal)
#       Si no existe, se usa la ruta canónica data/D_07/serie_dengue.rds,
#       que es donde los scripts 07 y 09 guardan la serie. El cambio
#       se avisa en consola.
#
#   M_dengue nacional mensual: promedio simple de las 12 regiones EMA
#   por mes, igual que en el script 17. M es un índice y no un conteo,
#   así que sumarlo entre regiones cambiaría la escala.
#
# SALIDAS
#   data/D_19/var_cvar_resultados.csv
#       Metodo, VaR_95, CVaR_95, VaR_99, CVaR_99, Kupiec_p,
#       Christoffersen_p
#   data/D_19/var_backtest.csv
#       fecha_sem, perdida_observada, var_95_pronosticado, excede
#   graficas/G_19/var_backtest.png
#
# NOTA SOBRE LA ALINEACIÓN
#   La serie semanal se une al panel mensual por
#   mes = lubridate::floor_date(fecha_sem, "month"), de modo que cada
#   semana hereda el M del mes al que pertenece. El filtro 2022-2025
#   se aplica después de la unión. La pérdida es un nivel (no un
#   rendimiento), así que filtrar antes de diferenciar es correcto;
#   el primer rendimiento disponible es el de la segunda semana de
#   2022, porque el panel no tiene M antes de 2022-01.
#
# NOTA SOBRE LA COMPARACIÓN ENTRE MÉTODOS
#   El VaR histórico es incondicional (percentil de toda la muestra),
#   mientras que el de Montecarlo es condicional a la última pérdida
#   observada. Por eso sus niveles no tienen que coincidir: el primero
#   describe la cola de todo el periodo y el segundo, la próxima
#   semana dado el nivel actual. La tabla los reporta juntos y el
#   comentario de la gráfica lo recuerda.
#
# NOTA SOBRE EL USO DE ::
#   Todas las funciones no base se llaman con paquete::funcion() para
#   evitar el enmascaramiento entre dplyr, stats, MASS, boot, rugarch
#   y ggplot2. No se adjuntan paquetes con library().
# ============================================================

# --- 0. Rutas, carpetas y paquetes ---------------------------------

dir.create(here::here("data", "D_19"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("graficas", "G_19"), recursive = TRUE, showWarnings = FALSE)

ruta_panel         <- here::here("data", "D_13", "panel_ema_den_chik.rds")
ruta_serie_series  <- here::here("data", "series", "serie_dengue.rds")
ruta_serie_d07     <- here::here("data", "D_07", "serie_dengue.rds")
ruta_resultados    <- here::here("data", "D_19", "var_cvar_resultados.csv")
ruta_backtest      <- here::here("data", "D_19", "var_backtest.csv")
ruta_png           <- here::here("graficas", "G_19", "var_backtest.png")

desde_anio <- 2022L
hasta_anio <- 2025L
semilla    <- 2025L

nivel_95   <- 0.95
nivel_99   <- 0.99
n_mc       <- 5000L
n_boot     <- 5000L
bloque     <- 8L
ventana    <- 52L
alpha_bt   <- 0.05

paquetes <- c(
  "rugarch", "boot", "MASS", "dplyr", "ggplot2",
  "lubridate", "readr", "here", "tibble"
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

# Ruta de la serie: se prefiere data/series y se cae a D_07, que es
# donde los scripts 07 y 09 guardan la serie nacional semanal.
if (file.exists(ruta_serie_series)) {
  ruta_serie <- ruta_serie_series
} else if (file.exists(ruta_serie_d07)) {
  ruta_serie <- ruta_serie_d07
  message(
    "No existe ", ruta_serie_series, "; se usa la ruta canónica ",
    ruta_serie_d07, "."
  )
} else {
  stop(
    "No se encontró la serie semanal de dengue. Se buscó en ",
    ruta_serie_series, " y en ", ruta_serie_d07, ". ",
    "Ejecute primero R/07_serietemporal_conimputacion.R.",
    call. = FALSE
  )
}

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel del script 13: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

# --- 1. Serie nacional semanal -------------------------------------

serie_nacional <- readRDS(ruta_serie)

requeridas_serie <- c("dpto", "fecha_sem", "casos")
faltan_serie <- setdiff(requeridas_serie, names(serie_nacional))
if (length(faltan_serie) > 0L) {
  stop(
    "Faltan columnas en la serie semanal: ",
    paste(faltan_serie, collapse = ", "), ". Columnas disponibles: ",
    paste(names(serie_nacional), collapse = ", "), ".",
    call. = FALSE
  )
}

serie_nacional <- serie_nacional |>
  dplyr::filter(dpto == "NACIONAL") |>
  dplyr::mutate(
    fecha_sem = as.Date(fecha_sem),
    casos = as.numeric(casos)
  ) |>
  dplyr::arrange(fecha_sem) |>
  dplyr::select(fecha_sem, casos)

if (nrow(serie_nacional) < 100L) {
  stop(
    "La serie nacional tiene solo ", nrow(serie_nacional),
    " semanas; se esperaban al menos 100.",
    call. = FALSE
  )
}
if (anyNA(serie_nacional$fecha_sem) || anyNA(serie_nacional$casos)) {
  stop("La serie nacional tiene fechas o casos faltantes.", call. = FALSE)
}
if (anyDuplicated(serie_nacional$fecha_sem)) {
  stop("La serie nacional tiene semanas duplicadas.", call. = FALSE)
}

saltos_serie <- as.numeric(diff(serie_nacional$fecha_sem), units = "days")
if (any(saltos_serie > 7)) {
  stop(
    "La serie nacional tiene ", sum(saltos_serie > 7),
    " saltos mayores a 7 días; revise las semanas faltantes.",
    call. = FALSE
  )
}

# --- 2. M_dengue nacional mensual y expansión semanal ---------------

panel <- readRDS(ruta_panel)

requeridas_panel <- c("mes", "M_dengue")
faltan_panel <- setdiff(requeridas_panel, names(panel))
if (length(faltan_panel) > 0L) {
  stop(
    "Faltan columnas en el panel: ",
    paste(faltan_panel, collapse = ", "), ". Columnas disponibles: ",
    paste(names(panel), collapse = ", "), ".",
    call. = FALSE
  )
}

# M_dengue es un índice (residuo de la regresión auxiliar del script
# 12), no un conteo: se promedia entre las 12 regiones EMA. Sumarlo
# produciría una escala distinta y una estacionalidad espuria.
mensual_miedo <- panel |>
  dplyr::mutate(mes = as.Date(mes)) |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    M_dengue = mean(M_dengue, na.rm = TRUE),
    n_regiones = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(mes)

if (nrow(mensual_miedo) == 0L) {
  stop("El panel no tiene meses con M_dengue.", call. = FALSE)
}
if (anyNA(mensual_miedo$M_dengue)) {
  stop(
    "El miedo mensual tiene meses sin M_dengue. Revise el script 13.",
    call. = FALSE
  )
}
if (dplyr::n_distinct(mensual_miedo$n_regiones) != 1L) {
  warning(
    "El número de regiones por mes no es constante: ",
    paste(sort(unique(mensual_miedo$n_regiones)), collapse = ", "), ".",
    call. = FALSE
  )
}

# El panel debe cubrir todos los meses del periodo de análisis; si
# falta un mes, el left_join lo delataría como NA más abajo, pero es
# mejor detenerse con un mensaje claro.
meses_esperados <- seq(
  as.Date(sprintf("%d-01-01", desde_anio)),
  as.Date(sprintf("%d-12-01", hasta_anio)),
  by = "1 month"
)
meses_faltantes <- setdiff(meses_esperados, mensual_miedo$mes)
if (length(meses_faltantes) > 0L) {
  stop(
    "Faltan meses en el panel: ",
    paste(format(meses_faltantes, "%Y-%m"), collapse = ", "), ".",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat(" SERIE NACIONAL SEMANAL Y MIEDO MENSUAL\n")
cat("============================================================\n")
cat(sprintf(
  "Semanas nacionales: %d | Periodo: %s a %s\n",
  nrow(serie_nacional),
  format(min(serie_nacional$fecha_sem), "%Y-%m-%d"),
  format(max(serie_nacional$fecha_sem), "%Y-%m-%d")
))
cat(sprintf(
  "Meses con M_dengue: %d | Regiones por mes: %s\n",
  nrow(mensual_miedo),
  paste(sort(unique(mensual_miedo$n_regiones)), collapse = ", ")
))
cat(sprintf(
  "M_dengue nacional: media %.3f | mínimo %.3f | máximo %.3f\n",
  mean(mensual_miedo$M_dengue),
  min(mensual_miedo$M_dengue),
  max(mensual_miedo$M_dengue)
))

# --- 3. Alineación semanal-mensual, pérdida y muestra 2022-2025 ----
#
# floor_date(fecha_sem, "month") asigna a cada semana el primer día
# del mes al que pertenece; el left_join le hereda el M de ese mes.

datos <- serie_nacional |>
  dplyr::mutate(mes = lubridate::floor_date(fecha_sem, "month")) |>
  dplyr::left_join(
    dplyr::select(mensual_miedo, mes, M_dengue),
    by = "mes"
  ) |>
  dplyr::filter(
    lubridate::year(fecha_sem) >= desde_anio,
    lubridate::year(fecha_sem) <= hasta_anio
  ) |>
  dplyr::mutate(
    perdida = casos * 1000000 + M_dengue * 10000
  ) |>
  dplyr::select(fecha_sem, mes, casos, M_dengue, perdida) |>
  dplyr::arrange(fecha_sem)

if (nrow(datos) < 150L) {
  stop(
    "La muestra 2022-2025 tiene solo ", nrow(datos),
    " semanas; se esperaban al menos 150.",
    call. = FALSE
  )
}
if (anyNA(datos$M_dengue) || anyNA(datos$perdida)) {
  stop(
    "Hay valores faltantes después de la alineación: ",
    sum(is.na(datos$M_dengue)), " en M_dengue y ",
    sum(is.na(datos$perdida)), " en la pérdida.",
    call. = FALSE
  )
}

saltos_datos <- as.numeric(diff(datos$fecha_sem), units = "days")
if (any(saltos_datos != 7)) {
  stop(
    "La muestra tiene ", sum(saltos_datos != 7),
    " saltos distintos de 7 días; la ventana móvil de 52 semanas ",
    "requiere una serie semanal completa.",
    call. = FALSE
  )
}
if (anyDuplicated(datos$fecha_sem)) {
  stop("Hay semanas duplicadas después de la alineación.", call. = FALSE)
}
if (any(!is.finite(datos$perdida))) {
  stop("La pérdida tiene valores no finitos.", call. = FALSE)
}
if (any(datos$perdida <= 0)) {
  stop(
    "Hay ", sum(datos$perdida <= 0),
    " semanas con pérdida no positiva; los rendimientos logarítmicos ",
    "del método Montecarlo requieren Loss > 0.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat(" MUESTRA ALINEADA Y PÉRDIDA\n")
cat("============================================================\n")
cat(sprintf(
  "Semanas: %d | Periodo: %s a %s | Meses: %d\n",
  nrow(datos),
  format(min(datos$fecha_sem), "%Y-%m-%d"),
  format(max(datos$fecha_sem), "%Y-%m-%d"),
  dplyr::n_distinct(datos$mes)
))
cat(sprintf(
  "Casos: %s en total | máximo semanal: %s (%s)\n",
  format(sum(datos$casos), big.mark = ".", decimal.mark = ","),
  format(max(datos$casos), big.mark = ".", decimal.mark = ","),
  format(datos$fecha_sem[which.max(datos$casos)], "%Y-%m-%d")
))
cat(sprintf(
  "Pérdida (miles de millones COP): media %.3f | mínimo %.3f | máximo %.3f\n",
  mean(datos$perdida) / 1e9,
  min(datos$perdida) / 1e9,
  max(datos$perdida) / 1e9
))

# --- 4. Utilidades de riesgo ---------------------------------------

# Percentil empírico tipo 7 (default de R) sin nombres, para que las
# tres metodologías usen la misma definición.
percentil <- function(x, p) {
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 7))
}

# CVaR: media de la cola derecha. Se usa >= para no perder las
# observaciones exactamente iguales al percentil (puede haber empates
# porque la pérdida está dominada por casos semanales enteros).
cvar_empirico <- function(x, umbral) {
  cola <- x[x >= umbral]
  if (length(cola) == 0L) {
    return(NA_real_)
  }
  mean(cola)
}

medidas_riesgo <- function(x) {
  q95 <- percentil(x, nivel_95)
  q99 <- percentil(x, nivel_99)
  c(
    VaR_95  = q95,
    CVaR_95 = cvar_empirico(x, q95),
    VaR_99  = q99,
    CVaR_99 = cvar_empirico(x, q99)
  )
}

formatear_p <- function(p) {
  if (is.na(p)) {
    return("= NA")
  }
  if (p < 2.2e-16) {
    return("< 2.2e-16")
  }
  sprintf("= %.4f", p)
}

# --- 5. Método 1: VaR y CVaR histórico ------------------------------

var_hist <- medidas_riesgo(datos$perdida)

if (any(!is.finite(var_hist))) {
  stop("El VaR histórico tiene valores no finitos.", call. = FALSE)
}

# --- 6. Método 2: VaR y CVaR Montecarlo -----------------------------
#
# Se ajusta una t-Student a los rendimientos logarítmicos de la
# pérdida y se simulan 5.000 choques semanales. Si el ajuste t falla
# o produce parámetros inválidos, se cae a una normal.

rendimientos <- diff(log(datos$perdida))
if (any(!is.finite(rendimientos))) {
  stop("Los rendimientos de la pérdida tienen valores no finitos.",
       call. = FALSE)
}

ajustar_rendimientos <- function(x) {
  ajuste_t <- tryCatch(
    suppressWarnings(MASS::fitdistr(x, "t")),
    error = function(e) NULL
  )

  if (!is.null(ajuste_t)) {
    est <- ajuste_t$estimate
    if (all(c("m", "s", "df") %in% names(est)) &&
        all(is.finite(est)) &&
        est[["s"]] > 0 &&
        est[["df"]] > 2) {
      return(list(
        tipo = "t-Student",
        m = unname(est[["m"]]),
        s = unname(est[["s"]]),
        df = unname(est[["df"]])
      ))
    }
  }

  message(
    "  El ajuste t-Student no fue válido; se usa una normal para ",
    "los rendimientos de la pérdida."
  )
  list(tipo = "normal", m = mean(x), s = stats::sd(x), df = NA_real_)
}

ajuste_mc <- ajustar_rendimientos(rendimientos)
if (!is.finite(ajuste_mc$s) || ajuste_mc$s <= 0) {
  stop("La escala ajustada de los rendimientos no es válida.",
       call. = FALSE)
}

set.seed(semilla)
if (ajuste_mc$tipo == "t-Student") {
  choques <- ajuste_mc$m + ajuste_mc$s * stats::rt(n_mc, df = ajuste_mc$df)
} else {
  choques <- stats::rnorm(n_mc, mean = ajuste_mc$m, sd = ajuste_mc$s)
}

ultima_perdida <- datos$perdida[nrow(datos)]
perdida_sim <- ultima_perdida * exp(choques)

if (length(perdida_sim) != n_mc || any(!is.finite(perdida_sim))) {
  stop("La simulación Montecarlo produjo valores no finitos.",
       call. = FALSE)
}
if (any(perdida_sim <= 0)) {
  stop("La simulación Montecarlo produjo pérdidas no positivas.",
       call. = FALSE)
}

var_mc <- medidas_riesgo(perdida_sim)

if (any(!is.finite(var_mc))) {
  stop("El VaR Montecarlo tiene valores no finitos.", call. = FALSE)
}

# --- 7. Método 3: VaR y CVaR bootstrap por bloques ------------------
#
# boot::tsboot remuestrea bloques fijos de 8 semanas para conservar
# la autocorrelación. En cada réplica se recalculan VaR y CVaR; el
# estimador bootstrap es la mediana de las réplicas y el intervalo
# del 95% son los percentiles 2,5 y 97,5 de esa distribución.

set.seed(semilla)
boot_out <- boot::tsboot(
  datos$perdida,
  statistic = medidas_riesgo,
  R = n_boot,
  l = bloque,
  sim = "fixed"
)

if (!is.matrix(boot_out$t) || ncol(boot_out$t) != 4L ||
    nrow(boot_out$t) != n_boot) {
  stop(
    "boot::tsboot devolvió una matriz de réplicas inesperada: ",
    paste(dim(boot_out$t), collapse = " x "), ".",
    call. = FALSE
  )
}
if (anyNA(boot_out$t) || any(!is.finite(boot_out$t))) {
  stop("La distribución bootstrap tiene valores no finitos.",
       call. = FALSE)
}

colnames(boot_out$t) <- c("VaR_95", "CVaR_95", "VaR_99", "CVaR_99")
var_boot <- apply(boot_out$t, 2, stats::median)
var_boot_ic <- apply(
  boot_out$t, 2, stats::quantile,
  probs = c(0.025, 0.975), names = FALSE
)

if (any(!is.finite(var_boot))) {
  stop("El VaR bootstrap tiene valores no finitos.", call. = FALSE)
}

# --- 8. Backtesting con ventana móvil de 52 semanas -----------------
#
# Para cada semana i = 53, ..., n, el VaR se estima con las 52
# pérdidas anteriores (i-52):(i-1) y se compara con la pérdida
# observada de la semana i. Es un pronóstico genuinamente fuera de
# muestra: la semana i nunca entra en su propia ventana.

n_semanas <- nrow(datos)
if (n_semanas <= ventana) {
  stop(
    "La muestra tiene ", n_semanas, " semanas y la ventana inicial es ",
    ventana, "; no quedan semanas para backtesting.",
    call. = FALSE
  )
}

indices_bt <- (ventana + 1L):n_semanas
var_pronosticado <- vapply(
  indices_bt,
  function(i) {
    percentil(datos$perdida[(i - ventana):(i - 1L)], nivel_95)
  },
  numeric(1)
)
perdida_observada <- datos$perdida[indices_bt]
fechas_bt <- datos$fecha_sem[indices_bt]

if (length(var_pronosticado) != length(indices_bt) ||
    length(perdida_observada) != length(indices_bt)) {
  stop("El backtesting no generó un pronóstico por semana.", call. = FALSE)
}
if (any(!is.finite(var_pronosticado)) ||
    any(!is.finite(perdida_observada))) {
  stop("El backtesting tiene valores no finitos.", call. = FALSE)
}

excede <- perdida_observada > var_pronosticado

backtest <- tibble::tibble(
  fecha_sem = fechas_bt,
  perdida_observada = perdida_observada,
  var_95_pronosticado = var_pronosticado,
  excede = excede
)

# rugarch::VaRTest espera "actual" y "VaR" con la convención de
# retornos: violación si actual < VaR. Al pasar -pérdida y -VaR, la
# violación queda pérdida_obs > VaR_pronosticado.
if (sum(excede) == 0L || sum(excede) == length(excede)) {
  warning(
    "El backtesting tiene ", sum(excede), " excedencias de ",
    length(excede), "; VaRTest no puede calcular los estadísticos. ",
    "Los p-valores se dejan en NA.",
    call. = FALSE
  )
  kupiec_p <- NA_real_
  christoffersen_p <- NA_real_
  prueba_var <- NULL
} else {
  prueba_var <- tryCatch(
    rugarch::VaRTest(
      alpha = alpha_bt,
      actual = -perdida_observada,
      VaR = -var_pronosticado
    ),
    error = function(e) e
  )

  if (inherits(prueba_var, "error")) {
    stop(
      "rugarch::VaRTest falló: ", conditionMessage(prueba_var), ".",
      call. = FALSE
    )
  }

  kupiec_p <- as.numeric(prueba_var$uc.LRp)
  christoffersen_p <- as.numeric(prueba_var$cc.LRp)

  if (!is.finite(kupiec_p)) {
    kupiec_p <- NA_real_
  }
  if (!is.finite(christoffersen_p)) {
    christoffersen_p <- NA_real_
  }

  # Coherencia entre el conteo propio y el de rugarch.
  if (as.integer(prueba_var$actual.exceed) != sum(excede)) {
    warning(
      "rugarch reporta ", prueba_var$actual.exceed,
      " excedencias y el conteo propio es ", sum(excede), ".",
      call. = FALSE
    )
  }
}

# --- 9. Tabla de resultados ----------------------------------------

resultados <- tibble::tibble(
  Metodo = c("Historico", "Montecarlo", "Bootstrap"),
  VaR_95 = c(
    unname(var_hist[["VaR_95"]]),
    unname(var_mc[["VaR_95"]]),
    unname(var_boot[["VaR_95"]])
  ),
  CVaR_95 = c(
    unname(var_hist[["CVaR_95"]]),
    unname(var_mc[["CVaR_95"]]),
    unname(var_boot[["CVaR_95"]])
  ),
  VaR_99 = c(
    unname(var_hist[["VaR_99"]]),
    unname(var_mc[["VaR_99"]]),
    unname(var_boot[["VaR_99"]])
  ),
  CVaR_99 = c(
    unname(var_hist[["CVaR_99"]]),
    unname(var_mc[["CVaR_99"]]),
    unname(var_boot[["CVaR_99"]])
  ),
  Kupiec_p = c(kupiec_p, NA_real_, NA_real_),
  Christoffersen_p = c(christoffersen_p, NA_real_, NA_real_)
)

if (any(!is.finite(as.matrix(resultados[, 2:5])))) {
  stop("La tabla de resultados tiene VaR o CVaR no finitos.",
       call. = FALSE)
}
if (any(resultados$CVaR_95 < resultados$VaR_95 - 1e-8) ||
    any(resultados$CVaR_99 < resultados$VaR_99 - 1e-8)) {
  warning(
    "Algún CVaR es menor que su VaR; revise el cálculo de la cola.",
    call. = FALSE
  )
}

# --- 10. Gráfica: pérdida observada vs. VaR pronosticado -----------

colores <- c(
  "Pérdida observada" = "#252525",
  "VaR 95% pronosticado" = "#B2182B"
)

p_backtest <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = datos,
    ggplot2::aes(
      x = fecha_sem,
      y = perdida / 1e9,
      colour = "Pérdida observada"
    ),
    linewidth = 0.55
  ) +
  ggplot2::geom_line(
    data = backtest,
    ggplot2::aes(
      x = fecha_sem,
      y = var_95_pronosticado / 1e9,
      colour = "VaR 95% pronosticado"
    ),
    linewidth = 0.75,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(backtest, excede),
    ggplot2::aes(x = fecha_sem, y = perdida_observada / 1e9),
    colour = "#B2182B",
    size = 1.1,
    shape = 16
  ) +
  ggplot2::scale_colour_manual(values = colores, name = NULL) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      override.aes = list(linetype = c("solid", "dashed"))
    )
  ) +
  ggplot2::scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y-%m",
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.1f", x),
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  ggplot2::labs(
    title = "Pérdida observada vs. VaR al 95% con ventana móvil",
    subtitle = sprintf(
      paste0(
        "Ventana de %d semanas | %d pronósticos | %d excedencias ",
        "(esperadas %.1f) | Kupiec p %s | Christoffersen p %s"
      ),
      ventana,
      nrow(backtest),
      sum(backtest$excede),
      alpha_bt * nrow(backtest),
      formatear_p(kupiec_p),
      formatear_p(christoffersen_p)
    ),
    x = NULL,
    y = "Pérdida semanal (miles de millones de COP)",
    caption = paste0(
      "Pérdida = casos x 1.000.000 + M_dengue x 10.000 (costos ficticios ",
      "de simulación). Punto rojo: la pérdida superó el VaR pronosticado. ",
      "El VaR histórico es incondicional; el de Montecarlo es condicional ",
      "a la última pérdida."
    )
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "left",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7.2, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_png,
  plot = p_backtest,
  width = 10.5,
  height = 5.4,
  dpi = 200,
  bg = "white"
)

# --- 11. Guardado --------------------------------------------------

readr::write_csv(resultados, ruta_resultados, na = "NA")
readr::write_csv(backtest, ruta_backtest, na = "NA")

rutas_salida <- c(ruta_resultados, ruta_backtest, ruta_png)
for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 12. Resumen en consola ----------------------------------------

cat("\n============================================================\n")
cat(" VaR Y CVaR DE LA PÉRDIDA POR DENGUE (miles de millones COP)\n")
cat("============================================================\n")
cat(sprintf(
  "Muestra histórica: %d semanas (%s a %s) | %d rendimientos\n",
  nrow(datos),
  format(min(datos$fecha_sem), "%Y-%m-%d"),
  format(max(datos$fecha_sem), "%Y-%m-%d"),
  length(rendimientos)
))
cat(sprintf(
  "Pérdida: casos x 1.000.000 + M_dengue x 10.000 | última = %.3f\n",
  ultima_perdida / 1e9
))

cat("\n--- VaR y CVaR por método ---\n")
for (i in seq_len(nrow(resultados))) {
  cat(sprintf(
    "%-11s | VaR 95%% = %7.3f | CVaR 95%% = %7.3f | VaR 99%% = %7.3f | CVaR 99%% = %7.3f\n",
    resultados$Metodo[i],
    resultados$VaR_95[i] / 1e9,
    resultados$CVaR_95[i] / 1e9,
    resultados$VaR_99[i] / 1e9,
    resultados$CVaR_99[i] / 1e9
  ))
}

cat(sprintf(
  "\nMontecarlo: %s(m = %.5f, s = %.5f, df = %s) | %d simulaciones\n",
  ajuste_mc$tipo,
  ajuste_mc$m,
  ajuste_mc$s,
  if (is.na(ajuste_mc$df)) "NA" else sprintf("%.2f", ajuste_mc$df),
  n_mc
))
cat(sprintf(
  "Bootstrap: %d réplicas de bloques de %d semanas | IC 95%% VaR 95%%: %.3f a %.3f\n",
  n_boot,
  bloque,
  var_boot_ic[1L, "VaR_95"] / 1e9,
  var_boot_ic[2L, "VaR_95"] / 1e9
))

cat("\n--- Backtesting del VaR histórico (ventana móvil) ---\n")
cat(sprintf(
  "Ventana inicial: %d semanas | Pronósticos: %d | Excedencias: %d (esperadas %.1f)\n",
  ventana,
  nrow(backtest),
  sum(backtest$excede),
  alpha_bt * nrow(backtest)
))
cat(sprintf(
  "Tasa de excedencia observada: %.2f %% (nominal %.1f %%)\n",
  100 * mean(backtest$excede),
  100 * alpha_bt
))
cat(sprintf("Kupiec (uc.LRp) %s", formatear_p(kupiec_p)))
cat(sprintf(
  " -> %s\n",
  if (is.na(kupiec_p)) {
    "no calculable"
  } else if (kupiec_p < 0.05) {
    "se rechaza la cobertura correcta al 5%"
  } else {
    "no se rechaza la cobertura correcta al 5%"
  }
))
cat(sprintf("Christoffersen (cc.LRp) %s", formatear_p(christoffersen_p)))
cat(sprintf(
  " -> %s\n",
  if (is.na(christoffersen_p)) {
    "no calculable"
  } else if (christoffersen_p < 0.05) {
    "se rechaza independencia y cobertura al 5%"
  } else {
    "no se rechaza independencia y cobertura al 5%"
  }
))

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

cat("\nNotas:\n")
cat("  - 1.000.000 COP por caso es un costo ficticio de simulación\n")
cat("    inspirado en la literatura de carga económica; 10.000 COP por\n")
cat("    punto de M_dengue es un factor de escala para el miedo.\n")
cat("  - M_dengue es el promedio mensual de las 12 regiones EMA y puede\n")
cat("    ser negativo; el término del miedo es pequeño frente al de casos.\n")
cat("  - El VaR de Montecarlo es condicional a la última pérdida\n")
cat("    observada (horizonte 1 semana); el histórico es incondicional.\n")
cat("  - Kupiec y Christoffersen solo se calculan para el VaR histórico\n")
cat("    rodante; las otras dos filas quedan en NA en la tabla.\n")
cat("  - Un p-valor bajo indica que el VaR rodante no cubre bien las\n")
cat("    pérdidas: en una epidemia creciente, la ventana de 52 semanas\n")
cat("    tarda en incorporar el aumento de casos.\n")

message("\n✅ Script 19 finalizado.")
