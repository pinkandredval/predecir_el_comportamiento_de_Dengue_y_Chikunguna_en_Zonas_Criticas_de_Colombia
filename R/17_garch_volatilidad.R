# ============================================================
# 17_garch_volatilidad.R
# Tabla 4 — Teoría del Riesgo: el miedo en la varianza (GARCH)
# ============================================================
#
# OBJETIVO
#   Contrastar si el miedo (M_dengue) explica la volatilidad
#   condicional del rendimiento epidémico semanal del dengue. La
#   pregunta de la Tabla 4 es si, más allá del nivel de casos, el
#   miedo agrega una prima de riesgo: un coeficiente vxreg1 > 0 y
#   significativo implica que el miedo aumenta la volatilidad.
#
# MODELO
#   Rendimiento epidémico semanal:
#
#       r_t = log(casos_t + 1) - log(casos_{t-1} + 1)
#           = diff(log(casos + 1))
#
#   Media (idéntica en los tres modelos):
#
#       r_t = mu + phi r_{t-1} + e_t,     e_t = sigma_t z_t,
#       z_t ~ N(0, 1)
#
#   Varianza (el miedo entra como regresor externo, vxreg1 = delta):
#
#     sGARCH:   sigma^2_t = omega + alpha e^2_{t-1}
#                           + beta sigma^2_{t-1} + delta M_t
#
#     GJR-GARCH: sigma^2_t = omega + alpha e^2_{t-1}
#                           + gamma e^2_{t-1} 1{e_{t-1} < 0}
#                           + beta sigma^2_{t-1} + delta M_t
#
#     EGARCH:   log sigma^2_t = omega
#                           + alpha (|z_{t-1}| - E|z|)
#                           + gamma z_{t-1}
#                           + beta log sigma^2_{t-1} + delta M_t
#
#   H0: delta = 0 (el miedo no explica la volatilidad).
#   H1: delta != 0. La lectura es direccional: delta > 0 aumenta la
#   volatilidad y delta < 0 la reduce. El signo se reporta tal como
#   se estima; no se fuerza la hipótesis de aumento.
#
# INSUMOS
#   data/series/serie_dengue.rds       (serie nacional semanal)
#       Si no existe, el script usa data/D_07/serie_dengue.rds, que
#       es la ruta canónica que dejan los scripts 07 y 09 de este
#       repositorio. El aviso se imprime en consola.
#   data/D_13/panel_ema_den_chik.rds   (panel región-mes EMA)
#
#   De la serie semanal se usa dpto == "NACIONAL". Del panel se
#   promedia M_dengue por mes (promedio simple de las 12 regiones
#   EMA), porque M es un índice y no un conteo.
#
# SALIDAS
#   data/D_17/garch_resultados.csv
#       Modelo, AIC, vxreg1_Estimate, vxreg1_p_value
#   data/D_17/serie_volatilidad.csv
#       fecha_sem, r_t, M_dengue, volatilidad, modelo
#   graficas/G_17/volatilidad_condicional.png
#
# NOTA SOBRE LA ALINEACIÓN (CRÍTICA)
#   La serie semanal tiene n semanas y diff() devuelve n - 1
#   rendimientos. Para que r_t y M_t tengan exactamente la misma
#   longitud se hace:
#
#     1. left_join semanal-mensual por
#        mes = lubridate::floor_date(fecha_sem, "month"), de modo que
#        cada semana hereda el M del mes al que pertenece.
#     2. Se calcula r_t = diff(log(casos + 1)) sobre la serie ordenada.
#     3. Se elimina la primera fila (la semana sin rendimiento) con
#        dplyr::slice(-1L). Así la fila i de la tabla queda con el
#        rendimiento r_t[i] y el miedo M_t[i] del mismo mes.
#
#   El filtro 2022-2025 se aplica DESPUÉS de calcular el rendimiento:
#   el primer rendimiento de 2022 usa la última semana de 2021 como
#   insumo (pre-muestra), igual que en los scripts 13 y 16.
#
# NOTA SOBRE EL AIC
#   rugarch::infocriteria() devuelve el AIC por observación. Para la
#   tabla se guarda el AIC total, 2k - 2 log L, que es la convención
#   del script 09 de este proyecto. Como los tres modelos usan la
#   misma muestra, la misma media y la misma distribución, el orden
#   por AIC total coincide con el orden por AIC por observación.
#
# NOTA SOBRE EL USO DE ::
#   Todas las funciones no base se llaman con paquete::funcion() para
#   evitar el enmascaramiento entre dplyr, stats, rugarch y ggplot2.
#   No se adjuntan paquetes con library().
# ============================================================

# --- 0. Rutas, carpetas y paquetes ---------------------------------

dir.create(here::here("data", "D_17"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("graficas", "G_17"), recursive = TRUE, showWarnings = FALSE)

ruta_serie_series <- here::here("data", "series", "serie_dengue.rds")
ruta_serie_d07    <- here::here("data", "D_07", "serie_dengue.rds")
ruta_panel        <- here::here("data", "D_13", "panel_ema_den_chik.rds")
ruta_resultados   <- here::here("data", "D_17", "garch_resultados.csv")
ruta_volatilidad  <- here::here("data", "D_17", "serie_volatilidad.csv")
ruta_png          <- here::here("graficas", "G_17", "volatilidad_condicional.png")

desde_anio <- 2022L
hasta_anio <- 2025L
semilla    <- 2025L

paquetes <- c("rugarch", "dplyr", "ggplot2", "lubridate", "readr", "here")
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

# Ruta de la serie: se prefiere la ruta del plan y se cae a D_07, que
# es donde el script 07 guarda la serie nacional semanal.
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
  dplyr::arrange(fecha_sem) |>
  dplyr::mutate(
    fecha_sem = as.Date(fecha_sem),
    casos = as.numeric(casos),
    mes = lubridate::floor_date(fecha_sem, "month")
  )

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

saltos <- as.numeric(diff(serie_nacional$fecha_sem), units = "days")
if (any(saltos > 7)) {
  stop(
    "La serie nacional tiene ", sum(saltos > 7),
    " saltos mayores a 7 días; revise las semanas faltantes.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat(" SERIE NACIONAL SEMANAL\n")
cat("============================================================\n")
cat(sprintf(
  "Semanas: %d | Periodo: %s a %s\n",
  nrow(serie_nacional),
  format(min(serie_nacional$fecha_sem), "%Y-%m-%d"),
  format(max(serie_nacional$fecha_sem), "%Y-%m-%d")
))
cat(sprintf(
  "Casos: %s en total | máximo semanal: %s (%s)\n",
  format(sum(serie_nacional$casos), big.mark = ".", decimal.mark = ","),
  format(max(serie_nacional$casos), big.mark = ".", decimal.mark = ","),
  format(
    serie_nacional$fecha_sem[which.max(serie_nacional$casos)],
    "%Y-%m-%d"
  )
))

# --- 2. Miedo mensual nacional -------------------------------------

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

cat("\n============================================================\n")
cat(" MIEDO MENSUAL NACIONAL\n")
cat("============================================================\n")
cat(sprintf(
  "Meses: %d | Periodo: %s a %s | Regiones por mes: %s\n",
  nrow(mensual_miedo),
  format(min(mensual_miedo$mes), "%Y-%m"),
  format(max(mensual_miedo$mes), "%Y-%m"),
  paste(sort(unique(mensual_miedo$n_regiones)), collapse = ", ")
))
cat(sprintf(
  "M_dengue: media %.3f | mínimo %.3f | máximo %.3f\n",
  mean(mensual_miedo$M_dengue),
  min(mensual_miedo$M_dengue),
  max(mensual_miedo$M_dengue)
))

# --- 3. Alineación semanal-mensual y rendimiento -------------------
#
# left_join por mes: cada semana hereda el M del mes de su fecha_sem.
# Luego diff() pierde la primera semana y se recorta esa misma fila
# para que M y r_t queden uno a uno.

serie_alineada <- serie_nacional |>
  dplyr::left_join(
    dplyr::select(mensual_miedo, mes, M_dengue),
    by = "mes"
  )

rendimiento <- diff(log(serie_alineada$casos + 1))

if (length(rendimiento) != nrow(serie_alineada) - 1L) {
  stop(
    "diff() devolvió ", length(rendimiento),
    " rendimientos para ", nrow(serie_alineada), " semanas.",
    call. = FALSE
  )
}

datos_modelo <- serie_alineada |>
  dplyr::slice(-1L) |>
  dplyr::mutate(r_t = as.numeric(rendimiento)) |>
  dplyr::filter(
    lubridate::year(fecha_sem) >= desde_anio,
    lubridate::year(fecha_sem) <= hasta_anio
  ) |>
  dplyr::select(fecha_sem, mes, casos, M_dengue, r_t) |>
  dplyr::arrange(fecha_sem)

n_esperado <- sum(
  lubridate::year(serie_alineada$fecha_sem[-1L]) >= desde_anio &
    lubridate::year(serie_alineada$fecha_sem[-1L]) <= hasta_anio
)
if (nrow(datos_modelo) != n_esperado) {
  warning(
    "La muestra tiene ", nrow(datos_modelo), " filas y se esperaban ",
    n_esperado, "; revise la alineación.",
    call. = FALSE
  )
}

if (nrow(datos_modelo) < 100L) {
  stop(
    "La muestra 2022-2025 tiene solo ", nrow(datos_modelo),
    " semanas; se esperaban al menos 100.",
    call. = FALSE
  )
}
if (anyNA(datos_modelo[, c("M_dengue", "r_t")])) {
  stop(
    "Hay valores faltantes en r_t o M_dengue después de la alineación: ",
    sum(is.na(datos_modelo$r_t)), " en r_t y ",
    sum(is.na(datos_modelo$M_dengue)), " en M_dengue.",
    call. = FALSE
  )
}
if (anyDuplicated(datos_modelo$fecha_sem)) {
  stop("Hay semanas duplicadas después de la alineación.", call. = FALSE)
}

saltos_modelo <- as.numeric(diff(datos_modelo$fecha_sem), units = "days")
if (any(saltos_modelo > 7)) {
  stop(
    "Hay ", sum(saltos_modelo > 7),
    " saltos semanales en la muestra del modelo.",
    call. = FALSE
  )
}

# Cada año debe aportar semanas y cada mes debe tener miedo. Si falta
# un mes completo, el left_join lo delataría como NA (ya comprobado).
meses_esperados <- seq(
  as.Date(sprintf("%d-01-01", desde_anio)),
  as.Date(sprintf("%d-12-01", hasta_anio)),
  by = "1 month"
)
meses_faltantes <- setdiff(meses_esperados, unique(datos_modelo$mes))
if (length(meses_faltantes) > 0L) {
  stop(
    "Faltan meses en la muestra: ",
    paste(format(meses_faltantes, "%Y-%m"), collapse = ", "), ".",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat(" ALINEACIÓN Y MUESTRA DEL MODELO\n")
cat("============================================================\n")
cat(sprintf(
  "Semanas con rendimiento: %d | Periodo: %s a %s\n",
  nrow(datos_modelo),
  format(min(datos_modelo$fecha_sem), "%Y-%m-%d"),
  format(max(datos_modelo$fecha_sem), "%Y-%m-%d")
))
cat(sprintf(
  "Meses representados: %d | M_dengue sin faltantes: %d de %d\n",
  dplyr::n_distinct(datos_modelo$mes),
  sum(!is.na(datos_modelo$M_dengue)), nrow(datos_modelo)
))
cat(sprintf(
  "r_t: media %.4f | desv. est. %.4f | mínimo %.4f | máximo %.4f\n",
  mean(datos_modelo$r_t), stats::sd(datos_modelo$r_t),
  min(datos_modelo$r_t), max(datos_modelo$r_t)
))
cat(sprintf(
  "Primera observación: %s (r_t = %.4f, M = %.3f)\n",
  format(datos_modelo$fecha_sem[1], "%Y-%m-%d"),
  datos_modelo$r_t[1], datos_modelo$M_dengue[1]
))
cat(sprintf(
  "Última observación:  %s (r_t = %.4f, M = %.3f)\n",
  format(datos_modelo$fecha_sem[nrow(datos_modelo)], "%Y-%m-%d"),
  datos_modelo$r_t[nrow(datos_modelo)],
  datos_modelo$M_dengue[nrow(datos_modelo)]
))

# --- 4. Especificaciones y ajuste ----------------------------------

especificaciones <- c(
  "sGARCH"    = "sGARCH",
  "GJR-GARCH" = "gjrGARCH",
  "EGARCH"    = "eGARCH"
)

ajustar_garch <- function(modelo, datos) {
  spec <- rugarch::ugarchspec(
    variance.model = list(
      model = modelo,
      garchOrder = c(1, 1),
      external.regressors = cbind(M_dengue = datos$M_dengue)
    ),
    mean.model = list(
      armaOrder = c(1, 0),
      include.mean = TRUE
    ),
    distribution.model = "norm"
  )

  # solver = "hybrid" es el solver principal (combina solnp y
  # nlminb). Si falla por completo, se reintenta con solnp para no
  # perder la tabla; el cambio se avisa en consola.
  fit <- tryCatch(
    rugarch::ugarchfit(spec = spec, data = datos$r_t, solver = "hybrid"),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    message(
      "  El solver híbrido falló para ", modelo,
      "; se reintenta con solver = \"solnp\"."
    )
    fit <- tryCatch(
      rugarch::ugarchfit(spec = spec, data = datos$r_t, solver = "solnp"),
      error = function(e) e
    )
  }

  if (inherits(fit, "error")) {
    stop(
      "No se pudo ajustar el modelo ", modelo, ": ",
      conditionMessage(fit), ".",
      call. = FALSE
    )
  }
  if (is.null(fit@fit$matcoef)) {
    stop(
      "El modelo ", modelo, " no devolvió la matriz de coeficientes.",
      call. = FALSE
    )
  }
  if (rugarch::convergence(fit) != 0L) {
    warning(
      "El modelo ", modelo, " terminó con código de convergencia ",
      rugarch::convergence(fit), ".",
      call. = FALSE
    )
  }

  fit
}

set.seed(semilla)

ajustes <- vector("list", length(especificaciones))
names(ajustes) <- names(especificaciones)

cat("\n============================================================\n")
cat(" AJUSTE DE LOS TRES MODELOS GARCH\n")
cat("============================================================\n")

for (etiqueta in names(especificaciones)) {
  message("Ajustando ", etiqueta, " ...")
  ajustes[[etiqueta]] <- ajustar_garch(
    especificaciones[[etiqueta]], datos_modelo
  )
  cat(sprintf(
    "  %-10s ajustado | convergencia = %d | T = %d\n",
    etiqueta, rugarch::convergence(ajustes[[etiqueta]]),
    length(datos_modelo$r_t)
  ))
}

# --- 5. Extracción de vxreg1 y AIC ---------------------------------

# vxreg1 es el nombre que rugarch le da al primer regresor externo de
# la ecuación de varianza. El contraste es H0: vxreg1 = 0; el p-valor
# es Pr(>|t|) de la matriz matcoef.
extraer_vxreg1 <- function(fit, etiqueta) {
  mat <- fit@fit$matcoef
  colnames(mat) <- trimws(colnames(mat))

  if (!"vxreg1" %in% rownames(mat)) {
    stop(
      "El modelo ", etiqueta, " no reporta vxreg1. Coeficientes: ",
      paste(rownames(mat), collapse = ", "), ".",
      call. = FALSE
    )
  }

  estimacion <- unname(mat["vxreg1", "Estimate"])
  p_valor <- unname(mat["vxreg1", "Pr(>|t|)"])
  error_estandar <- unname(mat["vxreg1", "Std. Error"])

  # AIC total = 2k - 2 log L. rugarch::infocriteria() devuelve el AIC
  # por observación, que se conserva para la impresión en consola.
  k <- length(fit@fit$coef)
  log_verosimilitud <- as.numeric(fit@fit$LLH)
  aic_total <- 2 * k - 2 * log_verosimilitud
  aic_observacion <- as.numeric(rugarch::infocriteria(fit)[1, 1])

  tibble::tibble(
    Modelo = etiqueta,
    AIC = aic_total,
    vxreg1_Estimate = estimacion,
    vxreg1_p_value = p_valor,
    vxreg1_Std_Error = error_estandar,
    AIC_por_observacion = aic_observacion,
    LogLik = log_verosimilitud,
    n_parametros = k,
    Convergencia = rugarch::convergence(fit)
  )
}

resultados <- dplyr::bind_rows(
  lapply(names(ajustes), function(etiqueta) {
    extraer_vxreg1(ajustes[[etiqueta]], etiqueta)
  })
)

if (any(!is.finite(resultados$AIC))) {
  stop("Algún modelo produjo un AIC no finito.", call. = FALSE)
}

# El AIC total debe ser coherente con el AIC por observación de
# rugarch multiplicado por T.
aic_implicito <- resultados$AIC_por_observacion * nrow(datos_modelo)
if (any(abs(aic_implicito - resultados$AIC) > 1e-4)) {
  warning(
    "El AIC total no coincide con infocriteria() * T; revise la ",
    "extracción de la verosimilitud.",
    call. = FALSE
  )
}

# --- 6. Mejor modelo por AIC (GJR gana el empate) ------------------

aics <- resultados$AIC
min_aic <- min(aics)
tolerancia_empate <- 1e-6
empatados <- resultados$Modelo[abs(aics - min_aic) <= tolerancia_empate]

mejor_modelo <- if ("GJR-GARCH" %in% empatados) {
  "GJR-GARCH"
} else {
  empatados[[1L]]
}

mejor_fit <- ajustes[[mejor_modelo]]
mejor_fila <- dplyr::filter(resultados, Modelo == mejor_modelo)

if (nrow(mejor_fila) != 1L) {
  stop("No se pudo identificar una única fila para el mejor modelo.",
       call. = FALSE)
}

# --- 7. Serie de volatilidad condicional y gráfico -----------------

volatilidad <- as.numeric(rugarch::sigma(mejor_fit))

if (length(volatilidad) != nrow(datos_modelo)) {
  stop(
    "sigma(fit) devolvió ", length(volatilidad),
    " valores para ", nrow(datos_modelo), " observaciones.",
    call. = FALSE
  )
}
if (anyNA(volatilidad) || any(!is.finite(volatilidad))) {
  stop("La volatilidad condicional tiene valores no finitos.",
       call. = FALSE)
}

serie_volatilidad <- tibble::tibble(
  fecha_sem = datos_modelo$fecha_sem,
  r_t = datos_modelo$r_t,
  M_dengue = datos_modelo$M_dengue,
  volatilidad = volatilidad,
  modelo = mejor_modelo
)

indice_maximo <- which.max(serie_volatilidad$volatilidad)
fecha_maximo <- serie_volatilidad$fecha_sem[indice_maximo]
valor_maximo <- serie_volatilidad$volatilidad[indice_maximo]

# Formato legible del p-valor para títulos y consola.
formatear_p <- function(p) {
  if (is.na(p)) {
    return("= NA")
  }
  if (p < 2.2e-16) {
    return("< 2.2e-16")
  }
  sprintf("= %.4f", p)
}

interpretar_miedo <- function(estimacion, p) {
  if (is.na(p)) {
    return("p-valor no disponible")
  }
  if (p >= 0.05) {
    return("no significativo al 5%")
  }
  if (estimacion > 0) {
    return("aumenta la volatilidad (p < 0.05)")
  }
  "reduce la volatilidad (p < 0.05)"
}

lectura_mejor <- interpretar_miedo(
  mejor_fila$vxreg1_Estimate, mejor_fila$vxreg1_p_value
)

p_volatilidad <- ggplot2::ggplot(
  serie_volatilidad,
  ggplot2::aes(x = fecha_sem, y = volatilidad)
) +
  ggplot2::geom_line(colour = "#B2182B", linewidth = 0.7) +
  ggplot2::geom_point(
    data = dplyr::filter(
      serie_volatilidad, volatilidad == valor_maximo
    ),
    colour = "#67000D", size = 1.8
  ) +
  ggplot2::scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y-%m",
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) sprintf("%.3f", x),
    expand = ggplot2::expansion(mult = c(0.02, 0.06))
  ) +
  ggplot2::labs(
    title = "Volatilidad condicional del dengue con el miedo en la varianza",
    subtitle = sprintf(
      paste0(
        "%s(1,1) seleccionado por AIC | r_t = \u0394log(casos + 1) | ",
        "2022-2025 | T = %d"
      ),
      mejor_modelo, nrow(datos_modelo)
    ),
    x = NULL,
    y = "Desviación estándar condicional de r_t",
    caption = sprintf(
      paste0(
        "vxreg1 = %.6f (%s); p %s: el miedo %s. ",
        "AIC = %.2f. Máximo: %.3f el %s."
      ),
      mejor_fila$vxreg1_Estimate,
      "coeficiente del miedo en la varianza",
      formatear_p(mejor_fila$vxreg1_p_value),
      lectura_mejor,
      mejor_fila$AIC,
      valor_maximo,
      format(fecha_maximo, "%Y-%m-%d")
    )
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7.5, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_png,
  plot = p_volatilidad,
  width = 10.5,
  height = 5.2,
  dpi = 200,
  bg = "white"
)

# --- 8. Guardado ---------------------------------------------------

# La tabla final conserva exactamente las cuatro columnas pedidas.
resultados_csv <- dplyr::select(
  resultados,
  Modelo, AIC, vxreg1_Estimate, vxreg1_p_value
)

readr::write_csv(resultados_csv, ruta_resultados)
readr::write_csv(serie_volatilidad, ruta_volatilidad)

rutas_salida <- c(ruta_resultados, ruta_volatilidad, ruta_png)
for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 9. Resumen en consola -----------------------------------------

cat("\n============================================================\n")
cat(" TABLA 4 — GARCH CON EL MIEDO EN LA VARIANZA\n")
cat("============================================================\n")
cat(sprintf(
  "Muestra: %d semanas (%s a %s) | %d meses | M_dengue mensual\n",
  nrow(datos_modelo),
  format(min(datos_modelo$fecha_sem), "%Y-%m-%d"),
  format(max(datos_modelo$fecha_sem), "%Y-%m-%d"),
  dplyr::n_distinct(datos_modelo$mes)
))
cat("Media: r_t = mu + phi r_{t-1} + e_t | distribución normal\n")
cat("Varianza: regresor externo M_dengue (vxreg1) en los tres modelos\n")

cat("\n--- COEFICIENTE DEL MIEDO (vxreg1) Y AIC ---\n")
for (i in seq_len(nrow(resultados))) {
  cat(sprintf(
    "%-10s | vxreg1 = % .6e | p %-12s | AIC = %12.2f | %s\n",
    resultados$Modelo[i],
    resultados$vxreg1_Estimate[i],
    formatear_p(resultados$vxreg1_p_value[i]),
    resultados$AIC[i],
    interpretar_miedo(
      resultados$vxreg1_Estimate[i],
      resultados$vxreg1_p_value[i]
    )
  ))
}

cat(sprintf(
  "\nAIC por observación (rugarch): %s\n",
  paste(
    sprintf("%s = %.4f", resultados$Modelo, resultados$AIC_por_observacion),
    collapse = " | "
  )
))
cat(sprintf(
  "Mejor modelo por AIC: %s (AIC total = %.2f)\n",
  mejor_modelo, mejor_fila$AIC
))
cat(sprintf(
  "vxreg1 del mejor modelo: %.6f (EE %.3e; p %s)\n",
  mejor_fila$vxreg1_Estimate, mejor_fila$vxreg1_Std_Error,
  formatear_p(mejor_fila$vxreg1_p_value)
))
cat(sprintf(
  "Lectura: el miedo %s.\n",
  lectura_mejor
))

cat("\n--- VOLATILIDAD CONDICIONAL ---\n")
cat(sprintf(
  "Rango: %.4f a %.4f | media %.4f | máximo el %s\n",
  min(serie_volatilidad$volatilidad),
  max(serie_volatilidad$volatilidad),
  mean(serie_volatilidad$volatilidad),
  format(fecha_maximo, "%Y-%m-%d")
))

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

cat("\nNotas:\n")
cat("  - r_t = diff(log(casos + 1)); se pierde la primera semana y se\n")
cat("    recorta esa fila para alinear M_t y r_t uno a uno.\n")
cat("  - M_dengue es el promedio mensual de las 12 regiones EMA.\n")
cat("  - El signo de vxreg1 importa: positivo aumenta la volatilidad;\n")
cat("    negativo la reduce. El p-valor es el de H0: vxreg1 = 0.\n")
cat("  - La significancia de vxreg1 no implica un efecto económico\n")
cat("    grande: compare su magnitud con la escala de sigma^2_t.\n")
cat("  - El AIC guardado es el total (2k - 2 log L); el de rugarch es\n")
cat("    por observación, pero el orden entre modelos es el mismo.\n")

message("\n✅ Script 17 finalizado.")
