# ============================================================
# 15_panel_descomposicion.R
# Panel con efectos fijos: incidencia, miedo excedente y
# ocupación hotelera. Participación de la atención (pi).
# ============================================================
#
# OBJETIVO (plan, §4.4 «Vía principal: panel con efectos fijos»)
#
#   Y_rt = a_r + g_t + s_{r,m(t)}
#          + sum_{k=1..3} b_k I_{r,t-k}
#          + sum_{k=1..3} theta_k M_{r,t-k}
#          + delta' C_{r,t-1} + e_rt
#
#   Y  ocupación hotelera (EMA, porcentaje)
#   I  incidencia (casos por 100.000 habitantes)
#   M  miedo excedente (residuo de la regresión auxiliar del script 12)
#   C  clima local (lluvia y temperatura), rezagado un mes
#   g_t efectos fijos de año-mes; s_{r,m(t)} estacionalidad por región
#
# El modelo usa solo rezagos (t-1 a t-3) para evitar la causalidad
# inversa: el turismo también puede traer casos y aumentar el reporte.
#
# La participación de la atención en el pico de 2024 se calcula con las
# contribuciones ajustadas de la página 7 del plan:
#
#   pi_atencion = |theta_hat * Delta_M|
#                 / (|b_hat * Delta_I| + |theta_hat * Delta_M|)
#
#   b_hat     = sum_k b_k       (efecto de la enfermedad)
#   theta_hat = sum_k theta_k   (efecto del miedo)
#   Delta_X   = media de X en el pico de 2024
#               menos media de X en los mismos meses de 2022-2023
#
# El pico de 2024 se define como el segundo trimestre (abril-junio),
# que el plan identifica como el pico histórico de dengue (más de
# 7.500 casos semanales). Los «mismos meses» son abril-junio de 2022
# y 2023.
#
# INSUMOS
#   data/D_13/panel_ema_den_chik.rds
#
# SALIDAS
#   data/D_15/modelos_panel.rds
#   data/D_15/tabla_coeficientes.csv
#   data/D_15/participacion_atencion.csv
#   graficas/G_15/coefplot_panel.png
#
# NOTA SOBRE LA INCIDENCIA
#   La Sección 4.4 define I como casos por 100.000 habitantes. El panel
#   ya trae I_dengue e I_chik (calculadas en el script 13 a partir de
#   casos_dengue / casos_chik y la población regional). Se usan esas
#   columnas, no los casos crudos, para que b_k sea el efecto de un
#   punto de incidencia y b_hat * Delta_I quede en puntos porcentuales
#   de ocupación, igual que theta_hat * Delta_M.
#
# NOTA SOBRE fixest
#   La vía principal usa fixest::feols() con vcov = "DK"
#   (Driscoll-Kraay) y, como contraste, errores agrupados por región.
#   Si fixest no está instalado, el script estima el mismo modelo con
#   lm() y calcula los dos tipos de errores estándar en base R. Esa
#   rama de respaldo permite reproducir las salidas en entornos donde
#   fixest no está disponible; cuando fixest está instalado, se usa
#   fixest y los objetos guardados son de clase fixest.
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# --- 0. Rutas y carpetas -------------------------------------------

dir.create(here("data", "D_15"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("graficas", "G_15"), recursive = TRUE, showWarnings = FALSE)

ruta_panel <- here("data", "D_13", "panel_ema_den_chik.rds")
ruta_modelos <- here("data", "D_15", "modelos_panel.rds")
ruta_tabla <- here("data", "D_15", "tabla_coeficientes.csv")
ruta_participacion <- here("data", "D_15", "participacion_atencion.csv")
ruta_coefplot <- here("graficas", "G_15", "coefplot_panel.png")

# --- 1. Motor de estimación ----------------------------------------
#
# fixest es la vía principal. El respaldo con lm() solo se activa si el
# paquete no está instalado; así el script no se detiene en un entorno
# sin fixest y las cifras se pueden contrastar.

motor <- if (requireNamespace("fixest", quietly = TRUE)) {
  "fixest"
} else {
  "lm"
}

if (motor == "fixest") {
  suppressPackageStartupMessages(library(fixest))
  message("Motor de estimación: fixest::feols().")
} else {
  message(
    "fixest no está instalado: se usa lm() con efectos fijos y errores ",
    "Driscoll-Kraay/agrupados calculados en base R. Instale fixest con ",
    "renv::install(\"fixest\") para la vía principal."
  )
}

# --- 2. Carga y validación del panel -------------------------------

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel del script 13: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

panel <- readRDS(ruta_panel)

requeridas <- c(
  "region_ema", "mes", "anio", "mes_num", "anio_mes", "Y",
  "I_dengue", "I_dengue_lag1", "I_dengue_lag2", "I_dengue_lag3",
  "M_dengue", "M_dengue_lag1", "M_dengue_lag2", "M_dengue_lag3",
  "I_chik", "I_chik_lag1", "I_chik_lag2", "I_chik_lag3",
  "M_chik", "M_chik_lag1", "M_chik_lag2", "M_chik_lag3",
  "lluvia_lag1", "temp_lag1"
)

faltantes <- setdiff(requeridas, names(panel))
if (length(faltantes) > 0) {
  stop(
    "Faltan columnas en el panel: ", paste(faltantes, collapse = ", "), ".",
    call. = FALSE
  )
}

panel <- panel |>
  mutate(
    region_ema = as.character(region_ema),
    mes = as.Date(mes),
    anio = as.integer(anio),
    mes_num = as.integer(mes_num),
    anio_mes = as.character(anio_mes),
    # Índice temporal 1..48 para panel.id de fixest (más seguro que Date).
    t = as.integer(factor(mes))
  ) |>
  arrange(region_ema, mes)

if (nrow(panel) != 576) {
  stop(
    "El panel debe tener 576 filas (12 regiones x 48 meses) y tiene ",
    nrow(panel), ".", call. = FALSE
  )
}

if (anyDuplicated(panel[, c("region_ema", "mes")])) {
  stop("Hay regiones-mes duplicadas en el panel.", call. = FALSE)
}

if (n_distinct(panel$region_ema) != 12) {
  stop("El panel debe tener 12 regiones EMA.", call. = FALSE)
}

if (n_distinct(panel$mes) != 48) {
  stop("El panel debe tener 48 meses (2022-01 a 2025-12).", call. = FALSE)
}

variables_modelo <- c(
  "Y",
  "I_dengue_lag1", "I_dengue_lag2", "I_dengue_lag3",
  "M_dengue_lag1", "M_dengue_lag2", "M_dengue_lag3",
  "I_chik_lag1", "I_chik_lag2", "I_chik_lag3",
  "M_chik_lag1", "M_chik_lag2", "M_chik_lag3",
  "lluvia_lag1", "temp_lag1"
)

na_modelo <- colSums(is.na(panel[, variables_modelo]))
if (any(na_modelo > 0)) {
  stop(
    "Hay valores faltantes en las variables del modelo: ",
    paste(names(na_modelo)[na_modelo > 0], na_modelo[na_modelo > 0],
          sep = " = ", collapse = ", "),
    ".", call. = FALSE
  )
}

periodo_panel <- sprintf(
  "%s a %s", format(min(panel$mes), "%Y-%m"),
  format(max(panel$mes), "%Y-%m")
)

message(
  "Panel cargado: ", nrow(panel), " filas | ",
  n_distinct(panel$region_ema), " regiones | ",
  n_distinct(panel$mes), " meses (", periodo_panel, ")."
)

# --- 3. Especificación del modelo ----------------------------------

# Términos que entran a la tabla y al coefplot, en el orden del plan:
# enfermedad (3 rezagos), miedo (3 rezagos) y clima (2 controles).
terminos_dengue <- c(
  "I_dengue_lag1", "I_dengue_lag2", "I_dengue_lag3",
  "M_dengue_lag1", "M_dengue_lag2", "M_dengue_lag3",
  "lluvia_lag1", "temp_lag1"
)

terminos_chik <- c(
  "I_chik_lag1", "I_chik_lag2", "I_chik_lag3",
  "M_chik_lag1", "M_chik_lag2", "M_chik_lag3",
  "lluvia_lag1", "temp_lag1"
)

# El coefplot se concentra en los canales de interés (enfermedad y
# miedo); los controles de clima quedan en la tabla, no en la figura.
terminos_plot_dengue <- terminos_dengue[1:6]
terminos_plot_chik <- terminos_chik[1:6]

# Fórmula fixest: los efectos fijos van después de |.
#   region_ema            efecto fijo de región
#   anio_mes              efecto fijo de año-mes (g_t)
#   region_ema^mes_num    estacionalidad propia de cada región (s_r,m)
fml_dengue_fixest <- Y ~
  I_dengue_lag1 + I_dengue_lag2 + I_dengue_lag3 +
  M_dengue_lag1 + M_dengue_lag2 + M_dengue_lag3 +
  lluvia_lag1 + temp_lag1 |
  region_ema + anio_mes + region_ema^mes_num

fml_chik_fixest <- Y ~
  I_chik_lag1 + I_chik_lag2 + I_chik_lag3 +
  M_chik_lag1 + M_chik_lag2 + M_chik_lag3 +
  lluvia_lag1 + temp_lag1 |
  region_ema + anio_mes + region_ema^mes_num

# Fórmula equivalente para lm(): los efectos fijos se expanden como
# factores. La interacción región x mes_num reproduce s_r,m; la región
# principal queda absorbida por la interacción y lm() la marca aliased.
fml_dengue_lm <- Y ~
  I_dengue_lag1 + I_dengue_lag2 + I_dengue_lag3 +
  M_dengue_lag1 + M_dengue_lag2 + M_dengue_lag3 +
  lluvia_lag1 + temp_lag1 +
  factor(region_ema) + factor(anio_mes) +
  factor(region_ema):factor(mes_num)

fml_chik_lm <- Y ~
  I_chik_lag1 + I_chik_lag2 + I_chik_lag3 +
  M_chik_lag1 + M_chik_lag2 + M_chik_lag3 +
  lluvia_lag1 + temp_lag1 +
  factor(region_ema) + factor(anio_mes) +
  factor(region_ema):factor(mes_num)

# --- 4. Errores estándar de respaldo (solo si no hay fixest) -------

# Driscoll-Kraay: Newey-West sobre los momentos agregados por mes.
# X es la matriz de diseño (sin columnas aliased), e los residuos y
# h_t = sum_i x_it e_it. La matriz de covarianzas es
#   (X'X)^-1 [ sum_t h_t h_t'
#              + sum_l w_l sum_t (h_t h_{t-l}' + h_{t-l} h_t') ] (X'X)^-1
# con pesos de Bartlett w_l = 1 - l/(L+1) y L = floor(4*(T/100)^(2/9)).
vcov_dk <- function(modelo, tiempo, lag = NULL) {
  X <- stats::model.matrix(modelo)
  e <- stats::residuals(modelo)

  # Columnas que lm() eliminó por colinealidad.
  keep <- !is.na(stats::coef(modelo))
  X <- X[, keep, drop = FALSE]

  n <- nrow(X)
  k <- ncol(X)
  XtX_inv <- solve(crossprod(X))

  t_idx <- match(tiempo, sort(unique(tiempo)))
  h <- rowsum(X * e, t_idx, reorder = TRUE)
  Tt <- nrow(h)

  if (is.null(lag)) {
    lag <- floor(4 * (Tt / 100)^(2 / 9))
  }

  S <- crossprod(h)
  if (lag > 0) {
    for (l in seq_len(lag)) {
      w <- 1 - l / (lag + 1)
      G <- crossprod(
        h[(l + 1):Tt, , drop = FALSE],
        h[1:(Tt - l), , drop = FALSE]
      )
      S <- S + w * (G + t(G))
    }
  }

  V <- XtX_inv %*% S %*% XtX_inv

  # Corrección de muestra pequeña, análoga a la de fixest.
  V * n / (n - k)
}

# Errores agrupados por región. Se usa sandwich si está instalado; si
# no, se calcula el clásico cluster-robust CR1.
vcov_cluster <- function(modelo, cluster) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    return(sandwich::vcovCL(modelo, cluster = cluster, type = "HC1"))
  }

  X <- stats::model.matrix(modelo)
  e <- stats::residuals(modelo)
  keep <- !is.na(stats::coef(modelo))
  X <- X[, keep, drop = FALSE]

  n <- nrow(X)
  k <- ncol(X)
  G <- length(unique(cluster))
  XtX_inv <- solve(crossprod(X))
  meat <- crossprod(rowsum(X * e, cluster, reorder = TRUE))

  V <- XtX_inv %*% meat %*% XtX_inv
  V * (G / (G - 1)) * ((n - 1) / (n - k))
}

# --- 5. Estimación -------------------------------------------------

# Estima un modelo con la fórmula fixest o con lm(), según el motor, y
# devuelve un objeto con coeficientes, covarianza, p-valores y bondad
# de ajuste. Los modelos fixest se guardan tal cual; los de respaldo
# llevan clase feols_fallback.
estimar_panel <- function(fml_fixest, fml_lm, vcov_tipo) {
  if (motor == "fixest") {
    if (vcov_tipo == "DK") {
      modelo <- feols(
        fml_fixest,
        data = panel,
        panel.id = ~ region_ema + t,
        vcov = "DK"
      )
    } else {
      modelo <- feols(
        fml_fixest,
        data = panel,
        panel.id = ~ region_ema + t,
        vcov = ~ region_ema
      )
    }

    return(list(
      modelo = modelo,
      clase = "fixest",
      n_obs = as.integer(nobs(modelo)),
      r2 = as.numeric(r2(modelo, "r2")),
      df = NULL
    ))
  }

  modelo <- stats::lm(fml_lm, data = panel)
  idx <- as.integer(rownames(stats::model.frame(modelo)))
  tiempo <- panel$mes[idx]
  cluster <- panel$region_ema[idx]

  if (vcov_tipo == "DK") {
    V <- vcov_dk(modelo, tiempo)
    df <- length(unique(tiempo)) - 1L
  } else {
    V <- vcov_cluster(modelo, cluster)
    df <- length(unique(cluster)) - 1L
  }

  cf <- stats::coef(modelo)
  se_keep <- sqrt(diag(V))
  names(se_keep) <- colnames(V)

  # cf incluye los coeficientes aliased (NA); se conserva su longitud
  # para que los nombres sigan alineados con la matriz de diseño.
  se <- rep(NA_real_, length(cf))
  names(se) <- names(cf)
  se[names(se_keep)] <- se_keep

  p <- rep(NA_real_, length(cf))
  names(p) <- names(cf)
  p[names(se_keep)] <- 2 * stats::pt(
    -abs(cf[names(se_keep)] / se_keep), df = df
  )

  fallback <- list(
    coefficients = cf,
    vcov = V,
    se = se,
    pvalue = p,
    nobs = stats::nobs(modelo),
    r2 = summary(modelo)$r.squared,
    df = df,
    formula = stats::formula(modelo),
    call = match.call()
  )
  class(fallback) <- c("feols_fallback", "list")

  list(
    modelo = fallback,
    clase = "lm_fallback",
    n_obs = as.integer(stats::nobs(modelo)),
    r2 = summary(modelo)$r.squared,
    df = df
  )
}

message("\nEstimando el modelo principal (dengue) y el placebo (chikungunya)...")

m_dengue_dk <- estimar_panel(fml_dengue_fixest, fml_dengue_lm, "DK")
m_dengue_cl <- estimar_panel(fml_dengue_fixest, fml_dengue_lm, "cluster")
m_chik_dk <- estimar_panel(fml_chik_fixest, fml_chik_lm, "DK")
m_chik_cl <- estimar_panel(fml_chik_fixest, fml_chik_lm, "cluster")

# --- 6. Extracción de coeficientes ---------------------------------

# Devuelve una fila por término con estimación, error estándar y
# p-valor del modelo recibido.
extraer_coefs <- function(fit, terminos, modelo) {
  if (motor == "fixest") {
    cf <- stats::coef(fit)
    V <- stats::vcov(fit)
    p <- fixest::pvalue(fit)
  } else {
    cf <- fit$coefficients
    V <- fit$vcov
    p <- fit$pvalue
  }

  se <- if (motor == "fixest") sqrt(diag(V)) else fit$se
  if (motor == "fixest") names(se) <- names(cf)

  faltan <- setdiff(terminos, names(cf))
  if (length(faltan) > 0) {
    stop(
      "El modelo ", modelo, " no estimó los términos: ",
      paste(faltan, collapse = ", "), ".",
      call. = FALSE
    )
  }

  tibble(
    modelo = modelo,
    termino = terminos,
    estimacion = as.numeric(cf[terminos]),
    se = as.numeric(se[terminos]),
    p = as.numeric(p[terminos])
  )
}

# Combina la estimación puntual con los dos tipos de error estándar.
construir_tabla <- function(fit_dk, fit_cl, terminos, modelo) {
  a <- extraer_coefs(fit_dk, terminos, modelo)
  b <- extraer_coefs(fit_cl, terminos, modelo)
  
  a |>
    left_join(
      b |> dplyr::select(modelo, termino, se_cluster = se, p_cluster = p),
      by = c("modelo", "termino")
    ) |>
    rename(se_dk = se, p_dk = p) |>
    mutate(
      canal = case_when(
        str_detect(termino, "^I_") ~ "Enfermedad",
        str_detect(termino, "^M_") ~ "Miedo",
        TRUE ~ "Clima"
      ),
      sig_dk = case_when(
        is.na(p_dk) ~ "",
        p_dk < 0.01 ~ "***",
        p_dk < 0.05 ~ "**",
        p_dk < 0.10 ~ "*",
        TRUE ~ ""
      ),
      sig_cluster = case_when(
        is.na(p_cluster) ~ "",
        p_cluster < 0.01 ~ "***",
        p_cluster < 0.05 ~ "**",
        p_cluster < 0.10 ~ "*",
        TRUE ~ ""
      )
    )
}

tabla_dengue <- construir_tabla(
  m_dengue_dk$modelo, m_dengue_cl$modelo, terminos_dengue,
  "Dengue (principal)"
)

tabla_chik <- construir_tabla(
  m_chik_dk$modelo, m_chik_cl$modelo, terminos_chik,
  "Chikungunya (placebo)"
)

tabla_coeficientes <- bind_rows(tabla_dengue, tabla_chik) |>
  mutate(
    modelo = factor(
      modelo,
      levels = c("Dengue (principal)", "Chikungunya (placebo)")
    ),
    canal = factor(canal, levels = c("Enfermedad", "Miedo", "Clima")),
    termino = factor(
      termino,
      levels = unique(c(terminos_dengue, terminos_chik))
    )
  ) |>
  arrange(modelo, canal, termino) |>
  mutate(
    modelo = as.character(modelo),
    canal = as.character(canal),
    termino = as.character(termino)
  )

# Metadatos del modelo en cada fila, para que la tabla sea autónoma.
tabla_coeficientes <- tabla_coeficientes |>
  mutate(
    n_obs = if_else(
      modelo == "Dengue (principal)",
      m_dengue_dk$n_obs, m_chik_dk$n_obs
    ),
    r2 = if_else(
      modelo == "Dengue (principal)",
      m_dengue_dk$r2, m_chik_dk$r2
    )
  ) |>
  dplyr::select(
    modelo, canal, termino, estimacion,
    se_dk, p_dk, sig_dk, se_cluster, p_cluster, sig_cluster,
    n_obs, r2
  )

write_csv(tabla_coeficientes, ruta_tabla)

# --- 7. Participación de la atención (pi) --------------------------

# Pico de 2024: segundo trimestre (abril, mayo y junio), el pico
# histórico de dengue que documenta el plan. La base son los mismos
# meses de 2022 y 2023.
meses_pico <- 4:6
anio_pico <- 2024L
anios_base <- c(2022L, 2023L)

calcular_participacion <- function(fit, enfermedad, modelo) {
  terminos_i <- paste0("I_", enfermedad, "_lag", 1:3)
  terminos_m <- paste0("M_", enfermedad, "_lag", 1:3)

  if (motor == "fixest") {
    cf <- stats::coef(fit)
  } else {
    cf <- fit$coefficients
  }

  b_hat <- sum(cf[terminos_i])
  theta_hat <- sum(cf[terminos_m])

  pico <- panel |>
    filter(anio == anio_pico, mes_num %in% meses_pico)

  base <- panel |>
    filter(anio %in% anios_base, mes_num %in% meses_pico)

  var_i <- paste0("I_", enfermedad)
  var_m <- paste0("M_", enfermedad)

  delta_i <- mean(pico[[var_i]]) - mean(base[[var_i]])
  delta_m <- mean(pico[[var_m]]) - mean(base[[var_m]])

  contrib_enfermedad <- b_hat * delta_i
  contrib_miedo <- theta_hat * delta_m
  denominador <- abs(contrib_enfermedad) + abs(contrib_miedo)

  pi_atencion <- if (is.finite(denominador) && denominador > 0) {
    abs(contrib_miedo) / denominador
  } else {
    NA_real_
  }

  tibble(
    modelo = modelo,
    enfermedad = enfermedad,
    b_hat = b_hat,
    theta_hat = theta_hat,
    delta_I = delta_i,
    delta_M = delta_m,
    contribucion_enfermedad = contrib_enfermedad,
    contribucion_miedo = contrib_miedo,
    pi_atencion = pi_atencion,
    pi_atencion_pct = 100 * pi_atencion,
    pico = "2024-Q2 (abril-junio)",
    base = "2022-Q2 y 2023-Q2 (abril-junio)",
    n_obs = if (motor == "fixest") as.integer(nobs(fit)) else fit$nobs,
    r2 = if (motor == "fixest") {
      as.numeric(r2(fit, "r2"))
    } else {
      fit$r2
    }
  )
}

participacion <- calcular_participacion(
  m_dengue_dk$modelo, "dengue", "Dengue (principal)"
)

write_csv(participacion, ruta_participacion)

# --- 8. Modelos estimados ------------------------------------------

modelos_panel <- list(
  dengue_dk = m_dengue_dk$modelo,
  dengue_cluster = m_dengue_cl$modelo,
  chik_dk = m_chik_dk$modelo,
  chik_cluster = m_chik_cl$modelo
)

attr(modelos_panel, "motor") <- motor
attr(modelos_panel, "periodo") <- periodo_panel
attr(modelos_panel, "formula_dengue") <- fml_dengue_fixest
attr(modelos_panel, "formula_chik") <- fml_chik_fixest
attr(modelos_panel, "nota") <- paste(
  "Modelo principal: Y ~ I(t-1..t-3) + M(t-1..t-3) + lluvia(t-1) +",
  "temp(t-1) | region + anio_mes + region^mes_num.",
  "dengue_dk y chik_dk usan Driscoll-Kraay;",
  "dengue_cluster y chik_cluster usan errores agrupados por región."
)

saveRDS(modelos_panel, ruta_modelos)

# --- 9. Coefplot ----------------------------------------------------

# Etiquetas legibles para el gráfico.
etiquetas_terminos <- c(
  I_dengue_lag1 = "I dengue (t-1)",
  I_dengue_lag2 = "I dengue (t-2)",
  I_dengue_lag3 = "I dengue (t-3)",
  M_dengue_lag1 = "M miedo (t-1)",
  M_dengue_lag2 = "M miedo (t-2)",
  M_dengue_lag3 = "M miedo (t-3)",
  I_chik_lag1 = "I chik (t-1)",
  I_chik_lag2 = "I chik (t-2)",
  I_chik_lag3 = "I chik (t-3)",
  M_chik_lag1 = "M miedo chik (t-1)",
  M_chik_lag2 = "M miedo chik (t-2)",
  M_chik_lag3 = "M miedo chik (t-3)",
  lluvia_lag1 = "Lluvia (t-1)",
  temp_lag1 = "Temperatura (t-1)"
)

# Vía preferida: coefplot de fixest, con dos paneles (dengue y placebo).
coefplot_ok <- FALSE

if (motor == "fixest") {
  coefplot_ok <- tryCatch(
    {
      grDevices::png(
        filename = ruta_coefplot,
        width = 1800, height = 1500, res = 180
      )
      op <- graphics::par(
        mfrow = c(2, 1),
        mar = c(4.2, 10.5, 3.2, 1.2)
      )

      coefplot(
        m_dengue_dk$modelo,
        keep = terminos_plot_dengue,
        main = "Dengue: enfermedad y miedo (rezagos 1-3)",
        xlab = "Coeficiente (puntos porcentuales de ocupación)",
        ci_level = 0.95,
        horiz = TRUE
      )
      coefplot(
        m_chik_dk$modelo,
        keep = terminos_plot_chik,
        main = "Chikungunya (placebo): enfermedad y miedo",
        xlab = "Coeficiente (puntos porcentuales de ocupación)",
        ci_level = 0.95,
        horiz = TRUE
      )

      graphics::par(op)
      grDevices::dev.off()
      TRUE
    },
    error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      message(
        "coefplot() de fixest no se pudo dibujar (", conditionMessage(e),
        "); se usa el gráfico ggplot equivalente."
      )
      FALSE
    }
  )
}

# Respaldo y también gráfico equivalente si coefplot falla: ggplot con
# intervalos al 95% calculados con los errores Driscoll-Kraay.
if (!coefplot_ok) {
  datos_plot <- bind_rows(
    extraer_coefs(
      m_dengue_dk$modelo, terminos_plot_dengue, "Dengue (principal)"
    ),
    extraer_coefs(
      m_chik_dk$modelo, terminos_plot_chik, "Chikungunya (placebo)"
    )
  ) |>
    mutate(
      canal = case_when(
        str_detect(termino, "^I_") ~ "Enfermedad",
        str_detect(termino, "^M_") ~ "Miedo",
        TRUE ~ "Clima"
      ),
      etiqueta = unname(etiquetas_terminos[termino]),
      inferior = estimacion - 1.96 * se,
      superior = estimacion + 1.96 * se
    )

  p_coefplot <- ggplot(
    datos_plot,
    aes(x = estimacion, y = etiqueta, colour = canal)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_pointrange(
      aes(xmin = inferior, xmax = superior),
      size = 0.45, linewidth = 0.55
    ) +
    facet_wrap(~modelo, scales = "free", ncol = 1) +
    scale_colour_manual(
      values = c(
        "Enfermedad" = "#08519C",
        "Miedo" = "#B2182B",
        "Clima" = "#6E6E6E"
      ),
      name = NULL
    ) +
    labs(
      title = "Coeficientes del panel con efectos fijos",
      subtitle = paste0(
        "Y = ocupación hotelera; I = incidencia por 100.000 hab.; ",
        "M = miedo excedente. Intervalos al 95% con errores Driscoll-Kraay."
      ),
      x = "Coeficiente (puntos porcentuales de ocupación)",
      y = NULL,
      caption = paste0(
        "Controles: lluvia y temperatura rezagadas un mes; ",
        "efectos fijos de región, año-mes y región x mes."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, size = 7, colour = "grey35")
    )

  ggsave(
    ruta_coefplot, p_coefplot,
    width = 10, height = 7.5, dpi = 180, bg = "white"
  )
}

if (!file.exists(ruta_coefplot) || file.info(ruta_coefplot)$size == 0) {
  stop("No se generó el coefplot: ", ruta_coefplot, call. = FALSE)
}

# --- 10. Resumen en consola ----------------------------------------

# Impresión compacta de la tabla de coeficientes cuando no hay fixest.
imprimir_tabla <- function(tabla) {
  for (i in seq_len(nrow(tabla))) {
    cat(sprintf(
      "  %-18s %10.5f  %8.5f  %7.4f  %8.5f  %7.4f\n",
      tabla$termino[i],
      tabla$estimacion[i],
      tabla$se_dk[i], tabla$p_dk[i],
      tabla$se_cluster[i], tabla$p_cluster[i]
    ))
  }
}

cat("\n============================================================\n")
cat(" PANEL CON EFECTOS FIJOS: DENGUE, MIEDO Y OCUPACIÓN\n")
cat("============================================================\n")
cat(sprintf(
  "Panel: %d filas | %d regiones | %d meses (%s)\n",
  nrow(panel), n_distinct(panel$region_ema),
  n_distinct(panel$mes), periodo_panel
))
cat(sprintf("Motor: %s\n", if (motor == "fixest") "fixest" else "lm (respaldo)"))
cat("Efectos fijos: región + año-mes + región x mes\n")
cat("Errores: Driscoll-Kraay (principal) y agrupados por región (contraste)\n")

cat("\n--- MODELO PRINCIPAL: DENGUE ---\n")
if (motor == "fixest") {
  print(summary(m_dengue_dk$modelo))
} else {
  cat(sprintf(
    "  %-18s %10s  %8s  %7s  %8s  %7s\n",
    "Término", "Estimación", "SE DK", "p DK", "SE cluster", "p cluster"
  ))
  imprimir_tabla(tabla_dengue)
  cat(sprintf(
    "\n  N = %d | R2 = %.4f\n",
    m_dengue_dk$n_obs, m_dengue_dk$r2
  ))
}

cat("\n--- MODELO PLACEBO: CHIKUNGUNYA ---\n")
if (motor == "fixest") {
  print(summary(m_chik_dk$modelo))
} else {
  cat(sprintf(
    "  %-18s %10s  %8s  %7s  %8s  %7s\n",
    "Término", "Estimación", "SE DK", "p DK", "SE cluster", "p cluster"
  ))
  imprimir_tabla(tabla_chik)
  cat(sprintf(
    "\n  N = %d | R2 = %.4f\n",
    m_chik_dk$n_obs, m_chik_dk$r2
  ))
}

cat("\n--- PARTICIPACIÓN DE LA ATENCIÓN (pi) ---\n")
cat(sprintf("Pico: %s; base: %s\n", participacion$pico, participacion$base))
cat(sprintf("b_hat (suma I)        = %10.6f\n", participacion$b_hat))
cat(sprintf("theta_hat (suma M)    = %10.6f\n", participacion$theta_hat))
cat(sprintf("Delta_I               = %10.6f\n", participacion$delta_I))
cat(sprintf("Delta_M               = %10.6f\n", participacion$delta_M))
cat(sprintf(
  "Contribución enfermedad = %10.6f puntos porcentuales\n",
  participacion$contribucion_enfermedad
))
cat(sprintf(
  "Contribución miedo      = %10.6f puntos porcentuales\n",
  participacion$contribucion_miedo
))
cat(sprintf(
  "\npi_atencion = %.6f  (%.2f%% de la caída atribuible al miedo)\n",
  participacion$pi_atencion, participacion$pi_atencion_pct
))

cat("\nSalidas:\n")
cat("  ", ruta_modelos, "\n", sep = "")
cat("  ", ruta_tabla, "\n", sep = "")
cat("  ", ruta_participacion, "\n", sep = "")
cat("  ", ruta_coefplot, "\n", sep = "")

cat("\nNotas:\n")
cat("  - La incidencia I se mide por 100.000 habitantes (I_dengue, I_chik).\n")
cat("  - El pico de 2024 son abril-junio; la base son abril-junio de 2022-2023.\n")
cat("  - La tabla incluye SE Driscoll-Kraay y agrupados por región.\n")

message("\n✅ Script 15 finalizado.")
