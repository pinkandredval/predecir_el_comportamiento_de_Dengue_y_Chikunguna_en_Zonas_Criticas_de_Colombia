# ============================================================
# 08_econometria.R
# Análisis econométrico de series temporales (Paper 2)
# ADF, Box-Cox, diferenciación, ACF/PACF, SARIMA, validación
#
# Correcciones Tabla 2 (#4, #5, #6):
#   #4 auto.arima(lambda = "auto", biasadj = TRUE) y MASE con
#      accuracy() en una ventana reservada (últimos 12 meses).
#   #5 Se sustituye Breusch-Pagan por McLeod-Li (Ljung-Box sobre
#      residuos al cuadrado), que detecta heterocedasticidad ARCH.
#   #6 Se elimina Durbin-Watson. Ljung-Box con
#      fitdf = p + q + P + Q y lag = 24 (mensual).
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(tseries)
  library(forecast)
  library(scales)
})

serie_dengue <- readRDS(here("data", "D_07", "serie_dengue.rds"))
serie_chik   <- readRDS(here("data", "D_07", "serie_chik.rds"))

dir.create(here("data", "D_08"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("graficas", "G_08"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Serie mensual nacional (rejilla completa, 72 meses)
# ============================================================

preparar_mensual <- function(serie) {
  serie %>%
    filter(dpto == "NACIONAL") %>%
    mutate(mes = floor_date(fecha_sem + 3, "month")) %>%  # jueves de la semana (ISO)
    group_by(mes) %>%
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") %>%
    complete(mes = seq(min(mes), max(mes), by = "month"),
             fill = list(casos = 0)) %>%
    arrange(mes)
}

mensual_dengue <- preparar_mensual(serie_dengue)
mensual_chik   <- preparar_mensual(serie_chik)

stopifnot(nrow(mensual_dengue) == 72L, nrow(mensual_chik) == 72L)

ts_dengue <- ts(mensual_dengue$casos,
                start = c(year(min(mensual_dengue$mes)),
                          month(min(mensual_dengue$mes))),
                frequency = 12)

ts_chik <- ts(mensual_chik$casos,
              start = c(year(min(mensual_chik$mes)),
                        month(min(mensual_chik$mes))),
              frequency = 12)

cat("\n===== SERIE MENSUAL =====\n")
cat("Dengue:      ", length(ts_dengue), "obs |",
    format(min(mensual_dengue$mes)), "→", format(max(mensual_dengue$mes)), "\n")
cat("Chikungunya: ", length(ts_chik),   "obs |",
    format(min(mensual_chik$mes)), "→", format(max(mensual_chik$mes)), "\n")

# ============================================================
# 2. Utilidades
# ============================================================

# MAPE calculado solo sobre meses con casos > 0. Así se evita el Inf
# que produce la división por cero en chikungunya.
mape_positivo <- function(actual, pred) {
  actual <- as.numeric(actual)
  pred   <- as.numeric(pred)
  ok     <- is.finite(actual) & is.finite(pred) & actual > 0
  if (!any(ok)) return(NA_real_)
  mean(abs((actual[ok] - pred[ok]) / actual[ok])) * 100
}

# auto.arima con lambda="auto" y biasadj=TRUE. Box-Cox exige datos
# estrictamente positivos; chikungunya tiene 4 meses en cero, así que
# en ese caso se ajusta sin transformación y se documenta.
ajustar_sarima <- function(y, nombre) {
  if (all(y > 0)) {
    cat("Box-Cox [", nombre, "]: lambda = 'auto', biasadj = TRUE\n")
    fit <- auto.arima(y, seasonal = TRUE, stepwise = TRUE,
                      approximation = FALSE,
                      lambda = "auto", biasadj = TRUE)
    cat("Lambda estimado:", round(fit$lambda, 4), "\n")
  } else {
    cat("Box-Cox [", nombre, "]: no aplicable (",
        sum(y == 0), "ceros). Se ajusta sin transformación;\n",
        "  lambda = 'auto' no está definido para datos no positivos.\n")
    fit <- auto.arima(y, seasonal = TRUE, stepwise = TRUE,
                      approximation = FALSE,
                      lambda = NULL, biasadj = TRUE)
  }
  fit
}

# ============================================================
# 3. Función: análisis econométrico completo
# ============================================================

analisis_econometrico <- function(ts_serie, nombre, color = "#08519C") {

  n <- length(ts_serie)

  cat("\n=========================================================\n")
  cat("ANÁLISIS ECONOMÉTRICO:", nombre, "\n")
  cat("=========================================================\n")

  # --- ADF serie original ---
  cat("\n--- ADF (serie original) ---\n")
  adf_orig <- adf.test(ts_serie)
  print(adf_orig)

  # --- Box-Cox descriptivo (solo si la serie es positiva) ---
  cat("\n--- Box-Cox ---\n")
  bc_lambda <- if (all(ts_serie > 0)) BoxCox.lambda(ts_serie) else NA_real_
  if (is.na(bc_lambda)) {
    cat("Lambda óptimo: no aplicable (la serie contiene ceros)\n")
    ts_bc <- ts_serie
  } else {
    cat("Lambda óptimo (Guerrero):", round(bc_lambda, 4), "\n")
    ts_bc <- BoxCox(ts_serie, bc_lambda)
  }

  # --- Diferenciación ordinaria ---
  ts_diff <- diff(ts_bc, differences = 1)
  cat("\n--- ADF (serie diferenciada) ---\n")
  adf_diff <- adf.test(ts_diff)
  print(adf_diff)

  # --- Diferenciación estacional ---
  ts_diff_seas <- diff(ts_diff, lag = 12)
  cat("\n--- ADF (diferenciada estacionalmente) ---\n")
  adf_seas <- adf.test(ts_diff_seas)
  print(adf_seas)

  # --- Gráficos ---
  p_season <- ggseasonplot(ts_bc) +
    labs(title = paste0("Season plot — ", nombre))
  p_subserie <- ggsubseriesplot(ts_bc) +
    labs(title = paste0("Subseries plot — ", nombre))

  ggsave(here("graficas", "G_08",
              paste0("seasonplot_", tolower(nombre), ".png")),
         p_season, width = 9, height = 5, dpi = 150, bg = "white")
  ggsave(here("graficas", "G_08",
              paste0("subseries_", tolower(nombre), ".png")),
         p_subserie, width = 9, height = 5, dpi = 150, bg = "white")

  # --- ACF/PACF ---
  png(here("graficas", "G_08",
           paste0("acf_pacf_", tolower(nombre), ".png")),
      width = 1200, height = 500, res = 150)
  par(mfrow = c(1, 2))
  acf(ts_diff, main = paste0("ACF — ", nombre, " (dif.)"))
  pacf(ts_diff, main = paste0("PACF — ", nombre, " (dif.)"))
  dev.off()

  # --- SARIMA automático con lambda="auto" y biasadj=TRUE ---
  cat("\n--- Ajuste SARIMA (auto.arima) ---\n")
  fit_sarima <- ajustar_sarima(ts_serie, nombre)
  print(summary(fit_sarima))

  cat("\nAIC:", AIC(fit_sarima), "\n")
  cat("BIC:", BIC(fit_sarima), "\n")

  # --- Órdenes para el fitdf del Ljung-Box ---
  p <- fit_sarima$arma[1]
  q <- fit_sarima$arma[2]
  P <- fit_sarima$arma[3]
  Q <- fit_sarima$arma[4]
  fitdf_lb <- as.integer(p + q + P + Q)
  cat("\nÓrdenes: p =", p, "| q =", q, "| P =", P, "| Q =", Q,
      "| fitdf = p+q+P+Q =", fitdf_lb, "\n")

  # --- Validación de residuos ---
  residuos <- residuals(fit_sarima)
  if (anyNA(residuos)) residuos <- residuos[is.finite(residuos)]

  cat("\n--- Jarque-Bera (normalidad) ---\n")
  jb <- jarque.bera.test(residuos)
  print(jb)

  # McLeod-Li: Ljung-Box sobre residuos al cuadrado. Sustituye a
  # Breusch-Pagan, que no detecta heterocedasticidad condicional (ARCH).
  cat("\n--- McLeod-Li (ARCH) sobre residuos^2, lag = 12 ---\n")
  mcleod_li <- Box.test(residuos^2, lag = 12, type = "Ljung-Box")
  print(mcleod_li)

  cat("\n--- McLeod-Li (ARCH) sobre residuos^2, lag = 24 ---\n")
  mcleod_li_24 <- Box.test(residuos^2, lag = 24, type = "Ljung-Box")
  print(mcleod_li_24)

  # Ljung-Box con fitdf = p + q + P + Q. Se elimina Durbin-Watson.
  cat("\n--- Ljung-Box (autocorrelación), lag = 24, fitdf =",
      fitdf_lb, "---\n")
  lb <- Box.test(residuos, lag = 24, type = "Ljung-Box", fitdf = fitdf_lb)
  print(lb)

  # --- Validación fuera de muestra: últimos 12 meses reservados ---
  h <- 12L
  cat("\n--- Validación con ventana reservada (h =", h, "meses) ---\n")

  if (n > h + 12L) {
    train <- window(ts_serie, end   = time(ts_serie)[n - h])
    test  <- window(ts_serie, start = time(ts_serie)[n - h + 1])

    fit_train <- ajustar_sarima(train, paste0(nombre, " [entrenamiento]"))
    fc_train  <- forecast(fit_train, h = h)
    acc_test  <- accuracy(fc_train, test)

    print(acc_test)
    cat("Nota: la columna MAPE de accuracy() puede ser Inf si el test",
        "tiene meses con 0 casos.\n")
    cat("La métrica que se reporta abajo excluye esos meses (MAPE positivo).\n")

    validacion <- list(
      h         = h,
      n_train   = length(train),
      n_test    = length(test),
      train_end = end(train),
      test_start = start(test),
      accuracy  = acc_test,
      mae       = unname(acc_test["Test set", "MAE"]),
      rmse      = unname(acc_test["Test set", "RMSE"]),
      mape      = mape_positivo(test, fc_train$mean),
      mase      = unname(acc_test["Test set", "MASE"])
    )

    cat("MASE (ventana reservada):", round(validacion$mase, 4), "\n")
    cat("MAPE (meses con casos > 0):",
        round(validacion$mape, 2), "%\n")
  } else {
    validacion <- list(h = h, mae = NA_real_, rmse = NA_real_,
                       mape = NA_real_, mase = NA_real_)
    cat("Serie demasiado corta para reservar", h, "meses.\n")
  }

  # --- Pronóstico 12 meses ---
  forecast_sarima <- forecast(fit_sarima, h = h)

  # La transformación Box-Cox se aplica dentro de auto.arima; con
  # biasadj = TRUE el pronóstico ya vuelve a la escala de casos.
  # Se acotan los límites negativos (chikungunya los produce).
  forecast_sarima$lower[forecast_sarima$lower < 0] <- 0
  forecast_sarima$upper[forecast_sarima$upper < 0] <- 0

  png(here("graficas", "G_08",
           paste0("forecast_", tolower(nombre), ".png")),
      width = 1200, height = 600, res = 150)
  plot(forecast_sarima,
       main = paste0("Pronóstico SARIMA 12 meses — ", nombre),
       xlab = "Año", ylab = "Casos")
  dev.off()

  # Serie transformada en la escala del modelo (para trazabilidad).
  lambda_fit <- fit_sarima$lambda
  serie_bc <- if (!is.null(lambda_fit) && all(ts_serie > 0)) {
    BoxCox(ts_serie, lambda_fit)
  } else {
    ts_serie
  }

  # --- Guardar resultados ---
  list(
    serie_original = ts_serie,
    serie_bc       = serie_bc,
    lambda_bc      = lambda_fit,
    fit            = fit_sarima,
    forecast       = forecast_sarima,
    validacion     = validacion,
    tests = list(adf = adf_orig, adf_diff = adf_diff,
                 adf_seas = adf_seas,
                 jb = jb,
                 mcleod_li = mcleod_li,
                 mcleod_li_24 = mcleod_li_24,
                 lb = lb)
  )
}

# ============================================================
# 4. Ejecutar para Dengue y Chikungunya
# ============================================================

res_econ_dengue <- analisis_econometrico(ts_dengue, "Dengue",      "#08519C")
res_econ_chik   <- analisis_econometrico(ts_chik,   "Chikungunya", "#B2182B")

# ============================================================
# 5. Guardar resultados
# ============================================================

saveRDS(list(dengue = res_econ_dengue, chik = res_econ_chik),
        here("data", "D_08", "econometria_sarima.rds"))

cat("\n✅ Análisis econométrico completado\n")
