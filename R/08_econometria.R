# ============================================================
# 08_econometria.R
# Análisis econométrico de series temporales (Paper 2)
# ADF, Box-Cox, diferenciación, ACF/PACF, SARIMA, validación
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(tseries)
  library(forecast)
  library(lmtest)
  library(MASS)
  library(scales)
})

serie_dengue <- readRDS(here("data", "series", "serie_dengue.rds"))
serie_chik   <- readRDS(here("data", "series", "serie_chik.rds"))

dir.create(here("graficas", "econometria"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Serie mensual nacional
# ============================================================

preparar_mensual <- function(serie) {
  serie %>%
    filter(dpto == "NACIONAL") %>%
    mutate(mes = floor_date(fecha_sem, "month")) %>%
    group_by(mes) %>%
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") %>%
    arrange(mes) %>%
    filter(!is.na(casos))
}

mensual_dengue <- preparar_mensual(serie_dengue)
mensual_chik   <- preparar_mensual(serie_chik)

ts_dengue <- ts(mensual_dengue$casos,
                start = c(year(min(mensual_dengue$mes)),
                          month(min(mensual_dengue$mes))),
                frequency = 12)

ts_chik <- ts(mensual_chik$casos,
              start = c(year(min(mensual_chik$mes)),
                        month(min(mensual_chik$mes))),
              frequency = 12)

cat("\n===== SERIE MENSUAL =====\n")
cat("Dengue:      ", length(ts_dengue), "obs\n")
cat("Chikungunya: ", length(ts_chik),   "obs\n")

# ============================================================
# 2. Función: análisis econométrico completo
# ============================================================

analisis_econometrico <- function(ts_serie, nombre, color = "#08519C") {
  
  cat("\n=========================================================\n")
  cat("ANÁLISIS ECONOMÉTRICO:", nombre, "\n")
  cat("=========================================================\n")
  
  # --- ADF serie original ---
  cat("\n--- ADF (serie original) ---\n")
  adf_orig <- adf.test(ts_serie)
  print(adf_orig)
  
  # --- Box-Cox ---
  cat("\n--- Box-Cox ---\n")
  bc_lambda <- BoxCox.lambda(ts_serie)
  cat("Lambda óptimo:", round(bc_lambda, 4), "\n")
  ts_bc <- BoxCox(ts_serie, bc_lambda)
  
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
  
  ggsave(here("graficas", "econometria",
              paste0("seasonplot_", tolower(nombre), ".png")),
         p_season, width = 9, height = 5, dpi = 150, bg = "white")
  ggsave(here("graficas", "econometria",
              paste0("subseries_", tolower(nombre), ".png")),
         p_subserie, width = 9, height = 5, dpi = 150, bg = "white")
  
  # --- ACF/PACF ---
  png(here("graficas", "econometria",
           paste0("acf_pacf_", tolower(nombre), ".png")),
      width = 1200, height = 500, res = 150)
  par(mfrow = c(1, 2))
  acf(ts_diff, main = paste0("ACF — ", nombre, " (dif.)"))
  pacf(ts_diff, main = paste0("PACF — ", nombre, " (dif.)"))
  dev.off()
  
  # --- SARIMA automático ---
  cat("\n--- Ajuste SARIMA (auto.arima) ---\n")
  fit_sarima <- auto.arima(ts_bc, seasonal = TRUE,
                           stepwise = TRUE, approximation = FALSE)
  print(summary(fit_sarima))
  
  cat("\nAIC:", AIC(fit_sarima), "\n")
  cat("BIC:", BIC(fit_sarima), "\n")
  
  # --- Validación de residuos ---
  residuos <- residuals(fit_sarima)
  
  cat("\n--- Jarque-Bera (normalidad) ---\n")
  jb <- jarque.bera.test(residuos)
  print(jb)
  
  cat("\n--- Breusch-Pagan (homocedasticidad) ---\n")
  bp <- bptest(residuos ~ fitted(fit_sarima))
  print(bp)
  
  cat("\n--- Durbin-Watson (autocorrelación) ---\n")
  dw <- dwtest(residuos ~ fitted(fit_sarima))
  print(dw)
  
  cat("\n--- Ljung-Box (autocorrelación global) ---\n")
  lb <- Box.test(residuos, lag = 12, type = "Ljung-Box")
  print(lb)
  
  # --- Pronóstico 12 meses ---
  forecast_sarima <- forecast(fit_sarima, h = 12)
  
  # Volver a escala original
  forecast_orig <- forecast_sarima
  forecast_orig$mean  <- InvBoxCox(forecast_sarima$mean,  bc_lambda)
  forecast_orig$lower <- InvBoxCox(forecast_sarima$lower, bc_lambda)
  forecast_orig$upper <- InvBoxCox(forecast_sarima$upper, bc_lambda)
  
  png(here("graficas", "econometria",
           paste0("forecast_", tolower(nombre), ".png")),
      width = 1200, height = 600, res = 150)
  plot(forecast_orig,
       main = paste0("Pronóstico SARIMA 12 meses — ", nombre),
       xlab = "Año", ylab = "Casos")
  dev.off()
  
  # --- Guardar resultados ---
  list(
    serie_original  = ts_serie,
    serie_bc        = ts_bc,
    lambda_bc       = bc_lambda,
    fit             = fit_sarima,
    forecast        = forecast_orig,
    tests = list(adf = adf_orig, adf_diff = adf_diff,
                 adf_seas = adf_seas,
                 jb = jb, bp = bp, dw = dw, lb = lb)
  )
}

# ============================================================
# 3. Ejecutar para Dengue y Chikungunya
# ============================================================

res_econ_dengue <- analisis_econometrico(ts_dengue, "Dengue",      "#08519C")
res_econ_chik   <- analisis_econometrico(ts_chik,   "Chikungunya", "#B2182B")

# ============================================================
# 4. Guardar resultados
# ============================================================

saveRDS(list(dengue = res_econ_dengue, chik = res_econ_chik),
        here("data", "series", "econometria_sarima.rds"))

cat("\n✅ Análisis econométrico completado\n")
