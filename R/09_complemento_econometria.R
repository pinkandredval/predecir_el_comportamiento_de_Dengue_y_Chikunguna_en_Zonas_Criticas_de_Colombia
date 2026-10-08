# ============================================================
# 09_econometria_avanzada.R
# Análisis econométrico avanzado para Dengue y Chikungunya
# Complementa al 08_econometria.R sin modificarlo.
#   - VAR + Causalidad de Granger + IRF + FEVD
#   - GARCH(1,1) para Dengue
#   - Suavizamiento Exponencial (ETS) para Chikungunya
#   - Test de Encompassing (SARIMA vs VAR)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(forecast)
  library(vars)
  library(rugarch)
  library(lmtest)
  library(tseries)
  library(scales)
})

# --- Cargar resultados del 08 (NO recalcular SARIMA) ---
res_econ <- readRDS(here("data", "series", "econometria_sarima.rds"))

# ============================================================
# 0. Preparar series mensuales para análisis multivariado
# ============================================================

serie_dengue <- readRDS(here("data", "series", "serie_dengue.rds"))
serie_chik   <- readRDS(here("data", "series", "serie_chik.rds"))

preparar_mensual <- function(serie, nombre) {
  serie %>%
    filter(dpto == "NACIONAL") %>%
    mutate(mes = floor_date(fecha_sem, "month")) %>%
    group_by(mes) %>%
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") %>%
    arrange(mes) %>%
    filter(!is.na(casos)) %>%
    mutate(evento = nombre)
}

mensual_dengue <- preparar_mensual(serie_dengue, "Dengue")
mensual_chik   <- preparar_mensual(serie_chik,   "Chikungunya")

# --- Alinear por fecha (dplyr::select explícito) ---
df_multi <- full_join(
  mensual_dengue %>% dplyr::select(mes, dengue = casos),
  mensual_chik   %>% dplyr::select(mes, chik   = casos),
  by = "mes"
) %>%
  arrange(mes) %>%
  filter(!is.na(dengue), !is.na(chik))

ts_multi <- ts(
  as.matrix(df_multi[, c("dengue", "chik")]),
  start = c(year(min(df_multi$mes)), month(min(df_multi$mes))),
  frequency = 12
)

cat("\n===== SERIES MENSUALES ALINEADAS =====\n")
cat("Observaciones:", nrow(df_multi), "\n")
cat("Rango:", format(min(df_multi$mes)), "→", format(max(df_multi$mes)), "\n")
print(head(df_multi, 3))

# ============================================================
# 1. VAR + Causalidad de Granger
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("1. MODELO VAR + CAUSALIDAD DE GRANGER\n")
cat("=========================================================\n")

lag_sel <- VARselect(ts_multi, lag.max = 8, type = "const")
cat("\n--- Selección de rezagos óptimos ---\n")
print(lag_sel$selection)

lag_optimo <- as.numeric(lag_sel$selection["AIC(n)"])
cat("\nRezago óptimo según AIC:", lag_optimo, "\n")

fit_var <- VAR(ts_multi, p = lag_optimo, type = "const")
cat("\n--- Resumen del VAR ---\n")
print(summary(fit_var))

cat("\n--- Estabilidad del VAR (raíces) ---\n")
roots_var <- roots(fit_var)
cat("Módulo máximo:", round(max(roots_var), 3), "\n")
cat("¿Estable?", all(roots_var < 1), "\n")

cat("\n--- Test de Portmanteau ---\n")
port_var <- serial.test(fit_var, lags.pt = 12, type = "PT.asymptotic")
print(port_var)

cat("\n--- Causalidad de Granger ---\n")
cat("H0: 'Dengue' NO causa 'Chikungunya'\n")
granger_dengue_a_chik <- causality(fit_var, cause = "dengue")
print(granger_dengue_a_chik$Granger)

cat("\nH0: 'Chikungunya' NO causa 'Dengue'\n")
granger_chik_a_dengue <- causality(fit_var, cause = "chik")
print(granger_chik_a_dengue$Granger)

# --- IRF ---
irf_var <- irf(fit_var, n.ahead = 12, boot = TRUE, ci = 0.95)
plot(irf_var)

# --- FEVD ---
fevd_var <- fevd(fit_var, n.ahead = 12)
plot(fevd_var)

# --- Pronóstico ---
forecast_var <- predict(fit_var, n.ahead = 12)
print(forecast_var)

# ============================================================
# 2. GARCH(1,1) para Dengue
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("2. GARCH(1,1) — Dengue\n")
cat("=========================================================\n")

ts_bc_dengue <- res_econ$dengue$serie_bc

spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(4, 0), include.mean = FALSE),
  distribution.model = "std"
)

fit_garch <- ugarchfit(spec = spec_garch, data = ts_bc_dengue)

cat("\n--- Resumen GARCH(1,1) ---\n")
print(fit_garch)

resid_std <- residuals(fit_garch, standardize = TRUE)

cat("\n--- Jarque-Bera (residuos estandarizados) ---\n")
print(jarque.bera.test(resid_std))

cat("\n--- Ljung-Box residuos ---\n")
print(Box.test(resid_std, lag = 12, type = "Ljung-Box"))

cat("\n--- Ljung-Box residuos^2 ---\n")
print(Box.test(resid_std^2, lag = 12, type = "Ljung-Box"))

forecast_garch <- ugarchforecast(fit_garch, n.ahead = 12)
plot(forecast_garch, which = 1)

aic_sarima <- AIC(res_econ$dengue$fit)
aic_garch  <- infocriteria(fit_garch)[1]

cat("\n--- Comparación AIC ---\n")
cat("SARIMA (08)     :", round(aic_sarima, 2), "\n")
cat("SARIMA+GARCH(09):", round(aic_garch, 2), "\n")
cat("Mejora:", round(aic_sarima - aic_garch, 2), "\n")

# ============================================================
# 3. ETS para Chikungunya
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("3. ETS — Chikungunya\n")
cat("=========================================================\n")

ts_chik <- res_econ$chik$serie_original

fit_ets <- ets(ts_chik)
cat("\n--- Resumen ETS ---\n")
print(summary(fit_ets))

fitted_ets <- fitted(fit_ets)
mape_ets <- mean(abs((ts_chik - fitted_ets) / ts_chik), na.rm = TRUE) * 100
cat("\nMAPE ETS:", round(mape_ets, 2), "%\n")

resid_ets <- residuals(fit_ets)

cat("\n--- Ljung-Box residuos ETS ---\n")
print(Box.test(resid_ets, lag = 12, type = "Ljung-Box"))

forecast_ets <- forecast(fit_ets, h = 12)
plot(forecast_ets, main = "Pronóstico ETS — Chikungunya",
     xlab = "Año", ylab = "Casos")

# ============================================================
# 4. Test de Encompassing
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("4. TEST DE ENCOMPASSING (SARIMA vs VAR)\n")
cat("=========================================================\n")

pred_sarima <- fitted(res_econ$dengue$fit)
pred_var    <- fitted(fit_var)[, "dengue"]

n_min  <- min(length(pred_sarima), length(pred_var))
y_real <- as.numeric(res_econ$dengue$serie_bc)[1:n_min]
pred_s <- pred_sarima[1:n_min]
pred_v <- pred_var[1:n_min]

df_enc <- tibble(y = y_real, pred_s = pred_s, pred_v = pred_v) %>% drop_na()

if (nrow(df_enc) > 10) {
  lm_enc <- lm(y ~ pred_s + pred_v, data = df_enc)
  cat("\n--- Regresión de Encompassing ---\n")
  print(summary(lm_enc))
  
  cat("\n--- Test de Wald (coeficiente de pred_v = 0) ---\n")
  wald_test <- lmtest::waldtest(lm_enc, . ~ . - pred_v)
  print(wald_test)
} else {
  cat("Datos insuficientes para encompassing\n")
}

# ============================================================
# 5. Resumen comparativo consolidado
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("5. RESUMEN COMPARATIVO FINAL\n")
cat("=========================================================\n")

# --- MAPE SARIMA (in-sample) ---
mape_sarima_dengue <- mean(
  abs((res_econ$dengue$serie_bc - fitted(res_econ$dengue$fit)) /
        res_econ$dengue$serie_bc),
  na.rm = TRUE
) * 100

mape_sarima_chik <- mean(
  abs((res_econ$chik$serie_bc - fitted(res_econ$chik$fit)) /
        res_econ$chik$serie_bc),
  na.rm = TRUE
) * 100

# --- MAPE SARIMA+GARCH (media igual al SARIMA) ---
mape_garch_dengue <- mape_sarima_dengue

# --- MAPE VAR por variable ---
fitted_var  <- fitted(fit_var)
real_dengue <- as.numeric(ts_multi[, "dengue"])
real_chik   <- as.numeric(ts_multi[, "chik"])
pred_var_d  <- as.numeric(fitted_var[, "dengue"])
pred_var_c  <- as.numeric(fitted_var[, "chik"])

mape_var_dengue <- mean(abs((real_dengue - pred_var_d) / real_dengue),
                        na.rm = TRUE) * 100
mape_var_chik   <- mean(abs((real_chik - pred_var_c) / real_chik),
                        na.rm = TRUE) * 100

cat("\n--- MAPEs calculados ---\n")
cat("SARIMA Dengue:       ", round(mape_sarima_dengue, 2), "%\n")
cat("SARIMA+GARCH Dengue: ", round(mape_garch_dengue, 2), "%\n")
cat("VAR Dengue:          ", round(mape_var_dengue, 2),   "%\n")
cat("SARIMA Chikungunya:  ", round(mape_sarima_chik, 2),  "%\n")
cat("VAR Chikungunya:     ", round(mape_var_chik, 2),    "%\n")
cat("ETS Chikungunya:     ", round(mape_ets, 2),          "%\n")

# --- Utilidad ---
util_sarima_dengue <- is.finite(mape_sarima_dengue) && mape_sarima_dengue < 50
util_sarima_chik   <- is.finite(mape_sarima_chik)   && mape_sarima_chik   < 50

resumen_avanzado <- tibble(
  evento = c("Dengue", "Dengue", "Dengue",
             "Chikungunya", "Chikungunya", "Chikungunya"),
  modelo = c("SARIMA", "SARIMA+GARCH", "VAR",
             "SARIMA", "ETS", "VAR"),
  mape   = c(round(mape_sarima_dengue, 2),
             round(mape_garch_dengue, 2),
             round(mape_var_dengue, 2),
             round(mape_sarima_chik, 2),
             round(mape_ets, 2),
             round(mape_var_chik, 2)),
  aic    = c(as.numeric(AIC(res_econ$dengue$fit))[1],
             as.numeric(aic_garch)[1],
             as.numeric(AIC(fit_var))[1],
             as.numeric(AIC(res_econ$chik$fit))[1],
             as.numeric(AIC(fit_ets))[1],
             as.numeric(AIC(fit_var))[1]),
  util   = c(util_sarima_dengue, TRUE, FALSE,
             util_sarima_chik,   TRUE, FALSE)
)

print(as.data.frame(resumen_avanzado), row.names = FALSE)

# ============================================================
# 6. Guardar resultados
# ============================================================

saveRDS(
  list(
    var            = fit_var,
    var_forecast   = forecast_var,
    var_lag        = lag_optimo,
    granger        = list(dengue_a_chik = granger_dengue_a_chik,
                          chik_a_dengue = granger_chik_a_dengue),
    irf            = irf_var,
    fevd           = fevd_var,
    garch_dengue   = fit_garch,
    garch_forecast = forecast_garch,
    ets_chik       = fit_ets,
    ets_forecast   = forecast_ets,
    mapes = list(sarima_dengue = mape_sarima_dengue,
                 garch_dengue  = mape_garch_dengue,
                 var_dengue    = mape_var_dengue,
                 sarima_chik   = mape_sarima_chik,
                 ets_chik      = mape_ets,
                 var_chik      = mape_var_chik),
    resumen        = resumen_avanzado
  ),
  here("data", "series", "econometria_avanzada.rds")
)

write.csv(resumen_avanzado,
          here("data", "series", "resumen_econometria_avanzada.csv"),
          row.names = FALSE)

cat("\n✅ Análisis econométrico avanzado completado\n")
