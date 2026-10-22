# ============================================================
# 09_econometria_avanzada.R
# Análisis econométrico avanzado para Dengue y Chikungunya
# Complementa al 08_econometria.R sin modificarlo.
#   - VAR + Causalidad de Granger + IRF + FEVD
#   - GARCH(1,1) mensual para la comparación de AIC
#   - GJR-GARCH semanal para el diagnóstico de volatilidad
#   - Suavizamiento Exponencial (ETS) para Chikungunya
#   - Test de Encompassing (SARIMA vs VAR), correctamente alineado
#
# Correcciones Tabla 2 (#2, #3, #7, #8, #9):
#   #2 AIC del GARCH = infocriteria(fit)[1] * length(r), para llevarlo
#      a la escala total del SARIMA (no AIC por observación).
#   #3 El VAR pierde los primeros p meses; y_real y pred_sarima se
#      recortan con el mismo índice (p+1):n. No se recicla nada.
#   #7/#9 Se retira el GARCH mensual si no estima ARCH (alpha1 no
#      significativo). Se añade un GJR-GARCH semanal con Nyblom y
#      sign bias.
#   #8 El AIC del VAR no se compara con SARIMA/ETS: es una
#      verosimilitud conjunta de dos ecuaciones.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(forecast)
  library(vars)
  library(rugarch)
  library(lmtest)
  library(tseries)
  library(scales)
})

# --- Cargar resultados del 08 (NO recalcular SARIMA) ---
stopifnot(file.exists(here("data", "series", "econometria_sarima.rds")))
res_econ <- readRDS(here("data", "series", "econometria_sarima.rds"))

serie_dengue <- readRDS(here("data", "series", "serie_dengue.rds"))
serie_chik   <- readRDS(here("data", "series", "serie_chik.rds"))

dir.create(here("graficas", "econometria_avanzada"),
           recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0. Preparar series mensuales para análisis multivariado
# ============================================================

preparar_mensual <- function(serie, nombre) {
  serie %>%
    filter(dpto == "NACIONAL") %>%
    mutate(mes = floor_date(fecha_sem, "month")) %>%
    group_by(mes) %>%
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") %>%
    complete(mes = seq(min(mes), max(mes), by = "month"),
             fill = list(casos = 0)) %>%
    arrange(mes) %>%
    mutate(evento = nombre)
}

mensual_dengue <- preparar_mensual(serie_dengue, "Dengue")
mensual_chik   <- preparar_mensual(serie_chik,   "Chikungunya")

stopifnot(nrow(mensual_dengue) == 72L, nrow(mensual_chik) == 72L)

# --- Alinear por fecha (dplyr::select explícito) ---
df_multi <- full_join(
  mensual_dengue %>% dplyr::select(mes, dengue = casos),
  mensual_chik   %>% dplyr::select(mes, chik   = casos),
  by = "mes"
) %>%
  arrange(mes) %>%
  filter(!is.na(dengue), !is.na(chik))

stopifnot(nrow(df_multi) == 72L)

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
# Utilidades de error
# ============================================================

# MAPE solo sobre valores observados > 0. Evita el Inf de chikungunya.
mape_positivo <- function(actual, pred) {
  actual <- as.numeric(actual)
  pred   <- as.numeric(pred)
  ok     <- is.finite(actual) & is.finite(pred) & actual > 0
  if (!any(ok)) return(NA_real_)
  mean(abs((actual[ok] - pred[ok]) / actual[ok])) * 100
}

# MASE in-sample para el VAR (no tiene ventana reservada en este script).
mase_insample <- function(actual, pred) {
  actual <- as.numeric(actual)
  pred   <- as.numeric(pred)
  ok     <- is.finite(actual) & is.finite(pred)
  if (sum(ok) < 2L) return(NA_real_)
  denom <- mean(abs(diff(actual[ok])))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  mean(abs(actual[ok] - pred[ok])) / denom
}

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
set.seed(123)
irf_var <- irf(fit_var, n.ahead = 12, boot = TRUE, ci = 0.95)

png(here("graficas", "econometria_avanzada", "var_irf.png"),
    width = 1200, height = 800, res = 150)
plot(irf_var)
dev.off()

# --- FEVD ---
fevd_var <- fevd(fit_var, n.ahead = 12)

png(here("graficas", "econometria_avanzada", "var_fevd.png"),
    width = 1200, height = 800, res = 150)
plot(fevd_var)
dev.off()

# --- Pronóstico ---
forecast_var <- predict(fit_var, n.ahead = 12)
print(forecast_var)

# --- Alineación del VAR (corrección #3) ---
# fitted(fit_var) tiene n - p filas: las primeras p observaciones reales
# se descartan. Se usa el mismo índice para y_real, SARIMA y VAR.
p_var    <- fit_var$p
n_var    <- nrow(df_multi)
idx_var  <- (p_var + 1):n_var

fitted_var <- fitted(fit_var)
pred_var_d <- as.numeric(fitted_var[, "dengue"])
pred_var_c <- as.numeric(fitted_var[, "chik"])
y_var_d    <- as.numeric(ts_multi[idx_var, "dengue"])
y_var_c    <- as.numeric(ts_multi[idx_var, "chik"])

stopifnot(length(pred_var_d) == length(idx_var),
          length(pred_var_c) == length(idx_var),
          length(y_var_d)    == length(idx_var),
          length(y_var_c)    == length(idx_var))

cat("\n--- Alineación del VAR ---\n")
cat("p =", p_var, "| observaciones reales descartadas:", p_var,
    "| pares alineados:", length(idx_var), "\n")

# ============================================================
# 2. GARCH: mensual (comparación AIC) y semanal (diagnóstico)
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("2. GARCH — Dengue\n")
cat("=========================================================\n")

# ------------------------------------------------------------
# 2.1 GARCH(1,1) mensual sobre la serie Box-Cox del SARIMA
# ------------------------------------------------------------
# Se usa la misma muestra (72), la misma transformación Box-Cox y una
# distribución normal, igual que el SARIMA, para que el AIC total sea
# comparable. Si alpha1 no es significativo, el GARCH mensual no estima
# ARCH y se retira como modelo (corrección #7); su AIC total se reporta
# solo para dejar constancia de la escala.
ts_bc_dengue <- res_econ$dengue$serie_bc
r_mensual    <- as.numeric(ts_bc_dengue)

set.seed(123)
spec_garch_mensual <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_garch_mensual <- ugarchfit(
  spec   = spec_garch_mensual,
  data   = r_mensual,
  solver = "hybrid"
)

stopifnot(convergence(fit_garch_mensual) == 0)

cat("\n--- GARCH(1,1) mensual (misma escala Box-Cox que el SARIMA) ---\n")
print(fit_garch_mensual)

mat_m <- fit_garch_mensual@fit$matcoef
colnames(mat_m) <- trimws(colnames(mat_m))
alpha_p_m  <- mat_m["alpha1", "Pr(>|t|)"]
gamma_p_m  <- NA_real_
hay_arch_m <- alpha_p_m < 0.05
persist_m  <- mat_m["alpha1", "Estimate"] + mat_m["beta1", "Estimate"]

# Corrección #2: AIC total = AIC por observación * T.
aic_garch_mensual <- as.numeric(infocriteria(fit_garch_mensual)[1]) *
  length(r_mensual)

cat("\nAIC GARCH mensual por observación:",
    round(infocriteria(fit_garch_mensual)[1], 4), "\n")
cat("AIC GARCH mensual total (x", length(r_mensual), "):",
    round(aic_garch_mensual, 2), "\n")
cat("¿Estima ARCH? alpha p =", round(alpha_p_m, 4),
    "→", if (hay_arch_m) "SÍ" else "NO", "\n")
if (!hay_arch_m) {
  cat("El GARCH mensual NO estima ARCH: se retira como modelo.\n")
  cat("Su AIC total se conserva solo para mostrar la escala correcta.\n")
}

ny_m <- nyblom(fit_garch_mensual)
sb_m <- signbias(fit_garch_mensual)

cat("\n--- Nyblom (estabilidad de parámetros, GARCH mensual) ---\n")
cat("Estadístico conjunto:", round(ny_m$JointStat, 3),
    "| crítico 5%:", round(ny_m$JointCritical["5%"], 3), "\n")
print(round(ny_m$IndividualStat, 3))

cat("\n--- Sign bias (GARCH mensual) ---\n")
print(sb_m)

resid_std_m <- residuals(fit_garch_mensual, standardize = TRUE)

cat("\n--- Ljung-Box residuos estandarizados (GARCH mensual) ---\n")
print(Box.test(resid_std_m, lag = 12, type = "Ljung-Box", fitdf = 2))
cat("\n--- Ljung-Box residuos^2 (GARCH mensual) ---\n")
print(Box.test(resid_std_m^2, lag = 12, type = "Ljung-Box"))

# ------------------------------------------------------------
# 2.2 GJR-GARCH semanal sobre log1p(casos)
# ------------------------------------------------------------
# La serie mensual tiene 72 observaciones: demasiado corta para
# identificar un GARCH. El diagnóstico de volatilidad se hace sobre
# las 314 semanas, con estacionalidad anual (armónicos) y AR(1).
sem_dengue <- serie_dengue %>%
  filter(dpto == "NACIONAL") %>%
  arrange(anio, semana)

stopifnot(nrow(sem_dengue) == 314L)

y_sem <- ts(log1p(sem_dengue$casos), frequency = 52)
r_semanal <- as.numeric(y_sem)
t_sem <- seq_along(y_sem)

X_sem <- cbind(
  sin1 = sin(2 * pi * t_sem / 52),
  cos1 = cos(2 * pi * t_sem / 52),
  sin2 = sin(4 * pi * t_sem / 52),
  cos2 = cos(4 * pi * t_sem / 52)
)

set.seed(123)
spec_garch_semanal <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE,
                            external.regressors = X_sem),
  distribution.model = "std"
)

fit_garch_semanal <- tryCatch(
  ugarchfit(spec = spec_garch_semanal, data = r_semanal, solver = "hybrid"),
  error = function(e) NULL
)

if (is.null(fit_garch_semanal) || convergence(fit_garch_semanal) != 0) {
  cat("\nGJR-GARCH semanal no converge; se intenta EGARCH semanal.\n")
  spec_garch_semanal <- ugarchspec(
    variance.model     = list(model = "eGARCH", garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE,
                              external.regressors = X_sem),
    distribution.model = "std"
  )
  fit_garch_semanal <- ugarchfit(
    spec   = spec_garch_semanal,
    data   = r_semanal,
    solver = "hybrid"
  )
}

stopifnot(!is.null(fit_garch_semanal),
          convergence(fit_garch_semanal) == 0)

cat("\n--- GARCH semanal (log1p, 314 semanas) ---\n")
print(fit_garch_semanal)

mat_w <- fit_garch_semanal@fit$matcoef
colnames(mat_w) <- trimws(colnames(mat_w))
alpha_p_w  <- mat_w["alpha1", "Pr(>|t|)"]
gamma_p_w  <- if ("gamma1" %in% rownames(mat_w)) {
  mat_w["gamma1", "Pr(>|t|)"]
} else {
  NA_real_
}
hay_arch_w <- (alpha_p_w < 0.05) ||
  (!is.na(gamma_p_w) && gamma_p_w < 0.05)
modelo_w   <- fit_garch_semanal@model$modeldesc$vmodel

persist_w <- if (identical(modelo_w, "gjrGARCH")) {
  mat_w["alpha1", "Estimate"] + mat_w["beta1", "Estimate"] +
    mat_w["gamma1", "Estimate"] / 2
} else {
  mat_w["alpha1", "Estimate"] + mat_w["beta1", "Estimate"]
}

# Corrección #2, aplicada también al GARCH semanal: AIC total.
aic_garch_semanal <- as.numeric(infocriteria(fit_garch_semanal)[1]) *
  length(r_semanal)

cat("\nModelo semanal:", modelo_w, "| T =", length(r_semanal), "\n")
cat("AIC semanal por observación:", round(infocriteria(fit_garch_semanal)[1], 4), "\n")
cat("AIC semanal total (x", length(r_semanal), "):",
    round(aic_garch_semanal, 2), "\n")
cat("Persistencia:", round(persist_w, 4), "\n")
cat("¿Estima ARCH? alpha p =", round(alpha_p_w, 4),
    "| gamma p =", round(gamma_p_w, 4),
    "→", if (hay_arch_w) "SÍ" else "NO", "\n")

ny_w <- nyblom(fit_garch_semanal)
sb_w <- signbias(fit_garch_semanal)

cat("\n--- Nyblom (GARCH semanal) ---\n")
cat("Estadístico conjunto:", round(ny_w$JointStat, 3),
    "| crítico 5%:", round(ny_w$JointCritical["5%"], 3), "\n")
print(round(ny_w$IndividualStat, 3))

cat("\n--- Sign bias (GARCH semanal) ---\n")
print(sb_w)

resid_std_w <- residuals(fit_garch_semanal, standardize = TRUE)

cat("\n--- Ljung-Box residuos estandarizados (GARCH semanal) ---\n")
print(Box.test(resid_std_w, lag = 24, type = "Ljung-Box", fitdf = 4))
cat("\n--- Ljung-Box residuos^2 (GARCH semanal) ---\n")
print(Box.test(resid_std_w^2, lag = 24, type = "Ljung-Box"))

# Pronóstico semanal (referencia; no se compara con el SARIMA mensual)
t_fut <- (length(y_sem) + 1):(length(y_sem) + 12)
X_fut <- cbind(
  sin1 = sin(2 * pi * t_fut / 52),
  cos1 = cos(2 * pi * t_fut / 52),
  sin2 = sin(4 * pi * t_fut / 52),
  cos2 = cos(4 * pi * t_fut / 52)
)
forecast_garch_semanal <- tryCatch(
  ugarchforecast(fit_garch_semanal, n.ahead = 12,
                 external.forecasts = list(mregfor = X_fut)),
  error = function(e) NULL
)

# ------------------------------------------------------------
# 2.3 Comparación de AIC en la misma escala
# ------------------------------------------------------------
aic_sarima <- as.numeric(AIC(res_econ$dengue$fit))

cat("\n--- Comparación de AIC en escala total ---\n")
cat("SARIMA (mensual)          :", round(aic_sarima, 2), "\n")
cat("GARCH mensual (total)     :", round(aic_garch_mensual, 2), "\n")
if (is.finite(aic_garch_mensual)) {
  cat("Diferencia GARCH - SARIMA :",
      round(aic_garch_mensual - aic_sarima, 2),
      "(positiva = GARCH peor)\n")
}
cat("GJR-GARCH (semanal, total):", round(aic_garch_semanal, 2), "\n")
cat("El AIC semanal NO es comparable con el SARIMA mensual:",
    "distinta frecuencia y transformación.\n")

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

h <- 12L
n_chik <- length(ts_chik)

train_ets <- window(ts_chik, end   = time(ts_chik)[n_chik - h])
test_ets  <- window(ts_chik, start = time(ts_chik)[n_chik - h + 1])

fit_ets_train <- ets(train_ets)
fc_ets_train  <- forecast(fit_ets_train, h = h)
acc_ets       <- accuracy(fc_ets_train, test_ets)

mape_ets <- mape_positivo(test_ets, fc_ets_train$mean)
mase_ets <- unname(acc_ets["Test set", "MASE"])

cat("\n--- ETS: validación con ventana reservada (h =", h, ") ---\n")
print(acc_ets)
cat("Nota: la columna MAPE de accuracy() puede ser Inf si el test",
    "tiene meses con 0 casos.\n")
cat("La métrica que se reporta abajo excluye esos meses (MAPE positivo).\n")
cat("MASE ETS (ventana reservada):", round(mase_ets, 4), "\n")
cat("MAPE ETS (meses con casos > 0):", round(mape_ets, 2), "%\n")

resid_ets <- residuals(fit_ets)

cat("\n--- Ljung-Box residuos ETS ---\n")
print(Box.test(resid_ets, lag = 24, type = "Ljung-Box"))

forecast_ets <- forecast(fit_ets, h = h)

png(here("graficas", "econometria_avanzada", "ets_chikungunya.png"),
    width = 1200, height = 600, res = 150)
plot(forecast_ets, main = "Pronóstico ETS — Chikungunya",
     xlab = "Año", ylab = "Casos")
dev.off()

# ============================================================
# 4. Test de Encompassing (alineado)
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("4. TEST DE ENCOMPASSING (SARIMA vs VAR)\n")
cat("=========================================================\n")

pred_sarima <- as.numeric(fitted(res_econ$dengue$fit))
y_real      <- as.numeric(res_econ$dengue$serie_original)[idx_var]
pred_s      <- pred_sarima[idx_var]
pred_v      <- pred_var_d

stopifnot(length(y_real) == length(pred_s),
          length(pred_s) == length(pred_v))

df_enc <- tibble(y = y_real, pred_s = pred_s, pred_v = pred_v) %>%
  drop_na()

cat("\nPares alineados:", nrow(df_enc), "\n")
cat("Primer mes usado:", format(min(df_multi$mes[idx_var])), "\n")

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

# --- Métricas del SARIMA: ventana reservada del 08 ---
mape_sarima_dengue <- res_econ$dengue$validacion$mape
mase_sarima_dengue <- res_econ$dengue$validacion$mase
mape_sarima_chik   <- res_econ$chik$validacion$mape
mase_sarima_chik   <- res_econ$chik$validacion$mase

# --- MAPE y MASE del VAR, con la alineación (p+1):n ---
mape_var_dengue <- mape_positivo(y_var_d, pred_var_d)
mape_var_chik   <- mape_positivo(y_var_c, pred_var_c)
mase_var_dengue <- mase_insample(y_var_d, pred_var_d)
mase_var_chik   <- mase_insample(y_var_c, pred_var_c)

# --- Utilidad: MASE < 1 frente al benchmark naive (regla única) ---
util_sarima_dengue <- is.finite(mase_sarima_dengue) && mase_sarima_dengue < 1
util_sarima_chik   <- is.finite(mase_sarima_chik)   && mase_sarima_chik   < 1
util_ets           <- is.finite(mase_ets)           && mase_ets < 1

cat("\n--- MAPEs (finitos; sin Inf) ---\n")
cat("SARIMA Dengue:       ", round(mape_sarima_dengue, 2), "%\n")
cat("VAR Dengue:          ", round(mape_var_dengue, 2),   "%\n")
cat("SARIMA Chikungunya:  ", round(mape_sarima_chik, 2),  "%\n")
cat("ETS Chikungunya:     ", round(mape_ets, 2),          "%\n")
cat("VAR Chikungunya:     ", round(mape_var_chik, 2),     "%\n")

cat("\n--- MASE (ventana reservada; VAR in-sample) ---\n")
cat("SARIMA Dengue:       ", round(mase_sarima_dengue, 4), "\n")
cat("VAR Dengue:          ", round(mase_var_dengue, 4),   "\n")
cat("SARIMA Chikungunya:  ", round(mase_sarima_chik, 4),  "\n")
cat("ETS Chikungunya:     ", round(mase_ets, 4),          "\n")
cat("VAR Chikungunya:     ", round(mase_var_chik, 4),     "\n")

# --- Tabla de modelos univariados y VAR ---
resumen_modelos <- tibble(
  evento   = c("Dengue", "Dengue", "Chikungunya", "Chikungunya", "Chikungunya"),
  modelo   = c("SARIMA", "VAR", "SARIMA", "ETS", "VAR"),
  escala   = "mensual",
  mape     = c(mape_sarima_dengue, mape_var_dengue,
               mape_sarima_chik,  mape_ets, mape_var_chik),
  mase     = c(mase_sarima_dengue, mase_var_dengue,
               mase_sarima_chik,  mase_ets, mase_var_chik),
  aic      = c(as.numeric(AIC(res_econ$dengue$fit)), NA_real_,
               as.numeric(AIC(res_econ$chik$fit)),  as.numeric(AIC(fit_ets)),
               NA_real_),
  comparable_sarima = c(TRUE, FALSE, TRUE, TRUE, FALSE),
  util     = c(util_sarima_dengue, NA, util_sarima_chik, util_ets, NA),
  nota     = c(
    "MAPE/MASE en ventana reservada; MAPE excluye meses con 0 casos",
    "MAPE/MASE in-sample alineados; AIC conjunta del VAR no comparable con univariados",
    "MAPE/MASE en ventana reservada; MAPE excluye meses con 0 casos",
    "MAPE/MASE en ventana reservada; MAPE excluye meses con 0 casos",
    "MAPE/MASE in-sample alineados; AIC conjunta del VAR no comparable con univariados"
  )
)

# --- Tabla de GARCH: AIC total y diagnósticos ---
nota_garch_m <- paste0(
  "AIC total = infocriteria*T (T=", length(r_mensual), "); ",
  "misma muestra y transformación que el SARIMA; ",
  if (hay_arch_m) "estima ARCH" else "NO estima ARCH (retirado)",
  "; Nyblom=", round(ny_m$JointStat, 3),
  " (crítico 5%=", round(ny_m$JointCritical["5%"], 3), ")",
  "; sign bias conjunto p=", round(sb_m["Joint Effect", "prob"], 3)
)

nota_garch_w <- paste0(
  "AIC total = infocriteria*T (T=", length(r_semanal), "); ",
  "NO comparable con el SARIMA mensual (frecuencia/transformación); ",
  if (hay_arch_w) "estima ARCH" else "alpha/gamma no significativos",
  "; Nyblom=", round(ny_w$JointStat, 3),
  " (crítico 5%=", round(ny_w$JointCritical["5%"], 3), ")",
  "; sign bias conjunto p=", round(sb_w["Joint Effect", "prob"], 3)
)

resumen_garch <- tibble(
  evento   = c("Dengue", "Dengue"),
  modelo   = c("GARCH mensual (sGARCH)", "GJR-GARCH semanal"),
  escala   = c("mensual (Box-Cox)", "semanal (log1p)"),
  mape     = NA_real_,
  mase     = NA_real_,
  aic      = c(aic_garch_mensual, aic_garch_semanal),
  comparable_sarima = c(TRUE, FALSE),
  util     = NA,
  nota     = c(nota_garch_m, nota_garch_w)
)

resumen_avanzado <- bind_rows(resumen_modelos, resumen_garch)

cat("\n--- Tabla resumen consolidada (métrica principal) ---\n")
print(
  as.data.frame(
    resumen_avanzado %>%
      dplyr::select(evento, modelo, escala, mape, mase, aic, util)
  ),
  row.names = FALSE
)

cat("\n--- Notas de comparabilidad ---\n")
print(
  as.data.frame(
    resumen_avanzado %>%
      dplyr::filter(nzchar(nota)) %>%
      dplyr::select(modelo, nota)
  ),
  row.names = FALSE
)

cat("\n--- Lectura de la comparación de AIC ---\n")
cat("AIC SARIMA (total):          ", round(aic_sarima, 2), "\n")
cat("AIC GARCH mensual (total):   ", round(aic_garch_mensual, 2), "\n")
cat("AIC GJR semanal (total):     ", round(aic_garch_semanal, 2),
    "(escala semanal; no comparable)\n")
cat("El AIC del VAR no se compara: es la verosimilitud conjunta de 2 ecuaciones.\n")

# ============================================================
# 6. Guardar resultados
# ============================================================

saveRDS(
  list(
    var               = fit_var,
    var_forecast      = forecast_var,
    var_lag           = lag_optimo,
    var_align         = list(p = p_var, idx = idx_var,
                             n = length(idx_var)),
    granger           = list(dengue_a_chik = granger_dengue_a_chik,
                             chik_a_dengue = granger_chik_a_dengue),
    irf               = irf_var,
    fevd              = fevd_var,
    garch_dengue      = fit_garch_mensual,
    garch_semanal_dengue = fit_garch_semanal,
    garch_forecast_semanal = forecast_garch_semanal,
    ets_chik          = fit_ets,
    ets_forecast      = forecast_ets,
    aic               = list(sarima = aic_sarima,
                             garch_mensual = aic_garch_mensual,
                             garch_semanal = aic_garch_semanal),
    garch_diagnostico = list(
      mensual = list(alpha_p = alpha_p_m, gamma_p = gamma_p_m,
                     hay_arch = hay_arch_m, persistencia = persist_m,
                     nyblom = ny_m, sign_bias = sb_m),
      semanal = list(alpha_p = alpha_p_w, gamma_p = gamma_p_w,
                     hay_arch = hay_arch_w, persistencia = persist_w,
                     nyblom = ny_w, sign_bias = sb_w)
    ),
    mapes = list(sarima_dengue = mape_sarima_dengue,
                 var_dengue    = mape_var_dengue,
                 sarima_chik   = mape_sarima_chik,
                 ets_chik      = mape_ets,
                 var_chik      = mape_var_chik),
    mases = list(sarima_dengue = mase_sarima_dengue,
                 var_dengue    = mase_var_dengue,
                 sarima_chik   = mase_sarima_chik,
                 ets_chik      = mase_ets,
                 var_chik      = mase_var_chik),
    resumen           = resumen_avanzado
  ),
  here("data", "series", "econometria_avanzada.rds")
)

write.csv(
  resumen_avanzado,
  here("data", "series", "resumen_econometria_avanzada.csv"),
  row.names = FALSE,
  na = ""
)

cat("\n✅ Análisis econométrico avanzado completado\n")

