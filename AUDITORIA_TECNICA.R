# ============================================================
# AUDITORIA_TECNICA.R
# Reproduce las verificaciones del informe AUDITORIA_TECNICA.md
# No modifica ningún archivo del proyecto: solo lee y reporta.
#
# Ejecutar desde la raíz del proyecto:
#   Rscript --vanilla AUDITORIA_TECNICA.R
# ============================================================

options(width = 100)
suppressPackageStartupMessages({
  library(forecast); library(vars); library(rugarch); library(tseries)
})

hr <- function(titulo) cat("\n", strrep("=", 70), "\n", titulo, "\n", strrep("=", 70), "\n", sep = "")

# ------------------------------------------------------------
hr("1. SINTAXIS")
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  ok <- tryCatch({ parse(f); "OK" }, error = function(e) conditionMessage(e))
  cat(sprintf("  %-38s %s\n", basename(f), ok))
}

# ------------------------------------------------------------
hr("2. INTEGRIDAD DE LA BASE")
d <- readRDS("data/D_02/dengue_prep.rds")
ch <- readRDS("data/D_02/chikungunya_prep.rds")
cat("Dengue:        ", nrow(d), "filas x", ncol(d), "columnas\n")
cat("Chikungunya:   ", nrow(ch), "filas x", ncol(ch), "columnas\n")
cat("Duplicados de 'consecutive' (dengue):", sum(duplicated(d$consecutive)), "\n")
cat("Filas totalmente duplicadas (dengue):", sum(duplicated(d)), "\n")

cat("\nCasos por año (dengue):\n"); print(table(d$anio))
cat("\ncon_fin (dengue):\n"); print(table(d$con_fin, useNA = "ifany"))
cat("fec_def NA (dengue):", sum(is.na(d$fec_def)), "de", nrow(d), "\n")

fl <- grep("^na_estructural_", names(d), value = TRUE)
cat("\nCobertura de los flags 'estructurales' (dengue):\n")
print(round(sapply(d[fl], function(x) sum(x, na.rm = TRUE)) / nrow(d) * 100, 2))

cat("\nDepartamentos en cod_dpto_o:", length(unique(d$cod_dpto_o)),
    "| incluye Bogotá (11)?", "11" %in% d$cod_dpto_o, "\n")
cat("Filas con 'BOGOT' en departamento_ocurrencia:",
    sum(grepl("BOGOT", d$departamento_ocurrencia, ignore.case = TRUE)), "\n")

# ------------------------------------------------------------
hr("3. COBERTURA TEMPORAL DE LAS SERIES")
sd <- readRDS("data/D_07/serie_dengue.rds")
sc <- readRDS("data/D_07/serie_chik.rds")
for (nm in c("dengue", "chik")) {
  s <- if (nm == "dengue") sd else sc
  nac <- subset(s, dpto == "NACIONAL")
  meses <- format(seq(min(as.Date(nac$fecha_sem)), max(as.Date(nac$fecha_sem)), by = "month"), "%Y-%m")
  presentes <- unique(format(as.Date(nac$fecha_sem), "%Y-%m"))
  cat(sprintf("%-12s semanas con casos: %3d | meses: %2d de %2d | meses sin casos: %s\n",
              nm, nrow(nac), length(presentes), length(meses),
              paste(setdiff(meses, presentes), collapse = ", ")))
}

# ------------------------------------------------------------
hr("4. MODELO SARIMA (08)")
e <- readRDS("data/D_08/econometria_sarima.rds")
for (nm in names(e)) {
  x <- e[[nm]]
  cat("\n--", nm, "--\n")
  cat("orden:", paste(arimaorder(x$fit), collapse = ","),
      "| lambda:", round(x$lambda_bc, 4), "| AIC:", round(AIC(x$fit), 2), "\n")
  cat("observaciones:", length(x$serie_bc),
      "| inicio:", paste(start(x$serie_bc), collapse = "-"),
      "| fin:", paste(end(x$serie_bc), collapse = "-"), "\n")
  cat("rango historia ($x):", round(range(x$forecast$x), 1),
      "| rango pronóstico ($mean):", round(range(x$forecast$mean), 0), "\n")
  cat("rango IC inferior:", round(range(x$forecast$lower), 1),
      "| IC negativo?", any(x$forecast$lower < 0), "\n")
}

# ------------------------------------------------------------
hr("5. GARCH (09)")
a <- readRDS("data/D_09/econometria_avanzada.rds")
g <- a$garch_dengue
mc <- g@fit$matcoef
ar <- mc[grep("^ar", rownames(mc)), 1]
cat("Matriz de coeficientes:\n"); print(round(mc, 5))
cat("\nSuma de coeficientes AR:", round(sum(ar), 4), "\n")
cat("Módulos de las raíces AR:", round(Mod(polyroot(c(1, -ar))), 4), "\n")
cat("Persistencia (alpha1 + beta1):", round(mc["alpha1", 1] + mc["beta1", 1], 4), "\n")
cat("Convergencia:", g@fit$convergence, "\n")
r <- as.numeric(residuals(g))
cat("ACF de residuos (lags 1-3):", round(acf(r, lag.max = 3, plot = FALSE)$acf[-1], 3), "\n")
cat("Ljung-Box residuos  p:", format.pval(Box.test(r, lag = 12, type = "Ljung-Box")$p.value, digits = 3), "\n")
cat("Ljung-Box residuos² p:", format.pval(Box.test(r^2, lag = 12, type = "Ljung-Box")$p.value, digits = 3), "\n")

cat("\nComparación de AIC:\n")
aic_sarima <- as.numeric(AIC(e$dengue$fit))
aic_garch_n <- as.numeric(infocriteria(g)[1])
cat("  SARIMA (total)          :", round(aic_sarima, 2), "\n")
cat("  GARCH  (por observación):", round(aic_garch_n, 4), "\n")
cat("  GARCH  (escala total)   :", round(aic_garch_n * g@model$modeldata$T, 2), "\n")
cat("  'Mejora' que reporta 09 :", round(aic_sarima - aic_garch_n, 2),
    "| diferencia real:", round(aic_garch_n * g@model$modeldata$T - aic_sarima, 2), "\n")

# ------------------------------------------------------------
hr("6. VAR, GRANGER Y ENCOMPASSING (09)")
cat("Rezago:", a$var$p, "| observaciones:", a$var$obs,
    "| módulo máximo de raíces:", round(max(roots(a$var)), 4), "\n")
cat("Granger chik -> dengue p:", round(a$granger$chik_a_dengue$Granger$p.value, 4), "\n")
cat("Granger dengue -> chik p:", round(a$granger$dengue_a_chik$Granger$p.value, 4), "\n")

pred_s <- as.numeric(fitted(e$dengue$fit))
pred_v <- as.numeric(fitted(a$var)[, "dengue"])
y <- as.numeric(e$dengue$serie_bc)
n_min <- min(length(pred_s), length(pred_v))

m1 <- lm(y ~ ps + pv, data = data.frame(y = y[1:n_min], ps = pred_s[1:n_min], pv = pred_v[1:n_min]))
p <- a$var$p
idx <- (p + 1):length(y)
m2 <- lm(y ~ ps + pv, data = data.frame(y = y[idx][1:length(pred_v)],
                                        ps = pred_s[idx][1:length(pred_v)],
                                        pv = pred_v))
cat("\nCoeficiente de pred_var:\n")
cat("  desalineado (como en 09):", round(coef(m1)["pv"], 5),
    "| p:", format.pval(summary(m1)$coefficients["pv", 4], digits = 3), "\n")
cat("  alineado correctamente  :", round(coef(m2)["pv"], 5),
    "| p:", format.pval(summary(m2)$coefficients["pv", 4], digits = 3), "\n")
cat("  cor(y, pred_var):", round(cor(y[1:n_min], pred_v), 3), "->",
    round(cor(y[idx][1:length(pred_v)], pred_v), 3), "\n")

# ------------------------------------------------------------
hr("7. ARTEFACTOS")
cat("Figuras con 0 bytes:\n")
figs <- list.files("graficas", pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
if (length(figs)) print(figs[file.info(figs)$size == 0]) else cat("  (ninguna)\n")
cat("\n¿09 guarda figuras? (ggsave/png en el script):",
    any(grepl("ggsave|png\\(", readLines("R/09_complemento_econometria.R"))), "\n")

cat("\nAuditoría finalizada.\n")

