# Auditoría técnica y metodológica
## Proyecto *Comportamiento del dengue y la chikungunya en Colombia (2020–2025)*

| | |
|---|---|
| **Fecha** | 21 de octubre de 2026 |
| **Commit auditado** | `9452c25` *(Actualización descriptiva y exploratoria)* |
| **Alcance** | 9 scripts R (`R/01`–`R/09`, 2.440 líneas), objetos `.rds` de `data/`, tablas de `data/series/` y figuras de `graficas/` |
| **Herramientas** | R 4.5.2 + `renv` del proyecto, `Rscript --vanilla`; verificación de datos crudos y de objetos ajustados |

---

## 1. Resumen ejecutivo

**Veredicto.** El pipeline de ingesta y depuración (`01`–`03`) es defendible y está bien estructurado. El capítulo de faltantes e imputación (`05`–`06`) es técnicamente correcto en su ejecución pero **desproporcionado y desconectado del análisis**. El capítulo econométrico (`08`–`09`) contiene errores de especificación y de comparación que **invalidan las conclusiones que hoy se leen en `resumen_econometria_avanzada.csv`**.

**No hay errores de sintaxis:** los 9 scripts pasan `parse()` sin errores. El problema es semántico y metodológico, no de compilación.

**Cómo leer el informe frente a los tres criterios del encargo:**

| Criterio | Sección | Contenido |
|---|---|---|
| Auditoría técnica (sintaxis, advertencias, malas prácticas) | §3 (H-01…H-08) y §5 (T-01…T-30) | errores que alteran resultados + tabla de 30 hallazgos de implementación |
| Validez metodológica (GARCH, estacionariedad, datos de riesgo) | §4 (M-01…M-12) | especificación, supuestos y comparaciones |
| Robustez y optimización | §7, §8, §9 | código corregido, pruebas exigibles y plan priorizado |

### 1.1 Lo que hoy no se puede sostener

1. **«El GARCH(1,1) mejora al SARIMA en 407,9 puntos de AIC».** Es un artefacto de unidades: el AIC del SARIMA es total (414,70) y el de `rugarch` es *por observación* (6,80). En la misma escala el GARCH vale 6,80 × 72 = **489,55**, es decir, es **74,85 puntos peor**.
2. **«El GARCH captura la volatilidad del dengue».** La ecuación de varianza no estima nada: α₁ = 7,4 × 10⁻²¹ (p = 1,000) y β₁ = 0,193 (p = 0,721). Además la media AR(4) es **explosiva** (Σφ = 1,0068; raíz de módulo 0,991) y los residuos siguen autocorrelacionados (ρ₁ = 0,67). El modelo no es un GARCH útil: es un AR(4) inestable con varianza casi constante.
3. **«SARIMA+GARCH tiene MAPE 5,63 %».** Esa cifra es una **copia literal** del MAPE del SARIMA (`mape_garch_dengue <- mape_sarima_dengue`, línea 242). Nunca se calculó.
4. **«El modelo de chikungunya pronostica».** La serie mensual tiene 4 meses sin ningún caso que `count()` elimina; `ts()` los trata como consecutivos y **desplaza todas las etiquetas temporales 1 a 4 meses**. El ARIMA(0,1,0) termina etiquetado en 2025-08 cuando el dato real termina en 2025-12, y su MAPE es `Inf`.
5. **«El VAR muestra relaciones dengue–chikungunya».** VAR en **niveles no estacionarios**, 64 observaciones, sin estacionalidad y con 18 parámetros por ecuación. Granger: p = 0,104 y p = 0,257 (no hay evidencia). El *test de encompassing* que lo «confirma» compara series **desfasadas exactamente 4 meses**.
6. **«Las figuras de `08` son válidas».** `forecast_dengue.png` dibuja la historia en escala Box-Cox (35–97) y el pronóstico en casos (6.742–12.237): la serie histórica aparece como **una línea plana en cero**. Ver §7.4.
7. **«Se hizo imputación múltiple».** m = 3, `maxit` = 2, `ntree` = 5 y `complete(imp, 1)`: se usa **una sola** imputación, nunca se hace *pooling* y la incertidumbre no se propaga. Además **solo 2 variables tenían faltantes reales** (estrato y `cod_ase`) y ninguna entra al análisis de `07`–`09`.

### 1.2 Fortalezas que conviene conservar

- Estructura modular `01 → 09` con `here()`, `renv.lock` y `.gitignore` correcto (no se versionan datos pesados ni la librería).
- **Marcado de NA estructurales** (`02`): idea correcta y poco frecuente en trabajos de pregrado.
- **Visualización de tres estados** (presente / no aplica / faltante real) en `05`: es un aporte metodológico real y defendible ante el docente.
- Cero duplicados y cero identificadores repetidos en 748.037 registros (verificado).
- Diagnóstico pre-imputación documentado (Little, coocurrencia): buena práctica, aunque mal aplicada (M-05).
- `09` no modifica `08`: buena práctica de trazabilidad.

---

## 2. Método de auditoría y evidencia

1. `parse()` de los 9 archivos → **9/9 OK**.
2. Inspección de los objetos guardados (dimensiones, tipos, tablas de frecuencia, órdenes ARIMA, parámetros GARCH, raíces, ACF de residuos, escalas de pronóstico, AIC).
3. Verificación de integridad de los datos crudos y procesados (duplicados, años, semanas, códigos DIVIPOLA, fechas).
4. Re-ejecución de los cálculos clave con alineación correcta para cuantificar el sesgo introducido.

### 2.1 Tabla de evidencia

| # | Verificación | Resultado |
|---|---|---|
| 1 | `parse()` de los 9 scripts | 9/9 sin errores de sintaxis |
| 2 | Dengue: filas / duplicados de `consecutive` | 748.037 / 0 |
| 3 | Chikungunya: filas | 419 |
| 4 | Semanas con al menos un caso (dengue / chik) | 314 / **197** de 314 |
| 5 | Meses con registro (dengue / chik) | 72 / **68** de 72 |
| 6 | Dengue: orden ARIMA y λ | (4,1,0)(1,0,0)[12], λ = 0,3293, AIC = 414,70 |
| 7 | Chik: orden ARIMA y fin de serie | (0,1,0), λ = 0,5379, fin etiquetado **2025-08** (real: 2025-12) |
| 8 | GARCH: α₁ / β₁ | 7,4 × 10⁻²¹ (p = 1,000) / 0,193 (p = 0,721) |
| 9 | GARCH: Σ coeficientes AR(4) | **1,0068** → raíz 0,991 (no estacionario) |
| 10 | GARCH: ACF de residuos (lag 1) | **0,673**; Ljung-Box p = 6,5 × 10⁻⁷ |
| 11 | AIC SARIMA vs GARCH en la misma escala | 414,70 vs **489,55** |
| 12 | VAR: rezago / observaciones / Granger | 4 / 64 / p = 0,104 y 0,257 |
| 13 | Encompassing: coeficiente de `pred_var` | 0,0002 (desalineado) → **0,0005** (alineado) |
| 14 | Variables con faltantes reales (dengue) | `estrato` (34.310) y `cod_ase` (17.525) |
| 15 | Flags «estructurales» | `fec_def` 100 %, `sem_ges` 99,5 %, `nom_grupo` 98,6 % |
| 16 | Departamentos en `cod_dpto_o` | 33; **falta Bogotá D.C. (11)** |
| 17 | `forecast_dengue.png` | historia 35–97 (escala BC) vs pronóstico 6.742–12.237 |
| 18 | `var_irf.png` | **0 bytes** (archivo vacío) |

---

## 3. Hallazgos críticos (invalidan resultados ya producidos)

### H-01 · El MAPE del «SARIMA+GARCH» nunca se calculó
**Archivo:** `R/09_complemento_econometria.R:242` · **Severidad: Crítica (integridad)**

```r
mape_garch_dengue <- mape_sarima_dengue   # ← copia, no cálculo
```

La tabla final (`resumen_econometria_avanzada.csv`) reporta `Dengue, SARIMA+GARCH, 5.63` y `Dengue, SARIMA, 5.63`: son el mismo número por asignación. Un GARCH no cambia la media condicional, así que el MAPE *podría* coincidir, pero entonces debe **calcularse** y decirse explícitamente; hoy es una cifra presentada como resultado sin haber sido estimada.

**Corrección.** Calcular el MAPE (mejor MAE/RMSE) con la media condicional del GARCH en escala original, o eliminar la fila y comparar solo la parte de varianza (VaR, log-score).

---

### H-02 · Comparación de AIC entre modelos no comparables (y en unidades distintas)
**Archivo:** `R/09_complemento_econometria.R:152–158` · **Severidad: Crítica**

```r
aic_sarima <- AIC(res_econ$dengue$fit)          # AIC TOTAL  = 414,70
aic_garch  <- infocriteria(fit_garch)[1]        # AIC POR OBS = 6,7993
cat("Mejora:", round(aic_sarima - aic_garch, 2))  # 407,90  ← artefacto
```

Tres problemas acumulados:
1. **Unidades:** `infocriteria()` devuelve AIC/n. Multiplicando por T = 72 → 489,55.
2. **Verosimilitudes distintas:** SARIMA gaussiano con d = 1 y D = 0; GARCH con t-Student y sin diferenciar. Los AIC no son comparables entre familias de distribución ni entre transformaciones distintas.
3. **Muestras distintas:** el SARIMA pierde una observación por diferenciación; el GARCH usa 72.

**Corrección.** No comparar AIC entre un SARIMA y un GARCH. Si se insiste, fijar la misma muestra, la misma transformación, la misma distribución y escalar ambos AIC a total:

```r
aic_total <- function(fit_rugarch) {
  as.numeric(infocriteria(fit_rugarch)[1]) * fit_rugarch@model$modeldata$T
}
```

Para comparar capacidad predictiva, usar evaluación **fuera de muestra** (MAE/RMSE/MASE + Diebold-Mariano), no AIC.

---

### H-03 · La especificación GARCH es inválida (y el resultado lo demuestra)
**Archivo:** `R/09_complemento_econometria.R:127–133` · **Severidad: Crítica**

```r
spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(4, 0), include.mean = FALSE),
  distribution.model = "std"
)
```

| Problema | Evidencia numérica | Consecuencia |
|---|---|---|
| `include.mean = FALSE` sobre una serie Box-Cox con media ≈ 60 | rango 35,4–96,8 | el AR(4) absorbe el nivel en lugar de modelarlo |
| AR(4) no justificado y sin estacionalidad | Σφ = **1,0068**; raíz 0,991 | **proceso explosivo**; `ugarchfit` no impone estacionariedad |
| Ecuación de varianza sin señal | α₁ = 7,4 × 10⁻²¹ (p = 1,000); β₁ = 0,193 (p = 0,721) | **no hay efecto ARCH estimado**: el «GARCH» es varianza constante |
| Residuos no blanqueados | ρ₁ = 0,673; Ljung-Box p = 6,5 × 10⁻⁷ | los errores estándar y los p-valores no son válidos |
| Desviación de residuos estandarizados | sd = 1,28 (debería ser ≈ 1) | el modelo subestima la escala de la varianza |
| 72 observaciones mensuales | AR(4) + ω + α + β + ν = 8 parámetros | potencia insuficiente para GARCH |

**Diagnóstico de fondo.** El GARCH se aplica a una serie **mensual** de conteos, agregada, fuertemente estacional y con una media que no es constante. GARCH supone una serie de media (casi) cero, sin autocorrelación y con agrupamiento de volatilidad; nada de eso se verificó. Al agregar a mensual se destruye justamente la variabilidad semanal donde podría existir el ARCH (314 semanas disponibles).

**Corrección mínima** (sobre la serie **semanal**, en `log1p`, con media estacional y sobre-dispersión):

```r
y  <- ts(log1p(casos_semanales), frequency = 52)
t  <- seq_along(y)
X  <- cbind(sin1 = sin(2*pi*t/52), cos1 = cos(2*pi*t/52),
            sin2 = sin(4*pi*t/52), cos2 = cos(4*pi*t/52))

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE,
                        external.regressors = X),
  distribution.model = "std"
)
fit <- ugarchfit(spec, data = y, solver = "hybrid")
stopifnot(convergence(fit) == 0)

# Diagnóstico obligatorio
persistencia <- fit@fit$matcoef["alpha1", 1] + fit@fit$matcoef["beta1", 1]
Box.test(residuals(fit, standardize = TRUE), lag = 24, type = "Ljung-Box", fitdf = 3)
Box.test(residuals(fit, standardize = TRUE)^2, lag = 24, type = "Ljung-Box")
```

Y, para datos de conteo, contrastar con un modelo específico de conteos (ver §7.7).

---

### H-04 · La serie de chikungunya está desalineada y le faltan ceros
**Archivos:** `R/08_econometria.R:26–34`, `R/09_complemento_econometria.R:32–41`, `R/07` · **Severidad: Crítica**

```r
preparar_mensual <- function(serie) {
  serie |> filter(dpto == "NACIONAL") |>
    mutate(mes = floor_date(fecha_sem, "month")) |>
    group_by(mes) |> summarise(casos = sum(casos, na.rm = TRUE)) |>
    arrange(mes) |> filter(!is.na(casos))      # ← elimina meses sin casos
}
```

Chikungunya tiene **4 meses sin ningún caso** (2020-10, 2025-06, 2025-08, 2025-09). `summarise()`/`count()` no generan filas para ellos y `ts()` asume 68 meses consecutivos desde 2020-01: **desde 2020-10 cada observación queda etiquetada un mes antes**, y el desfase crece hasta 4 meses al final. Verificado: `end(serie_bc)` = **2025-08** cuando el último dato real es **2025-12**.

Consecuencias: el ADF (p = 0,01), el λ de Box-Cox, el ARIMA(0,1,0), el ETS y el VAR de chikungunya se estiman sobre una serie con el calendario corrido y con los ceros suprimidos. Todo el capítulo de chikungunya debe rehacerse.

**Corrección.** Rejilla completa y relleno con ceros (el cero es un dato, no un faltante):

```r
completar_mensual <- function(serie, evento) {
  serie |>
    filter(dpto == "NACIONAL") |>
    mutate(mes = floor_date(fecha_sem, "month")) |>
    group_by(mes) |>
    summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") |>
    complete(mes = seq(min(mes), max(mes), by = "month"),
             fill = list(casos = 0)) |>
    arrange(mes) |>
    mutate(evento = evento)
}
```

Esto también obliga a cambiar el tratamiento estadístico: con ceros, `BoxCox.lambda()` emite *warning* (Guerrero exige datos estrictamente positivos) y `ETS(M,...)` puede fallar. Es la señal de que el modelo correcto no es gaussiano-transformado sino de conteos (§7.7).

---

### H-05 · El test de encompassing compara series desfasadas 4 meses
**Archivo:** `R/09_complemento_econometria.R:197–214` · **Severidad: Crítica**

```r
pred_sarima <- fitted(res_econ$dengue$fit)     # 72 valores
pred_var    <- fitted(fit_var)[, "dengue"]     # 64 valores (VAR p = 4)
y_real      <- as.numeric(res_econ$dengue$serie_bc)[1:n_min]
```

El VAR pierde las **primeras 4** observaciones (p = 4). Al tomar `[1:64]` se empareja 2020-01…2025-04 del SARIMA con 2020-05…2025-04 del VAR: **desfase exacto de 4 meses**.

Evidencia de la magnitud del sesgo (recalculado):

| | Coef. `pred_var` | R² | cor(y, pred_var) |
|---|---|---|---|
| Como está escrito | 0,0002 (p = 0,002) | 0,949 | 0,636 |
| Alineado correctamente | **0,0005** (p < 0,001) | 0,955 | **0,933** |

Además, usar valores ajustados *in-sample* de dos modelos estimados sobre la misma muestra no es un test de encompassing válido (regresores generados; sobreajuste). El estándar es evaluación fuera de muestra + Diebold-Mariano.

```r
p   <- fit_var$p
idx <- (p + 1):length(y_bc)
n   <- length(pred_var)
dat <- data.frame(y  = y_bc[idx][1:n],
                  ps = pred_sarima[idx][1:n],
                  pv = pred_var)
summary(lm(y ~ ps + pv, data = dat))
```

---

### H-06 · La tabla resumen tiene columnas hardcodeadas
**Archivo:** `R/09_complemento_econometria.R:265–287` · **Severidad: Crítica (integridad)**

```r
util = c(util_sarima_dengue, TRUE, FALSE,
         util_sarima_chik,   TRUE, FALSE)
```

Resultado en el CSV: **Chikungunya / ETS → MAPE 73,76 % → `util = TRUE`**, cuando la regla del propio script es MAPE < 50 %. Y **Chikungunya / SARIMA → MAPE `Inf` → `util = FALSE`**, aunque su AIC (280,75) es mejor que el del ETS (471,96). La columna no mide utilidad: está escrita a mano.

**Corrección.** Calcular la utilidad con una regla única y explícita (p. ej. MASE < 1 frente al *seasonal naive*) y aplicarla a las 6 filas.

---

### H-07 · El MAPE se calcula en la escala Box-Cox (y por eso da `Inf`)
**Archivo:** `R/09_complemento_econometria.R:229–239` · **Severidad: Crítica**

```r
mape_sarima_dengue <- mean(abs((serie_bc - fitted(fit)) / serie_bc)) * 100
```

- Un 5,63 % **no es un error porcentual en casos**: es el error en la escala λ = 0,329.
- En chikungunya `serie_bc` contiene **0** (los meses con exactamente 1 caso dan BoxCox(1) = 0) → división por cero → `MAPE = Inf`, que aparece en la tabla final.
- MAPE no es una métrica válida para conteos con ceros.

**Corrección.** MAE, RMSE y MASE en la escala original de casos:

```r
accuracy(fc, test)   # forecast::accuracy ya devuelve MAE/RMSE/MAPE/MASE
```

---

### H-08 · El pronóstico mezcla dos escalas y la figura queda inservible
**Archivo:** `R/08_econometria.R:140–151` · **Severidad: Crítica (figura y objeto guardado)**

```r
forecast_orig$mean  <- InvBoxCox(forecast_sarima$mean,  bc_lambda)
forecast_orig$lower <- InvBoxCox(forecast_sarima$lower, bc_lambda)
forecast_orig$upper <- InvBoxCox(forecast_sarima$upper, bc_lambda)
plot(forecast_orig, ...)
```

`$x` (la historia que dibuja `plot.forecast`) **no se transforma**. Verificado: `range($x)` = 35,4–96,8 (escala Box-Cox) frente a `range($mean)` = 6.742–12.237 (casos). La figura `graficas/econometria/forecast_dengue.png` muestra la serie histórica como una línea plana en 0 y un abanico de pronóstico enorme. **Esa figura no puede ir al informe.**

Adicionalmente, la back-transformación de los límites gaussianos **no produce un intervalo predictivo válido** para conteos: el IC95 % de noviembre de 2026 es **[63, 42.883]** con media 7.311, cuando el mínimo histórico mensual es 2.226. En chikungunya los límites inferiores son **negativos** (hasta −24,8).

**Corrección.** Transformar el objeto completo (no solo `mean/lower/upper`) y validar que los límites tengan sentido:

```r
fc <- forecast(fit_sarima, h = 12)
fc$x      <- InvBoxCox(fc$x,      lambda)
fc$fitted <- InvBoxCox(fc$fitted, lambda)
fc$mean   <- InvBoxCox(fc$mean,   lambda)
fc$lower  <- InvBoxCox(fc$lower,  lambda)
fc$upper  <- InvBoxCox(fc$upper,  lambda)
fc$lower[fc$lower < 0] <- 0
fc$upper[fc$upper < 0] <- 0
```

> Nótese que la media sigue sin corrección de sesgo (ver §7.4): el pronóstico puntual es, en rigor, una mediana aproximada.

---

## 4. Validez metodológica (hallazgos altos)

### M-01 · Estrategia de estacionariedad incoherente
`08:65–84` aplica ADF a la serie original (p = 0,409), a la diferenciada (p = 0,012) y a la **diferenciada estacionalmente** (p = 0,066), pero el modelo final es **(4,1,0)(1,0,0)**, es decir, **D = 0**: el resultado del test estacional (que sugiere D = 1 al 5 %) se ignora. Además:
- solo se usa ADF, que tiene baja potencia frente a raíces cercanas a la unidad y es sensible a la especificación de tendencia;
- no se usa KPSS ni Phillips-Perron para contrastar;
- la prueba correcta de raíz unitaria **estacional** no es `adf.test(diff(...,12))` sino OCSB, Canova-Hansen o HEGY;
- `adf.test` usa `k` por defecto y no se reporta.

**Corrección.** Triplete ADF + KPSS + PP sobre niveles y primeras diferencias; `forecast::nsdiffs()` con OCSB y `urca::ur.za`/`hegy` para estacionalidad; documentar la decisión d/D con la tabla de tests.

### M-02 · Box-Cox sin justificación y con problemas de interpretación
λ = 0,329 (dengue) y 0,538 (chik) estimados por Guerrero. Para conteos, la transformación de Box-Cox no estabiliza la varianza Poisson (la raíz cuadrada o `log1p` son las estándar) y λ ≈ 1/3 no tiene una interpretación natural en este contexto: la transformación no está motivada por la estructura de varianza de los datos, sino por la conveniencia del ajuste. La estimación con ceros emite el *warning* «Guerrero's method ... strictly positive data» (verificado).

**Corrección.** Usar `log1p` (interpretable, admite ceros, equivale a λ → 0) o, mejor, un modelo de conteos con función de enlace log. Reportar el análisis de sensibilidad a λ.

### M-03 · VAR en niveles no estacionarios y sin estacionalidad
`09:75–102`:
- las dos series entran **en niveles** (conteos crudos, ADF p = 0,409 en dengue) → riesgo de **regresión espuria** y de **causalidad de Granger espuria**;
- 64 observaciones efectivas con p = 4 → 9 coeficientes por ecuación (18 en total, más las varianzas);
- sin dummies estacionales: los residuos conservan estacionalidad de 12 meses y el portmanteau la detectaría;
- `VARselect` con AIC sobre 64 obs sobreajusta (elige el máximo razonable);
- el resultado (Granger p = 0,104 y 0,257) **no permite afirmar nada** y, aun así, se presenta como hallazgo.

**Corrección.** VAR sobre `log1p` (o tasas por 100.000 hab.), con dummies mensuales o VAR estacional, contraste de estacionariedad y de estabilidad, y **CFC pre-blanqueada** antes de interpretar Granger. Con 72 meses, el VAR bivariado es marginal: mejor un modelo de rezagos distribuidos (ARDL) con corrección de errores.

### M-04 · Diagnósticos de residuos inadecuados
`08:120–134` usa Jarque-Bera (aceptable), **Breusch-Pagan** y **Durbin-Watson**:
- DW está sesgado hacia 2 cuando hay retardos de la endógena: **no es válido** aquí;
- BP asume homocedasticidad condicional y **no detecta ARCH**, que es justamente el objeto del proyecto;
- Ljung-Box con `lag = 12` sin `fitdf` y sin lag estacional (24) para datos mensuales.

**Corrección.** `Box.test(..., lag = 24, fitdf = p + q + P + Q)`, `FinTS::ArchTest` (o Ljung-Box sobre residuos²), `forecast::checkresiduals`, y pruebas de no linealidad (BDS) si se quiere justificar GARCH.

### M-05 · El test de Little se aplica sobre datos mal preparados
`05:459–478` corre `naniar::mcar_test` sobre 71 columnas en las que **la mayor parte de los NA son estructurales por construcción**: `fec_def` 100 %, `sem_ges` 99,5 %, `nom_grupo` 98,6 %, `fec_hos` 59,6 %. Un test de MCAR sobre esa matriz **rechaza por diseño** (p = 0), no porque el mecanismo sea MAR. Además se ejecuta sobre una submuestra aleatoria de 5.000 filas (el p-valor ya no corresponde a la base) e incluye una columna 100 % vacía (`fec_def`) que el test no puede utilizar.

Los **faltantes reales son solo dos**: `estrato` (4,59 %) y `cod_ase` (2,34 %). El test debía correrse sobre esas dos variables (con predictores relevantes) y no sobre el universo.

### M-06 · Imputación múltiple por debajo del estándar
`06:295–304`:

```r
res_dengue <- imputar_evento(..., m = 3, maxit = 2, ntree = 5)
```

- `maxit = 2` **no permite convergencia** de las cadenas;
- `m = 3` es insuficiente (regla práctica: m ≥ 20 cuando el porcentaje de faltantes es alto; mínimo 5 con faltantes bajos);
- `ntree = 5` es un bosque degenerado (el valor por defecto es 500);
- el propio diagnóstico (`diagnostico_imputacion.csv`) recomienda «m = 5, maxit = 3» para dengue: **lo documentado y lo ejecutado no coinciden**;
- nunca se revisa la convergencia (`plot(imp)`), solo se imprime `imp$loggedEvents`.

### M-07 · Se descarta la imputación múltiple al primer paso
`06:244` → `df_imputado <- complete(imp, 1)`. Se guardan los objetos `mids` «por si se quiere hacer pooling», pero **nunca se usan**. Toda la inferencia de `07`–`09` trabaja con una sola imputación, ignorando la incertidumbre de imputación. Para un trabajo con capítulo de imputación, esto es una contradicción interna.

### M-08 · El modelo de imputación no tiene predictores útiles
`06:118–136` selecciona **solo** las variables a imputar (`estrato`, `cod_ase`) y descarta `semana`, `departamento`, `edad`, `tip_cas`, régimen de salud. Es decir: se imputa el estrato sin usar la información contextual que lo predice. Un modelo de imputación debe incluir las variables del análisis (o al menos un subconjunto) para ser válido bajo MAR.

### M-09 · Semana epidemiológica aproximada con el 1 de enero
`02:39` y `07:49–51`:

```r
fecha_epi <- make_date(anio) + weeks(semana - 1)
```

La semana epidemiológica (MMWR/PAHO, la que usa el INS) empieza **domingo** y la semana 1 es la que contiene el **4 de enero**. Anclar en el 1 de enero introduce un sesgo sistemático. Verificado: `fecha_epi − fec_not` tiene **mediana −7 días** y **41 % de los registros con |diferencia| > 7 días**. Como la serie mensual se construye con `floor_date()`, ese sesgo desplaza la fase estacional y mueve casos entre meses.

**Corrección.**

```r
# Semana 1 = la que contiene el 4 de enero (MMWR/PAHO)
ancla  <- as.Date(paste0(anio, "-01-04"))
inicio <- ancla - as.integer(format(ancla, "%w"))   # %w: 0 = domingo
fecha_epi <- inicio + (semana - 1) * 7
```

Validar contra el calendario oficial de semanas del INS y documentar la convención.

### M-10 · La serie «NACIONAL» no incluye Bogotá
`04:343`, `07:26`. En dengue, `cod_dpto_o` tiene 33 códigos y **no aparece el 11 (Bogotá D.C.)**; la columna `departamento_ocurrencia` **no tiene ni una fila de Bogotá** (verificado: 0 registros). En chikungunya sí aparece el 11 (2 casos). Es decir:

- la serie rotulada `NACIONAL` es **Colombia sin Bogotá** (y sin el peso epidemiológico de la capital);
- las coberturas de los dos eventos **no son comparables**, lo que sesga el VAR y cualquier comparación dengue–chikungunya.

Hay que verificar el filtro de descarga en el portal del INS o declararlo como limitación explícita.

### M-11 · `con_fin` es constante y `fec_def` está 100 % vacía
Verificado en dengue: `con_fin` solo toma el valor **«1»** en los 748.037 registros y `fec_def` es **NA en el 100 %**. Con la regla del script (`is.na(fec_def) & con_fin == "1"`) **todos** los registros quedan marcados como «estructurales»: el 100 % de la tabla sale en la columna `na_estructural_fec_def`. Eso no es un faltante estructural: es una **variable sin información** (o un extracto que excluye las defunciones). Cualquier análisis de letalidad/gravedad es hoy imposible, y la auditoría de `03` no lo advierte porque el flag la «justifica».

Además hay una **contradicción entre comentario y código**: el comentario de `02:58` dice «fec_def: solo aplica si con_fin = 2» y el código usa `con_fin == "1"`. Verificar el diccionario del INS (1 = vivo, 2 = muerto) y alinear ambos.

### M-12 · No hay evaluación predictiva fuera de muestra
Todo el capítulo de pronóstico (`08:137`, `09:113`) ajusta y pronostica sobre la misma muestra. El MAPE es *in-sample*, el GARCH se compara *in-sample* y no hay *backtesting*. Para un trabajo cuyo objetivo es «predecir el comportamiento», esto es el vacío metodológico más serio después de H-03.

**Corrección.** Partición temporal (p. ej. 2020–2024 entrenamiento, 2025 prueba), origen móvil (*rolling origin*), `accuracy()` contra `snaive()` y `meanf()`, y Diebold-Mariano:

```r
dm.test(e1, e2, h = 12, power = 2)
```

---

## 5. Hallazgos técnicos y de robustez

| ID | Archivo:línea | Hallazgo | Corrección |
|---|---|---|---|
| T-01 | `02:41` | `select(-any_of(c("gru_pob","cbmete")))`: la columna real es **`cbmte`** (así viene en el xlsx). `cbmete` no existe → no se elimina nada y `cbmte` sobrevive (y luego se marca 100 % estructural). | Corregir el nombre a `cbmte` o eliminarlo explícitamente. |
| T-02 | `02:113–116` | `na_estructural_nom_grupo = is.na(nom_grupo)` marca el **98,6 %** de la base como estructural. No es estructural: es faltante real masivo o columna inútil. | Condicionar a los indicadores `gp_*` o eliminar la columna. |
| T-03 | `02:76–116` | `if (tiene(...)) ... else FALSE` dentro de `mutate()`: funciona por reciclado, pero mezcla escalar y vector y devolvería `NA` si el predicado fuese `NA`; `sum()` sobre esos flags (`02:128`) devolvería `NA` (hoy no ocurre). | `sum(..., na.rm = TRUE)` y `tidyr::replace_na()`. |
| T-04 | `02:16–20` | `parsear_fecha()` no reporta cuántas fechas no pudo convertir (hoy 0, pero es silencioso). | Contar `NA` antes/después y `warning()` si aumenta. |
| T-05 | `05:171–195` | `patron_faltantes()` grafica **todos** los patrones únicos (no un top-N) y la altura se calcula con `nrow(unique(...))` (`05:294`): figuras ilegibles y potencialmente gigantes. | Limitar a los 20–30 patrones más frecuentes (o al 95 % acumulado). |
| T-06 | `05:282` | `names(clasificar_celdas(df_evt))` dentro del bucle de categorías: recalcula una matriz de 748.037 × 71 caracteres **10 veces por evento**. | Calcular `estados` una vez fuera del bucle. |
| T-07 | `05:137–140` | `do.call(paste0, lapply(...))` es frágil; además se materializa el patrón para 748k filas. | `do.call(paste, c(mat, sep = ""))` o `apply(..., 1, paste, collapse="")`. |
| T-08 | `05:519` | Correlaciones con `NaN` si alguna columna es constante; `which()` ignora `NA` silenciosamente. | `cor(..., use = "pairwise.complete.obs")` + filtro de varianza 0. |
| T-09 | `06:257–263` | `n_na_antes` se calcula sobre `df_imp`, donde los estructurales **ya** se pusieron a `NA`: el denominador de `pct_imputado` incluye no-aplica. | Calcular antes de enmascarar y reportar por separado. |
| T-10 | `06:384–410` | Compara `na_orig` (solo faltantes reales) con `na_final` (todos los NA, incluidos los estructurales restaurados): las cifras no son comparables. | Separar reales y estructurales en el reporte. |
| T-11 | `06` (nombre de archivo) | `06_imputación.R` lleva tilde: genera problemas de portabilidad y de *encoding* (en `source()` desde otra *locale*, en `list.files()` con regex y en sistemas de archivos no UTF-8). | Renombrar a `06_imputacion.R`. |
| T-12 | `07:26–31` | `intersect(...)[1]` selecciona en silencio la primera columna disponible: si mañana se agrega otra candidata, cambia la prioridad sin aviso. El `if (is.na(col_dpto))` posterior sí maneja el caso vacío. | Hacer explícita la prioridad y registrar la columna usada. |
| T-13 | `07:43–44` | `grupo_edad` con `ifelse` anidado y solo dos grupos (Joven/Adulto), distinto de los grupos de `04`; los `NA` se descartan y la suma de las dos series no coincide con la nacional. | `cut()` con los mismos cortes de `04` y una categoría «Sin dato». |
| T-14 | `08:100–106` | `png(...)` + `par(mfrow)` + `dev.off()`: si `acf()` falla, el dispositivo queda abierto. | `on.exit(dev.off(), add = TRUE)` o `ggsave`. |
| T-15 | `08:120–134` | Ljung-Box sin `fitdf` y con `lag = 12` en serie mensual; DW inválido (ver M-04). | `lag = 24`, `fitdf = 5`, `checkresiduals()`. |
| T-16 | `09:105–110, 149–150, 184–186` | `plot()` sin dispositivo: en `Rscript` van a `Rplots.pdf`. Las cuatro figuras de `graficas/econometria_avanzada/` **no las produce el script** (no hay `ggsave`/`png` en `09`) y `var_irf.png` pesa **0 bytes**. | Estandarizar `ggsave(here(...))` y regenerar. |
| T-17 | `09:105` | `irf(..., boot = TRUE)` sin `set.seed()`: los intervalos de confianza no son reproducibles. | `set.seed(123)` antes. |
| T-18 | `09:23` | `readRDS()` sin verificar existencia: si `08` no se corrió, el error es críptico. | `stopifnot(file.exists(...))`. |
| T-19 | `04` | Ningún gráfico se guarda (`ggplot` sin `ggsave`): en `Rscript` se pierden y `graficas/` no tiene nada de `04`. | `ggsave()` para cada figura. |
| T-20 | `04:419–430` | `outliers_edad` aplica el criterio IQR a la edad: marca como «atípicos» a los mayores por construcción. No aporta. | Usar el criterio solo si se justifica; o reportar percentiles. |
| T-21 | `04:12–13` | `select <- dplyr::select; filter <- dplyr::filter` en el entorno global: parche que oculta el conflicto (y en `09` se resuelve de otra forma con `dplyr::select`). | `conflicted::conflict_prefer()` o prefijo explícito. |
| T-22 | `01:24–33` | `map(read_xlsx)` sin `col_types`: lectura por inferencia de 6 archivos × ~70 columnas de texto. En memoria el objeto ronda 1–2 GB (y `bind_rows` duplica temporalmente). | `col_types = "text"` y `col_select` con las columnas necesarias; o `arrow`/`data.table`. |
| T-23 | `01:25` | `set_names(str_extract(basename, "\\d{4}"))`: si un archivo no trae año, el id es `NA`. | `stopifnot(!anyNA(...))`. |
| T-24 | `08:39–47` | `ts()` con `start = c(año, mes)` asume serie mensual completa (ver H-04). | Rejilla completa antes de `ts()`. |
| T-25 | `09:46–53` | `full_join` + `filter(!is.na(dengue), !is.na(chik))` elimina el mes completo de **ambos** eventos si falta en uno: acorta la muestra del VAR y rompe la comparación con el SARIMA (que usa 72). | Rejilla completa con ceros antes del *join*. |
| T-26 | `03:26–30` | `summarise(across(everything(), ~ class(.)[1]))` pierde las clases múltiples (`POSIXct` + `POSIXt`). | `sapply(df, function(x) paste(class(x), collapse = "/"))`. |
| T-27 | `README:15, 59–65` | Fecha de descarga sin completar; sección «Cómo replicar» duplicada; texto residual `bhfbh`. | Limpiar. |
| T-28 | `renv` | La activación automática de `renv` no terminó en 20 min en ejecución no interactiva (tuve que usar `Rscript --vanilla` con `R_LIBS_USER`). Puede ser el entorno, pero conviene verificar `renv::status()`. | Documentar `Rscript --vanilla` para ejecución por lotes y re-snapshot. |
| T-29 | raíz del repo | No hay script maestro (`00_run_all.R`) ni `set.seed` global; el orden depende de ejecutar 01→09 a mano. | Script maestro + semilla. |
| T-30 | raíz del repo | `.RDataTmp` (11,8 MB) sin ignorar; `output/` vacío. | Añadir a `.gitignore`. |

---

## 6. Revisión archivo por archivo

**`01_carga_datos.R`** — Correcto y simple. Falta `col_types`, verificación de columnas mínimas y control de que cada archivo aporte su año. El `bind_rows(.id = "anio")` sobrescribe el `ANO` original con el año del nombre del archivo: hoy coinciden (verificado: 748.037/748.037), pero es una coincidencia frágil.

**`02_preparacion.R`** — La idea de marcar NA estructurales es la mejor del pipeline. Tiene tres defectos: el typo `cbmete` (T-01), la regla de `nom_grupo` (T-02) y la contradicción `con_fin` (M-11). Además `select(-any_of("gru_pob"))` elimina una columna que el resto del código sigue nombrando (`05`, `06`), lo que obliga a `intersect()` defensivos.

**`03_auditoria.R`** — Cumple su función, pero **no detecta** los dos problemas que sí aparecen en la base: `con_fin` constante y `fec_def` 100 % vacía (la columna vacía se lista, pero no se advierte que el flag estructural la está «justificando»). Añadir un chequeo de columnas con un único valor y de columnas 100 % NA cruzado con sus flags.

**`04_descrip_explo.R`** — Bien escrito, pero no guarda ninguna figura (T-19) y el análisis de atípicos en edad no aporta (T-20). La selección de la columna de departamento funciona (verificado: `cod_dpto_o` es character de 2 dígitos y `nombrar_dpto` cubre los 33 códigos), pero deja fuera a Bogotá (M-10).

**`05_faltantes.R`** — El diagnóstico de tres estados es valioso. Problemas: el test de Little mal aplicado (M-05), gráficos con todos los patrones (T-05), recálculo de la matriz de estados dentro del bucle (T-06). La lógica de decisión del método (`05:534–556`) es razonable pero **no se conecta programáticamente** con `06`: si el diagnóstico cambiara, `06` seguiría usando rf/pmm a mano.

**`06_imputación.R`** — Ejecuta bien lo que se propone, pero el planteamiento es desproporcionado: **1.081 líneas para imputar 2 variables** (`estrato`, `cod_ase`) que no entran en ningún modelo posterior. Y la configuración (m = 3, `maxit` = 2, `ntree` = 5) contradice el propio diagnóstico. Con los datos actuales, el capítulo podría reducirse a una tabla y una nota, o eliminarse.

**`07_serietemporal_conimputacion.R`** — El núcleo del problema de alineación (H-04). Además: `grupo_edad` binario (T-13), `fecha_sem` con anclaje incorrecto (M-09) y, si el código de departamento faltara, `count(dpto)` lo convertiría en una categoría `NA` más. Lo bueno: guarda CSV/RDS y una figura nacional.

**`08_econometria.R`** — Buena estructura de función, pero: serie incompleta (H-04), estacionariedad incoherente (M-01), diagnósticos inadecuados (M-04), pronóstico con escalas mezcladas (H-08), ACF/PACF sobre una diferenciación que el modelo final no usa (`08:104`, d = 1 vs D = 0), y ningún *split* de validación (M-12).

**`09_complemento_econometria.R`** — Es el script con más problemas críticos: AIC mal comparado (H-02), GARCH inválido (H-03), encompassing desalineado (H-05), tabla hardcodeada (H-06), MAPE en escala equivocada (H-07), VAR en niveles (M-03), figuras no reproducibles (T-16) y MAPE copiado (H-01). La separación «sin modificar `08`» es una buena práctica que conviene mantener al rehacerlo.

---

## 7. Correcciones prioritarias (código listo para usar)

### 7.1 Rejilla temporal completa (bloquea 07–09)

```r
completar_mensual <- function(serie, evento) {
  stopifnot(all(c("dpto", "fecha_sem", "casos") %in% names(serie)))
  serie |>
    dplyr::filter(dpto == "NACIONAL") |>
    dplyr::mutate(mes = lubridate::floor_date(fecha_sem, "month")) |>
    dplyr::group_by(mes) |>
    dplyr::summarise(casos = sum(casos, na.rm = TRUE), .groups = "drop") |>
    tidyr::complete(mes = seq(min(mes), max(mes), by = "month"),
                    fill = list(casos = 0)) |>
    dplyr::arrange(mes) |>
    dplyr::mutate(evento = evento)
}
```

### 7.2 Semana epidemiológica correcta

```r
fecha_epi_mmwr <- function(anio, semana) {
  ancla  <- as.Date(paste0(anio, "-01-04"))          # la semana 1 contiene el 4 de enero
  inicio <- ancla - as.integer(format(ancla, "%w"))  # domingo de esa semana
  inicio + (semana - 1) * 7
}
```

### 7.3 Comparación honesta de modelos

```r
# Nunca comparar AIC entre familias distintas. Si se hace, misma muestra y misma escala:
aic_total <- function(fit) as.numeric(infocriteria(fit)[1]) * fit@model$modeldata$T

# La comparación válida es fuera de muestra:
h <- 12; n <- length(y)
train <- window(y, end = time(y)[n - h])
test  <- window(y, start = time(y)[n - h + 1])

fc_sarima <- forecast(auto.arima(train), h = h)
fc_naive  <- snaive(train, h = h)
accuracy(fc_sarima, test)
accuracy(fc_naive,  test)
dm.test(test - fc_sarima$mean, test - fc_naive$mean, h = h, power = 2)
```

### 7.4 Figura de pronóstico en una sola escala

```r
fc <- forecast(fit_sarima, h = 12)
fc$x      <- InvBoxCox(fc$x,      lambda)
fc$fitted <- InvBoxCox(fc$fitted, lambda)
fc$mean   <- InvBoxCox(fc$mean,   lambda)
fc$lower  <- InvBoxCox(fc$lower,  lambda)
fc$upper  <- InvBoxCox(fc$upper,  lambda)
fc$lower[fc$lower < 0] <- 0        # pmax() corrompe los atributos ts
fc$upper[fc$upper < 0] <- 0
```

> La media back-transformada **no es** E[y]: falta la corrección de sesgo. `InvBoxCox(..., biasadj = TRUE)` exige `fvar` (la varianza del pronóstico, que `forecast.Arima` no siempre conserva); si no se dispone de ella, hay que declarar que el pronóstico puntual es la mediana aproximada, no la media.

### 7.5 VAR estacionario y correctamente alineado

```r
Y <- cbind(dengue = log1p(dengue_mensual$casos),
           chik   = log1p(chik_mensual$casos))
Y <- ts(Y, start = c(2020, 1), frequency = 12)

# estacionalidad determinista
mes <- factor(cycle(Y[, 1]))
X   <- model.matrix(~ mes)

VARselect(Y, lag.max = 6, type = "const", exogen = X[, -1])
fit_var <- VAR(Y, p = 2, type = "const", exogen = X[, -1])

# alinear antes de cualquier comparación
p    <- fit_var$p
idx  <- (p + 1):nrow(Y)
pred <- fitted(fit_var)[, "dengue"]
y    <- as.numeric(Y[idx, "dengue"])
```

### 7.6 Imputación múltiple bien hecha (si se conserva el capítulo)

```r
metodos <- make.method(datos)
metodos[!colSums(is.na(datos))] <- ""      # no imputar columnas completas

imp <- mice(datos, m = 20, maxit = 10, method = metodos,
            seed = 123, printFlag = FALSE)
plot(imp)                                   # convergencia: obligatorio
dens <- densityplot(imp); dens              # comparación de densidades

# analizar en cada imputación y combinar
resultados <- lapply(1:20, function(i) {
  d <- complete(imp, i)
  # ... modelo ...
})
```

### 7.7 Modelo de conteos (alternativa recomendada al GARCH gaussiano)

```r
# renv::install("tscount")
fit_ing <- tscount::tsglm(
  y,
  model = list(past_obs = 1, past_mean = 1,
               external = X_estacional),
  distr = "nbinom",        # overdispersión típica de datos epidemiológicos
  link  = "log"
)
summary(fit_ing)
```

---

## 8. Pruebas adicionales exigibles para el estándar del curso

| Ámbito | Prueba | Para qué |
|---|---|---|
| Estacionariedad | ADF + KPSS + Phillips-Perron; OCSB/Canova-Hansen; HEGY | decidir d y D con evidencia, no por defecto |
| Estacionalidad | `nsdiffs`, seasonplot, subseries, densidad espectral | justificar la parte estacional |
| Heterocedasticidad | Ljung-Box sobre residuos², `FinTS::ArchTest` | **antes** de imponer un GARCH |
| Especificación GARCH | Nyblom (estabilidad), sign bias, weighted LB, `convergence()` | validez de los parámetros |
| No linealidad | BDS sobre residuos del modelo lineal | justificar el cambio de modelo |
| Quiebres | `strucchange::breakpoints`, `tsoutliers` | 2020 (COVID) y 2024 (epidemia) son candidatos |
| Sobredispersión | `AER::dispersiontest` (Poisson vs. cuasi-Poisson) | elegir familia de conteos |
| Ceros | proporción de ceros vs. esperada; test de inflación | decidir NB vs. ZIP |
| Pronóstico | `accuracy()` fuera de muestra, `snaive` como *benchmark*, Diebold-Mariano, Mincer-Zarnowitz | validez predictiva real |
| Riesgo | VaR/ES del cuantil condicional, Kupiec POF, Christoffersen (independencia y cobertura) | conecta GARCH con «teoría de riesgo» |
| Series bivariadas | CCF pre-blanqueada, ARDL/cointegración (Engle-Granger o Johansen) | sustituir el Granger en niveles |

**Sobre la conexión con teoría de riesgo.** El proyecto se presenta como «econometría y teoría de riesgo». El puente natural es definir un **VaR epidemiológico**: el cuantil α de la distribución condicional de casos semanales, `VaR_α(t) = F_t⁻¹(α)`, y su *expected shortfall*. Eso permite usar exactamente las herramientas de riesgo (Kupiec, Christoffersen) sobre las alertas tempranas, que es un aporte mucho más defendible que comparar AIC entre un SARIMA y un GARCH.

---

## 9. Plan de trabajo priorizado

| Prioridad | Acción | Archivos | Esfuerzo |
|---|---|---|---|
| **P0** | Rehacer la rejilla temporal (ceros + semana MMWR) | `07`, `08`, `09` | 2 h |
| **P0** | Eliminar la fila «SARIMA+GARCH» o calcularla; corregir la columna `util` | `09` | 1 h |
| **P0** | Corregir la comparación de AIC y el encompassing (o retirarlos) | `09` | 1 h |
| **P0** | Re-especificar el GARCH (intercepto, estacionalidad, sin AR explosivo, serie semanal) | `09` | 3 h |
| **P0** | Regenerar figuras de `08` (escala) y `09` (guardar) | `08`, `09` | 1 h |
| **P0** | Verificar `con_fin`/`fec_def` y la exclusión de Bogotá con el diccionario y la fuente | `02`, `04` | 2 h |
| **P1** | `log1p` en lugar de Box-Cox; intervalos válidos o modelo de conteos | `08` | 2 h |
| **P1** | Estacionariedad completa y VAR en diferencias/log con dummies | `08`, `09` | 3 h |
| **P1** | Imputación: m = 20, `maxit` = 10, `pool()`, o retirar el capítulo | `05`, `06` | 2 h |
| **P1** | Validación fuera de muestra (2025) + Diebold-Mariano | `08`, `09` | 3 h |
| **P1** | Dummies COVID (2020–2021) y epidemia 2024 | `08`, `09` | 1 h |
| **P2** | INGARCH/NB; tasas por 100.000; VaR/ES con *backtesting* | nuevo `10` | 6 h |
| **P2** | Renombrar `06`, script maestro, semilla, *tests*, README | varios | 2 h |

**Regla de oro para la entrega:** ninguna cifra de la tabla resumen debe provenir de una asignación manual ni de una escala distinta a la del modelo comparado.

---

## 10. Checklist antes de entregar

- [ ] `resumen_econometria_avanzada.csv` no contiene `Inf`, ni valores copiados, ni `util` escrito a mano.
- [ ] Todos los AIC comparados están en la misma escala y con la misma muestra.
- [ ] La serie mensual tiene 72 meses por evento (con ceros) y `start/end` correctos.
- [ ] Las semanas epidemiológicas coinciden con el calendario del INS (validar 10 fechas al azar).
- [ ] El GARCH converge y α + β < 1. Si α y β no son significativos, se concluye que **no hay efecto ARCH** y no se presenta el GARCH como una mejora.
- [ ] Los intervalos de pronóstico son positivos y contienen la media.
- [ ] Cada figura del informe se regenera con un `ggsave()` del repositorio (ningún archivo de 0 bytes).
- [ ] El VAR es estacionario (raíces < 1) y sus residuos pasan el portmanteau.
- [ ] La imputación se hace con m ≥ 20 y se reporta el *pooling*, o se justifica su ausencia.
- [ ] `con_fin` y `fec_def` están verificados contra el diccionario del INS.
- [ ] Se declara explícitamente que la serie «nacional» excluye Bogotá (o se corrige la descarga).
- [ ] Hay un *split* temporal y un *benchmark* naive en la evaluación.

---

## Anexo A · Cómo reproducir las verificaciones

Con el proyecto abierto, desde la raíz y sin activar `renv` (evita el bloqueo de carga que documenta T-28):

```powershell
$env:R_LIBS_USER = Join-Path (Get-Location) "renv\library\windows\R-4.5\x86_64-w64-mingw32"
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" --vanilla "AUDITORIA_TECNICA.R"
```

Si `renv` se activa con normalidad (RStudio), basta `source("AUDITORIA_TECNICA.R")`.

Las verificaciones de la tabla §2.1 se obtienen de:

```r
parse("R/09_complemento_econometria.R")                 # sintaxis
d <- readRDS("data/processed/dengue_prep.rds"); dim(d)  # 748037 x 80
table(d$anio); table(d$con_fin); sum(is.na(d$fec_def))
a <- readRDS("data/series/econometria_avanzada.rds")
a$garch_dengue@fit$matcoef                              # α, β
Mod(polyroot(c(1, -a$garch_dengue@fit$matcoef[1:4, 1])))  # raíces AR
acf(residuals(a$garch_dengue), lag.max = 12, plot = FALSE)
```

---

## Anexo B · Nota sobre el alcance

El repositorio es **R** (no hay scripts de Python). La auditoría cubre los 9 scripts y los objetos generados. Si existe una libreta o *script* en Python con los mismos modelos, conviene auditarla por separado: los problemas de especificación (H-03, M-03) y de comparación (H-02) suelen repetirse, pero las trampas de alineación y de escalas cambian.

**Sobre el veredicto global.** El trabajo tiene una base de datos seria, un pipeline reproducible y una idea metodológica valiosa (los tres estados de faltantes). Lo que falta no es más código: es **una decisión explícita sobre qué se está modelando** —conteos semanales de notificación con estacionalidad y ceros— y la disciplina de comparar solo lo comparable. Con las correcciones P0, el capítulo econométrico pasa de «no defendible» a «defendible»; con las P1, a «publicable como trabajo de curso».
