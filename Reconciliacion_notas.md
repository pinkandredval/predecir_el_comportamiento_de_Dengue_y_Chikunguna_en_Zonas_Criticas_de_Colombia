# Reconciliación: auditoría ↔ Plan «Dengue, miedo y riesgo» (USTA 2026-II)

**Fecha:** 21 de octubre de 2026 · **Plan:** versión del 8 de octubre de 2026 (26 pp.) · **Auditoría:** `AUDITORIA_TECNICA.md`, commit `9452c25`.

Este documento no cambia el plan ni el diseño del estudio. Clasifica los 50 hallazgos de la auditoría frente a las 12 correcciones de la Tabla 2, señala lo que la auditoría no cubrió y corrige tres cosas de la propia auditoría.

---

## 0. Verificaciones hechas para reconciliar

| Verificación | Resultado | Fuente |
|---|---|---|
| Microdatos versionados | **Sí**: `data/processed/dengue_crudo.rds` (37,4 MB) y `chikungunya_crudo.rds` desde el commit `05079be`; `.RDataTmp` también | `git ls-files`, `git cat-file -s` |
| Código «01» | **1.600 casos, todos del exterior** (`cod_pais_o` ≠ 170: Venezuela 1.057, Brasil 190, Perú 77, México 43, R. Dominicana 31…) | `dengue_prep.rds` |
| `edad` según `uni_med` | 735.521 en años; **12.219 en meses**; **289 en días**; 8 en horas. 79 registros en días caen en «Adulto (18+)» | `dengue_prep.rds` |
| `ajuste` | 0: 400.794 · 3: 269.522 · 5: 31.549 · 7: 46.172. **No hay códigos 6 ni D** | `dengue_prep.rds` |
| 2025-W53 | 1.928 casos en la semana 53 de 2025, dentro de las «314 semanas» | `dengue_prep.rds` |

---

## 1. Correcciones a mi auditoría

| # | Lo que dije | Lo que hay | Efecto |
|---|---|---|---|
| 1 | «`.gitignore` correcto: no se versionan datos pesados» | Los dos `_crudo.rds` (37,4 MB) y `.RDataTmp` **están versionados**. Mi comprobación con `git ls-files \| Select-String` devolvió líneas vacías y la leí como «sin coincidencias» | El plan (§8) tiene razón. Es el P0 más urgente y mi auditoría lo omitió |
| 2 | «`nombrar_dpto` cubre los 33 códigos» | Cubre 32: el código «01» no está en el `case_when` y cae en «Cod. 01»; Bogotá (11) sí está mapeada pero no aparece en los datos | En `04` se grafica un pseudo-departamento con 1.600 casos. Corregir el rótulo o filtrar |
| 3 | «La serie NACIONAL excluye Bogotá» (hueco de cobertura) | Bogotá **no tiene dengue autóctono** (2.600 m, sin vector). El plan la usa con `I = 0` en el panel y estima con y sin ella como robustez | Retiro el hallazgo (M-10). Solo queda la precisión de rotular «sin casos autóctonos» |

---

## 2. Clasificación de los 50 hallazgos

La tricotomía pedida no cubre un caso frecuente: hallazgos que **no están en la Tabla 2 pero sí en el plan** (otra sección o anexo). Los separo como **(a′)** para no llamarlos «nuevos».

| Clasificación | N.º | Cuáles |
|---|---|---|
| **(a)** Coincide con la Tabla 2 | **14** | H-01…H-07 · M-04 · M-06 · M-07 · M-12 · T-15 · T-24 · T-25 |
| **(a′)** No está en la Tabla 2, pero sí en el plan | **8** | M-01, M-02, M-03 (script 14 y §4.5) · M-09 (§4.2, `ini_sin`) · T-17, T-29 (§11, semillas) · T-27 (MANIFIESTO.md) · T-30 (§8) |
| **(b)** Nuevo, no está en el plan | **27** | H-08 · M-05 · M-08 · M-11 · T-01…T-14, T-16, T-18…T-23, T-26, T-28 |
| **(c)** Contradice el plan | **1** | M-10 (y el plan tiene razón) |

N-01 a N-04 son hallazgos **nuevos de esta reconciliación**, no de la auditoría: por eso no entran en los 50.

Los 30 hallazgos T se agrupan en la tabla de una página por familia para que quepa; el detalle individual está en `AUDITORIA_TECNICA.md` §5.

---

## 3. Cobertura de las 12 correcciones de la Tabla 2

| # | Corrección | ¿La cubrió mi auditoría? | Detalle |
|---|---|---|---|
| 1 | Ceros en el calendario de 314 semanas; rehacer ADF, SARIMA, ETS, VAR y Granger | **Sí** (H-04) | Mismo diagnóstico y misma causa raíz |
| 2 | Retirar la mejora de AIC; GARCH sobre residuos o serie semanal | **Sí** (H-01, H-02, H-03) | Añado α₁ ≈ 0 y la raíz explosiva como evidencia extra; no cambia la corrección |
| 3 | Alinear el VAR y el encompassing | **Sí** (H-05) | Cuantifiqué el sesgo: 0,0002 → 0,0005 |
| 4 | MAPE en casos, ventana reservada, MASE, columna «útil» | **Sí** (H-06, H-07, M-12) | |
| 5 | ARCH-LM o McLeod-Li en vez de Breusch-Pagan | **Sí** (M-04) | La premisa de que la omití es incorrecta |
| 6 | Breusch-Godfrey o Ljung-Box con `fitdf` en vez de DW | **Sí** (M-04, T-15) | Añadir `P + Q` al `fitdf` si el modelo tiene parte estacional |
| 7 | INGARCH Poisson o binomial negativa | **Parcial** | Recomendé INGARCH (§7.7), pero **no vi la contradicción** «ADF estacionaria vs ARIMA(0,1,0)». La alternativa estaba, el diagnóstico no |
| 8 | AIC del VAR no comparable con univariados | **No** | Mi crítica de AIC se quedó en SARIMA vs GARCH |
| 9 | GARCH con 72 meses; Nyblom y sign bias; GJR/EGARCH | **Parcial** | Cuestioné las 72 observaciones y propuse la serie semanal; **no** calculé Nyblom ni sign bias |
| 10 | Tasas por 100.000; excluir el código «01» | **No** | Mi hallazgo sobre Bogotá es otro y estaba mal planteado |
| 11 | `m` y `maxit`; `pool()` | **Sí** (M-06, M-07) | |
| 12 | Inmunidad de rebaño como hipótesis | **No** | Es una afirmación del paper, no del código; fuera del alcance |

**Cobertura: 7 completas** (#1, #2, #3, #4, #5, #6, #11), **2 parciales** (#7, #9) y **3 no cubiertas** (#8, #10, #12).

**¿Omisión o solución indirecta?**

- **#7:** omisión en el diagnóstico; mi §7.7 da la alternativa (solución indirecta, sin el argumento). Sugerencia de implementación: hacer el INGARCH sobre la serie **semanal** (314 obs, 117 ceros) y no sobre 68 meses; con 419 casos, la mensual es demasiado corta para identificar `past_obs` y `past_mean`.
- **#8:** omisión sin solución indirecta. El 1.520,32 es la verosimilitud conjunta de dos ecuaciones; se puede comparar entre VAR de distinto `p`, nunca con un AIC univariado.
- **#10:** omisión sin solución indirecta. Confirmado con `cod_pais_o`: los 1.600 casos del «01» son importados. Excluirlos de las tasas es obligatorio (el denominador es población colombiana); en conteos solo si se rotulan como importados.
- **#12:** no aplica; el plan acierta al bajarla a hipótesis.

**Dónde el plan mejora mi propuesta de código:**

- **#1:** el `completar()` de B.1 (rejilla semanal, 314 semanas) es superior a mi `completar_mensual()` (solo rellena meses), porque conserva la granularidad semanal que necesitan el GARCH de #9 y la verificación de §11. Adoptar B.1.
- **#4 y H-08:** `auto.arima(lambda = "auto", biasadj = TRUE)` resuelve de raíz el problema de escalas de `plot(forecast)`; mi bloque manual de `InvBoxCox` queda obsoleto y no debe usarse.
- **#3:** mi alineación es equivalente a la del plan; B.1 la implementa recortando la serie real por `p` meses.

---

## 4. P0 del plan que mi auditoría no mencionó

### §8 · Microdatos en GitHub

- **Confirmado el diagnóstico del plan:** están versionados `dengue_crudo.rds` (37,4 MB), `chikungunya_crudo.rds` y `.RDataTmp` (11,8 MB), este último desde el mismo primer commit.
- Lo que **no** está versionado: `dengue_prep.rds`, `chikungunya_prep.rds`, los CSV de diagnóstico y las series. El `.gitignore` ya excluye `data/processed/`, pero no destrackea lo ya añadido.
- El `git rm --cached` del día 1 es correcto. Dos matices:
  1. `rm --cached` saca los blobs del árbol, **no del historial**. Si el grupo decide no reescribir, el dato sigue accesible por el commit `05079be`. GitHub además conserva vistas cacheadas de commits antiguos un tiempo; el purgado requiere soporte.
  2. `.RDataTmp` puede contener el mismo entorno con microdatos; conviene inspeccionarlo antes de decidir si basta con sacarlo del árbol.
- No hay credenciales en riesgo; el problema es de datos personales (fecha de nacimiento, municipio, etnia) en un repo público.

### §10 · Google Trends departamental sin probar

- No se puede resolver desde la auditoría de código: `gtrendsR` no está en el `renv` y este entorno no tiene red. Es una pregunta empírica del día 1, no un hallazgo.
- Prueba mínima (una llamada, tres chequeos): `gtrends("brote dengue", geo = "CO-VAC", time = "2020-01-01 2024-12-31")` debe devolver (i) `interest_over_time` no vacío, (ii) granularidad semanal, (iii) `interest_by_region` con departamentos para `geo = "CO"`.
- El plan B de la Tabla 8 (vía nacional como principal) es correcto y ya está decidido; no requiere cambios.
- Riesgo real: 429 de Google y el reescalado de las dos ventanas de 5 años. El empalme de B.8 (razón de medias en el tramo común) es el método correcto, pero el tramo común debe ser ≥ 12 meses y el ancla debe reportarse.

### #5 y #6

La premisa de que mi auditoría no las mencionó es incorrecta: están en **M-04** y **T-15**, con la misma corrección que propone el plan (ARCH-LM o McLeod-Li; Breusch-Godfrey o `checkresiduals`). Único matiz: en un SARIMA, `fitdf` debe incluir también `P + Q`.

### #7 · Chikungunya: INGARCH

- Mi auditoría recomendó INGARCH como alternativa, pero no formuló la contradicción que sí ve el plan (ADF estacionaria y, a la vez, ARIMA(0,1,0)).
- Implementación: `tscount::tsglm` sobre la serie **semanal** (314 obs), `distr = "nbinom"`, `link = "log"`, `past_obs = 1`, `past_mean = 1`; reportar ACF de residuos de Pearson y sobredispersión. Con 68 meses y 419 casos, el mensual no identifica el modelo.

### #8 · AIC del VAR

- No lo cubrí. El 1.520,32 es la verosimilitud conjunta del sistema (2 ecuaciones), no un AIC univariado: comparable entre VAR con distinto `p`, no contra SARIMA/ETS.
- La corrección del plan (evaluar por pronóstico fuera de muestra) es la adecuada. Si se quiere un criterio de información, que sea solo entre VAR.

### #10 · Tasas y código «01»

- No lo cubrí. Verificado: los 1.600 casos del «01» son del exterior (Venezuela 1.057, Brasil 190, Perú 77, México 43, R. Dominicana 31, EE. UU. 12 y otros). La hipótesis del plan («lo más probable es que sean del exterior») queda confirmada, no probable.
- Consecuencia doble: excluirlos de las tasas y corregir el rótulo «Cod. 01» que hoy aparece en el ranking de `04`.
- Las tasas necesitan las proyecciones del DANE por departamento (hoy «por buscar») y los pesos poblacionales de las 12 regiones EMA.

### #12 · Inmunidad de rebaño

- No lo cubrí porque es una afirmación del paper, no del código. Coincido con el plan: con cuatro serotipos, seis años y datos agregados no se identifica; presentarlo como hipótesis es lo correcto. Sin acción de código.

---

## 5. VaR epidemiológico: no sustituye, y probablemente no haga falta

El plan ya tiene VaR y CVaR **económicos** (pérdida = casos × costo por caso + pérdida asociada a la atención), con tres métodos de simulación y backtesting de Kupiec y Christoffersen (Tabla 4, B.5, script 19). Eso cubre el vínculo con teoría del riesgo.

Un «VaR epidemiológico» (cuantil de la incidencia condicional) **no debe reemplazar nada**. Además, como módulo nuevo, duplicaría dos cosas que el plan ya hace: el **probit de alerta** (P(zona de brote), script 21) y la maquinaria de backtesting del VaR económico.

Como mucho, cabe como **subproducto opcional** (P2): reportar el percentil 95 de la distribución condicional de incidencia como umbral de alerta adicional, con los mismos Kupiec/Christoffersen, dentro del script 19 o 21. Si no añade una cifra al paper, no vale la pena.

---

## 6. Observaciones a los esqueletos (no cambian el diseño)

- **B.7:** `uni_med %in % c("4","5")` no compila (`%in%`); `entrena` y `prueba` no están definidos (se define `en_entrena`); falta `TRUE ~ NA_real_` en el `case_when`. Además, la corrección de `uni_med` solo se aplica al scorecard: `04` y la serie por edad de `07` siguen leyendo meses y días como años (12.508 registros en meses o días, más 8 en horas).
- **B.1:** `expand_grid(dpto = unique(serie$dpto), ...)` no crea los departamentos que nunca tuvieron casos: en chikungunya solo 31 de los 33 códigos tienen algún caso. Para el panel regional hay que partir de los 33 códigos DIVIPOLA. Para la serie nacional no afecta. Mantener el `stopifnot(nrow(...) == 72)`.
- **B.1:** el calendario de 314 semanas incluye 2025-W53 (1.928 casos); el calendario MMWR asigna 52 semanas a 2025. Verificar con el INS: el efecto en la serie mensual es nulo (cae en diciembre), pero la verificación de §11 («sus 314 semanas») quedaría mal si el INS usa 313.
- **B.4:** `diff(log(casos_sem + 1))` pierde una semana; `M_sem_rezagada` debe recortarse para alinear longitudes antes de `ugarchspec`.
- **B.5:** verificar qué devuelve `fitted()` sobre un objeto `uGARCHsim` antes de usarlo como serie simulada; si devuelve la media condicional, el VaR debe calcularse sobre la serie simulada (`seriesSim`).

---

## 7. Orden sugerido para el día 1

1. `git rm --cached` de los dos `_crudo.rds` y `.RDataTmp`; `.gitignore`; push. Decidir si se reescribe el historial.
2. Corregir `07` (ceros, «01», `uni_med`) y regenerar `08`/`09` con las correcciones 1–4.
3. Prueba de Google Trends (una llamada) y decisión del plan B.
4. Corregir la sintaxis de B.7 antes de que alguien lo copie.
