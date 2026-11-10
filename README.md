# Enfermedad y atención: impacto económico del dengue y la chikungunya en 12 regiones EMA de Colombia (2022–2025)

**Proyecto de investigación reproducible · Econometría y Teoría del Riesgo · Universidad Santo Tomás (USTA)**

> Un proyecto de investigación reproducible que separa dos canales —la enfermedad y el miedo— y traduce el segundo a riesgo financiero y de crédito. Todo el análisis se reconstruye desde los datos crudos con **22 scripts de R** y un entorno fijado con **`renv`**.

---

## Ficha del proyecto

| Campo | Detalle |
|---|---|
| **Título** | *Enfermedad y atención: impacto económico del dengue y la chikungunya en 12 regiones EMA de Colombia (2022–2025)* |
| **Pregunta** | ¿Cuánto del efecto económico de una epidemia corresponde a la enfermedad y cuánto a la percepción de riesgo? |
| **Materias** | Econometría · Teoría del Riesgo |
| **Universidad** | Universidad Santo Tomás (USTA), Bogotá D. C. |
| **Periodo de datos** | 2020–2025 (panel región-mes 2022–2025) |
| **Unidad de análisis** | 12 regiones EMA × 48 meses = **576 observaciones** |
| **Stack** | R 4.5.2 · renv · Google Trends · NASA POWER · API Socrata |
| **Salidas** | `data/D_01` … `data/D_22` · `graficas/G_01` … `graficas/G_22` · `paper/` |

La cadena lógica del proyecto es una sola línea:

**casos de dengue → gente que busca información (miedo) → caída del turismo → deterioro de la cartera → pérdida en pesos.**

Este README explica cómo reconstruir esa cadena de principio a fin, sin atajos y sin datos inventados.

---

## Hallazgos principales

Tres cifras resumen el aporte del proyecto. **El 76,47 %** del efecto estimado sobre la ocupación hotelera en el pico de 2024-T2 se atribuye al canal del miedo, frente al 23,53 % de la incidencia: el miedo pesa más que el virus. Un probit agregado región-mes clasifica las alertas epidémicas con **AUC = 0,8267**, muy por encima del scorecard individual de hospitalización (0,5691), lo que demuestra que la señal predictiva está en el agregado y no en la ficha del paciente. Y en la cartera de crédito, la pérdida esperada agregada pasa de **1.435,2** a **2.220,3** miles de millones de COP en los trimestres epidémicos: **+54,71 %**.

> **Lectura rigurosa:** el +54,71 % es un cociente de promedios agregados, no un efecto causal. Comparado dentro de cada región, el aumento es de **+7,46 %**. El paper reporta ambas cifras y el repositorio permite verificarlas en `data/D_22/perdida_esperada.csv`.

**El miedo es medible, aunque su medición sea subjetiva.** No se observa directamente: se aproxima con un índice de atención construido a partir de Google Trends y depurado de la incidencia y de la estacionalidad. Esa es la variable `M` del proyecto. Siguiendo una recomendación del profesor, se declara de forma explícita que se trata de una **medición subjetiva (un proxy)**, no de una observación directa del miedo; sus resultados deben leerse con esa cautela.

---

## Estructura del repositorio

```text
proyecto_dengue_chikungunya/
│
├── data/                                  # Datos y resultados intermedios
│   ├── raw/                               # ⇦ DATOS ORIGINALES (no versionados)
│   │   ├── dengue/                        # SIVIGILA 210 · Datos_2020_210.xlsx … Datos_2025_210.xlsx
│   │   ├── chikungunya/                   # SIVIGILA 217 · Datos_2020_217.xlsx … Datos_2025_217.xlsx
│   │   ├── EMA/                           # DANE · anex-EMA-ago2026.xlsx
│   │   └── proyecciones_dane_2018-2042/   # DANE · PPED-AreaMun-2018-2042_VP.xlsx
│   │
│   └── D_01/ … D_22/                      # Salidas de cada script (una carpeta por paso)
│       ├── D_01/                          # Bases crudas consolidadas
│       ├── D_06/                          # Bases imputadas
│       ├── D_12/                          # Google Trends y miedo excedente
│       ├── D_13/                          # Panel EMA + clima (caché NASA POWER)
│       └── D_22/                          # Pérdida esperada de la cartera
│
├── R/                                     # ⇦ PIPELINE REPRODUCIBLE (22 scripts)
│   ├── 01_carga_datos.R
│   ├── 02_preparacion.R
│   ├── …
│   ├── 21_logit_probit.R
│   └── 22_riesgo_credito.R
│
├── graficas/                              # ⇦ RESULTADOS GRÁFICOS MODULARES
│   ├── G_01/ … G_22/                      # Una carpeta G_xx por script
│   ├── G_15/coefplot_panel.png
│   ├── G_19/var_backtest.png
│   ├── G_21/probit_alerta.png
│   └── G_22/el_estres.png
│
├── paper/                                 # ⇦ Informe final y poster
│   ├── Dengue__miedo_y_economía__el_precio_de_la_alerta_en_Colombia__2022__2025_ (1).pdf
│   └── Informe_Dengue_Miedo_Economia_USTA_v2.pdf
│
├── renv/                                  # Entorno reproducible (librería local)
├── renv.lock                              # Versiones exactas de más de 230 paquetes
├── .Rprofile                              # Activa renv al abrir el proyecto
├── proyecto_dengue_chikungunya.Rproj      # Proyecto de RStudio
└── README.md                              # Este archivo
```

**Regla de oro:** `data/raw/` contiene lo que se descarga una vez; `data/D_xx/` contiene lo que produce el pipeline; `graficas/G_xx/` contiene las figuras. Nada de eso se mezcla a mano.

---

## Requisitos

- **R 4.5.2** (la versión registrada en `renv.lock`).
- **RStudio** (recomendado) o cualquier IDE que respete el directorio del proyecto.
- **Git** para clonar el repositorio.
- **Conexión a internet** para los scripts **12**, **13** y **22**, que descargan datos por API.
- **Rtools** en Windows si algún paquete necesita compilación durante `renv::restore()`.
- Se recomiendan **8 GB de RAM o más**: la imputación múltiple del script 06 y las simulaciones del 19 son las etapas más exigentes.

---

## Clonar y restaurar el entorno

### 1. Clona el repositorio

```bash
git clone https://github.com/pinkandredval/predecir_el_comportamiento_de_Dengue_y_Chikunguna_en_Zonas_Criticas_de_Colombia.git proyecto_dengue_chikungunya
cd proyecto_dengue_chikungunya
```

> Usa `git clone`, no la descarga en ZIP: así conservas el historial, el `renv.lock` y la estructura exacta del proyecto.

### 2. Abre el proyecto y restaura `renv`

Abre `proyecto_dengue_chikungunya.Rproj` en RStudio. El archivo `.Rprofile` activa `renv` automáticamente; en la consola verás algo como:

```text
- Project 'proyecto_dengue_chikungunya' loaded. [renv 1.2.4]
```

Ahora restaura las versiones exactas de los paquetes:

```r
renv::restore()
```

- Responde **`Y`** cuando pregunte si deseas continuar.
- La primera restauración descarga más de 230 paquetes y puede tardar varios minutos.
- Si se interrumpe, vuelve a ejecutar `renv::restore()`: es idempotente y retoma lo que falta.

Verifica que el entorno quedó consistente:

```r
renv::status()
```

La respuesta esperada es:

```text
No issues found -- the project is in a consistent state.
```

> Si `renv::status()` lista paquetes descargados pero aún no usados, no es un error: el entorno está completo. Si `restore()` falla, revisa primero que estés en R 4.5.2 y que Rtools esté instalado.

---

## Fuentes de datos

El proyecto separa **datos originales** (`data/raw/`) de **datos derivados** (`data/D_xx/`). Los originales no se versionan porque pesan cientos de MB y porque deben descargarse de su fuente oficial. Hay dos formas de obtenerlos.

### 1. Descarga manual (la haces tú)

Debes descargar **cuatro fuentes** y dejarlas exactamente en las rutas indicadas. Sin estos archivos, los scripts 01, 10 y 13 no pueden correr.

| Fuente | Qué descargar | Enlace oficial | Guardar en |
|---|---|---|---|
| **SIVIGILA — Dengue (evento 210)** | Un archivo por año, 2020–2025 (6 XLSX) | [Portal SIVIGILA — Buscador de microdatos](https://portalsivigila.ins.gov.co/Paginas/Buscador.aspx) | `data/raw/dengue/` |
| **SIVIGILA — Chikungunya (evento 217)** | Un archivo por año, 2020–2025 (6 XLSX) | [Portal SIVIGILA — Buscador de microdatos](https://portalsivigila.ins.gov.co/Paginas/Buscador.aspx) | `data/raw/chikungunya/` |
| **DANE — Proyecciones de población** | En la página del DANE, sección **«Proyecciones y retroproyecciones de población municipal 2018–2042»**: descargar «Serie municipal de población por área» | [Proyecciones de población (DANE)](https://www.dane.gov.co/index.php/estadisticas-por-tema/demografia-y-poblacion/proyecciones-de-poblacion) | `data/raw/proyecciones_dane_2018-2042/` |
| **DANE — Encuesta Mensual de Alojamiento (EMA)** | En la página del DANE, en **«Información [mes]»**: abrir **«Anexos»** y descargar el XLSX del mes de referencia | [Encuesta Mensual de Alojamiento (DANE)](https://www.dane.gov.co/index.php/estadisticas-por-tema/servicios/encuesta-mensual-de-alojamiento-ema) | `data/raw/EMA/` |

**Nombres de archivo y rutas:**

- **Dengue:** `Datos_2020_210.xlsx`, `Datos_2021_210.xlsx`, …, `Datos_2025_210.xlsx` en `data/raw/dengue/`.
- **Chikungunya:** `Datos_2020_217.xlsx`, `Datos_2021_217.xlsx`, …, `Datos_2025_217.xlsx` en `data/raw/chikungunya/`.
- **Proyecciones DANE:** `PPED-AreaMun-2018-2042_VP.xlsx` en `data/raw/proyecciones_dane_2018-2042/`. El nombre está escrito en el script 10, no lo cambies.
- **EMA:** `anex-EMA-ago2026.xlsx` en `data/raw/EMA/`. El nombre está escrito en el script 13. Si el DANE publica un mes más reciente, renómbralo a este nombre o actualiza `ruta_ema` en `R/13_panel_ema.R`.

> En SIVIGILA el nombre no es obligatorio (el script 01 lee todos los `.xlsx` de la carpeta), pero debe haber **un archivo por año**. En DANE y EMA sí es obligatorio: las rutas están escritas en el código.

> **Verificación rápida antes de ejecutar.** Desde la raíz del proyecto:
>
> ```r
> list.files("data/raw", pattern = "\\.xlsx$", recursive = TRUE)
> ```
>
> Deben aparecer **14 archivos**: 6 de dengue, 6 de chikungunya, 1 de proyecciones y 1 de EMA. El script 01 lee **todos** los `.xlsx` de las carpetas `dengue/` y `chikungunya/`, así que no dejes duplicados ni archivos temporales de Excel (`~$…`).

### 2. Descarga automática (la hace el código R)

**No necesitas descargar a mano Google Trends, los datos climáticos de NASA POWER ni la cartera de la Superfinanciera: se descargan automáticamente al ejecutar los scripts 12, 13 y 22.** El código los obtiene por API y los deja en las carpetas del proyecto.

| Fuente | Script | Qué descarga | Dónde queda |
|---|---|---|---|
| **Google Trends** | `12_indice_atencion.R` | Índice de atención para 33 códigos DIVIPOLA (32 departamentos + Bogotá D. C.), vía `gtrendsR` | `data/D_12/trends_*.csv` y `data/D_12/miedo_excedente*.csv` |
| **NASA POWER** | `13_panel_ema.R` | Precipitación (`PRECTOTCORR`) y temperatura (`T2M`) mensuales de las 12 capitales, vía API o `nasapower` | `data/D_13/clima/clima_*.csv` (caché) |
| **Superfinanciera — Socrata** | `22_riesgo_credito.R` | Saldos de cartera por municipio, dataset `u2wk-tfe3` | `data/D_22/perdida_esperada.csv` y `data/D_22/var_granger_credito.csv` |

**Notas importantes de estas descargas:**

- **Google Trends (script 12):** no requiere API key, pero hace una consulta por departamento y por enfermedad (dengue y chikungunya), con una pausa de 5 segundos entre cada una; en total son 66 consultas y puede tardar 10 minutos o más. Si Google responde con error 429, espera y vuelve a ejecutar el script. Los índices de Trends se reescalan con el tiempo, así que pequeñas diferencias entre descargas son esperables: los CSV de `data/D_12/` son la foto de la muestra que se usó.
- **NASA POWER (script 13):** usa `nasapower` si está instalado; si no, consulta la misma API mensual con `httr2`. No requiere key. La caché de `data/D_13/clima/` evita repetir la descarga; bórrala o pon `forzar_descarga_clima <- TRUE` en el script si necesitas refrescarla.
- **Socrata (script 22):** descarga aproximadamente 94.000 filas en páginas de 50.000, con timeout de 180 segundos y tres reintentos. Si la API falla, el script lanza un **warning** y genera datos simulados de respaldo para no romper el pipeline. Para resultados reales, verifica que la consola **no** muestre ese warning y vuelve a ejecutar.

---

## Orden de ejecución

**Los 22 scripts deben ejecutarse del 01 al 22, en orden secuencial y sin saltos.** Cada script consume las salidas del anterior: el 13 necesita el panel del 10, el 15 necesita el panel del 13 y el 22 necesita la cartera y el panel ya construidos. Alterar el orden produce errores de archivos faltantes o, peor, resultados silenciosamente inconsistentes.

### Opción A — Ejecución completa en un solo bucle

Desde la raíz del proyecto, con `renv` restaurado:

```r
scripts <- sort(list.files("R", pattern = "\\.R$", full.names = TRUE))

for (s in scripts) {
  message("\n>>> Ejecutando: ", basename(s))
  source(s, echo = FALSE)
}
```

### Opción B — Script por script

Abre cada archivo en RStudio y pulsa **Source**. El orden es estrictamente:

```text
01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11
   → 12 → 13 → 14 → 15 → 16 → 17 → 18 → 19 → 20 → 21 → 22
```

**Recomendaciones:**

- Ejecuta siempre desde la raíz del proyecto (o con el `.Rproj` abierto): todos los scripts usan `here()` para resolver rutas.
- La primera corrida completa puede tardar varios minutos; la imputación (06) y las simulaciones (19) son las etapas más lentas.
- Si un script falla, corrige la causa y retoma desde ese punto: las salidas anteriores quedan guardadas en `data/D_xx/`.
- No borres `data/D_10/serie_dengue_tasas.csv` ni `data/D_10/serie_chik_tasas.csv`: son la entrada de los scripts 11, 12 y 13.

---

## Scripts del pipeline

| Número | Nombre del archivo | Descripción |
|---|---|---|
| 01 | `01_carga_datos.R` | Lee y consolida los seis archivos anuales de dengue (210) y chikungunya (217) del SIVIGILA; guarda las bases crudas en `D_01`. |
| 02 | `02_preparacion.R` | Normaliza nombres, fechas y tipos; marca faltantes estructurales («no aplica»); genera las bases preparadas en `D_02`. |
| 03 | `03_auditoria.R` | Audita la calidad de las bases preparadas: dimensiones, duplicados, consistencia y valores extremos; reporte en consola. |
| 04 | `04_descrip_explo.R` | Análisis descriptivo y exploratorio; produce las figuras de edad, sexo, hospitalización y series en `G_04`. |
| 05 | `05_faltantes.R` | Diagnóstico de faltantes (MCAR/MAR), coocurrencia y mapas; deja el diagnóstico pre-imputación en `D_05` y `G_05`. |
| 06 | `06_imputación.R` | Imputación múltiple con `mice`: bosque aleatorio para dengue y `pmm` para chikungunya, respetando los NA estructurales; bases completas en `D_06`. |
| 07 | `07_serietemporal_conimputacion.R` | Construye series semanales MMWR (33 códigos DIVIPOLA × 313 semanas) y mensuales, nacionales, departamentales y municipales, en `D_07`. |
| 08 | `08_econometria.R` | ADF, Box-Cox, diferenciación, ACF/PACF, SARIMA y validación fuera de muestra (MASE a 12 meses), con figuras en `G_08`. |
| 09 | `09_complemento_econometria.R` | Complemento avanzado: VAR, Granger, IRF, FEVD, GARCH, ETS y encompassing SARIMA vs. VAR; salidas en `D_09` y `G_09`. |
| 10 | `10_tasas_incidencia.R` | Calcula tasas por 100.000 habitantes con las proyecciones DANE para las 12 regiones EMA y genera rankings en `D_10` y `G_10`. |
| 11 | `11_canal_endemico_tecnico.R` | Canal endémico por el método de Bortman y análisis técnico (bandas, tendencia, ruptura) para dengue y chikungunya, con figuras en `G_11`. |
| 12 | `12_indice_atencion.R` | Descarga Google Trends, construye el índice de atención y calcula el miedo excedente depurado de incidencia y estacionalidad, en `D_12`. |
| 13 | `13_panel_ema.R` | Arma el panel región-mes (12 × 48): ocupación EMA, incidencia, miedo y clima de NASA POWER; salidas en `D_13`. |
| 14 | `14_estacionariedad.R` | Pruebas de raíz unitaria y estacionariedad del panel: Box-Cox, ADF, KPSS y ACF/PACF; 41 pruebas en `D_14` y `G_14`. |
| 15 | `15_panel_descomposicion.R` | Panel de efectos fijos con errores de Driscoll-Kraay y agrupados; descomposición π del pico 2024-T2 (76,47 % miedo) en `D_15` y `G_15`. |
| 16 | `16_var_its_granger.R` | VAR nacional, FEVD, causalidad de Granger y serie interrumpida de chikungunya 2025; resultados en `D_16` y `G_16`. |
| 17 | `17_garch_volatilidad.R` | Modela la volatilidad condicional del rendimiento epidémico semanal con GARCH y GJR-GARCH, incluyendo el miedo como regresor; salidas en `D_17` y `G_17`. |
| 18 | `18_capm_eventos.R` | CAPM regional con validación leave-one-out y estudio de eventos alrededor de los picos de miedo; resultados en `D_18` y `G_18`. |
| 19 | `19_var_cvar_simulation.R` | VaR y CVaR al 95 % y 99 % por método histórico, Monte Carlo y bootstrap por bloques, con backtesting de Kupiec y Christoffersen; salidas en `D_19` y `G_19`. |
| 20 | `20_markowitz.R` | Frontera eficiente de Markowitz con las 12 regiones como activos, comparación normal vs. epidémico y mapa de isomedia/isovarianza; salidas en `D_20` y `G_20`. |
| 21 | `21_logit_probit.R` | Scorecard individual de hospitalización (logit/probit) y probit agregado de alerta temprana; curvas ROC, efectos marginales y métricas en `D_21` y `G_21`. |
| 22 | `22_riesgo_credito.R` | Descarga la cartera de la Superfinanciera vía Socrata y estima la pérdida esperada `EL = PD × LGD × EAD`, con contraste epidémico y Granger dengue → morosidad; salidas en `D_22` y `G_22`. |

---

## Resultados y trazabilidad

- **Datos intermedios:** `data/D_01` … `data/D_22`. La carpeta `D_xx` corresponde siempre al script `xx` que la produce. Así se sabe de dónde salió cada tabla sin abrir el código.
- **Figuras:** `graficas/G_01` … `graficas/G_22`. Cada script que produce gráficos los guarda en su propia carpeta modular **`G_xx`**; no hay una carpeta única donde se mezclen resultados. En total se generan más de **100 figuras**.
- **Cifras estrella:** `data/D_15/participacion_atencion.csv` (76,47 %), `data/D_21/probit_alerta_metricas.csv` (AUC 0,8267) y `data/D_22/perdida_esperada.csv` (+54,71 %).
- **Trazabilidad:** cada cifra del paper puede rastrearse hasta un `.csv` de `data/D_xx` y hasta la línea del script que la produce. Esa es la prueba de que el resultado no es un resumen manual.
- **Informe final:** `paper/Informe_Dengue_Miedo_Economia_USTA_v2.pdf`.

---

## Notas de reproducibilidad

- **Entorno fijado:** `renv.lock` registra R 4.5.2 y las versiones exactas de más de 230 paquetes. `renv::restore()` reconstruye ese entorno en cualquier máquina.
- **Datos crudos fuera del control de versiones:** `data/raw/**` está ignorado por Git, igual que los `.rds` y los CSV pesados generados por los scripts. El repositorio guarda la estructura (`.gitkeep`) y el código; los datos se reconstruyen.
- **Rutas portables:** todos los scripts usan `here()`, nunca rutas absolutas. Por eso deben ejecutarse desde la raíz del proyecto.
- **Semillas fijas:** las rutinas aleatorias (imputación, partición train/test, simulaciones, bootstrap) fijan `set.seed()`, de modo que dos corridas sobre los mismos insumos producen los mismos resultados.
- **Dependencias externas:** Google Trends puede cambiar sus índices o limitar consultas; NASA POWER y Socrata pueden estar temporalmente fuera de servicio. El código lo maneja con caché, reintentos y advertencias explícitas; revisa siempre la consola.
- **Verificación de consistencia:** después de restaurar, `renv::status()` debe decir `No issues found`. Después de correr el pipeline, las carpetas `D_xx` y `G_xx` deben tener los archivos que el script reporta en consola.

---

## Equipo

**Integrantes:**

- **Kevin Leonardo Chaparro Reyes**
- **Valentina Muñoz Palma**
- **Paula Margarita Triana Ancinez**

**Materias:** Econometría y Teoría del Riesgo

**Universidad:** Universidad Santo Tomás (USTA) — Bogotá D. C.

**Facultades:** Facultad de Estadística · Facultad de Ciencias Económicas y Administrativas

**Docentes:** Jhonier Rangel · William Rincón

**Periodo:** 2026-II

---

## Cierre

Replicar este proyecto no es solo ejecutar 22 scripts: es reconstruir una cadena completa de evidencia, desde un registro del SIVIGILA hasta una pérdida esperada en pesos. Cada carpeta `D_xx` y cada figura `G_xx` existen para que cualquier persona pueda auditar el camino. Si el pipeline corre limpio y `renv::status()` no encuentra problemas, el resultado está listo para defenderse.

**El miedo es una variable económica medible; la econometría la mide y la teoría del riesgo le pone precio.**
