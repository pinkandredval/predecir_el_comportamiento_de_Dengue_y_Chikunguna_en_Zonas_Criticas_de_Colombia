# ============================================================
# 14_estacionariedad.R
# Pruebas de raíz unitaria y estacionariedad del panel EMA
# Box-Cox, ADF, KPSS y ACF/PACF (nacional y 12 regiones)
# ============================================================
#
# OBJETIVO
#   Antes de estimar el panel con efectos fijos (script 15), verificar
#   el orden de integración de las series del panel:
#
#     Y  ocupación hotelera (EMA, porcentaje)
#     I  incidencia de dengue y chikungunya (casos)
#     M  miedo excedente (residuo de la regresión auxiliar, script 12)
#
#   El script hace tres cosas:
#     1. Box-Cox: lambda óptimo (forecast::BoxCox.lambda) de las series
#        estrictamente positivas. Si la serie tiene ceros o negativos,
#        la transformación no está definida y se documenta el motivo.
#     2. ADF (H0: raíz unitaria) y KPSS (H0: estacionariedad en nivel)
#        sobre la serie nacional y sobre las 12 regiones EMA en Y,
#        casos_dengue y M_dengue. Los p-valores y estadísticos se
#        guardan en una tabla ordenada.
#     3. ACF y PACF de las series nacionales (en nivel y, cuando el ADF
#        no rechaza, en primera diferencia) y de las regionales.
#
# INSUMOS
#   data/D_13/panel_ema_den_chik.rds   (576 filas: 12 regiones x 48 meses)
#
# SALIDAS
#   data/D_14/pruebas_estacionariedad.csv
#   data/D_14/pruebas_estacionariedad.rds
#   graficas/G_14/acf_pacf_nacional.png
#   graficas/G_14/acf_pacf_regiones.png
#
# CRITERIOS (nivel de significación alfa = 0,05)
#   - ADF:  p < 0,05  => se rechaza la raíz unitaria (estacionaria).
#   - KPSS: p < 0,05  => se rechaza la estacionariedad en nivel.
#   - requiere_diferenciacion se define con el ADF (p >= 0,05), tal
#     como pide el plan; el KPSS entra como contraste, no como criterio.
#   - Las pruebas se aplican a la serie en nivel, sin transformar y sin
#     ajuste estacional: el Box-Cox es descriptivo (lambda), no entra
#     en las pruebas.
#
# NOTA SOBRE LA AGREGACIÓN NACIONAL
#   Y y M son porcentajes/residuos: se promedian por mes. Los casos de
#   dengue y chikungunya se suman. Así la serie nacional conserva la
#   escala de las variables del panel.
#
# NOTA SOBRE CHIKUNGUNYA
#   A nivel regional sus series son casi todas ceros (446 de 576
#   observaciones) y en dos regiones son constantes. Por eso las pruebas
#   regionales se limitan a Y, casos_dengue y M_dengue; chikungunya se
#   prueba a nivel nacional y aparece documentada en la tabla.
# ============================================================

suppressPackageStartupMessages({
  library(tseries)
  library(forecast)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(here)
})

# --- 0. Constantes y rutas -----------------------------------------

dir.create(here("data", "D_14"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("graficas", "G_14"), recursive = TRUE, showWarnings = FALSE)

ruta_panel   <- here("data", "D_13", "panel_ema_den_chik.rds")
ruta_csv     <- here("data", "D_14", "pruebas_estacionariedad.csv")
ruta_rds     <- here("data", "D_14", "pruebas_estacionariedad.rds")
ruta_acf_nac <- here("graficas", "G_14", "acf_pacf_nacional.png")
ruta_acf_reg <- here("graficas", "G_14", "acf_pacf_regiones.png")

alfa <- 0.05

# Rezagos que se grafican en cada ACF/PACF. La serie nacional tiene 48
# observaciones (4 años) y las regionales también; 18 y 12 rezagos
# cubren un ciclo anual y medio sin saturar los paneles.
lag_max_nacional <- 18L
lag_max_regiones <- 12L

# Orden canónico de las 12 regiones EMA (el mismo del script 13). Se usa
# para ordenar la tabla y las filas del gráfico regional.
regiones_ema <- c(
  "Bogotá D.C.",
  "Cartagena",
  "Región Central",
  "Región Costa Caribe",
  "Región Eje Cafetero",
  "Antioquia",
  "San Andrés y Providencia",
  "Región Pacífico",
  "Región Santanderes",
  "Región Llanos Orinoquía",
  "Golfo Morrosquillo y Sabana",
  "Región Amazonía"
)

# Etiquetas legibles para la tabla y los gráficos.
etiquetas_variable <- c(
  Y            = "Ocupación hotelera (%)",
  casos_dengue = "Casos de dengue (suma mensual)",
  casos_chik   = "Casos de chikungunya (suma mensual)",
  M_dengue     = "Miedo excedente de dengue (residuo)",
  M_chik       = "Miedo excedente de chikungunya (residuo)"
)

# Orden de las variables en la tabla (Y primero, luego casos, luego M).
niveles_variable <- c("Y", "casos_dengue", "casos_chik", "M_dengue", "M_chik")

# --- 1. Carga y validación del panel -------------------------------

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel del script 13: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

panel <- readRDS(ruta_panel)

requeridas <- c("region_ema", "mes", "Y", "casos_dengue", "M_dengue")
faltantes <- setdiff(requeridas, names(panel))
if (length(faltantes) > 0) {
  stop(
    "Faltan columnas en el panel: ", paste(faltantes, collapse = ", "), ".",
    call. = FALSE
  )
}

# Chikungunya puede no estar si el script 13 no encontró la fuente: en
# ese caso la columna existe pero con NA. Se garantiza que exista para
# que la agregación nacional no falle.
for (col in c("casos_chik", "M_chik")) {
  if (!col %in% names(panel)) panel[[col]] <- NA_real_
}

panel <- panel |>
  mutate(region_ema = as.character(region_ema), mes = as.Date(mes)) |>
  arrange(region_ema, mes)

if (nrow(panel) != 576) {
  stop(
    "El panel debe tener 576 filas (12 regiones x 48 meses) y tiene ",
    nrow(panel), ".", call. = FALSE
  )
}

if (dplyr::n_distinct(panel$region_ema) != 12) {
  stop("El panel debe tener 12 regiones EMA.", call. = FALSE)
}

if (dplyr::n_distinct(panel$mes) != 48) {
  stop("El panel debe tener 48 meses (2022-01 a 2025-12).", call. = FALSE)
}

periodo_panel <- sprintf(
  "%s a %s",
  format(min(panel$mes), "%Y-%m"),
  format(max(panel$mes), "%Y-%m")
)

sobrantes <- setdiff(unique(panel$region_ema), regiones_ema)
if (length(sobrantes) > 0) {
  warning(
    "Regiones fuera de la lista EMA: ", paste(sobrantes, collapse = ", "),
    ". Se excluyen de las pruebas regionales.",
    call. = FALSE
  )
}

# tsibble: se usa para validar que la rejilla región-mes esté completa y
# regular. El resto del script trabaja sobre el data.frame, de modo que
# si tsibble no está instalado el script sigue funcionando.
#
# El índice se declara como yearmonth: con un Date, tsibble infiere un
# intervalo diario (los meses tienen 28-31 días) y has_gaps() marcaría
# huecos que no existen.
if (requireNamespace("tsibble", quietly = TRUE)) {
  panel_ts <- panel |>
    mutate(mes_ym = tsibble::yearmonth(mes)) |>
    tsibble::as_tsibble(key = region_ema, index = mes_ym, regular = TRUE)

  huecos <- tsibble::has_gaps(panel_ts)
  if (any(huecos$.gaps)) {
    warning(
      "El panel tiene meses faltantes en al menos una región. ",
      "Las pruebas se calculan sobre las observaciones disponibles.",
      call. = FALSE
    )
  }
} else {
  panel_ts <- panel |> group_by(region_ema)
  message(
    "tsibble no está instalado: se trabaja con un data.frame agrupado ",
    "por region_ema."
  )
}

n_na <- colSums(is.na(panel[, c("Y", "casos_dengue", "M_dengue")]))
if (any(n_na > 0)) {
  warning(
    "Hay valores faltantes en las variables del panel: ",
    paste(names(n_na)[n_na > 0], n_na[n_na > 0], sep = " = ",
          collapse = ", "),
    ". Las pruebas usan las observaciones completas.",
    call. = FALSE
  )
}

message(
  "Panel cargado: ", nrow(panel), " filas | ",
  dplyr::n_distinct(panel$region_ema), " regiones | ",
  dplyr::n_distinct(panel$mes), " meses (", periodo_panel, ")."
)

# --- 2. Agregación nacional ----------------------------------------
#
# Y y M se promedian porque son un porcentaje y un residuo; los casos
# se suman. El resultado son cinco series mensuales de 48 observaciones.

nacional <- panel |>
  group_by(mes) |>
  summarise(
    Y            = mean(Y, na.rm = TRUE),
    casos_dengue = sum(casos_dengue, na.rm = TRUE),
    casos_chik   = sum(casos_chik, na.rm = TRUE),
    M_dengue     = mean(M_dengue, na.rm = TRUE),
    M_chik       = mean(M_chik, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(mes)

if (nrow(nacional) != 48) {
  stop(
    "La serie nacional debe tener 48 meses y tiene ", nrow(nacional), ".",
    call. = FALSE
  )
}

# mean() de una columna toda NA devuelve NaN. NaN es un valor ausente
# para las pruebas, pero se normaliza a NA para que la tabla y el CSV no
# arrastren "NaN".
nacional <- nacional |>
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

message(
  "Serie nacional: ", nrow(nacional), " meses. ",
  "Dengue: ", format(sum(nacional$casos_dengue, na.rm = TRUE),
                     big.mark = ".", decimal.mark = ",", scientific = FALSE),
  " casos | Chikungunya: ",
  format(sum(nacional$casos_chik, na.rm = TRUE),
         big.mark = ".", decimal.mark = ",", scientific = FALSE),
  " casos."
)

# --- 3. Utilidades --------------------------------------------------

# Etiqueta corta de la serie para los paneles de los gráficos.
etiqueta_serie <- function(v) {
  unname(c(
    Y            = "Y: ocupación hotelera (%)",
    casos_dengue = "Casos de dengue",
    casos_chik   = "Casos de chikungunya",
    M_dengue     = "M: miedo excedente (dengue)",
    M_chik       = "M: miedo excedente (chik.)"
  )[v])
}

# Prueba una serie y devuelve una fila con Box-Cox, ADF, KPSS y el
# diagnóstico conjunto.
#
#   x            vector numérico de la serie (nivel)
#   nivel        "Nacional" o "Regional"
#   region_ema   nombre de la región (NA en las filas nacionales)
#   variable     Y, casos_dengue, casos_chik, M_dengue o M_chik
#   etiqueta     texto largo de la variable
#   periodo      periodo cubierto ("2022-01 a 2025-12")
probar_serie <- function(x, nivel, region_ema, variable, etiqueta, periodo) {

  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)

  n_ceros     <- if (n > 0) sum(x == 0) else NA_integer_
  n_negativos <- if (n > 0) sum(x < 0) else NA_integer_
  minimo      <- if (n > 0) min(x) else NA_real_
  maximo      <- if (n > 0) max(x) else NA_real_

  # Una serie constante no tiene varianza: ADF y KPSS no están definidos
  # (devuelven NaN) y diferenciarla no aporta nada.
  constante <- n > 0 && length(unique(x)) == 1L
  es_miedo  <- str_detect(variable, "^M_")

  # --- Box-Cox -----------------------------------------------------
  # BoxCox.lambda() exige datos estrictamente positivos. Se documenta
  # el motivo cuando no aplica en lugar de devolver un lambda vacío.
  boxcox_aplica <- FALSE
  boxcox_lambda <- NA_real_
  boxcox_nota   <- NA_character_

  if (n < 8) {
    boxcox_nota <- sprintf("No aplica: solo %d observaciones", n)
  } else if (constante) {
    boxcox_nota <- "No aplica: serie constante"
  } else if (n_negativos > 0) {
    boxcox_nota <- sprintf(
      "No aplica: %d valores negativos%s",
      n_negativos,
      if (es_miedo) " (M es un residuo con signo)" else ""
    )
  } else if (n_ceros > 0) {
    boxcox_nota <- sprintf("No aplica: %d ceros en la serie", n_ceros)
  } else {
    boxcox_lambda <- tryCatch(
      suppressWarnings(forecast::BoxCox.lambda(x, method = "guerrero")),
      error = function(e) NA_real_
    )
    if (is.finite(boxcox_lambda)) {
      boxcox_aplica <- TRUE
      boxcox_nota <- "lambda óptimo de Box-Cox (Guerrero) sobre la serie en nivel"
    } else {
      boxcox_nota <- "No se pudo estimar lambda (BoxCox.lambda falló)"
    }
  }

  # --- ADF (H0: raíz unitaria) -------------------------------------
  adf_estadistico <- NA_real_
  adf_p           <- NA_real_
  adf_rezagos     <- NA_integer_
  adf_conclusion  <- NA_character_

  if (constante) {
    adf_conclusion <- "No evaluable: serie constante"
  } else if (n < 10) {
    adf_conclusion <- "No evaluable: menos de 10 observaciones"
  } else {
    adf <- tryCatch(
      suppressWarnings(tseries::adf.test(x)),
      error = function(e) e
    )
    if (!inherits(adf, "error") && is.finite(adf$p.value)) {
      adf_estadistico <- unname(adf$statistic)
      adf_p           <- adf$p.value
      adf_rezagos     <- as.integer(unname(adf$parameter))
      adf_conclusion  <- if (adf_p < alfa) {
        "Rechaza raíz unitaria (estacionaria)"
      } else {
        "No rechaza raíz unitaria (no estacionaria)"
      }
    } else if (inherits(adf, "error")) {
      adf_conclusion <- paste0("No evaluable: ", conditionMessage(adf))
    } else {
      adf_conclusion <- "No evaluable: ADF sin p-valor"
    }
  }

  # --- KPSS (H0: estacionariedad en nivel) -------------------------
  kpss_estadistico <- NA_real_
  kpss_p           <- NA_real_
  kpss_conclusion  <- NA_character_

  if (constante) {
    kpss_conclusion <- "No evaluable: serie constante"
  } else if (n < 10) {
    kpss_conclusion <- "No evaluable: menos de 10 observaciones"
  } else {
    kpss <- tryCatch(
      suppressWarnings(tseries::kpss.test(x, null = "Level")),
      error = function(e) e
    )
    if (!inherits(kpss, "error") && is.finite(kpss$p.value)) {
      kpss_estadistico <- unname(kpss$statistic)
      kpss_p           <- kpss$p.value
      kpss_conclusion  <- if (kpss_p < alfa) {
        "Rechaza estacionariedad en nivel (no estacionaria)"
      } else {
        "No rechaza estacionariedad en nivel"
      }
    } else if (inherits(kpss, "error")) {
      kpss_conclusion <- paste0("No evaluable: ", conditionMessage(kpss))
    } else {
      kpss_conclusion <- "No evaluable: KPSS sin p-valor"
    }
  }

  # --- Diagnóstico conjunto ----------------------------------------
  # ADF y KPSS tienen hipótesis nulas opuestas: cuando coinciden el
  # diagnóstico es claro; cuando discrepan se deja constancia.
  adf_rechaza  <- is.finite(adf_p)  && adf_p  < alfa
  kpss_rechaza <- is.finite(kpss_p) && kpss_p < alfa

  diagnostico <- dplyr::case_when(
    constante ~ "No evaluable: serie constante",
    !is.finite(adf_p) || !is.finite(kpss_p) ~ "No evaluable: prueba no disponible",
    adf_rechaza & !kpss_rechaza ~ "Estacionaria (ADF y KPSS coinciden)",
    !adf_rechaza & kpss_rechaza ~
      "No estacionaria: raíz unitaria (ADF y KPSS coinciden)",
    adf_rechaza & kpss_rechaza ~
      "Mixta: estacionaria con tendencia o cambio estructural",
    TRUE ~ "Mixta: evidencia insuficiente"
  )

  # El criterio del plan es el ADF. NA = no evaluable (serie constante o
  # prueba no disponible).
  requiere_diferenciacion <- if (is.finite(adf_p)) adf_p >= alfa else NA

  # --- Notas de la fila --------------------------------------------
  notas <- character(0)
  if (constante) {
    notas <- c(notas, "Serie constante: ADF y KPSS no están definidos")
  }
  if (is.finite(adf_p) && adf_p <= 0.01) {
    notas <- c(notas, "ADF: p en el mínimo tabulado (0,01)")
  }
  if (is.finite(kpss_p) && kpss_p >= 0.10) {
    notas <- c(notas, "KPSS: p en el máximo tabulado (0,10)")
  }
  if (adf_rechaza && kpss_rechaza) {
    notas <- c(notas, "ADF y KPSS discrepan: revisar tendencia o quiebre")
  }

  tibble(
    nivel                = nivel,
    region_ema           = region_ema,
    variable             = variable,
    variable_etiqueta    = etiqueta,
    periodo              = periodo,
    n_obs                = n,
    minimo               = minimo,
    maximo               = maximo,
    n_ceros              = n_ceros,
    n_negativos          = n_negativos,
    boxcox_aplica        = boxcox_aplica,
    boxcox_lambda        = boxcox_lambda,
    boxcox_nota          = boxcox_nota,
    adf_estadistico      = adf_estadistico,
    adf_p                = adf_p,
    adf_rezagos          = adf_rezagos,
    adf_conclusion       = adf_conclusion,
    kpss_estadistico     = kpss_estadistico,
    kpss_p               = kpss_p,
    kpss_conclusion      = kpss_conclusion,
    diagnostico          = diagnostico,
    requiere_diferenciacion = requiere_diferenciacion,
    nota = if (length(notas) > 0) {
      paste(notas, collapse = "; ")
    } else {
      NA_character_
    }
  )
}

# Calcula ACF y PACF hasta lag_max y devuelve una tabla larga.
# Devuelve NULL si la serie es constante o demasiado corta, porque acf()
# no produce valores utilizables.
calcular_acf_pacf <- function(x, lag_max) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 8 || length(unique(x)) == 1L) return(NULL)

  lag_max <- min(lag_max, length(x) - 2L)
  a <- stats::acf(x, lag.max = lag_max, plot = FALSE)
  p <- stats::pacf(x, lag.max = lag_max, plot = FALSE)

  # Banda de ruido blanco al 95%: ±1,96/sqrt(n).
  banda <- stats::qnorm(1 - alfa / 2) / sqrt(a$n.used)

  tibble(
    lag   = c(as.integer(a$lag), as.integer(p$lag)),
    valor = c(as.numeric(a$acf), as.numeric(p$acf)),
    tipo  = rep(c("ACF", "PACF"), c(length(a$lag), length(p$lag))),
    banda = banda
  ) |>
    filter(lag > 0, is.finite(valor))
}

# --- 4. Pruebas: serie nacional ------------------------------------

pruebas_nacional <- bind_rows(lapply(niveles_variable, function(v) {
  probar_serie(
    x          = nacional[[v]],
    nivel      = "Nacional",
    region_ema = NA_character_,
    variable   = v,
    etiqueta   = unname(etiquetas_variable[v]),
    periodo    = periodo_panel
  )
}))

# --- 5. Pruebas: 12 regiones x 3 variables -------------------------
#
# Regiones: las 12 de la EMA. Variables: Y, casos_dengue y M_dengue
# (chikungunya regional es casi todo ceros; ver nota de cabecera).

variables_regionales <- c("Y", "casos_dengue", "M_dengue")

pruebas_regionales <- bind_rows(lapply(regiones_ema, function(r) {
  d <- panel |> filter(region_ema == r) |> arrange(mes)

  bind_rows(lapply(variables_regionales, function(v) {
    probar_serie(
      x          = d[[v]],
      nivel      = "Regional",
      region_ema = r,
      variable   = v,
      etiqueta   = unname(etiquetas_variable[v]),
      periodo    = periodo_panel
    )
  }))
}))

# --- 6. Tabla ordenada y guardado ----------------------------------

pruebas <- bind_rows(pruebas_nacional, pruebas_regionales) |>
  mutate(
    nivel    = factor(nivel, levels = c("Nacional", "Regional")),
    variable = factor(variable, levels = niveles_variable),
    region_ema = factor(region_ema, levels = regiones_ema)
  ) |>
  arrange(nivel, variable, region_ema) |>
  mutate(
    nivel      = as.character(nivel),
    variable   = as.character(variable),
    region_ema = as.character(region_ema)
  ) |>
  mutate(
    across(c(minimo, maximo, boxcox_lambda), ~ round(.x, 4)),
    across(c(adf_estadistico, kpss_estadistico), ~ round(.x, 4)),
    across(c(adf_p, kpss_p), ~ round(.x, 6))
  ) |>
  # dplyr::select explícito: MASS (y otros paquetes) exportan un select()
  # que, si queda adjunto después de dplyr, rompe este paso.
  dplyr::select(
    nivel, region_ema, variable, variable_etiqueta, periodo,
    n_obs, minimo, maximo, n_ceros, n_negativos,
    boxcox_aplica, boxcox_lambda, boxcox_nota,
    adf_estadistico, adf_p, adf_rezagos, adf_conclusion,
    kpss_estadistico, kpss_p, kpss_conclusion,
    diagnostico, requiere_diferenciacion, nota
  )

if (nrow(pruebas) != 5 + 12 * 3) {
  stop(
    "La tabla debe tener 41 filas (5 nacionales + 36 regionales) y tiene ",
    nrow(pruebas), ".", call. = FALSE
  )
}

if (anyDuplicated(pruebas[, c("nivel", "region_ema", "variable")])) {
  stop("Hay series duplicadas en la tabla de pruebas.", call. = FALSE)
}

saveRDS(pruebas, ruta_rds)
write_csv(pruebas, ruta_csv)

message(
  "\nTabla de pruebas guardada en:\n  ", ruta_csv, "\n  ", ruta_rds
)

# --- 7. Gráficos ACF / PACF ----------------------------------------
#
# Nacional: una fila por variable, dos columnas (ACF y PACF). Si el ADF
# no rechaza la raíz unitaria, se agregan las barras de la primera
# diferencia de esa misma serie.

series_nacionales <- niveles_variable
niveles_serie_nac <- unname(c(
  Y            = "Y: ocupación hotelera (%)",
  casos_dengue = "Casos de dengue",
  casos_chik   = "Casos de chikungunya",
  M_dengue     = "M: miedo excedente (dengue)",
  M_chik       = "M: miedo excedente (chik.)"
)[series_nacionales])

acf_nacional <- bind_rows(lapply(series_nacionales, function(v) {
  x <- as.numeric(nacional[[v]])
  x <- x[is.finite(x)]

  base <- calcular_acf_pacf(x, lag_max_nacional)
  if (is.null(base)) return(NULL)

  base$serie         <- etiqueta_serie(v)
  base$transformacion <- "Nivel"

  fila <- pruebas |> filter(nivel == "Nacional", variable == v)
  if (nrow(fila) == 1 && isTRUE(fila$requiere_diferenciacion)) {
    dif <- calcular_acf_pacf(diff(x), lag_max_nacional)
    if (!is.null(dif)) {
      dif$serie          <- etiqueta_serie(v)
      dif$transformacion <- "1ª diferencia"
      base <- bind_rows(base, dif)
    }
  }

  base
}))

acf_nacional <- acf_nacional |>
  mutate(
    serie = factor(serie, levels = niveles_serie_nac),
    transformacion = factor(transformacion, levels = c("Nivel", "1ª diferencia")),
    tipo = factor(tipo, levels = c("ACF", "PACF"))
  )

# Una fila por banda y signo, para dibujar ±1,96/sqrt(n) en cada panel.
bandas_nac <- acf_nacional |>
  distinct(serie, tipo, transformacion, banda) |>
  tidyr::expand_grid(signo = c(1, -1)) |>
  mutate(yintercept = signo * banda)

# Si alguna serie nacional fuera constante no tendría ACF/PACF: se rotula
# el panel para que la fila no desaparezca en silencio.
constantes_nac <- pruebas |>
  filter(nivel == "Nacional", str_detect(diagnostico, "constante")) |>
  tidyr::expand_grid(tipo = c("ACF", "PACF")) |>
  transmute(
    serie = factor(etiqueta_serie(variable), levels = niveles_serie_nac),
    tipo = factor(tipo, levels = c("ACF", "PACF")),
    lag = lag_max_nacional / 2,
    valor = 0,
    etiqueta = if_else(
      n_ceros == n_obs,
      "Serie constante\n(todos ceros)",
      "Serie constante\n(sin variación)"
    )
  )

colores_transformacion <- c("Nivel" = "#08519C", "1ª diferencia" = "#B2182B")

p_acf_nac <- ggplot(
  acf_nacional,
  aes(x = lag, y = valor, fill = transformacion)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = NA) +
  geom_hline(
    data = bandas_nac,
    aes(yintercept = yintercept, colour = transformacion),
    linetype = "dashed", linewidth = 0.25, inherit.aes = FALSE
  ) +
  geom_text(
    data = constantes_nac,
    aes(x = lag, y = valor, label = etiqueta),
    size = 2.6, colour = "grey35", lineheight = 0.9, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = colores_transformacion, name = NULL) +
  scale_colour_manual(values = colores_transformacion, guide = "none") +
  scale_x_continuous(breaks = seq(0, lag_max_nacional, by = 3)) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_grid(serie ~ tipo, drop = FALSE) +
  labs(
    title = "ACF y PACF de las series nacionales (2022-2025)",
    subtitle = paste0(
      "Nacional: Y y M promediadas por mes; casos sumados. ",
      "Barras azules: serie en nivel; barras rojas: primera diferencia ",
      "(solo si el ADF no rechaza la raíz unitaria)."
    ),
    x = "Rezago (meses)",
    y = "Autocorrelación",
    caption = paste0(
      "Bandas punteadas: ±1,96/sqrt(n) al 95% (azul: nivel; rojo: diferencia). ",
      "48 observaciones mensuales por serie (47 en la primera diferencia)."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(face = "bold", size = 8),
    strip.text.y = element_text(face = "bold", size = 7.5, angle = 0, hjust = 0),
    plot.title.position = "plot",
    plot.caption = element_text(hjust = 0, size = 7, colour = "grey35"),
    plot.caption.position = "plot",
    panel.spacing = grid::unit(0.35, "lines")
  )

ggsave(
  ruta_acf_nac, p_acf_nac,
  width = 10, height = 9.5, dpi = 150, bg = "white"
)

# Regional: una fila por región (12) y una columna por variable (3).
# En cada panel van la ACF (barras) y la PACF (línea y puntos), ambas en
# nivel; la tabla indica qué series requieren diferenciación.
etiquetas_facet <- c(
  Y            = "Y: ocupación\nhotelera (%)",
  casos_dengue = "Casos de\ndengue",
  M_dengue     = "M: miedo\nexcedente"
)

acf_regiones <- bind_rows(lapply(regiones_ema, function(r) {
  d <- panel |> filter(region_ema == r)

  bind_rows(lapply(variables_regionales, function(v) {
    a <- calcular_acf_pacf(d[[v]], lag_max_regiones)
    if (is.null(a)) return(NULL)
    a$region_ema <- r
    a$variable   <- v
    a
  }))
}))

acf_regiones <- acf_regiones |>
  mutate(
    region_ema = factor(region_ema, levels = regiones_ema),
    variable   = factor(variable, levels = variables_regionales),
    tipo       = factor(tipo, levels = c("ACF", "PACF"))
  )

bandas_reg <- acf_regiones |>
  distinct(region_ema, variable, banda) |>
  tidyr::expand_grid(signo = c(1, -1)) |>
  mutate(yintercept = signo * banda)

# Las series constantes no producen ACF/PACF: se rotulan dentro del panel
# para que no quede vacío.
constantes_reg <- pruebas |>
  filter(nivel == "Regional", str_detect(diagnostico, "constante")) |>
  transmute(
    region_ema = factor(region_ema, levels = regiones_ema),
    variable   = factor(variable, levels = variables_regionales),
    lag        = lag_max_regiones / 2,
    valor      = 0,
    etiqueta   = if_else(
      n_ceros == n_obs,
      "Serie constante\n(todos ceros)",
      "Serie constante\n(sin variación)"
    )
  )

n_dif_reg <- sum(
  pruebas$nivel == "Regional" & pruebas$requiere_diferenciacion %in% TRUE
)
n_reg <- sum(pruebas$nivel == "Regional")
n_reg_eval <- sum(
  pruebas$nivel == "Regional" & !is.na(pruebas$requiere_diferenciacion)
)
n_reg_noeval <- n_reg - n_reg_eval

caption_reg <- paste0(
  "Bandas punteadas: ±1,96/sqrt(n) al 95%. ",
  n_dif_reg, " de ", n_reg_eval, " series regionales evaluables requieren ",
  "diferenciación según el ADF",
  if (n_reg_noeval > 0) {
    paste0(
      "; ", n_reg_noeval, " no ",
      if (n_reg_noeval == 1) "es evaluable" else "son evaluables",
      " (serie constante)."
    )
  } else {
    "."
  },
  " Ver pruebas_estacionariedad.csv."
)

p_acf_reg <- ggplot() +
  geom_col(
    data = filter(acf_regiones, tipo == "ACF"),
    aes(x = lag, y = valor),
    fill = "#08519C", width = 0.7, alpha = 0.85
  ) +
  geom_hline(
    data = bandas_reg,
    aes(yintercept = yintercept),
    linetype = "dotted", colour = "grey30", linewidth = 0.25
  ) +
  geom_line(
    data = filter(acf_regiones, tipo == "PACF"),
    aes(x = lag, y = valor),
    colour = "#B2182B", linewidth = 0.45
  ) +
  geom_point(
    data = filter(acf_regiones, tipo == "PACF"),
    aes(x = lag, y = valor),
    colour = "#B2182B", size = 0.9
  ) +
  geom_text(
    data = constantes_reg,
    aes(x = lag, y = valor, label = etiqueta),
    size = 2.6, colour = "grey35", lineheight = 0.9
  ) +
  scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_grid(
    region_ema ~ variable,
    labeller = labeller(
      variable = etiquetas_facet,
      region_ema = label_wrap_gen(18)
    ),
    drop = FALSE
  ) +
  labs(
    title = "ACF y PACF por región EMA (series en nivel, 2022-2025)",
    subtitle = paste0(
      "Barras azules: ACF. Línea y puntos rojos: PACF. ",
      "Cada panel usa los 48 meses de la región."
    ),
    x = "Rezago (meses)",
    y = "Autocorrelación",
    caption = caption_reg
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(size = 7, face = "bold"),
    strip.text.y = element_text(size = 6.5, face = "bold", angle = 0, hjust = 0),
    panel.spacing = grid::unit(0.35, "lines"),
    plot.title.position = "plot",
    plot.caption = element_text(hjust = 0, size = 6.5, colour = "grey35"),
    plot.caption.position = "plot"
  )

ggsave(
  ruta_acf_reg, p_acf_reg,
  width = 10, height = 21, dpi = 150, bg = "white"
)

for (ruta in c(ruta_acf_nac, ruta_acf_reg)) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó el gráfico: ", ruta, call. = FALSE)
  }
}

message("Gráficos guardados en:\n  ", ruta_acf_nac, "\n  ", ruta_acf_reg)

# --- 8. Resumen en consola -----------------------------------------

# sprintf() cuenta bytes en Windows y desalinea los nombres con tilde
# (Bogotá, Región...). pad() rellena contando caracteres visibles.
pad <- function(x, ancho) {
  x <- as.character(x)
  paste0(x, strrep(" ", pmax(0, ancho - nchar(x, type = "width"))))
}

fmt_p <- function(p, ancho = 8) {
  ifelse(is.na(p), formatC("NA", width = ancho),
         formatC(p, format = "f", digits = 4, width = ancho))
}

celda_p <- function(adf_p, kpss_p) {
  a <- if (length(adf_p) == 0 || is.na(adf_p)) "  NA " else sprintf("%.3f", adf_p)
  k <- if (length(kpss_p) == 0 || is.na(kpss_p)) "  NA " else sprintf("%.3f", kpss_p)
  paste0(a, " / ", k)
}

cat("\n============================================================\n")
cat(" PRUEBAS DE ESTACIONARIEDAD (Script 14)\n")
cat("============================================================\n")
cat(sprintf(
  "Panel: %d filas | %d regiones | %d meses (%s)\n",
  nrow(panel), dplyr::n_distinct(panel$region_ema),
  dplyr::n_distinct(panel$mes), periodo_panel
))
cat(sprintf(
  "Series evaluadas: %d (%d nacionales + %d regionales)\n",
  nrow(pruebas), sum(pruebas$nivel == "Nacional"),
  sum(pruebas$nivel == "Regional")
))
cat("ADF  H0: raíz unitaria    -> p < 0,05: estacionaria\n")
cat("KPSS H0: estacionariedad  -> p < 0,05: no estacionaria\n")

# --- 8.1 Tabla nacional --------------------------------------------

cat("\n--- Series nacionales ---\n")
cat(sprintf(
  "%-15s %8s %8s %5s  %s\n",
  "Variable", "ADF p", "KPSS p", "n", "Diagnóstico"
))

nac_tab <- pruebas |> filter(nivel == "Nacional")
for (i in seq_len(nrow(nac_tab))) {
  cat(sprintf(
    "%-15s %s %s %5d  %s\n",
    nac_tab$variable[i],
    fmt_p(nac_tab$adf_p[i]),
    fmt_p(nac_tab$kpss_p[i]),
    nac_tab$n_obs[i],
    nac_tab$diagnostico[i]
  ))
}

# --- 8.2 p-valores regionales --------------------------------------

cat("\n--- Series regionales: p-valor ADF / KPSS ---\n")
cat(pad("Región", 30), " ", pad("Y", 15), " ",
    pad("Dengue (casos)", 15), " ", pad("M (miedo exc.)", 15), "\n", sep = "")

for (r in regiones_ema) {
  fila <- pruebas |> filter(nivel == "Regional", region_ema == r)
  celdas <- vapply(variables_regionales, function(v) {
    f <- fila[fila$variable == v, ]
    celda_p(f$adf_p, f$kpss_p)
  }, character(1))
  cat(pad(r, 30), " ", pad(celdas[1], 15), " ",
      pad(celdas[2], 15), " ", pad(celdas[3], 15), "\n", sep = "")
}

# --- 8.3 Conteo por diagnóstico ------------------------------------

cat("\n--- Diagnóstico por serie ---\n")
resumen_diag <- pruebas |>
  dplyr::count(diagnostico, name = "n") |>
  arrange(desc(n))

for (i in seq_len(nrow(resumen_diag))) {
  cat(sprintf("  %-62s %2d\n", resumen_diag$diagnostico[i],
              resumen_diag$n[i]))
}
cat(sprintf(
  "  %-62s %2d\n", "Total", sum(resumen_diag$n)
))

# --- 8.4 Series estacionarias y series que requieren diferencia ----

cat("\n--- Series estacionarias según el ADF (p < 0,05) ---\n")
for (v in niveles_variable) {
  sub <- pruebas |>
    filter(variable == v, requiere_diferenciacion %in% FALSE)
  if (nrow(sub) == 0) next
  quien <- ifelse(sub$nivel == "Nacional", "Nacional", sub$region_ema)
  cat(sprintf("  %-14s %d: %s\n", v, nrow(sub), paste(quien, collapse = ", ")))
}

n_dif <- sum(pruebas$requiere_diferenciacion %in% TRUE)
n_eval <- sum(!is.na(pruebas$requiere_diferenciacion))

cat("\n--- Series que requieren diferenciación (ADF p >= 0,05) ---\n")
cat(sprintf("Total: %d de %d series evaluables\n", n_dif, n_eval))

for (v in niveles_variable) {
  sub <- pruebas |>
    filter(variable == v, requiere_diferenciacion %in% TRUE)
  if (nrow(sub) == 0) next
  quien <- ifelse(sub$nivel == "Nacional", "Nacional", sub$region_ema)
  cat(sprintf("  %-14s %d: %s\n", v, nrow(sub), paste(quien, collapse = ", ")))
}

# --- 8.5 Series no evaluables --------------------------------------

no_eval <- pruebas |> filter(is.na(requiere_diferenciacion))
if (nrow(no_eval) > 0) {
  cat("\n--- Series no evaluables ---\n")
  for (i in seq_len(nrow(no_eval))) {
    quien <- if (no_eval$nivel[i] == "Nacional") {
      "Nacional"
    } else {
      no_eval$region_ema[i]
    }
    cat("  ", pad(quien, 30), " ", pad(no_eval$variable[i], 14), " ",
        no_eval$diagnostico[i], "\n", sep = "")
  }
}

# --- 8.6 Box-Cox ----------------------------------------------------

cat("\n--- Box-Cox (lambda óptimo, método Guerrero) ---\n")
bc <- pruebas |> filter(boxcox_aplica)
cat(sprintf(
  "Series con lambda estimado: %d de %d\n", nrow(bc), nrow(pruebas)
))

for (i in seq_len(nrow(bc))) {
  quien <- if (bc$nivel[i] == "Nacional") "Nacional" else bc$region_ema[i]
  cat("  ", pad(quien, 30), " ", pad(bc$variable[i], 14),
      sprintf("lambda = %7.4f", bc$boxcox_lambda[i]), "\n", sep = "")
}

# El motivo textual incluye el número de ceros o negativos de cada serie;
# para el resumen se agrupa por tipo de motivo.
motivos <- pruebas |>
  filter(!boxcox_aplica) |>
  mutate(
    motivo = dplyr::case_when(
      str_detect(boxcox_nota, "serie constante") ~ "Serie constante",
      str_detect(boxcox_nota, "valores negativos") ~
        "Valores negativos (M es un residuo con signo)",
      str_detect(boxcox_nota, "ceros") ~ "Contiene ceros (chikungunya o dengue)",
      TRUE ~ boxcox_nota
    )
  ) |>
  dplyr::count(motivo, name = "n") |>
  arrange(desc(n))

cat("\n  Series sin Box-Cox (motivo):\n")
for (i in seq_len(nrow(motivos))) {
  cat("    ", pad(motivos$motivo[i], 48), sprintf("%2d", motivos$n[i]),
      "\n", sep = "")
}

# --- 8.7 Cierre -----------------------------------------------------

cat("\nSalidas:\n")
cat("  ", ruta_csv, "\n", sep = "")
cat("  ", ruta_rds, "\n", sep = "")
cat("  ", ruta_acf_nac, "\n", sep = "")
cat("  ", ruta_acf_reg, "\n", sep = "")

cat("\nNotas:\n")
cat("  - Las pruebas se aplican a la serie en nivel, sin transformar y ",
    "sin ajuste estacional.\n", sep = "")
cat("  - Con 48 observaciones el ADF tiene poca potencia: el KPSS se usa ",
    "como contraste.\n", sep = "")
cat("  - Las series con ceros o valores negativos no admiten Box-Cox; ",
    "el motivo queda en la columna boxcox_nota.\n", sep = "")
if (nrow(no_eval) > 0) {
  cat("  - ", nrow(no_eval), " serie(s) sin ADF evaluable",
      " (constante o sin observaciones suficientes); el motivo está en la ",
      "tabla y en la sección anterior.\n", sep = "")
}

message("\n✅ Script 14 finalizado.")

