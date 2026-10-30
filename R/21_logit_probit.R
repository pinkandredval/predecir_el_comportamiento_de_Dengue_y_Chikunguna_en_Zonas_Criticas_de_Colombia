# ============================================================
# 21_logit_probit.R
# Teoría del Riesgo — Tabla 4:
#   scorecard individual de severidad (logit/probit) y
#   probit de alerta temprana con el panel EMA
# ============================================================
#
# OBJETIVO
#   Construir dos ejercicios de clasificación binaria:
#
#   1. SCORECARD DE SEVERIDAD (nivel individual)
#      Estimar la probabilidad de que un caso de dengue requiera
#      hospitalización (hosp = 1) en función de la edad, el sexo, el
#      estrato, el tipo de régimen de salud y el área de residencia.
#      Se ajustan un logit y un probit con la misma especificación, se
#      evalúan fuera de muestra (30% de prueba) con AUC, Gini y KS, y se
#      dibujan la curva ROC y los efectos marginales del logit.
#
#   2. PROBIT DE ALERTA TEMPRANA (nivel panel región-mes)
#      Definir una alerta cuando los casos de dengue de una región EMA
#      superan su percentil 90 regional y modelar esa alerta con el
#      índice de miedo M_dengue, el rezago de casos y la lluvia. El
#      rezago se calcula dentro de cada región para no mezclar series.
#
# PREPARACIÓN DEL SCORECARD
#   pac_hos es un código de texto del SIVIGILA: "1" = hospitalizado,
#   "2" = no hospitalizado. Se convierte a numérico y se crea
#
#       hosp = 1 si pac_hos == "1", 0 si no, NA si pac_hos es NA.
#
#   edad viene acompañada de uni_med, la unidad de medida:
#       1 = años, 2 = meses, 3 = días, 4 = horas.
#   Se homologa todo a años:
#
#       edad_a = edad                     si uni_med = 1
#              = edad / 12                si uni_med = 2
#              = edad / 365.25            si uni_med = 3
#              = edad / (365.25 * 24)     si uni_med = 4
#
#   Si uni_med no existe o viene vacío, se asume que edad ya está en
#   años. Los códigos desconocidos quedan en NA y se excluyen del ajuste
#   (en esta base solo aparecen 1, 2, 3 y 4).
#
#   estrato se convierte en factor con niveles "1" a "6" y "Sin dato".
#   Se asignan a "Sin dato" los NA, las cadenas vacías, el texto
#   "SIN DATO" y cualquier código fuera de 1-6. sexo, tip_ss y area se
#   convierten en factores.
#
#   No se imputan valores: el modelo se estima con las filas completas.
#   En la base de dengue preparada no hay NA en hosp, edad_a, sexo,
#   tip_ss ni area; estrato aporta 34.310 NA que pasan a "Sin dato".
#
# PARTICIÓN ENTRE ENTRENAMIENTO Y PRUEBA
#   set.seed(123) y muestra aleatoria simple: 70% entrena, 30% prueba.
#   La semilla fija la partición para que las métricas sean
#   reproducibles. No se estratifica por hosp: con 748.037 casos y una
#   prevalencia cercana al 40%, ambas muestras contienen las dos clases
#   con probabilidad prácticamente uno.
#
# MODELOS DEL SCORECARD
#   logit  : glm(hosp ~ edad_a + sexo + estrato + tip_ss + area,
#                family = binomial("logit"),  data = entrena)
#   probit : la misma fórmula con family = binomial("probit").
#
#   El logit se usa como modelo principal del scorecard (es el que se
#   evalúa y el que alimenta los efectos marginales). El probit se
#   conserva como contraste de la misma especificación.
#
# EVALUACIÓN
#   AUC     : pROC::roc(prueba$hosp, predict(logit, prueba,
#                                            type = "response")).
#   Gini    : 2 * AUC - 1.
#   KS      : max(abs(ecdf(p[hosp == 1])(p) - ecdf(p[hosp == 0])(p))),
#             con p = probabilidad predicha en prueba.
#
#   Si pROC no está instalado, se usa una implementación interna
#   equivalente: AUC por rangos (Mann-Whitney, empates a 0,5) y curva
#   ROC por umbrales. La ruta principal del script sigue siendo
#   pROC::roc(); la alternativa solo evita que el script se detenga en
#   un entorno sin ese paquete. El CSV registra la fuente usada.
#
#   Los efectos marginales se obtienen con
#   marginaleffects::avg_slopes(logit). Si marginaleffects falla o no
#   está instalado, se emite un warning y se usa broom::tidy(logit),
#   que devuelve coeficientes en log-odds (no efectos marginales). La
#   gráfica se rotula de acuerdo con la vía que se haya usado.
#
# PANEL DE ALERTA TEMPRANA
#   El panel D_13 tiene 12 regiones EMA × 48 meses (2022-01 a 2025-12).
#   Se ordena por region_ema y mes y, dentro de cada región:
#
#       alerta    = 1 si casos_dengue > quantile(casos_dengue, 0.90)
#       lag_casos = dplyr::lag(casos_dengue, 1)
#
#   El rezago se calcula dentro del group_by(region_ema) para que la
#   última observación de una región no se convierta en el "mes
#   anterior" de la siguiente. El primer mes de cada región queda en NA
#   y se excluye del modelo: se pierden 12 filas y quedan 564.
#
#   La columna de lluvia puede llamarse lluvia o lluvia_mm; el script
#   detecta la que exista y la homologa a lluvia. En D_13 la columna es
#   lluvia_mm.
#
#   Modelo:
#       glm(alerta ~ M_dengue + lag_casos + lluvia,
#           family = binomial("probit"), data = datos_alerta)
#
#   Se evalúa con AUC y se grafica la probabilidad predicha contra
#   M_dengue y contra el tiempo (promedio regional y banda 10-90).
#
# INSUMOS
#   data/D_02/dengue_prep.rds
#       Base individual preparada de dengue. Si no existe, se busca
#       data/imputed/dengue_imputado.rds y, como respaldo adicional,
#       data/D_06/dengue_imputado.rds.
#   data/D_13/panel_ema_den_chik.rds
#       Panel región-mes EMA con casos_dengue, M_dengue y lluvia_mm.
#
# SALIDAS
#   data/D_21/scorecard_metricas.csv
#       modelo, auc, ks, gini, n_entrena, n_prueba,
#       prevalencia_prueba, aic, fuente_roc.
#   data/D_21/probit_alerta_metricas.csv
#       modelo, auc, n_obs, n_regiones, n_alertas, prevalencia, aic,
#       fuente_roc.
#   data/D_21/probit_alerta_predicciones.csv
#       region_ema, mes, casos_dengue, lag_casos, M_dengue, lluvia,
#       alerta, prob_alerta (564 filas).
#   graficas/G_21/curva_roc.png
#   graficas/G_21/efectos_marginales.png
#   graficas/G_21/probit_alerta.png
#
# ADVERTENCIAS DE INTERPRETACIÓN
#   - El AUC del scorecard es moderado (alrededor de 0,57): la edad, el
#     estrato y el régimen aportan señal, pero la hospitalización no se
#     explica bien con estas covariables. Un AUC bajo no es un error del
#     script; es el resultado del poder predictivo de la especificación.
#   - La alerta del panel es un evento raro por construcción (10% de los
#     meses de cada región). El AUC se calcula sobre 55 alertas y 509
#     no alertas; es una evaluación in-sample del probit, no una
#     validación temporal.
#   - Bogotá D.C. tiene casos_dengue = 0 en los 48 meses del panel (así
#     viene de D_10/D_13). Por la regla > p90 no genera ninguna alerta;
#     se conserva en el modelo porque la regla se aplica por región.
#   - Las asociaciones son descriptivas, no causales.
#
# REGLAS DEL GUION
#   - Todas las funciones no base se llaman con paquete::funcion().
#   - Todas las rutas usan here::here().
#   - No se ejecuta git add ni git commit.
# ============================================================

# --- 0. Paquetes, rutas y carpetas ---------------------------------

paquetes <- c(
  "dplyr", "ggplot2", "readr", "tibble", "here", "broom", "gridExtra"
)

faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]
if (length(faltantes) > 0L) {
  stop(
    "Faltan paquetes: ", paste(faltantes, collapse = ", "), ". ",
    "Instálelos con renv::install(c(",
    paste0("\"", faltantes, "\"", collapse = ", "), ")).",
    call. = FALSE
  )
}

# pROC y marginaleffects son la vía principal, pero el script puede
# seguir funcionando sin ellos. Se avisa una sola vez al inicio.
tiene_pROC <- requireNamespace("pROC", quietly = TRUE)
tiene_marginaleffects <- requireNamespace("marginaleffects", quietly = TRUE)

if (!tiene_pROC) {
  message(
    "pROC no está instalado: se usará una ROC interna equivalente ",
    "(AUC por rangos). Instale pROC con renv::install(\"pROC\") para ",
    "usar la ruta canónica."
  )
}

if (!tiene_marginaleffects) {
  message(
    "marginaleffects no está instalado: se usará broom::tidy() como ",
    "respaldo para los efectos marginales."
  )
}

dir.create(here::here("data", "D_21"), recursive = TRUE, showWarnings = FALSE)
dir.create(
  here::here("graficas", "G_21"),
  recursive = TRUE, showWarnings = FALSE
)

rutas_scorecard <- c(
  here::here("data", "D_02", "dengue_prep.rds"),
  here::here("data", "imputed", "dengue_imputado.rds"),
  here::here("data", "D_06", "dengue_imputado.rds")
)

rutas_scorecard_existentes <- rutas_scorecard[file.exists(rutas_scorecard)]
if (length(rutas_scorecard_existentes) == 0L) {
  stop(
    "No se encontró la base individual de dengue. Se buscó en:\n  ",
    paste(rutas_scorecard, collapse = "\n  "), "\n",
    "Ejecute primero R/02_preparacion.R (o R/06_imputación.R).",
    call. = FALSE
  )
}

ruta_scorecard <- rutas_scorecard_existentes[1L]
ruta_panel <- here::here("data", "D_13", "panel_ema_den_chik.rds")

ruta_metricas_scorecard <- here::here(
  "data", "D_21", "scorecard_metricas.csv"
)
ruta_metricas_alerta <- here::here(
  "data", "D_21", "probit_alerta_metricas.csv"
)
ruta_predicciones_alerta <- here::here(
  "data", "D_21", "probit_alerta_predicciones.csv"
)
ruta_png_roc <- here::here("graficas", "G_21", "curva_roc.png")
ruta_png_efectos <- here::here(
  "graficas", "G_21", "efectos_marginales.png"
)
ruta_png_alerta <- here::here("graficas", "G_21", "probit_alerta.png")

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel de alerta temprana: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

# --- 1. Utilidades -------------------------------------------------

# Convierte edad a años según el código de uni_med. Los códigos son los
# del SIVIGILA: 1 = años, 2 = meses, 3 = días, 4 = horas. También acepta
# las etiquetas en texto por si la base viene con nombres en lugar de
# códigos. Un código desconocido produce NA (no se inventa la unidad).
convertir_edad_anios <- function(edad, uni_med) {
  edad <- suppressWarnings(as.numeric(as.character(edad)))
  unidad <- tolower(trimws(as.character(uni_med)))

  factor_anios <- rep(NA_real_, length(unidad))

  es_anios <- unidad %in% c("1", "año", "años", "ano", "anos") |
    (!is.na(unidad) & grepl("^a(ñ|n)o", unidad))
  es_meses <- unidad %in% c("2", "mes", "meses") |
    (!is.na(unidad) & grepl("^mes", unidad))
  es_dias <- unidad %in% c("3", "día", "días", "dia", "dias") |
    (!is.na(unidad) & grepl("^d[ií]a", unidad))
  es_horas <- unidad %in% c("4", "hora", "horas") |
    (!is.na(unidad) & grepl("^hora", unidad))
  sin_unidad <- is.na(unidad) | unidad == ""

  factor_anios[es_anios] <- 1
  factor_anios[es_meses] <- 1 / 12
  factor_anios[es_dias] <- 1 / 365.25
  factor_anios[es_horas] <- 1 / (365.25 * 24)
  factor_anios[sin_unidad] <- 1

  edad * factor_anios
}

# Estrato: "Sin dato" para NA, vacíos, el texto SIN DATO y códigos fuera
# de 1-6. Se fija el orden de los niveles para que la categoría de
# referencia sea el estrato 1 y "Sin dato" quede al final.
normalizar_estrato <- function(x) {
  x <- trimws(as.character(x))
  validos <- as.character(1:6)

  x[is.na(x) | x == "" | toupper(x) == "SIN DATO"] <- "Sin dato"
  x[!(x %in% validos)] <- "Sin dato"

  factor(x, levels = c(validos, "Sin dato"))
}

# ROC interna, usada solo si pROC no está disponible. Devuelve la misma
# información que se extrae de un objeto pROC::roc: specificities,
# sensitivities y auc. La orientación se fija en "auto": si el AUC
# crudo es menor que 0,5 se invierte el predictor, igual que pROC.
#
# El AUC se calcula por rangos (Mann-Whitney), que con empates a 0,5 es
# algebraicamente igual al área bajo la curva ROC empírica. No se usa
# outer() porque la muestra de prueba tiene 224.412 filas y una matriz
# de pares ocuparía decenas de GB.
roc_interna <- function(response, predictor) {
  y <- as.integer(response)
  p <- as.numeric(predictor)

  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]

  if (length(y) == 0L) {
    stop("No hay observaciones completas para calcular la ROC.",
         call. = FALSE)
  }
  if (!all(y %in% c(0L, 1L))) {
    stop("La respuesta de la ROC debe ser 0/1.", call. = FALSE)
  }

  n1 <- as.numeric(sum(y == 1L))
  n0 <- as.numeric(sum(y == 0L))
  if (n1 == 0 || n0 == 0) {
    stop(
      "Se necesitan las dos clases (0 y 1) para calcular la ROC.",
      call. = FALSE
    )
  }

  rangos <- rank(p)
  auc <- (sum(rangos[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)

  if (auc < 0.5) {
    p <- -p
    auc <- 1 - auc
  }

  orden <- order(p, decreasing = TRUE)
  p_orden <- p[orden]
  y_orden <- y[orden]

  # Un punto por valor único del predictor: todos los empates entran
  # juntos en el mismo umbral.
  indices <- which(!duplicated(p_orden))
  tp <- cumsum(y_orden == 1L)[indices]
  fp <- cumsum(y_orden == 0L)[indices]

  sensibilidad <- tp / n1
  especificidad <- 1 - fp / n0

  # Extremos de la curva: (0, 0) y (1, 1) en (1 - especificidad,
  # sensibilidad).
  especificidad <- c(1, especificidad, 0)
  sensibilidad <- c(0, sensibilidad, 1)

  orden_curva <- order(1 - especificidad, sensibilidad)

  list(
    specificities = especificidad[orden_curva],
    sensitivities = sensibilidad[orden_curva],
    auc = auc,
    fuente = "interna"
  )
}

# Calcula la ROC por la vía principal (pROC) y, si el paquete no está
# instalado o falla, por la vía interna. El argumento usar_pROC permite
# probar ambas rutas de forma explícita.
calcular_roc <- function(response, predictor, usar_pROC = tiene_pROC) {
  if (usar_pROC) {
    roc <- tryCatch(
      pROC::roc(
        response = response,
        predictor = predictor,
        quiet = TRUE,
        percent = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(roc)) {
      return(list(
        specificities = as.numeric(roc$specificities),
        sensitivities = as.numeric(roc$sensitivities),
        auc = as.numeric(pROC::auc(roc)),
        fuente = "pROC"
      ))
    }

    warning(
      "pROC::roc() falló; se usa la ROC interna equivalente.",
      call. = FALSE
    )
  }

  roc_interna(response, predictor)
}

# KS: máxima distancia vertical entre las ECDF de casos y no casos,
# evaluadas en las probabilidades predichas. Es la fórmula del guion.
calcular_ks <- function(probabilidad, respuesta) {
  y <- as.integer(respuesta)
  p <- as.numeric(probabilidad)

  if (length(p) != length(y)) {
    stop("probabilidad y respuesta deben tener la misma longitud.",
         call. = FALSE)
  }
  if (!all(y %in% c(0L, 1L))) {
    stop("La respuesta del KS debe ser 0/1.", call. = FALSE)
  }
  if (!any(y == 1L) || !any(y == 0L)) {
    stop("Se necesitan las dos clases para calcular el KS.",
         call. = FALSE)
  }

  max(abs(
    stats::ecdf(p[y == 1L])(p) - stats::ecdf(p[y == 0L])(p)
  ))
}

# Efectos marginales del logit. Vía principal: avg_slopes(). Respaldo:
# broom::tidy(), que devuelve coeficientes (log-odds), no efectos
# marginales; la etiqueta del eje y el subtítulo lo hacen explícito.
calcular_efectos <- function(modelo, usar_marginaleffects = TRUE) {
  if (usar_marginaleffects) {
    efectos <- tryCatch(
      marginaleffects::avg_slopes(modelo),
      error = function(e) {
        warning(
          "marginaleffects::avg_slopes() falló: ", conditionMessage(e),
          ". Se usa broom::tidy().",
          call. = FALSE
        )
        NULL
      }
    )

    if (!is.null(efectos) && nrow(efectos) > 0L) {
      efectos <- as.data.frame(efectos)

      if (!"term" %in% names(efectos)) {
        stop(
          "avg_slopes() no devolvió la columna term; no se puede ",
          "construir la gráfica.",
          call. = FALSE
        )
      }

      if ("contrast" %in% names(efectos)) {
        efectos$etiqueta <- ifelse(
          is.na(efectos$contrast) | efectos$contrast == "",
          efectos$term,
          paste0(efectos$term, ": ", efectos$contrast)
        )
      } else {
        efectos$etiqueta <- efectos$term
      }

      if (!"conf.low" %in% names(efectos)) {
        efectos$conf.low <- NA_real_
      }
      if (!"conf.high" %in% names(efectos)) {
        efectos$conf.high <- NA_real_
      }

      return(list(
        datos = efectos,
        etiqueta_x = "Efecto marginal promedio",
        fuente = "marginaleffects::avg_slopes()",
        es_ame = TRUE
      ))
    }

    if (is.null(efectos)) {
      # El warning ya se emitió en el tryCatch.
    } else {
      warning(
        "avg_slopes() devolvió una tabla vacía; se usa broom::tidy().",
        call. = FALSE
      )
    }
  } else {
    warning(
      "El paquete marginaleffects no está instalado; se usa ",
      "broom::tidy().",
      call. = FALSE
    )
  }

  efectos <- broom::tidy(modelo, conf.int = TRUE)
  efectos <- efectos[efectos$term != "(Intercept)", , drop = FALSE]

  if (!"conf.low" %in% names(efectos)) {
    efectos$conf.low <- NA_real_
  }
  if (!"conf.high" %in% names(efectos)) {
    efectos$conf.high <- NA_real_
  }

  efectos$etiqueta <- efectos$term

  list(
    datos = efectos,
    etiqueta_x = "Coeficiente (log-odds)",
    fuente = "broom::tidy()",
    es_ame = FALSE
  )
}

# --- 2. Lectura y preparación del scorecard ------------------------

scorecard <- readRDS(ruta_scorecard)

if (!inherits(scorecard, "data.frame")) {
  stop(
    "El objeto leído no es un data.frame: ", ruta_scorecard, ".",
    call. = FALSE
  )
}

requeridas_scorecard <- c("pac_hos", "edad", "sexo", "tip_ss", "area")
faltan_scorecard <- setdiff(requeridas_scorecard, names(scorecard))

if (length(faltan_scorecard) > 0L) {
  stop(
    "Faltan columnas en el scorecard: ",
    paste(faltan_scorecard, collapse = ", "), ". Columnas disponibles: ",
    paste(names(scorecard), collapse = ", "), ".",
    call. = FALSE
  )
}

# Si la columna no existe, se crea para poder aplicar la regla del
# enunciado sin detener el script. uni_med ausente equivale a edad en
# años; estrato ausente equivale a "Sin dato" en todas las filas.
if (!"uni_med" %in% names(scorecard)) {
  scorecard$uni_med <- NA_character_
}
if (!"estrato" %in% names(scorecard)) {
  scorecard$estrato <- "Sin dato"
}

datos_scorecard <- scorecard |>
  dplyr::transmute(
    edad = suppressWarnings(as.numeric(as.character(edad))),
    uni_med = as.character(uni_med),
    sexo = as.character(sexo),
    estrato = as.character(estrato),
    tip_ss = as.character(tip_ss),
    area = as.character(area),
    pac_hos_num = suppressWarnings(as.numeric(as.character(pac_hos)))
  ) |>
  dplyr::mutate(
    # hosp: 1 solo si el código original es "1"; 0 para cualquier otro
    # código observado; NA si pac_hos es NA.
    hosp = ifelse(
      is.na(pac_hos_num), NA_integer_,
      ifelse(pac_hos_num == 1, 1L, 0L)
    ),
    edad_a = convertir_edad_anios(edad, uni_med),
    estrato = normalizar_estrato(estrato),
    sexo = factor(sexo),
    tip_ss = factor(tip_ss),
    area = factor(area)
  )

# Códigos de uni_med desconocidos: se informan, no se imputan.
unidad_desconocida <- unique(
  datos_scorecard$uni_med[
    !is.na(datos_scorecard$uni_med) &
      trimws(datos_scorecard$uni_med) != "" &
      is.na(datos_scorecard$edad_a) &
      !is.na(datos_scorecard$edad)
  ]
)
if (length(unidad_desconocida) > 0L) {
  warning(
    "Hay unidades de edad no reconocidas en uni_med: ",
    paste(unidad_desconocida, collapse = ", "),
    ". Esas filas quedan con edad_a = NA y se excluyen del modelo.",
    call. = FALSE
  )
}

if (anyNA(datos_scorecard$hosp)) {
  message(
    "Se excluyen ", sum(is.na(datos_scorecard$hosp)),
    " filas con hosp = NA (pac_hos sin dato)."
  )
}
if (anyNA(datos_scorecard$edad_a)) {
  message(
    "Se excluyen ", sum(is.na(datos_scorecard$edad_a)),
    " filas con edad_a = NA."
  )
}
if (anyNA(datos_scorecard$sexo) ||
    anyNA(datos_scorecard$tip_ss) ||
    anyNA(datos_scorecard$area)) {
  message(
    "Se excluyen filas con sexo, tip_ss o area = NA."
  )
}

datos_modelo <- datos_scorecard |>
  dplyr::filter(
    !is.na(hosp),
    !is.na(edad_a),
    !is.na(sexo),
    !is.na(tip_ss),
    !is.na(area)
  )

if (nrow(datos_modelo) == 0L) {
  stop("No quedaron filas completas para el scorecard.", call. = FALSE)
}

if (length(unique(datos_modelo$hosp)) < 2L) {
  stop(
    "El scorecard quedó con una sola clase en hosp; no se puede ",
    "estimar un modelo binomial.",
    call. = FALSE
  )
}

if (any(datos_modelo$edad_a < 0)) {
  warning(
    "Hay edades negativas en edad_a; revise uni_med.",
    call. = FALSE
  )
}

n_hosp <- sum(datos_modelo$hosp == 1L)
n_no_hosp <- sum(datos_modelo$hosp == 0L)

cat("\n============================================================\n")
cat(" SCORECARD DE SEVERIDAD — PREPARACIÓN\n")
cat("============================================================\n")
cat(sprintf(
  "Base: %s\nFilas originales: %s | Filas completas: %s\n",
  ruta_scorecard,
  format(nrow(scorecard), big.mark = ".", decimal.mark = ","),
  format(nrow(datos_modelo), big.mark = ".", decimal.mark = ",")
))
cat(sprintf(
  "Hospitalizados: %s (%.2f%%) | No hospitalizados: %s (%.2f%%)\n",
  format(n_hosp, big.mark = ".", decimal.mark = ","),
  100 * n_hosp / nrow(datos_modelo),
  format(n_no_hosp, big.mark = ".", decimal.mark = ","),
  100 * n_no_hosp / nrow(datos_modelo)
))
cat(sprintf(
  "edad_a (años): mínimo %.4f | mediana %.2f | media %.2f | máximo %.2f\n",
  min(datos_modelo$edad_a),
  stats::median(datos_modelo$edad_a),
  mean(datos_modelo$edad_a),
  max(datos_modelo$edad_a)
))
cat(sprintf(
  "Estrato: %d categorías | 'Sin dato': %s filas\n",
  nlevels(datos_modelo$estrato),
  format(
    sum(datos_modelo$estrato == "Sin dato"),
    big.mark = ".", decimal.mark = ","
  )
))

# --- 3. Partición 70/30 --------------------------------------------

set.seed(123)

n_scorecard <- nrow(datos_modelo)
n_entrena <- floor(0.70 * n_scorecard)

indices_entrena <- sample.int(
  n = n_scorecard,
  size = n_entrena,
  replace = FALSE
)

entrena <- datos_modelo[indices_entrena, , drop = FALSE]
prueba <- datos_modelo[-indices_entrena, , drop = FALSE]

if (nrow(entrena) + nrow(prueba) != n_scorecard) {
  stop("La partición 70/30 no conserva todas las filas.", call. = FALSE)
}
if (length(unique(entrena$hosp)) < 2L ||
    length(unique(prueba$hosp)) < 2L) {
  stop(
    "Alguna de las muestras de entrenamiento o prueba quedó con una ",
    "sola clase de hosp.",
    call. = FALSE
  )
}

cat(sprintf(
  "\nPartición (set.seed(123)): entrena %s (%.2f%%) | prueba %s (%.2f%%)\n",
  format(nrow(entrena), big.mark = ".", decimal.mark = ","),
  100 * nrow(entrena) / n_scorecard,
  format(nrow(prueba), big.mark = ".", decimal.mark = ","),
  100 * nrow(prueba) / n_scorecard
))

# --- 4. Modelos logit y probit -------------------------------------

formula_scorecard <- hosp ~ edad_a + sexo + estrato + tip_ss + area

logit <- stats::glm(
  formula = formula_scorecard,
  family = stats::binomial("logit"),
  data = entrena
)

probit <- stats::glm(
  formula = formula_scorecard,
  family = stats::binomial("probit"),
  data = entrena
)

if (!isTRUE(logit$converged)) {
  stop("El logit del scorecard no convergió.", call. = FALSE)
}
if (!isTRUE(probit$converged)) {
  stop("El probit del scorecard no convergió.", call. = FALSE)
}

if (anyNA(stats::coef(logit)) || anyNA(stats::coef(probit))) {
  warning(
    "Hay coeficientes no estimados (aliased) en el scorecard; revise ",
    "la colinealidad de los factores.",
    call. = FALSE
  )
}

cat(sprintf(
  "\nLogit:  AIC = %.2f | coeficientes = %d\n",
  stats::AIC(logit), length(stats::coef(logit))
))
cat(sprintf(
  "Probit: AIC = %.2f | coeficientes = %d\n",
  stats::AIC(probit), length(stats::coef(probit))
))

# --- 5. Evaluación fuera de muestra: AUC, Gini y KS ----------------

p_logit <- as.numeric(
  stats::predict(logit, newdata = prueba, type = "response")
)
p_probit <- as.numeric(
  stats::predict(probit, newdata = prueba, type = "response")
)

if (length(p_logit) != nrow(prueba) ||
    length(p_probit) != nrow(prueba)) {
  stop(
    "Las predicciones no cubren todas las filas de prueba.",
    call. = FALSE
  )
}
if (anyNA(p_logit) || anyNA(p_probit) ||
    any(!is.finite(p_logit)) || any(!is.finite(p_probit))) {
  stop(
    "Hay probabilidades predichas no finitas en la muestra de prueba.",
    call. = FALSE
  )
}
if (any(p_logit < 0 | p_logit > 1) ||
    any(p_probit < 0 | p_probit > 1)) {
  stop(
    "Las probabilidades predichas están fuera de [0, 1].",
    call. = FALSE
  )
}

# Equivale a:
#   pROC::roc(prueba$hosp, predict(logit, prueba, type = "response"))
# calcular_roc() usa pROC::roc() cuando está disponible y, si no, la
# implementación interna equivalente.
roc_logit <- calcular_roc(prueba$hosp, p_logit)
roc_probit <- calcular_roc(prueba$hosp, p_probit)

auc_logit <- roc_logit$auc
auc_probit <- roc_probit$auc
gini_logit <- 2 * auc_logit - 1
gini_probit <- 2 * auc_probit - 1
ks_logit <- calcular_ks(p_logit, prueba$hosp)
ks_probit <- calcular_ks(p_probit, prueba$hosp)

if (auc_logit < 0.5 || auc_probit < 0.5) {
  warning(
    "Algún AUC quedó por debajo de 0,5; revise la orientación de las ",
    "probabilidades.",
    call. = FALSE
  )
}

cat("\n--- Métricas fuera de muestra (30% de prueba) ---\n")
cat(sprintf(
  "Logit : AUC = %.4f | Gini = %.4f | KS = %.4f\n",
  auc_logit, gini_logit, ks_logit
))
cat(sprintf(
  "Probit: AUC = %.4f | Gini = %.4f | KS = %.4f\n",
  auc_probit, gini_probit, ks_probit
))
cat(sprintf(
  "Fuente de la ROC: %s\n",
  if (identical(roc_logit$fuente, "pROC")) {
    "pROC::roc()"
  } else {
    "implementación interna (pROC no disponible)"
  }
))

# --- 6. Curva ROC --------------------------------------------------

datos_roc <- dplyr::bind_rows(
  tibble::tibble(
    modelo = "Logit",
    fpr = 1 - roc_logit$specificities,
    tpr = roc_logit$sensitivities
  ),
  tibble::tibble(
    modelo = "Probit",
    fpr = 1 - roc_probit$specificities,
    tpr = roc_probit$sensitivities
  )
) |>
  dplyr::arrange(modelo, fpr, tpr)

if (anyNA(datos_roc$fpr) || anyNA(datos_roc$tpr) ||
    max(datos_roc$fpr) > 1.001 || max(datos_roc$tpr) > 1.001 ||
    min(datos_roc$fpr) < -0.001 || min(datos_roc$tpr) < -0.001) {
  stop(
    "La curva ROC no está en la escala 0-1; revise specificities y ",
    "sensitivities.",
    call. = FALSE
  )
}

etiquetas_roc <- c(
  Logit = sprintf("Logit (AUC = %.4f)", auc_logit),
  Probit = sprintf("Probit (AUC = %.4f)", auc_probit)
)

caption_roc <- sprintf(
  paste0(
    "Muestra de prueba: %s casos (%.1f%% hospitalizados). ",
    "AUC logit = %.4f; AUC probit = %.4f. ",
    "Gini y KS del logit: %.4f y %.4f. Curva calculada con %s."
  ),
  format(nrow(prueba), big.mark = ".", decimal.mark = ","),
  100 * mean(prueba$hosp),
  auc_logit,
  auc_probit,
  gini_logit,
  ks_logit,
  if (identical(roc_logit$fuente, "pROC")) {
    "pROC::roc()"
  } else {
    "una implementación interna equivalente a pROC"
  }
)
caption_roc <- paste(strwrap(caption_roc, width = 130), collapse = "\n")

p_roc <- ggplot2::ggplot(
  datos_roc,
  ggplot2::aes(x = fpr, y = tpr, colour = modelo, linetype = modelo)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey55",
    linewidth = 0.35,
    linetype = "dotted"
  ) +
  ggplot2::geom_path(linewidth = 0.90) +
  ggplot2::scale_colour_manual(
    values = c("Logit" = "#2166AC", "Probit" = "#B2182B"),
    labels = etiquetas_roc,
    name = NULL
  ) +
  ggplot2::scale_linetype_manual(
    values = c("Logit" = "solid", "Probit" = "dashed"),
    labels = etiquetas_roc,
    name = NULL
  ) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  ggplot2::labs(
    title = "Curva ROC del scorecard de hospitalización",
    subtitle = sprintf(
      paste0(
        "Logit y probit con edad_a, sexo, estrato, tip_ss y area | ",
        "prueba: %s casos"
      ),
      format(nrow(prueba), big.mark = ".", decimal.mark = ",")
    ),
    x = "Tasa de falsos positivos (1 - especificidad)",
    y = "Tasa de verdaderos positivos (sensibilidad)",
    caption = caption_roc
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "left",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_png_roc,
  plot = p_roc,
  width = 7.8,
  height = 6.8,
  dpi = 200,
  bg = "white"
)

# --- 7. Efectos marginales del logit -------------------------------

# Vía principal: marginaleffects::avg_slopes(logit). Si falla, la
# función emite un warning y devuelve broom::tidy(logit).
efectos_logit <- calcular_efectos(
  logit,
  usar_marginaleffects = tiene_marginaleffects
)

datos_efectos <- efectos_logit$datos

if (nrow(datos_efectos) == 0L) {
  stop("La tabla de efectos marginales quedó vacía.", call. = FALSE)
}

if (!"estimate" %in% names(datos_efectos)) {
  stop(
    "La tabla de efectos no tiene columna estimate; revise ",
    "marginaleffects/broom.",
    call. = FALSE
  )
}

if (efectos_logit$es_ame) {
  titulo_efectos <- "Efectos marginales promedio del logit de hospitalización"
  subtitulo_efectos <- paste0(
    "marginaleffects::avg_slopes() | intervalo de confianza al 95% | ",
    "efectos sobre la probabilidad de hospitalización"
  )
} else {
  titulo_efectos <- "Coeficientes del logit de hospitalización (respaldo)"
  subtitulo_efectos <- paste0(
    "broom::tidy() | coeficientes en log-odds, no efectos marginales | ",
    "intervalo de confianza al 95%"
  )
}

caption_efectos <- paste0(
  "Modelo logit: hosp ~ edad_a + sexo + estrato + tip_ss + area, ",
  "estimado con ",
  format(nrow(entrena), big.mark = ".", decimal.mark = ","),
  " casos de entrenamiento. ",
  if (efectos_logit$es_ame) {
    "avg_slopes() promedia los efectos sobre la probabilidad predicha."
  } else {
    paste0(
      "marginaleffects no está disponible o falló; se grafican los ",
      "coeficientes de broom::tidy(). Un coeficiente positivo aumenta ",
      "el log-odds de hospitalización."
    )
  }
)
caption_efectos <- paste(strwrap(caption_efectos, width = 130), collapse = "\n")

p_efectos <- ggplot2::ggplot(
  datos_efectos,
  ggplot2::aes(
    x = estimate,
    y = stats::reorder(etiqueta, estimate)
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = "grey55",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = conf.low, xmax = conf.high),
    orientation = "y",
    width = 0.15,
    colour = "#2166AC",
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(colour = "#2166AC", size = 2.0) +
  ggplot2::labs(
    title = titulo_efectos,
    subtitle = subtitulo_efectos,
    x = efectos_logit$etiqueta_x,
    y = NULL,
    caption = caption_efectos
  ) +
  ggplot2::theme_minimal(base_size = 9.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_png_efectos,
  plot = p_efectos,
  width = 9.2,
  height = 6.8,
  dpi = 200,
  bg = "white"
)

# --- 8. Métricas del scorecard a CSV -------------------------------

metricas_scorecard <- tibble::tibble(
  modelo = c("Logit", "Probit"),
  auc = c(auc_logit, auc_probit),
  ks = c(ks_logit, ks_probit),
  gini = c(gini_logit, gini_probit),
  n_entrena = nrow(entrena),
  n_prueba = nrow(prueba),
  prevalencia_prueba = mean(prueba$hosp),
  aic = c(stats::AIC(logit), stats::AIC(probit)),
  fuente_roc = c(roc_logit$fuente, roc_probit$fuente)
)

readr::write_csv(
  metricas_scorecard,
  ruta_metricas_scorecard,
  na = "NA"
)

# --- 9. Panel de alerta temprana: preparación ----------------------

panel_crudo <- readRDS(ruta_panel)

if (!inherits(panel_crudo, "data.frame")) {
  stop(
    "El panel leído no es un data.frame: ", ruta_panel, ".",
    call. = FALSE
  )
}

requeridas_panel <- c("region_ema", "mes", "casos_dengue", "M_dengue")
faltan_panel <- setdiff(requeridas_panel, names(panel_crudo))

if (length(faltan_panel) > 0L) {
  stop(
    "Faltan columnas en el panel: ",
    paste(faltan_panel, collapse = ", "), ". Columnas disponibles: ",
    paste(names(panel_crudo), collapse = ", "), ".",
    call. = FALSE
  )
}

# La columna de lluvia puede llamarse lluvia o lluvia_mm. Se prefiere
# lluvia si existe; en D_13 la columna es lluvia_mm.
candidatas_lluvia <- c("lluvia", "lluvia_mm")
disponibles_lluvia <- candidatas_lluvia[
  candidatas_lluvia %in% names(panel_crudo)
]

if (length(disponibles_lluvia) == 0L) {
  stop(
    "No existe una columna de lluvia (lluvia o lluvia_mm) en el panel. ",
    "Columnas disponibles: ", paste(names(panel_crudo), collapse = ", "),
    ".",
    call. = FALSE
  )
}

col_lluvia <- disponibles_lluvia[1L]
message("Columna de lluvia detectada en el panel: ", col_lluvia)

panel_crudo$lluvia <- as.numeric(panel_crudo[[col_lluvia]])

panel <- panel_crudo |>
  dplyr::transmute(
    region_ema = as.character(region_ema),
    mes = as.Date(mes),
    casos_dengue = as.numeric(casos_dengue),
    M_dengue = as.numeric(M_dengue),
    lluvia = as.numeric(lluvia)
  ) |>
  dplyr::arrange(region_ema, mes)

if (anyNA(panel$region_ema) || anyNA(panel$mes)) {
  stop("Hay regiones o fechas faltantes en el panel.", call. = FALSE)
}
if (anyNA(panel$casos_dengue) || any(panel$casos_dengue < 0)) {
  stop(
    "casos_dengue tiene valores faltantes o negativos.",
    call. = FALSE
  )
}
if (anyNA(panel$M_dengue)) {
  stop("M_dengue tiene valores faltantes.", call. = FALSE)
}
if (anyNA(panel$lluvia)) {
  stop("La columna de lluvia tiene valores faltantes.", call. = FALSE)
}
if (anyDuplicated(panel[, c("region_ema", "mes")])) {
  stop("Hay pares región-mes duplicados en el panel.", call. = FALSE)
}

regiones <- sort(unique(panel$region_ema))
meses <- sort(unique(panel$mes))

if (length(regiones) != 12L) {
  warning(
    "Se esperaban 12 regiones EMA y hay ", length(regiones), ": ",
    paste(regiones, collapse = ", "), ".",
    call. = FALSE
  )
}

conteo_meses <- table(panel$region_ema)
if (length(unique(as.integer(conteo_meses))) != 1L) {
  warning(
    "El panel no está balanceado: las regiones tienen distinto número ",
    "de meses. El rezago se calcula igual dentro de cada región.",
    call. = FALSE
  )
}

# Alerta y rezago dentro de cada región. El arrange() anterior garantiza
# que lag() tome el mes inmediatamente anterior y no una fila de otra
# región.
panel <- panel |>
  dplyr::group_by(region_ema) |>
  dplyr::mutate(
    alerta = ifelse(
      casos_dengue > as.numeric(
        stats::quantile(
          casos_dengue,
          probs = 0.90,
          na.rm = TRUE
        )
      ),
      1L,
      0L
    ),
    lag_casos = dplyr::lag(casos_dengue, n = 1L)
  ) |>
  dplyr::ungroup()

if (anyNA(panel$alerta)) {
  warning(
    "Hay valores NA en alerta; se excluyen del modelo.",
    call. = FALSE
  )
}

alertas_region <- panel |>
  dplyr::group_by(region_ema) |>
  dplyr::summarise(
    n_meses = dplyr::n(),
    n_alertas = sum(alerta, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(region_ema)

regiones_sin_alerta <- alertas_region$region_ema[
  alertas_region$n_alertas == 0L
]

datos_alerta <- panel |>
  dplyr::filter(
    !is.na(alerta),
    !is.na(lag_casos),
    !is.na(M_dengue),
    !is.na(lluvia)
  )

if (nrow(datos_alerta) == 0L) {
  stop(
    "El panel quedó sin filas completas para el probit de alerta.",
    call. = FALSE
  )
}
if (length(unique(datos_alerta$alerta)) < 2L) {
  stop(
    "La variable alerta tiene una sola clase; no se puede estimar el ",
    "probit.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat(" PANEL DE ALERTA TEMPRANA — PREPARACIÓN\n")
cat("============================================================\n")
cat(sprintf(
  "Panel: %d filas | %d regiones | %s a %s\n",
  nrow(panel), length(regiones),
  format(min(meses), "%Y-%m"), format(max(meses), "%Y-%m")
))
cat(sprintf(
  "Alerta = 1 si casos_dengue > p90 regional | alertas: %d (%.2f%%)\n",
  sum(panel$alerta, na.rm = TRUE),
  100 * mean(panel$alerta, na.rm = TRUE)
))
cat(sprintf(
  paste0(
    "Modelo: %d filas (se pierde el primer mes de cada región por el ",
    "rezago: 12 filas)\n"
  ),
  nrow(datos_alerta)
))

cat("\nAlertas por región:\n")
for (i in seq_len(nrow(alertas_region))) {
  cat(sprintf(
    "  %-34s %2d alertas en %2d meses\n",
    alertas_region$region_ema[i],
    alertas_region$n_alertas[i],
    alertas_region$n_meses[i]
  ))
}

if (length(regiones_sin_alerta) > 0L) {
  cat(sprintf(
    "\nRegiones sin alertas con la regla > p90: %s\n",
    paste(regiones_sin_alerta, collapse = ", ")
  ))
  if ("Bogotá D.C." %in% regiones_sin_alerta) {
    cat(
      "  Bogotá D.C. tiene casos_dengue = 0 en todos los meses del ",
      "panel; la regla no puede generar alertas allí.\n",
      sep = ""
    )
  }
}

# --- 10. Probit de alerta temprana --------------------------------

probit_alerta <- stats::glm(
  alerta ~ M_dengue + lag_casos + lluvia,
  family = stats::binomial("probit"),
  data = datos_alerta
)

if (!isTRUE(probit_alerta$converged)) {
  stop("El probit de alerta temprana no convergió.", call. = FALSE)
}

coef_alerta <- stats::coef(probit_alerta)
if (!all(c("M_dengue", "lag_casos", "lluvia") %in% names(coef_alerta))) {
  stop(
    "El probit de alerta no tiene los tres coeficientes esperados. ",
    "Coeficientes: ", paste(names(coef_alerta), collapse = ", "), ".",
    call. = FALSE
  )
}

prob_alerta <- as.numeric(stats::fitted(probit_alerta))

if (length(prob_alerta) != nrow(datos_alerta) ||
    anyNA(prob_alerta) ||
    any(!is.finite(prob_alerta))) {
  stop(
    "Las probabilidades predichas del probit de alerta no son válidas.",
    call. = FALSE
  )
}

datos_alerta$prob_alerta <- prob_alerta

roc_alerta <- calcular_roc(datos_alerta$alerta, prob_alerta)
auc_alerta <- roc_alerta$auc

if (auc_alerta < 0.5) {
  warning(
    "El AUC del probit de alerta quedó por debajo de 0,5.",
    call. = FALSE
  )
}

cat("\n--- Coeficientes del probit de alerta ---\n")
print(stats::coef(summary(probit_alerta)))
cat(sprintf(
  "\nAUC del probit de alerta temprana: %.4f\n",
  auc_alerta
))
cat(sprintf(
  "Fuente de la ROC: %s | n = %d | prevalencia = %.4f\n",
  if (identical(roc_alerta$fuente, "pROC")) {
    "pROC::roc()"
  } else {
    "implementación interna"
  },
  nrow(datos_alerta),
  mean(datos_alerta$alerta)
))

# --- 11. Gráfica del probit de alerta -----------------------------

resumen_mes <- datos_alerta |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    prob_media = mean(prob_alerta),
    prob_p10 = as.numeric(
      stats::quantile(prob_alerta, probs = 0.10, na.rm = TRUE)
    ),
    prob_p90 = as.numeric(
      stats::quantile(prob_alerta, probs = 0.90, na.rm = TRUE)
    ),
    n_alertas = sum(alerta),
    .groups = "drop"
  ) |>
  dplyr::arrange(mes)

p_miedo <- ggplot2::ggplot(
  datos_alerta,
  ggplot2::aes(
    x = M_dengue,
    y = prob_alerta,
    colour = factor(alerta)
  )
) +
  ggplot2::geom_hline(
    yintercept = mean(datos_alerta$alerta),
    colour = "grey55",
    linewidth = 0.30,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(alpha = 0.65, size = 1.5) +
  ggplot2::scale_colour_manual(
    values = c("0" = "#2166AC", "1" = "#B2182B"),
    labels = c("No alerta", "Alerta"),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Probabilidad de alerta vs. M_dengue",
    subtitle = paste0(
      "Cada punto es una región-mes; la línea punteada marca la ",
      "prevalencia de alerta."
    ),
    x = "M_dengue (miedo excedente)",
    y = "Probabilidad predicha de alerta"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "left",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title.position = "plot"
  )

p_tiempo <- ggplot2::ggplot(
  resumen_mes,
  ggplot2::aes(x = mes, y = prob_media)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = prob_p10, ymax = prob_p90),
    fill = "#2166AC",
    alpha = 0.15,
    colour = NA
  ) +
  ggplot2::geom_line(colour = "#2166AC", linewidth = 0.80) +
  ggplot2::geom_point(
    data = dplyr::filter(resumen_mes, n_alertas > 0L),
    ggplot2::aes(x = mes, y = prob_media),
    inherit.aes = FALSE,
    colour = "#B2182B",
    size = 2.0
  ) +
  ggplot2::scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y-%m"
  ) +
  ggplot2::labs(
    title = "Probabilidad de alerta en el tiempo",
    subtitle = paste0(
      "Línea: promedio de las 12 regiones; banda: percentiles 10-90. ",
      "Puntos rojos: meses con al menos una alerta observada."
    ),
    x = "Mes",
    y = "Probabilidad predicha (media regional)"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

caption_alerta <- sprintf(
  paste0(
    "Probit: alerta ~ M_dengue + lag_casos + lluvia | %d regiones, ",
    "%d filas tras el rezago | alertas = %d (%.2f%%) | AUC = %.4f. ",
    "La alerta se define por región como casos_dengue > p90; ",
    "el rezago se calcula dentro de cada región."
  ),
  length(regiones),
  nrow(datos_alerta),
  sum(datos_alerta$alerta),
  100 * mean(datos_alerta$alerta),
  auc_alerta
)
caption_alerta <- paste(strwrap(caption_alerta, width = 150), collapse = "\n")

nota_alerta <- grid::textGrob(
  caption_alerta,
  gp = grid::gpar(fontsize = 7, col = "grey35"),
  hjust = 0,
  x = 0.01
)

# ggplotGrob() necesita un dispositivo gráfico. En Rscript se abriría
# Rplots.pdf; se usa un dispositivo PDF nulo para no dejar archivos.
grDevices::pdf(NULL)
grafico_alerta <- gridExtra::arrangeGrob(
  ggplot2::ggplotGrob(p_miedo),
  ggplot2::ggplotGrob(p_tiempo),
  ncol = 2,
  bottom = nota_alerta
)
invisible(grDevices::dev.off())

ggplot2::ggsave(
  filename = ruta_png_alerta,
  plot = grafico_alerta,
  width = 12,
  height = 5.8,
  dpi = 200,
  bg = "white"
)

# --- 12. Guardado de métricas y predicciones del panel ------------

metricas_alerta <- tibble::tibble(
  modelo = "Probit alerta temprana",
  auc = auc_alerta,
  n_obs = nrow(datos_alerta),
  n_regiones = length(unique(datos_alerta$region_ema)),
  n_alertas = sum(datos_alerta$alerta),
  prevalencia = mean(datos_alerta$alerta),
  aic = stats::AIC(probit_alerta),
  fuente_roc = roc_alerta$fuente
)

readr::write_csv(
  metricas_alerta,
  ruta_metricas_alerta,
  na = "NA"
)

predicciones_alerta <- datos_alerta |>
  dplyr::transmute(
    region_ema = region_ema,
    mes = mes,
    casos_dengue = casos_dengue,
    lag_casos = lag_casos,
    M_dengue = round(M_dengue, 6),
    lluvia = round(lluvia, 2),
    alerta = alerta,
    prob_alerta = round(prob_alerta, 6)
  ) |>
  dplyr::arrange(region_ema, mes)

readr::write_csv(
  predicciones_alerta,
  ruta_predicciones_alerta,
  na = "NA"
)

# --- 13. Validación de salidas -------------------------------------

rutas_salida <- c(
  ruta_metricas_scorecard,
  ruta_metricas_alerta,
  ruta_predicciones_alerta,
  ruta_png_roc,
  ruta_png_efectos,
  ruta_png_alerta
)

for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 14. Resumen en consola ----------------------------------------

cat("\n============================================================\n")
cat(" 21 — TABLA 4: SCORECARD LOGIT/PROBIT Y ALERTA TEMPRANA\n")
cat("============================================================\n")
cat(sprintf(
  "Scorecard: %s casos | %s entrena | %s prueba\n",
  format(nrow(datos_modelo), big.mark = ".", decimal.mark = ","),
  format(nrow(entrena), big.mark = ".", decimal.mark = ","),
  format(nrow(prueba), big.mark = ".", decimal.mark = ",")
))
cat(sprintf(
  "Scorecard Logit : AUC = %.4f | KS = %.4f | Gini = %.4f\n",
  auc_logit, ks_logit, gini_logit
))
cat(sprintf(
  "Scorecard Probit: AUC = %.4f | KS = %.4f | Gini = %.4f\n",
  auc_probit, ks_probit, gini_probit
))
cat(sprintf(
  "Probit de alerta temprana: AUC = %.4f\n",
  auc_alerta
))

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

cat("\nNotas:\n")
if (!tiene_pROC) {
  cat(
    "  - pROC no está instalado: la ROC y el AUC se calcularon con una\n",
    "    implementación interna equivalente (AUC por rangos).\n",
    sep = ""
  )
}
if (!tiene_marginaleffects) {
  cat(
    "  - marginaleffects no está instalado: los efectos marginales se\n",
    "    sustituyeron por los coeficientes de broom::tidy().\n",
    sep = ""
  )
}
if (length(regiones_sin_alerta) > 0L) {
  cat(
    "  - Regiones sin alertas con la regla > p90: ",
    paste(regiones_sin_alerta, collapse = ", "), "\n",
    sep = ""
  )
}
cat(
  "  - El AUC del scorecard es moderado: la especificación aporta\n",
  "    señal, pero no agota los determinantes de la hospitalización.\n",
  sep = ""
)
cat(
  "  - El AUC del panel es in-sample; no es una validación temporal.\n",
  sep = ""
)

message("\n✅ Script 21 finalizado.")
