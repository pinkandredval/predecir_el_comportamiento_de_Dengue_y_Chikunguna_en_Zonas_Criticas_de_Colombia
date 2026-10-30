# ============================================================
# 16_var_its_granger.R
# Vía de contraste del plan (§4.5): VAR nacional, FEVD,
# causalidad de Granger y serie interrumpida (ITS) de chikungunya
# ============================================================
#
# OBJETIVO
#   Contrastar, con un VAR nacional y una serie interrumpida, la
#   hipótesis del plan: el miedo (M) desplaza la ocupación hotelera (Y)
#   más allá de lo que explican la incidencia (I) y el clima.
#
#   (1) VAR(p) nacional en niveles:
#
#         [I_t, M_t, Y_t]' = c + sum_{k=1..p} A_k [I, M, Y]'_{t-k}
#                             + B C_t + e_t
#
#       con C_t = (lluvia_t, temp_t)' como exógeno contemporáneo.
#       Orden de Cholesky: I -> M -> Y (incidencia, miedo, economía).
#       Ese orden es el de las columnas de y en vars::VAR(): con él,
#       la FEVD y la IRF ortogonalizadas suponen que Y no afecta
#       contemporáneamente a I ni a M, M no afecta a I, y sí al revés.
#
#   (2) FEVD a 6 y 12 meses: fracción de la varianza del error de
#       pronóstico de la ocupación que explica un choque de miedo.
#
#   (3) Causalidad de Granger: M -> Y e I -> M.
#
#   (4) ITS de chikungunya en 2025:
#
#         Y_t = a + b M_chik_t + c D_t + d (M_chik_t x D_t)
#               + f1 lluvia_t + f2 temp_t + u_t
#
#       D_t = 1{t >= 2025-01-01}. El coeficiente d mide el cambio de
#       la pendiente Y-M_chik después del quiebre. En 2025 la
#       chikungunya no registra un brote local (el miedo es
#       internacional), de modo que un d significativo es consistente
#       con un canal de miedo puro (miedo sin enfermedad local).
#
# INSUMOS
#   data/D_13/panel_ema_den_chik.rds   (576 filas: 12 regiones x 48 meses)
#
# SALIDAS
#   data/D_16/modelo_var.rds
#   data/D_16/granger_resultados.csv
#   data/D_16/fevd_its_resultados.csv
#   graficas/G_16/fevd_dengue.png
#   graficas/G_16/irf_dengue.png
#   graficas/G_16/its_chikungunya.png
#
# NOTA SOBRE LA AGREGACIÓN NACIONAL
#   Los casos se SUMAN (son conteos) y el miedo, la ocupación (porcentaje),
#   la lluvia y la temperatura se PROMEDIAN. Sumar M o Y produciría series
#   en una escala distinta a la del panel y con una estacionalidad espuria.
#
# NOTA SOBRE EL NIVEL DE LAS SERIES
#   El plan especifica el VAR en niveles. El script 14 documenta el
#   orden de integración de las series; si alguna no es estacionaria,
#   la FEVD y la IRF deben leerse como descripción del sistema estimado
#   y no como impulso-respuesta estructurales. Las participaciones de
#   la FEVD son invariantes a la escala de las variables, de modo que
#   la comparación I vs. M vs. Y no depende de que los casos estén en
#   unidades o en miles.
#
# NOTA SOBRE vars::causality()
#   vars::causality(x, cause = "M_dengue") contrasta la exclusión
#   conjunta de los rezagos de M_dengue en TODAS las demás ecuaciones,
#   no solo en la de ocupación. Para responder la pregunta del plan
#   (¿M causa a Y? ¿I causa a M?) el script complementa ese contraste
#   con la prueba F de exclusión por ecuación, calculada sobre la misma
#   VAR (es la prueba de Granger estándar ecuación por ecuación). Ambas
#   quedan en granger_resultados.csv, identificadas por la columna
#   `prueba`.
#
# NOTA SOBRE EL USO DE ::
#   Todas las funciones no base se llaman con paquete::funcion() para
#   evitar el enmascaramiento entre dplyr, stats, vars y ggplot2. No se
#   adjuntan paquetes con library().
# ============================================================

# --- 0. Rutas y carpetas -------------------------------------------

dir.create(here::here("data", "D_16"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("graficas", "G_16"), recursive = TRUE, showWarnings = FALSE)

ruta_panel    <- here::here("data", "D_13", "panel_ema_den_chik.rds")
ruta_modelo   <- here::here("data", "D_16", "modelo_var.rds")
ruta_granger  <- here::here("data", "D_16", "granger_resultados.csv")
ruta_fevd_its <- here::here("data", "D_16", "fevd_its_resultados.csv")
ruta_fevd_png <- here::here("graficas", "G_16", "fevd_dengue.png")
ruta_irf_png  <- here::here("graficas", "G_16", "irf_dengue.png")
ruta_its_png  <- here::here("graficas", "G_16", "its_chikungunya.png")

# Constantes del modelo
rezago_maximo   <- 3L
horizonte_fevd  <- 12L
horizontes_clave <- c(6L, 12L)
fecha_quiebre   <- as.Date("2025-01-01")
semilla         <- 2025L

if (!requireNamespace("vars", quietly = TRUE)) {
  stop("Falta el paquete vars. Instálelo con renv::install(\"vars\").",
       call. = FALSE)
}

# --- 1. Carga y verificación previa del panel ----------------------

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel del script 13: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

panel <- readRDS(ruta_panel)

cat("\n============================================================\n")
cat(" VERIFICACIÓN PREVIA DEL PANEL\n")
cat("============================================================\n")
cat("names(panel):\n")
print(names(panel))
cat(sprintf("\nClase: %s | Filas: %d | Columnas: %d\n",
            paste(class(panel), collapse = "/"), nrow(panel), ncol(panel)))

# Los nombres del panel pueden variar entre corridas del script 13
# (por ejemplo, Y vs. ocupacion, lluvia_mm vs. lluvia, temp_c vs. temp).
# Se detectan aquí y se documenta el mapeo efectivo.
buscar_columna <- function(candidatos, panel, etiqueta) {
  presentes <- candidatos[candidatos %in% names(panel)]
  if (length(presentes) == 0L) {
    stop(
      "No se encontró la columna de ", etiqueta, ". ",
      "Candidatas: ", paste(candidatos, collapse = ", "), ". ",
      "Columnas disponibles: ", paste(names(panel), collapse = ", "), ".",
      call. = FALSE
    )
  }
  presentes[[1L]]
}

col_casos_dengue <- buscar_columna(c("casos_dengue"), panel, "casos de dengue")
col_casos_chik   <- buscar_columna(c("casos_chik"), panel, "casos de chikungunya")
col_m_dengue     <- buscar_columna(c("M_dengue"), panel, "miedo de dengue")
col_m_chik       <- buscar_columna(c("M_chik"), panel, "miedo de chikungunya")
col_ocupacion    <- buscar_columna(
  c("Y", "ocupacion", "ocupacion_hotelera"), panel, "ocupación hotelera"
)
col_lluvia <- buscar_columna(
  c("lluvia", "lluvia_mm", "precipitacion"), panel, "lluvia"
)
col_temp <- buscar_columna(
  c("temp", "temp_c", "temperatura"), panel, "temperatura"
)

if (!"mes" %in% names(panel)) {
  stop("El panel no tiene columna `mes`.", call. = FALSE)
}

cat("\nMapeo de columnas usado:\n")
cat(sprintf("  %-14s <- %s\n", "casos_dengue", col_casos_dengue))
cat(sprintf("  %-14s <- %s\n", "casos_chik", col_casos_chik))
cat(sprintf("  %-14s <- %s\n", "M_dengue", col_m_dengue))
cat(sprintf("  %-14s <- %s\n", "M_chik", col_m_chik))
cat(sprintf("  %-14s <- %s\n", "ocupacion", col_ocupacion))
cat(sprintf("  %-14s <- %s\n", "lluvia", col_lluvia))
cat(sprintf("  %-14s <- %s\n", "temp", col_temp))

# --- 2. Agregación nacional ----------------------------------------

# Se construye una tabla con nombres canónicos y luego se agrega por
# mes: suma para los conteos, media para índices y porcentajes.
nacional_base <- panel |>
  dplyr::transmute(
    mes = as.Date(.data[["mes"]]),
    casos_dengue = as.numeric(.data[[col_casos_dengue]]),
    casos_chik = as.numeric(.data[[col_casos_chik]]),
    M_dengue = as.numeric(.data[[col_m_dengue]]),
    M_chik = as.numeric(.data[[col_m_chik]]),
    ocupacion = as.numeric(.data[[col_ocupacion]]),
    lluvia = as.numeric(.data[[col_lluvia]]),
    temp = as.numeric(.data[[col_temp]])
  ) |>
  dplyr::mutate(anio = as.integer(format(mes, "%Y"))) |>
  dplyr::filter(anio >= 2022L, anio <= 2025L)

# Número de regiones por mes: se usa para verificar que la agregación
# promedia el mismo número de unidades en todos los meses.
regiones_por_mes <- nacional_base |>
  dplyr::count(mes, name = "n_regiones")

nacional <- nacional_base |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    casos_dengue = sum(casos_dengue, na.rm = TRUE),
    casos_chik = sum(casos_chik, na.rm = TRUE),
    M_dengue = mean(M_dengue, na.rm = TRUE),
    M_chik = mean(M_chik, na.rm = TRUE),
    ocupacion = mean(ocupacion, na.rm = TRUE),
    lluvia = mean(lluvia, na.rm = TRUE),
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(anio = as.integer(format(mes, "%Y"))) |>
  dplyr::arrange(mes)

# Validaciones de la serie nacional
if (nrow(nacional) != 48L) {
  stop(
    "La serie nacional debe tener 48 meses (2022-01 a 2025-12) y tiene ",
    nrow(nacional), ".", call. = FALSE
  )
}

if (anyNA(nacional)) {
  stop("Hay valores faltantes en la serie nacional.", call. = FALSE)
}

if (any(diff(nacional$mes) > 32)) {
  stop("Hay meses faltantes en la serie nacional.", call. = FALSE)
}

if (dplyr::n_distinct(regiones_por_mes$n_regiones) != 1L) {
  stop(
    "Las regiones no están completas en todos los meses: ",
    paste(unique(regiones_por_mes$n_regiones), collapse = ", "), ".",
    call. = FALSE
  )
}

n_regiones <- unique(regiones_por_mes$n_regiones)
periodo <- sprintf(
  "%s a %s", format(min(nacional$mes), "%Y-%m"),
  format(max(nacional$mes), "%Y-%m")
)

cat("\n--- AGREGACIÓN NACIONAL ---\n")
cat(sprintf(
  "Meses: %d (%s) | Regiones por mes: %d\n",
  nrow(nacional), periodo, n_regiones
))
cat(sprintf(
  "Casos de dengue: %s en total (máx. mensual %s en %s)\n",
  format(sum(nacional$casos_dengue), big.mark = ".", decimal.mark = ","),
  format(max(nacional$casos_dengue), big.mark = ".", decimal.mark = ","),
  format(nacional$mes[which.max(nacional$casos_dengue)], "%Y-%m")
))
cat(sprintf(
  "Casos de chikungunya: %s en total; en 2025: %s\n",
  format(sum(nacional$casos_chik), big.mark = ".", decimal.mark = ","),
  format(sum(nacional$casos_chik[nacional$anio == 2025L]),
         big.mark = ".", decimal.mark = ",")
))
cat(sprintf(
  "Ocupación media: %.2f%% (min %.2f, max %.2f)\n",
  mean(nacional$ocupacion), min(nacional$ocupacion), max(nacional$ocupacion)
))
cat(sprintf(
  "M_dengue media: %.3f | M_chik media: %.3f\n",
  mean(nacional$M_dengue), mean(nacional$M_chik)
))

# --- 3. Modelo VAR nacional ----------------------------------------

# Orden de Cholesky: I -> M -> Y. El orden de estas columnas es el que
# vars::VAR() usa para la descomposición de Cholesky de la matriz de
# covarianzas de los residuos; por eso no se reordenan después.
variables_var <- c("casos_dengue", "M_dengue", "ocupacion")
orden_cholesky <- c("I: casos_dengue", "M: M_dengue", "Y: ocupacion")

serie_var <- as.matrix(nacional[, variables_var])
colnames(serie_var) <- variables_var

exogen_var <- cbind(lluvia = nacional$lluvia, temp = nacional$temp)

if (!identical(colnames(serie_var), variables_var)) {
  stop("El orden de las columnas del VAR no es I -> M -> Y.", call. = FALSE)
}

cat("\n--- SELECCIÓN DE REZAGOS ---\n")

seleccion_rezagos <- vars::VARselect(
  serie_var,
  lag.max = rezago_maximo,
  type = "const",
  exogen = exogen_var
)

cat("Criterios de selección (vars::VARselect, lag.max = 3):\n")
print(seleccion_rezagos$selection)

criterios <- seleccion_rezagos$selection
p_optimo <- as.integer(criterios[["AIC(n)"]])

if (is.na(p_optimo)) {
  p_optimo <- as.integer(criterios[["SC(n)"]])
}
if (is.na(p_optimo) || p_optimo < 1L) {
  p_optimo <- 1L
  message(
    "VARselect no devolvió un rezago válido (AIC/SC); se usa p = 1 para ",
    "poder contrastar la causalidad de Granger."
  )
}

cat(sprintf(
  "\nRezago elegido por AIC: p = %d (HQ = %s, SC = %s, FPE = %s)\n",
  p_optimo,
  criterios[["HQ(n)"]], criterios[["SC(n)"]], criterios[["FPE(n)"]]
))

modelo_var <- vars::VAR(
  serie_var,
  p = p_optimo,
  type = "const",
  exogen = exogen_var
)

if (!inherits(modelo_var, "varest")) {
  stop("vars::VAR() no devolvió un objeto varest.", call. = FALSE)
}
if (modelo_var$K != 3L) {
  stop("La VAR debe tener 3 variables (I, M, Y).", call. = FALSE)
}

cat(sprintf(
  "VAR estimada: %d observaciones efectivas | %d ecuaciones | p = %d | exógenas: lluvia, temp\n",
  modelo_var$obs, modelo_var$K, modelo_var$p
))

# --- 4. FEVD a 6 y 12 meses ----------------------------------------

fevd_var <- vars::fevd(modelo_var, n.ahead = horizonte_fevd)

if (!"ocupacion" %in% names(fevd_var)) {
  stop(
    "La FEVD no contiene la respuesta `ocupacion`. Nombres: ",
    paste(names(fevd_var), collapse = ", "), ".",
    call. = FALSE
  )
}

fevd_ocupacion <- fevd_var[["ocupacion"]]
if (!is.matrix(fevd_ocupacion) || !"M_dengue" %in% colnames(fevd_ocupacion)) {
  stop("La FEVD de ocupación no tiene la columna M_dengue.", call. = FALSE)
}

# Fracción (0-1) de la varianza de ocupación explicada por el choque de
# miedo, en los horizontes clave. Es la lectura directa de la FEVD
# ortogonalizada: la columna del choque y la fila del horizonte.
fevd_miedo <- tibble::tibble(
  horizonte = horizontes_clave,
  participacion = as.numeric(fevd_ocupacion[horizontes_clave, "M_dengue"]),
  participacion_pct = 100 * participacion
)

# Serie completa de la FEVD de ocupación, para trazabilidad.
fevd_ocupacion_df <- tibble::tibble(
  horizonte = seq_len(nrow(fevd_ocupacion)),
  casos_dengue = as.numeric(fevd_ocupacion[, "casos_dengue"]),
  M_dengue = as.numeric(fevd_ocupacion[, "M_dengue"]),
  ocupacion = as.numeric(fevd_ocupacion[, "ocupacion"])
)

# Comprobación: las participaciones de cada fila deben sumar 1.
sumas_fevd <- rowSums(fevd_ocupacion)
if (any(abs(sumas_fevd - 1) > 1e-6)) {
  warning("Las participaciones de la FEVD no suman 1 en algún horizonte.",
          call. = FALSE)
}

cat("\n--- FEVD DE LA OCUPACIÓN (n.ahead = 12) ---\n")
print(round(fevd_ocupacion_df, 4))

cat(sprintf(
  "\nVarianza de ocupación explicada por M_dengue: %.2f%% a %d meses | %.2f%% a %d meses\n",
  fevd_miedo$participacion_pct[1], fevd_miedo$horizonte[1],
  fevd_miedo$participacion_pct[2], fevd_miedo$horizonte[2]
))

# --- 5. Causalidad de Granger --------------------------------------

# Prueba F de exclusión por ecuación: es la prueba de Granger estándar
# (los rezagos de `causa` no entran en la ecuación de `ecuacion`). Se
# calcula sobre la misma VAR estimada, comparando la ecuación sin
# restringir con la ecuación que excluye los p rezagos de la causa.
probar_granger_ecuacion <- function(modelo, ecuacion, causa) {
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop(
      "Se requiere lmtest para el contraste de Granger por ecuación. ",
      "Instálelo con renv::install(\"lmtest\").",
      call. = FALSE
    )
  }

  ajuste <- modelo$varresult[[ecuacion]]
  if (is.null(ajuste)) {
    stop("La VAR no tiene ecuación para ", ecuacion, ".", call. = FALSE)
  }

  marco <- stats::model.frame(ajuste)
  terminos <- names(marco)[-1L]
  rezagos_causa <- paste0(causa, ".l", seq_len(modelo$p))

  faltan <- setdiff(rezagos_causa, terminos)
  if (length(faltan) > 0L) {
    stop(
      "La ecuación ", ecuacion, " no contiene los rezagos: ",
      paste(faltan, collapse = ", "), ".",
      call. = FALSE
    )
  }

  conservar <- setdiff(terminos, rezagos_causa)
  fml_restringida <- stats::reformulate(
    conservar,
    response = names(marco)[1L],
    intercept = FALSE
  )
  restringido <- stats::lm(fml_restringida, data = marco)

  tabla <- as.data.frame(lmtest::waldtest(ajuste, restringido, test = "F"))

  tibble::tibble(
    prueba = "Granger por ecuación (Wald F)",
    causa = causa,
    efecto = ecuacion,
    estadistico = as.numeric(tabla[2L, "F"]),
    gl1 = abs(as.numeric(tabla[2L, "Df"])),
    gl2 = as.numeric(tabla[1L, "Res.Df"]),
    p_valor = as.numeric(tabla[2L, "Pr(>F)"])
  )
}

# Contraste conjunto que ofrece vars::causality(): exclusión de los
# rezagos de `causa` en todas las demás ecuaciones a la vez.
probar_granger_sistema <- function(modelo, causa) {
  objeto <- vars::causality(modelo, cause = causa)$Granger

  # htest$parameter guarda df1 (y df2 si la prueba es F); se normaliza a
  # escalares para que sprintf() no recicle la línea impresa.
  gl <- unname(objeto$parameter)

  tibble::tibble(
    prueba = "Granger sistema (vars::causality)",
    causa = causa,
    efecto = paste(setdiff(colnames(modelo$y), causa), collapse = " + "),
    estadistico = unname(objeto$statistic)[1L],
    gl1 = gl[1L],
    gl2 = if (length(gl) >= 2L) gl[2L] else NA_real_,
    p_valor = unname(objeto$p.value)[1L]
  )
}

granger_m_a_y <- probar_granger_ecuacion(modelo_var, "ocupacion", "M_dengue")
granger_i_a_m <- probar_granger_ecuacion(modelo_var, "M_dengue", "casos_dengue")

granger_sistema_m <- probar_granger_sistema(modelo_var, "M_dengue")
granger_sistema_i <- probar_granger_sistema(modelo_var, "casos_dengue")

granger_resultados <- dplyr::bind_rows(
  granger_m_a_y,
  granger_i_a_m,
  granger_sistema_m,
  granger_sistema_i
)

cat("\n--- CAUSALIDAD DE GRANGER ---\n")
cat(sprintf(
  "M_dengue -> ocupación        : F(%d, %d) = %.3f, p = %.4f\n",
  granger_m_a_y$gl1, granger_m_a_y$gl2,
  granger_m_a_y$estadistico, granger_m_a_y$p_valor
))
cat(sprintf(
  "casos_dengue -> M_dengue     : F(%d, %d) = %.3f, p = %.4f\n",
  granger_i_a_m$gl1, granger_i_a_m$gl2,
  granger_i_a_m$estadistico, granger_i_a_m$p_valor
))
cat(sprintf(
  "Sistema, M_dengue -> {%s}: F = %.3f, p = %.4f\n",
  granger_sistema_m$efecto, granger_sistema_m$estadistico,
  granger_sistema_m$p_valor
))
cat(sprintf(
  "Sistema, casos_dengue -> {%s}: F = %.3f, p = %.4f\n",
  granger_sistema_i$efecto, granger_sistema_i$estadistico,
  granger_sistema_i$p_valor
))

# --- 6. Gráfico de la FEVD -----------------------------------------

# Etiquetas legibles de las variables del VAR.
etiquetas_var <- c(
  casos_dengue = "I: casos de dengue",
  M_dengue = "M: miedo dengue",
  ocupacion = "Y: ocupación hotelera"
)
etiquetas_respuesta <- c(
  casos_dengue = "Respuesta: casos de dengue",
  M_dengue = "Respuesta: miedo dengue",
  ocupacion = "Respuesta: ocupación hotelera"
)

colores_var <- c(
  "I: casos de dengue" = "#08519C",
  "M: miedo dengue" = "#B2182B",
  "Y: ocupación hotelera" = "#1B7837"
)

fevd_largo <- dplyr::bind_rows(lapply(names(fevd_var), function(respuesta) {
  matriz <- as.data.frame(fevd_var[[respuesta]])
  matriz$horizonte <- seq_len(nrow(matriz))
  matriz$respuesta <- respuesta
  tidyr::pivot_longer(
    matriz,
    cols = -c(horizonte, respuesta),
    names_to = "choque",
    values_to = "participacion"
  )
})) |>
  dplyr::mutate(
    respuesta = factor(respuesta, levels = names(etiquetas_respuesta),
                       labels = unname(etiquetas_respuesta)),
    choque = factor(choque, levels = names(etiquetas_var),
                    labels = unname(etiquetas_var))
  )

p_fevd <- ggplot2::ggplot(
  fevd_largo,
  ggplot2::aes(x = horizonte, y = participacion, fill = choque)
) +
  ggplot2::geom_col(width = 0.82, colour = NA) +
  ggplot2::geom_vline(
    xintercept = horizontes_clave[1], linetype = "dashed",
    colour = "grey20", linewidth = 0.3
  ) +
  ggplot2::facet_wrap(~respuesta, ncol = 3) +
  ggplot2::scale_fill_manual(values = colores_var, name = NULL) +
  ggplot2::scale_x_continuous(
    breaks = c(1, 3, 6, 9, 12), limits = c(0.3, 12.7)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Descomposición de la varianza del error de pronóstico (FEVD)",
    subtitle = paste0(
      "VAR nacional con orden de Cholesky I \u2192 M \u2192 Y; ",
      "clima exógeno (lluvia y temperatura). La línea punteada marca 6 meses."
    ),
    x = "Horizonte (meses)",
    y = "Participación en la varianza",
    caption = sprintf(
      paste0(
        "Un choque de M_dengue explica el %.2f%% de la varianza de ",
        "ocupación a %d meses y el %.2f%% a %d meses. ",
        "Las barras de cada horizonte suman 100%%."
      ),
      fevd_miedo$participacion_pct[1], fevd_miedo$horizonte[1],
      fevd_miedo$participacion_pct[2], fevd_miedo$horizonte[2]
    )
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 8.5),
    axis.text.x = ggplot2::element_text(size = 7.5),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_fevd_png, plot = p_fevd,
  width = 10.5, height = 4.8, dpi = 200, bg = "white"
)

# --- 7. Gráfico de la IRF ------------------------------------------

# IRF ortogonalizadas (Cholesky) de los dos canales de interés:
# I -> Y (enfermedad) y M -> Y (miedo), con sus efectos sobre M.
impulsos <- c("casos_dengue", "M_dengue")
respuestas <- c("M_dengue", "ocupacion")

irf_var <- tryCatch(
  vars::irf(
    modelo_var,
    impulse = impulsos,
    response = respuestas,
    n.ahead = horizonte_fevd,
    ortho = TRUE,
    boot = TRUE,
    runs = 500,
    ci = 0.95,
    seed = semilla,
    progress = FALSE
  ),
  error = function(e) {
    message(
      "No se pudo bootstrapear la IRF (", conditionMessage(e),
      "); se grafica la IRF puntual sin bandas."
    )
    vars::irf(
      modelo_var,
      impulse = impulsos,
      response = respuestas,
      n.ahead = horizonte_fevd,
      ortho = TRUE,
      boot = FALSE,
      progress = FALSE
    )
  }
)

hay_bandas <- !is.null(irf_var$Lower) && !is.null(irf_var$Upper)

construir_irf_df <- function(objeto, con_bandas) {
  filas <- list()
  for (impulso in names(objeto$irf)) {
    central <- objeto$irf[[impulso]]
    inferior <- if (con_bandas) objeto$Lower[[impulso]] else NULL
    superior <- if (con_bandas) objeto$Upper[[impulso]] else NULL
    horizontes <- seq_len(nrow(central)) - 1L

    # Los vectores se calculan FUERA de tibble(): dentro de tibble() una
    # columna ya creada (p. ej. `respuesta`) tapa a la variable del mismo
    # nombre y se recicla a la longitud del horizonte.
    for (resp in colnames(central)) {
      valores_irf <- as.numeric(central[, resp])
      valores_inferior <- if (con_bandas) as.numeric(inferior[, resp]) else NA_real_
      valores_superior <- if (con_bandas) as.numeric(superior[, resp]) else NA_real_

      filas[[length(filas) + 1L]] <- tibble::tibble(
        choque = impulso,
        respuesta = resp,
        horizonte = horizontes,
        irf = valores_irf,
        inferior = valores_inferior,
        superior = valores_superior
      )
    }
  }
  dplyr::bind_rows(filas)
}

irf_largo <- construir_irf_df(irf_var, hay_bandas) |>
  dplyr::mutate(
    choque = factor(choque, levels = names(etiquetas_var),
                    labels = unname(etiquetas_var)),
    respuesta = factor(respuesta, levels = c("M_dengue", "ocupacion"),
                       labels = c("Respuesta: miedo (M)",
                                  "Respuesta: ocupación (Y)"))
  )

p_irf <- ggplot2::ggplot(
  irf_largo,
  ggplot2::aes(x = horizonte, y = irf)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.3) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = inferior, ymax = superior),
    fill = "#2166AC", alpha = 0.16, colour = NA, na.rm = TRUE
  ) +
  ggplot2::geom_line(colour = "#08519C", linewidth = 0.65, na.rm = TRUE) +
  ggplot2::geom_point(colour = "#08519C", size = 1.1, na.rm = TRUE) +
  ggplot2::facet_grid(respuesta ~ choque, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
  ggplot2::labs(
    title = "Funciones impulso-respuesta ortogonalizadas",
    subtitle = paste0(
      "VAR nacional, orden de Cholesky I \u2192 M \u2192 Y. ",
      if (hay_bandas) {
        "Bandas al 95% por bootstrap (500 réplicas)."
      } else {
        "IRF puntual (bootstrap no disponible)."
      }
    ),
    x = "Horizonte (meses)",
    y = "Respuesta (unidades de la variable)",
    caption = paste0(
      "Choques de un error estándar. La fila superior muestra el efecto ",
      "sobre el miedo; la inferior, el efecto sobre la ocupación hotelera."
    )
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_irf_png, plot = p_irf,
  width = 9, height = 6.4, dpi = 200, bg = "white"
)

# --- 8. Serie interrumpida (ITS) de chikungunya --------------------

# Dummy de quiebre: 1 desde el 1 de enero de 2025.
nacional <- nacional |>
  dplyr::mutate(
    dummy_2025 = as.integer(mes >= fecha_quiebre)
  )

if (sum(nacional$dummy_2025) == 0L) {
  stop("La dummy de 2025 no tiene observaciones en 1.", call. = FALSE)
}

cat("\n--- SERIE INTERRUMPIDA: CHIKUNGUNYA 2025 ---\n")
cat(sprintf(
  "Meses antes del quiebre: %d | meses después: %d\n",
  sum(nacional$dummy_2025 == 0L), sum(nacional$dummy_2025 == 1L)
))
cat(sprintf(
  "Casos de chikungunya notificados en 2025: %s (2022-2024: %s)\n",
  format(sum(nacional$casos_chik[nacional$anio == 2025L]),
         big.mark = ".", decimal.mark = ","),
  format(sum(nacional$casos_chik[nacional$anio < 2025L]),
         big.mark = ".", decimal.mark = ",")
))
cat(sprintf(
  "M_chik medio: 2022-2024 = %.2f | 2025 = %.2f\n",
  mean(nacional$M_chik[nacional$anio < 2025L]),
  mean(nacional$M_chik[nacional$anio == 2025L])
))

its <- stats::lm(
  ocupacion ~ M_chik * dummy_2025 + lluvia + temp,
  data = nacional
)

its_resumen <- summary(its)

its_tabla <- tibble::tibble(
  termino = rownames(its_resumen$coefficients),
  estimacion = as.numeric(its_resumen$coefficients[, 1]),
  error_estandar = as.numeric(its_resumen$coefficients[, 2]),
  estadistico_t = as.numeric(its_resumen$coefficients[, 3]),
  p_valor = as.numeric(its_resumen$coefficients[, 4])
)

termino_interaccion <- "M_chik:dummy_2025"
if (!termino_interaccion %in% its_tabla$termino) {
  stop(
    "La ITS no estimó el término ", termino_interaccion, ". Términos: ",
    paste(its_tabla$termino, collapse = ", "), ".",
    call. = FALSE
  )
}

its_interaccion <- dplyr::filter(its_tabla, termino == termino_interaccion)

cat("\nCoeficientes de la ITS:\n")
print(
  dplyr::mutate(
    its_tabla,
    dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4))
  )
)
cat(sprintf(
  "\nInteracción %s: estimación = %.4f, EE = %.4f, t = %.3f, p = %.4f\n",
  termino_interaccion, its_interaccion$estimacion,
  its_interaccion$error_estandar, its_interaccion$estadistico_t,
  its_interaccion$p_valor
))

lectura_significacion <- if (its_interaccion$p_valor < 0.01) {
  "significativa al 1%"
} else if (its_interaccion$p_valor < 0.05) {
  "significativa al 5%"
} else if (its_interaccion$p_valor < 0.10) {
  "significativa solo al 10% (no al 5%)"
} else {
  "no significativa a los niveles habituales (1%, 5% y 10%)"
}
cat(sprintf("Lectura: la interacción es %s.\n", lectura_significacion))

# Ajuste y contrafactual sin quiebre (dummy_2025 = 0), para el gráfico.
nacional_sin_quiebre <- dplyr::mutate(nacional, dummy_2025 = 0L)

nacional <- nacional |>
  dplyr::mutate(
    its_ajuste = as.numeric(stats::fitted(its)),
    its_contrafactual = as.numeric(
      stats::predict(its, newdata = nacional_sin_quiebre)
    )
  )

# --- 9. Gráfico de la serie interrumpida ---------------------------

panel_y <- "Ocupación hotelera (%)"
panel_m <- "Miedo chikungunya (M_chik)"

its_largo <- dplyr::bind_rows(
  tibble::tibble(
    variable = panel_y,
    mes = nacional$mes,
    observado = nacional$ocupacion,
    ajuste = nacional$its_ajuste,
    contrafactual = nacional$its_contrafactual
  ),
  tibble::tibble(
    variable = panel_m,
    mes = nacional$mes,
    observado = nacional$M_chik,
    ajuste = NA_real_,
    contrafactual = NA_real_
  )
) |>
  dplyr::mutate(
    variable = factor(variable, levels = c(panel_y, panel_m))
  )

p_its <- ggplot2::ggplot(its_largo, ggplot2::aes(x = mes)) +
  ggplot2::annotate(
    "rect",
    xmin = fecha_quiebre, xmax = max(its_largo$mes),
    ymin = -Inf, ymax = Inf, fill = "grey90", colour = NA
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = observado), colour = "grey20", linewidth = 0.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = observado), colour = "grey20", size = 0.9
  ) +
  ggplot2::geom_line(
    data = dplyr::filter(its_largo, variable == panel_y),
    ggplot2::aes(y = ajuste), colour = "#08519C", linewidth = 0.75
  ) +
  ggplot2::geom_line(
    data = dplyr::filter(its_largo, variable == panel_y),
    ggplot2::aes(y = contrafactual),
    colour = "#B2182B", linetype = "dashed", linewidth = 0.7
  ) +
  ggplot2::geom_vline(
    xintercept = fecha_quiebre, linetype = "dashed", colour = "grey25",
    linewidth = 0.4
  ) +
  ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y") +
  ggplot2::scale_x_date(
    date_breaks = "6 months", date_labels = "%Y-%m",
    expand = ggplot2::expansion(mult = c(0.01, 0.03))
  ) +
  ggplot2::labs(
    title = "Serie interrumpida de chikungunya: quiebre de enero de 2025",
    subtitle = sprintf(
      paste0(
        "Y ~ M_chik * D_2025 + lluvia + temp. Interacción = %.3f ",
        "(EE %.3f; p = %.3f): %s."
      ),
      its_interaccion$estimacion, its_interaccion$error_estandar,
      its_interaccion$p_valor, lectura_significacion
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Negro: observado. Azul: ajuste del modelo ITS. ",
      "Rojo punteado: contrafactual sin quiebre (D_2025 = 0). ",
      "Franja gris: 2025. ", sprintf("N = %d.", stats::nobs(its))
    )
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    axis.text.x = ggplot2::element_text(size = 7.5),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_its_png, plot = p_its,
  width = 10, height = 6.6, dpi = 200, bg = "white"
)

# --- 10. Guardado de resultados ------------------------------------

# El objeto guardado es la varest, con metadatos en atributos para que
# el RDS sea autónomo (orden, rezago y período no dependen del script).
attr(modelo_var, "seleccion_rezagos") <- criterios
attr(modelo_var, "criterio") <- "AIC"
attr(modelo_var, "rezago") <- p_optimo
attr(modelo_var, "orden_cholesky") <- orden_cholesky
attr(modelo_var, "exogenas") <- c("lluvia", "temp")
attr(modelo_var, "periodo") <- periodo
attr(modelo_var, "fecha_estimacion") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

saveRDS(modelo_var, ruta_modelo)

readr::write_csv(granger_resultados, ruta_granger)

# Un solo CSV con la FEVD (participación por horizonte) y la ITS
# (coeficientes con error estándar y p-valor). Las columnas que no
# aplican a un bloque quedan en NA.
fevd_its_resultados <- dplyr::bind_rows(
  fevd_miedo |>
    dplyr::transmute(
      bloque = "FEVD",
      respuesta = "ocupacion",
      choque = "M_dengue",
      horizonte = horizonte,
      termino = NA_character_,
      estimacion = NA_real_,
      error_estandar = NA_real_,
      estadistico_t = NA_real_,
      p_valor = NA_real_,
      participacion_pct = participacion_pct,
      nota = "Fracción de la varianza de ocupación explicada por M_dengue"
    ),
  its_tabla |>
    dplyr::transmute(
      bloque = "ITS",
      respuesta = "ocupacion",
      choque = NA_character_,
      horizonte = NA_integer_,
      termino = termino,
      estimacion = estimacion,
      error_estandar = error_estandar,
      estadistico_t = estadistico_t,
      p_valor = p_valor,
      participacion_pct = NA_real_,
      nota = "ITS chikungunya 2025: Y ~ M_chik * dummy_2025 + lluvia + temp"
    )
)

readr::write_csv(fevd_its_resultados, ruta_fevd_its)

# Comprobación de que las salidas existen y no están vacías.
rutas_salida <- c(
  ruta_modelo, ruta_granger, ruta_fevd_its,
  ruta_fevd_png, ruta_irf_png, ruta_its_png
)
for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 11. Resumen en consola ----------------------------------------

cat("\n============================================================\n")
cat(" VAR NACIONAL, FEVD, GRANGER E ITS DE CHIKUNGUNYA\n")
cat("============================================================\n")
cat(sprintf("Período: %s | %d meses | %d regiones agregadas\n",
            periodo, nrow(nacional), n_regiones))
cat(sprintf("VAR: [I, M, Y] = [casos_dengue, M_dengue, ocupacion] | p = %d\n",
            p_optimo))
cat("Cholesky: I -> M -> Y | exógenas: lluvia, temp\n")

cat("\n--- RESULTADO PRINCIPAL ---\n")
cat(sprintf(
  "Varianza de ocupación explicada por M_dengue: %.2f%% a %d meses | %.2f%% a %d meses\n",
  fevd_miedo$participacion_pct[1], fevd_miedo$horizonte[1],
  fevd_miedo$participacion_pct[2], fevd_miedo$horizonte[2]
))
cat(sprintf(
  "Granger M_dengue -> ocupación: F(%d, %d) = %.3f, p = %.4f\n",
  granger_m_a_y$gl1, granger_m_a_y$gl2,
  granger_m_a_y$estadistico, granger_m_a_y$p_valor
))
cat(sprintf(
  "Granger casos_dengue -> M_dengue: F(%d, %d) = %.3f, p = %.4f\n",
  granger_i_a_m$gl1, granger_i_a_m$gl2,
  granger_i_a_m$estadistico, granger_i_a_m$p_valor
))
cat(sprintf(
  "ITS, interacción %s: %.4f (EE %.4f), p = %.4f\n",
  termino_interaccion, its_interaccion$estimacion,
  its_interaccion$error_estandar, its_interaccion$p_valor
))

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

cat("\nNotas:\n")
cat("  - Los casos se suman y los índices (M, Y, lluvia, temp) se promedian.\n")
cat("  - La FEVD y la IRF son ortogonalizadas con el orden I -> M -> Y.\n")
cat("  - vars::causality() contrasta la exclusión conjunta en todas las\n")
cat("    ecuaciones; la prueba por ecuación aísla M -> Y e I -> M.\n")
cat("  - La ITS usa la dummy de 2025 y su interacción con M_chik.\n")

message("\n✅ Script 16 finalizado.")
