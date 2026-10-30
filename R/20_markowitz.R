# ============================================================
# 20_markowitz.R
# Teoría del Riesgo — Tabla 4: frontera eficiente de Markowitz
# ============================================================
#
# OBJETIVO
#   Tratar las 12 regiones EMA como activos y la ocupación hotelera
#   mensual (Y) como su rendimiento, para construir la frontera
#   eficiente de Markowitz: para cada nivel de rendimiento objetivo se
#   minimiza la varianza del portafolio con quadprog::solve.QP, sujeto a
#   que los pesos sumen 1 y sean no negativos (sin ventas en corto).
#
#   Además se compara la frontera de los meses normales con la de los
#   meses epidémicos y se dibuja el mapa de isomedia e isovarianza para
#   tres regiones.
#
# RENDIMIENTOS
#
#   r_rt = log(Y_rt) - log(Y_r,t-1) = diff(log(Y_rt)),  por región r.
#
#   Y es la ocupación hotelera de la EMA en porcentaje (26,9 a 78,1 en
#   el panel), siempre positiva, así que el logaritmo está bien definido.
#   El panel 2022-2025 tiene 48 meses por región: se obtienen 47
#   rendimientos mensuales por región y 564 en total. Al diferenciar se
#   pierde enero de 2022, que no tiene mes anterior dentro de la muestra.
#
#   mu    = media muestral de los 47 rendimientos de cada región (12x1).
#   Sigma = covarianza muestral de esos rendimientos (12x12).
#
# FRONTERA EFICIENTE
#
#   Para un rendimiento objetivo t se resuelve
#
#       min_w  w' Sigma w
#       s.a.   sum(w) = 1,  mu'w = t,  w >= 0.
#
#   En quadprog::solve.QP la función objetivo es 1/2 b'Db - d'b, por eso
#   se pasa Dmat = 2*Sigma y dvec = 0. Las dos primeras columnas de Amat
#   son las igualdades (suma de pesos y rendimiento); las 12 columnas de
#   la identidad imponen w >= 0.
#
#   Se calculan 50 niveles de rendimiento entre min(mu) y max(mu). Cada
#   llamada a solve.QP va dentro de tryCatch: si una iteración falla
#   (por ejemplo, en el extremo exacto del intervalo), devuelve NA y el
#   bucle continúa. Por eso la frontera guardada puede tener hasta 50
#   filas por muestra, algunas con NA.
#
#   El portafolio de mínima varianza (MVP) es el caso sin restricción de
#   rendimiento: solo sum(w) = 1 y w >= 0. Se estima con la muestra
#   completa y es el que se guarda en pesos_portafolio.csv.
#
# MESES NORMALES VS. EPIDÉMICOS
#
#   M_dengue es un índice región-mes (residuo de la regresión auxiliar
#   del script 12), no un conteo. Para separar meses se promedia entre
#   las 12 regiones y se obtiene un único calendario epidémico nacional
#   de 48 meses. El umbral es el percentil 75 de ese promedio mensual
#   (tipo 7, el default de R):
#
#       epidémico  si  M_dengue_nacional_mes > P75
#       normal     si  M_dengue_nacional_mes <= P75
#
#   Se clasifica por mes y no por celda región-mes para conservar una
#   matriz completa T x 12 en cada submuestra. Si se descartaran celdas
#   sueltas, la covarianza entre regiones quedaría calculada con meses
#   distintos para cada par y dejaría de ser una matriz de portafolio.
#
#   Con 48 meses, el percentil 75 deja 12 meses epidémicos y 36
#   normales. Como enero de 2022 no tiene rendimiento, la submuestra
#   epidémica tiene 12 retornos (el mes epidémico no es enero) y la
#   normal, 35.
#
#   La submuestra epidémica tiene 12 retornos para 12 activos: por
#   construcción Sigma es singular (rango 11). Se regulariza la diagonal
#   con Sigma + diag(1e-8) para que solve.QP pueda factorizarla. La
#   regularización es mínima y no cambia la interpretación económica:
#   solo permite invertir numéricamente una matriz que, con 12
#   observaciones, no tiene rango completo.
#
# ISOMEDIA E ISOVARIANZA (PASO 4 DEL PROFESOR)
#
#   Para tres regiones (Cartagena, Bogotá D.C. y Antioquia) se recorre
#   la malla w1, w2 en [0, 1] con paso 0,0025 y se conservan los puntos
#   con w1 + w2 <= 1, porque w3 = 1 - w1 - w2. En cada punto:
#
#       media(w)    = w1*mu1 + w2*mu2 + w3*mu3
#       varianza(w) = w' Sigma w.
#
#   Se dibujan las curvas de nivel de ambas superficies. metR no está en
#   el renv.lock del proyecto, así que se usa ggplot2::geom_contour_filled(),
#   que produce el mismo mapa sin añadir dependencias.
#
# INSUMOS
#   data/D_13/panel_ema_den_chik.rds
#
# SALIDAS
#   data/D_20/pesos_portafolio.csv
#       region_ema, peso, mu_mensual (12 filas, ordenadas por peso;
#       mu_mensual en tanto por uno).
#   data/D_20/frontera_datos.csv
#       muestra, rendimiento_objetivo, rendimiento_esperado, varianza,
#       desviacion_estandar (hasta 100 filas: 50 por muestra;
#       rendimientos y desviación en tanto por uno).
#   graficas/G_20/frontera_contraste.png
#   graficas/G_20/isomedia_isovarianza.png
#
# ADVERTENCIAS DE INTERPRETACIÓN
#   - 47 rendimientos mensuales son una muestra corta: la frontera es
#     in-sample y no incorpora incertidumbre de estimación.
#   - No hay activo libre de riesgo, así que no se calcula la línea del
#     mercado de capitales ni el portafolio tangente.
#   - El tramo de la frontera con rendimiento menor que el del MVP es la
#     rama ineficiente: para un riesgo mayor ofrece menos rendimiento.
#     Se conserva porque el ejercicio pide 50 niveles entre min(mu) y
#     max(mu); en la gráfica el rombo marca el MVP de cada muestra.
#   - La frontera epidémica se apoya en 12 retornos y una covarianza
#     regularizada: describe esos 12 meses, no un régimen estable.
#
# NOTA SOBRE EL USO DE ::
#   Todas las funciones no base se llaman con paquete::funcion() para
#   evitar el enmascaramiento entre dplyr, tidyr, quadprog, ggplot2,
#   gridExtra y stats. No se adjuntan paquetes con library().
# ============================================================

# --- 0. Rutas, carpetas y paquetes ---------------------------------

dir.create(here::here("data", "D_20"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("graficas", "G_20"), recursive = TRUE, showWarnings = FALSE)

ruta_panel         <- here::here("data", "D_13", "panel_ema_den_chik.rds")
ruta_pesos         <- here::here("data", "D_20", "pesos_portafolio.csv")
ruta_frontera      <- here::here("data", "D_20", "frontera_datos.csv")
ruta_png_frontera  <- here::here("graficas", "G_20", "frontera_contraste.png")
ruta_png_iso       <- here::here("graficas", "G_20", "isomedia_isovarianza.png")

desde_anio <- 2022L
hasta_anio <- 2025L

n_activos           <- 12L
n_niveles           <- 50L
cuantil_epidemico   <- 0.75
regularizacion      <- 1e-8
paso_malla          <- 0.0025
tolerancia_pesos    <- 1e-8

regiones_iso <- c("Cartagena", "Bogotá D.C.", "Antioquia")

paquetes <- c(
  "quadprog", "dplyr", "tidyr", "ggplot2", "readr",
  "tibble", "here", "gridExtra"
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

if (!file.exists(ruta_panel)) {
  stop(
    "No existe el panel del script 13: ", ruta_panel, ". ",
    "Ejecute primero R/13_panel_ema.R.",
    call. = FALSE
  )
}

# --- 1. Carga y validación del panel -------------------------------

panel <- readRDS(ruta_panel)

if (!inherits(panel, "data.frame")) {
  stop(
    "El panel no es un data.frame: ", paste(class(panel), collapse = ", "), ".",
    call. = FALSE
  )
}

# La columna de economía puede llamarse Y u ocupacion (con o sin tilde).
# Se detecta sin distinguir mayúsculas ni acentos y se homologa a Y.
detectar_columna_y <- function(datos) {
  nombres <- names(datos)
  normalizados <- tolower(iconv(nombres, to = "ASCII//TRANSLIT"))

  if ("Y" %in% nombres) {
    return("Y")
  }

  candidatas <- which(normalizados %in% c("y", "ocupacion"))
  if (length(candidatas) == 0L) {
    return(NA_character_)
  }

  nombres[candidatas[1L]]
}

col_y <- detectar_columna_y(panel)
if (is.na(col_y)) {
  stop(
    "No se encontró la columna de economía (Y u ocupacion). ",
    "Columnas disponibles: ", paste(names(panel), collapse = ", "), ".",
    call. = FALSE
  )
}

if (col_y != "Y") {
  names(panel)[names(panel) == col_y] <- "Y"
  message("Se renombró la columna '", col_y, "' a 'Y'.")
}

requeridas <- c("region_ema", "mes", "Y", "M_dengue")
faltan_panel <- setdiff(requeridas, names(panel))
if (length(faltan_panel) > 0L) {
  stop(
    "Faltan columnas en el panel: ", paste(faltan_panel, collapse = ", "), ". ",
    "Columnas disponibles: ", paste(names(panel), collapse = ", "), ".",
    call. = FALSE
  )
}

panel <- panel |>
  dplyr::mutate(
    mes = as.Date(mes),
    Y = as.numeric(Y),
    M_dengue = as.numeric(M_dengue)
  ) |>
  dplyr::filter(
    mes >= as.Date(sprintf("%d-01-01", desde_anio)),
    mes <= as.Date(sprintf("%d-12-01", hasta_anio))
  ) |>
  dplyr::arrange(region_ema, mes)

if (nrow(panel) == 0L) {
  stop("El panel no tiene observaciones entre 2022 y 2025.", call. = FALSE)
}

regiones <- sort(unique(panel$region_ema))

if (length(regiones) != n_activos) {
  stop(
    "Se esperaban ", n_activos, " regiones EMA y hay ", length(regiones),
    ": ", paste(regiones, collapse = ", "), ".",
    call. = FALSE
  )
}

meses <- sort(unique(panel$mes))
if (length(meses) != 48L) {
  warning(
    "Se esperaban 48 meses (2022-01 a 2025-12) y hay ", length(meses),
    ". La frontera se calcula con los meses disponibles.",
    call. = FALSE
  )
}

if (anyDuplicated(panel[, c("region_ema", "mes")])) {
  stop("Hay regiones-mes duplicadas en el panel.", call. = FALSE)
}

conteo_meses <- table(panel$region_ema)
if (length(unique(as.integer(conteo_meses))) != 1L) {
  stop(
    "El panel no está balanceado: las regiones tienen distinto número de ",
    "meses (", paste(names(conteo_meses), conteo_meses, sep = " = ",
                     collapse = "; "), ").",
    call. = FALSE
  )
}

if (anyNA(panel$Y) || any(!is.finite(panel$Y))) {
  stop(
    "Y tiene ", sum(is.na(panel$Y) | !is.finite(panel$Y)),
    " valores faltantes o no finitos; no se pueden calcular logaritmos.",
    call. = FALSE
  )
}

if (any(panel$Y <= 0)) {
  stop(
    "Hay ", sum(panel$Y <= 0),
    " observaciones de Y no positivas; diff(log(Y)) no está definido.",
    call. = FALSE
  )
}

if (anyNA(panel$M_dengue) || any(!is.finite(panel$M_dengue))) {
  stop(
    "M_dengue tiene valores faltantes o no finitos; no se pueden separar ",
    "los meses epidémicos.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat(" PANEL DE ENTRADA\n")
cat("============================================================\n")
cat(sprintf(
  "Filas: %d | Regiones: %d | Meses: %d | Periodo: %s a %s\n",
  nrow(panel), length(regiones), length(meses),
  format(min(meses), "%Y-%m"), format(max(meses), "%Y-%m")
))
cat(sprintf(
  "Y (ocupación): mínimo %.2f %% | media %.2f %% | máximo %.2f %%\n",
  min(panel$Y), mean(panel$Y), max(panel$Y)
))

# --- 2. Rendimientos logarítmicos, mu y Sigma -----------------------
#
# diff(log(Y)) dentro de cada región usa el orden temporal del panel.
# El primer mes de cada región queda en NA y se elimina.

retornos <- panel |>
  dplyr::group_by(region_ema) |>
  dplyr::mutate(r = c(NA_real_, diff(log(Y)))) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(r), is.finite(r))

n_ret_region <- table(retornos$region_ema)
if (length(unique(as.integer(n_ret_region))) != 1L) {
  stop(
    "Las regiones no tienen el mismo número de rendimientos: ",
    paste(names(n_ret_region), n_ret_region, sep = " = ", collapse = "; "),
    ".",
    call. = FALSE
  )
}

ret_ancho <- retornos |>
  dplyr::select(mes, region_ema, r) |>
  tidyr::pivot_wider(names_from = region_ema, values_from = r) |>
  dplyr::arrange(mes)

if (!all(regiones %in% names(ret_ancho))) {
  stop(
    "La matriz de rendimientos no contiene todas las regiones. ",
    "Faltan: ", paste(setdiff(regiones, names(ret_ancho)), collapse = ", "),
    ".",
    call. = FALSE
  )
}

matriz_ret <- as.matrix(ret_ancho[, regiones])
rownames(matriz_ret) <- format(ret_ancho$mes, "%Y-%m")

if (anyNA(matriz_ret)) {
  stop(
    "La matriz de rendimientos tiene ", sum(is.na(matriz_ret)),
    " celdas vacías; alguna región no tiene rendimiento en algún mes.",
    call. = FALSE
  )
}

if (nrow(matriz_ret) <= n_activos) {
  stop(
    "Hay ", nrow(matriz_ret), " retornos para ", n_activos,
    " activos; la covarianza no tendría rango completo.",
    call. = FALSE
  )
}

mu <- colMeans(matriz_ret)
Sigma <- stats::cov(matriz_ret)

if (any(!is.finite(mu)) || any(!is.finite(Sigma))) {
  stop("mu o Sigma tienen valores no finitos.", call. = FALSE)
}

# Simetría exacta: cov() la garantiza, pero conviene fijarla para que
# quadprog no reciba una matriz con diferencias de redondeo.
Sigma <- (Sigma + t(Sigma)) / 2

cat("\n============================================================\n")
cat(" RENDIMIENTOS MENSUALES DE LA OCUPACIÓN HOTELERA\n")
cat("============================================================\n")
cat(sprintf(
  "Retornos: %d meses x %d regiones = %d observaciones | %s a %s\n",
  nrow(matriz_ret), length(regiones), length(matriz_ret),
  format(min(ret_ancho$mes), "%Y-%m"), format(max(ret_ancho$mes), "%Y-%m")
))
cat(sprintf(
  "mu mensual: mínimo %+.4f %% (%s) | máximo %+.4f %% (%s)\n",
  100 * min(mu), names(mu)[which.min(mu)],
  100 * max(mu), names(mu)[which.max(mu)]
))

# --- 3. Utilidades de optimización ---------------------------------
#
# resolver_qp() devuelve NA (no un error) cuando solve.QP falla, de modo
# que el bucle de la frontera pueda continuar con el siguiente nivel.

resolver_qp <- function(Sigma, mu, objetivo = NULL,
                        tolerancia = tolerancia_pesos) {
  n <- nrow(Sigma)
  Dmat <- 2 * Sigma
  dvec <- rep(0, n)

  if (is.null(objetivo)) {
    # MVP: solo suma de pesos y no negatividad.
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, rep(0, n))
    meq <- 1L
  } else {
    # Frontera: suma de pesos, rendimiento objetivo y no negatividad.
    Amat <- cbind(rep(1, n), mu, diag(n))
    bvec <- c(1, objetivo, rep(0, n))
    meq <- 2L
  }

  solucion <- tryCatch(
    quadprog::solve.QP(
      Dmat = Dmat,
      dvec = dvec,
      Amat = Amat,
      bvec = bvec,
      meq = meq
    ),
    error = function(e) NULL
  )

  if (is.null(solucion)) {
    return(rep(NA_real_, n))
  }

  w <- as.numeric(solucion$solution)
  if (length(w) != n || any(!is.finite(w))) {
    return(rep(NA_real_, n))
  }

  # solve.QP respeta w >= 0; se toleran negativos de redondeo y se
  # renormaliza para que la suma sea exactamente 1. Un negativo grande
  # sería señal de una solución inválida y se descarta.
  if (any(w < -tolerancia)) {
    return(rep(NA_real_, n))
  }
  w[w < 0] <- 0
  if (sum(w) <= 0) {
    return(rep(NA_real_, n))
  }
  w <- w / sum(w)
  names(w) <- names(mu)
  w
}

varianza_portafolio <- function(w, Sigma) {
  if (anyNA(w)) {
    return(NA_real_)
  }
  as.numeric(crossprod(w, Sigma %*% w))
}

frontera_eficiente <- function(mu, Sigma, niveles = n_niveles) {
  n <- length(mu)
  objetivos <- seq(min(mu), max(mu), length.out = niveles)
  pesos <- matrix(
    NA_real_,
    nrow = length(objetivos),
    ncol = n,
    dimnames = list(NULL, names(mu))
  )

  for (i in seq_along(objetivos)) {
    w <- resolver_qp(Sigma, mu, objetivo = objetivos[i])
    if (anyNA(w)) {
      next
    }
    pesos[i, ] <- w
  }

  rendimiento <- as.numeric(pesos %*% mu)
  varianza <- rowSums((pesos %*% Sigma) * pesos)

  tibble::tibble(
    rendimiento_objetivo = objetivos,
    rendimiento_esperado = rendimiento,
    varianza = varianza,
    desviacion_estandar = sqrt(varianza)
  )
}

# Regulariza la diagonal solo si Sigma no es definida positiva. Con 12
# retornos y 12 activos la matriz epidémica es singular por construcción.
asegurar_definida <- function(Sigma, etiqueta, refuerzo = regularizacion) {
  Sigma <- (Sigma + t(Sigma)) / 2
  n <- nrow(Sigma)
  rango <- qr(Sigma)$rank

  es_definida_positiva <- !is.null(
    tryCatch(chol(Sigma), error = function(e) NULL)
  )

  if (rango < n || !es_definida_positiva || any(!is.finite(Sigma))) {
    message(
      "  Sigma de ", etiqueta, " no es definida positiva (rango ",
      rango, " de ", n, "); se añade ", format(refuerzo),
      " a la diagonal."
    )
    Sigma <- Sigma + diag(refuerzo, n)
  }

  if (is.null(tryCatch(chol(Sigma), error = function(e) NULL))) {
    warning(
      "Sigma de ", etiqueta, " sigue sin ser definida positiva tras la ",
      "regularización; solve.QP podría devolver NA en algunos niveles.",
      call. = FALSE
    )
  }

  Sigma
}

# --- 4. Portafolio de mínima varianza (muestra completa) ------------

w_mvp <- resolver_qp(Sigma, mu, objetivo = NULL)

if (anyNA(w_mvp)) {
  stop(
    "No se pudo resolver el portafolio de mínima varianza de la muestra ",
    "completa. Revise la matriz Sigma.",
    call. = FALSE
  )
}

ret_mvp <- sum(w_mvp * mu)
var_mvp <- varianza_portafolio(w_mvp, Sigma)
sd_mvp <- sqrt(var_mvp)

if (abs(sum(w_mvp) - 1) > 1e-6 || any(w_mvp < -1e-10)) {
  stop("El portafolio de mínima varianza no cumple las restricciones.",
       call. = FALSE)
}

# Comprobación económica: el MVP no puede tener más varianza que
# cualquier activo individual.
if (var_mvp > min(diag(Sigma)) + 1e-10) {
  warning(
    "La varianza del MVP (", sprintf("%.8f", var_mvp),
    ") supera la del activo de menor varianza (",
    sprintf("%.8f", min(diag(Sigma))),
    "); revise la optimización.",
    call. = FALSE
  )
}

pesos_tbl <- tibble::tibble(
  region_ema = regiones,
  peso = as.numeric(w_mvp[regiones]),
  mu_mensual = as.numeric(mu[regiones])
) |>
  dplyr::arrange(dplyr::desc(peso))

if (abs(sum(pesos_tbl$peso) - 1) > 1e-6) {
  stop("Los pesos guardados no suman 1.", call. = FALSE)
}

# --- 5. Meses normales y epidémicos --------------------------------
#
# M_dengue se promedia entre regiones para obtener un calendario único.
# Así todas las regiones comparten los mismos meses en cada submuestra.

m_mensual <- panel |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    M_dengue = mean(M_dengue, na.rm = TRUE),
    n_regiones = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(mes)

if (anyNA(m_mensual$M_dengue) || any(!is.finite(m_mensual$M_dengue))) {
  stop("El promedio mensual de M_dengue tiene valores no finitos.",
       call. = FALSE)
}

if (dplyr::n_distinct(m_mensual$n_regiones) != 1L) {
  warning(
    "El número de regiones por mes no es constante: ",
    paste(sort(unique(m_mensual$n_regiones)), collapse = ", "), ".",
    call. = FALSE
  )
}

umbral_epi <- stats::quantile(
  m_mensual$M_dengue,
  probs = cuantil_epidemico,
  names = FALSE,
  type = 7
)

meses_epi <- m_mensual$mes[m_mensual$M_dengue > umbral_epi]
meses_norm <- m_mensual$mes[m_mensual$M_dengue <= umbral_epi]

if (length(meses_epi) == 0L || length(meses_norm) == 0L) {
  stop(
    "La separación epidémica quedó vacía: ", length(meses_epi),
    " meses epidémicos y ", length(meses_norm), " normales.",
    call. = FALSE
  )
}

matriz_de_meses <- function(ancho, meses, regiones) {
  filas <- ancho$mes %in% meses
  if (!any(filas)) {
    return(NULL)
  }
  as.matrix(ancho[filas, regiones])
}

mat_normal <- matriz_de_meses(ret_ancho, meses_norm, regiones)
mat_epi <- matriz_de_meses(ret_ancho, meses_epi, regiones)

if (is.null(mat_normal) || is.null(mat_epi)) {
  stop("Alguna submuestra quedó sin retornos.", call. = FALSE)
}
if (anyNA(mat_normal) || anyNA(mat_epi)) {
  stop("Las submuestras normal o epidémica tienen celdas vacías.",
       call. = FALSE)
}
if (nrow(mat_normal) <= n_activos || nrow(mat_epi) <= 1L) {
  stop(
    "Submuestras insuficientes: normal ", nrow(mat_normal),
    " retornos, epidémica ", nrow(mat_epi), " retornos.",
    call. = FALSE
  )
}

mu_normal <- colMeans(mat_normal)
Sigma_normal <- stats::cov(mat_normal)
mu_epi <- colMeans(mat_epi)
Sigma_epi <- stats::cov(mat_epi)

rango_sigma_normal <- qr(Sigma_normal)$rank
rango_sigma_epi <- qr(Sigma_epi)$rank

Sigma_normal <- asegurar_definida(Sigma_normal, "los meses normales")
Sigma_epi <- asegurar_definida(Sigma_epi, "los meses epidémicos")

frontera_normal <- frontera_eficiente(mu_normal, Sigma_normal) |>
  dplyr::mutate(muestra = "Normal")

frontera_epi <- frontera_eficiente(mu_epi, Sigma_epi) |>
  dplyr::mutate(muestra = "Epidémica")

fronteras <- dplyr::bind_rows(frontera_normal, frontera_epi) |>
  dplyr::mutate(
    muestra = factor(muestra, levels = c("Normal", "Epidémica"))
  )

n_validos_normal <- sum(!is.na(frontera_normal$varianza))
n_validos_epi <- sum(!is.na(frontera_epi$varianza))

if (n_validos_normal < 2L || n_validos_epi < 2L) {
  warning(
    "Fronteras incompletas: normal ", n_validos_normal, "/", n_niveles,
    " puntos, epidémica ", n_validos_epi, "/", n_niveles,
    ". La gráfica omite los niveles sin solución.",
    call. = FALSE
  )
}

# MVP de cada submuestra: se marca en la gráfica de contraste.
w_mvp_normal <- resolver_qp(Sigma_normal, mu_normal)
w_mvp_epi <- resolver_qp(Sigma_epi, mu_epi)

ret_mvp_normal <- if (anyNA(w_mvp_normal)) NA_real_ else
  sum(w_mvp_normal * mu_normal)
sd_mvp_normal <- if (anyNA(w_mvp_normal)) NA_real_ else
  sqrt(varianza_portafolio(w_mvp_normal, Sigma_normal))

ret_mvp_epi <- if (anyNA(w_mvp_epi)) NA_real_ else
  sum(w_mvp_epi * mu_epi)
sd_mvp_epi <- if (anyNA(w_mvp_epi)) NA_real_ else
  sqrt(varianza_portafolio(w_mvp_epi, Sigma_epi))

mvp_puntos <- tibble::tibble(
  muestra = factor(c("Normal", "Epidémica"), levels = c("Normal", "Epidémica")),
  rendimiento = c(ret_mvp_normal, ret_mvp_epi),
  desviacion = c(sd_mvp_normal, sd_mvp_epi)
)

cat("\n============================================================\n")
cat(" MESES NORMALES Y EPIDÉMICOS\n")
cat("============================================================\n")
cat(sprintf(
  "M_dengue nacional: media %.3f | P75 = %.3f\n",
  mean(m_mensual$M_dengue), umbral_epi
))
cat(sprintf(
  "Meses epidémicos: %d de %d | meses normales: %d\n",
  length(meses_epi), nrow(m_mensual), length(meses_norm)
))
cat(sprintf(
  "Retornos: normal %d | epidémica %d (enero de 2022 no tiene retorno)\n",
  nrow(mat_normal), nrow(mat_epi)
))
cat(sprintf(
  "Rango de Sigma: normal %d/%d | epidémica %d/%d\n",
  rango_sigma_normal, n_activos, rango_sigma_epi, n_activos
))

# --- 6. Gráfica de contraste de fronteras --------------------------
#
# Riesgo (desviación estándar) en el eje x y rendimiento esperado en el
# eje y. El rombo marca el MVP de cada submuestra. La rama con
# rendimiento menor que el del MVP es ineficiente y se dibuja tenue.

mvp_retornos <- tibble::tibble(
  muestra = factor(c("Normal", "Epidémica"), levels = c("Normal", "Epidémica")),
  ret_mvp = c(ret_mvp_normal, ret_mvp_epi)
)

fronteras_grafico <- fronteras |>
  dplyr::left_join(mvp_retornos, by = "muestra") |>
  dplyr::filter(
    !is.na(desviacion_estandar),
    !is.na(rendimiento_esperado)
  ) |>
  dplyr::mutate(
    eficiente = !is.na(ret_mvp) & rendimiento_esperado >= ret_mvp - 1e-10
  )

fronteras_eficientes <- dplyr::filter(fronteras_grafico, eficiente)
fronteras_ineficientes <- dplyr::filter(fronteras_grafico, !eficiente)

caption_frontera <- paste0(
  "Rendimientos logarítmicos mensuales de la ocupación hotelera ",
  "(diff(log(Y))). Frontera calculada con quadprog::solve.QP ",
  "(sum(w) = 1, w >= 0) para 50 niveles entre min(mu) y max(mu). ",
  "Meses epidémicos: M_dengue nacional > percentil 75. La Sigma ",
  "epidémica se regulariza con +1e-8 en la diagonal porque tiene 12 ",
  "retornos para 12 activos. Rombo: portafolio de mínima varianza de ",
  "cada muestra; las líneas tenues son la rama ineficiente ",
  "(rendimiento menor que el del MVP). Muestra corta: in-sample."
)
caption_frontera <- paste(
  strwrap(caption_frontera, width = 160),
  collapse = "\n"
)

p_frontera <- ggplot2::ggplot(
  fronteras_grafico,
  ggplot2::aes(
    x = 100 * desviacion_estandar,
    y = 100 * rendimiento_esperado,
    colour = muestra
  )
) +
  ggplot2::geom_path(
    data = fronteras_ineficientes,
    linewidth = 0.65,
    alpha = 0.30
  ) +
  ggplot2::geom_point(
    data = fronteras_ineficientes,
    size = 0.9,
    alpha = 0.25
  ) +
  ggplot2::geom_path(
    data = fronteras_eficientes,
    linewidth = 0.90
  ) +
  ggplot2::geom_point(
    data = fronteras_eficientes,
    size = 1.0,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    data = mvp_puntos,
    ggplot2::aes(
      x = 100 * desviacion,
      y = 100 * rendimiento,
      colour = muestra
    ),
    size = 3.4,
    shape = 18,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c("Normal" = "#1B7837", "Epidémica" = "#B2182B"),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Frontera eficiente: meses normales vs. meses epidémicos",
    subtitle = sprintf(
      paste0(
        "12 regiones EMA | %d retornos mensuales (%s a %s) | ",
        "normal: %d retornos | epidémica: %d retornos"
      ),
      nrow(matriz_ret),
      format(min(ret_ancho$mes), "%Y-%m"),
      format(max(ret_ancho$mes), "%Y-%m"),
      nrow(mat_normal),
      nrow(mat_epi)
    ),
    x = "Riesgo: desviación estándar mensual (%)",
    y = "Rendimiento esperado mensual (%)",
    caption = caption_frontera
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "left",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 8),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 7, colour = "grey35"
    )
  )

ggplot2::ggsave(
  filename = ruta_png_frontera,
  plot = p_frontera,
  width = 10.5,
  height = 6.6,
  dpi = 200,
  bg = "white"
)

# --- 7. Isomedia e isovarianza para tres regiones -------------------

faltan_iso <- setdiff(regiones_iso, regiones)
if (length(faltan_iso) > 0L) {
  stop(
    "No están en region_ema: ", paste(faltan_iso, collapse = ", "), ". ",
    "Regiones disponibles: ", paste(regiones, collapse = ", "), ".",
    call. = FALSE
  )
}

malla <- expand.grid(
  w1 = seq(0, 1, by = paso_malla),
  w2 = seq(0, 1, by = paso_malla)
)

# w3 = 1 - w1 - w2: solo existe el triángulo w1 + w2 <= 1.
malla <- malla[malla$w1 + malla$w2 <= 1 + 1e-12, , drop = FALSE]
malla$w3 <- 1 - malla$w1 - malla$w2

pesos_malla <- as.matrix(malla[, c("w1", "w2", "w3")])
mu_iso <- mu[regiones_iso]
Sigma_iso <- Sigma[regiones_iso, regiones_iso]

malla$media <- as.numeric(pesos_malla %*% mu_iso)
malla$varianza <- as.numeric(
  rowSums((pesos_malla %*% Sigma_iso) * pesos_malla)
)

malla_larga <- malla |>
  dplyr::select(w1, w2, media, varianza) |>
  tidyr::pivot_longer(
    cols = c(media, varianza),
    names_to = "medida",
    values_to = "valor"
  )

datos_media <- dplyr::filter(malla_larga, medida == "media")
datos_var <- dplyr::filter(malla_larga, medida == "varianza")

# La base y el lado izquierdo los aportan los ejes; se dibuja solo la
# hipotenusa w1 + w2 = 1 para cerrar el simplex.
borde_diagonal <- tibble::tibble(
  w1 = c(1, 0),
  w2 = c(0, 1)
)

grafico_contorno <- function(datos, titulo, subtitulo, nombre_fill) {
  ggplot2::ggplot(datos, ggplot2::aes(x = w1, y = w2, z = valor)) +
    ggplot2::geom_contour_filled(
      bins = 9,
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_path(
      data = borde_diagonal,
      ggplot2::aes(x = w1, y = w2),
      inherit.aes = FALSE,
      colour = "grey35",
      linewidth = 0.45
    ) +
    ggplot2::scale_fill_viridis_d(option = "viridis", name = nombre_fill) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = titulo,
      subtitle = subtitulo,
      x = "w1: peso de Cartagena",
      y = "w2: peso de Bogotá D.C."
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 7.5),
      legend.text = ggplot2::element_text(size = 6.5),
      legend.key.width = grid::unit(0.9, "lines"),
      legend.key.height = grid::unit(0.55, "lines"),
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.subtitle = ggplot2::element_text(size = 7.8, colour = "grey30"),
      plot.title.position = "plot"
    )
}

p_iso_media <- grafico_contorno(
  datos = datos_media,
  titulo = "Isomedia: rendimiento esperado",
  subtitulo = paste0(
    "Rendimiento esperado (mu'w) de la muestra completa; ",
    "w3 = 1 - w1 - w2 es Antioquia y solo se grafica w1 + w2 <= 1."
  ),
  nombre_fill = "Rendimiento\nmensual"
)

p_iso_var <- grafico_contorno(
  datos = datos_var,
  titulo = "Isovarianza: varianza del portafolio",
  subtitulo = paste0(
    "Curvas de nivel de w' Sigma w con mu y Sigma de la muestra completa ",
    "(47 rendimientos)."
  ),
  nombre_fill = "Varianza\nmensual"
)

nota_iso <- grid::textGrob(
  paste0(
    "Malla de paso ", sprintf("%.4f", paso_malla), " en w1 y w2: ",
    format(nrow(malla), big.mark = ".", decimal.mark = ","),
    " portafolios factibles. Contornos con ggplot2::geom_contour_filled(); ",
    "metR no está instalado en el proyecto."
  ),
  gp = grid::gpar(fontsize = 7, col = "grey35"),
  hjust = 0,
  x = 0.01
)

# ggplotGrob() necesita un dispositivo gráfico. En Rscript se abriría
# Rplots.pdf; se usa un dispositivo PDF nulo para no dejar archivos.
grDevices::pdf(NULL)
grafico_iso <- gridExtra::arrangeGrob(
  ggplot2::ggplotGrob(p_iso_media),
  ggplot2::ggplotGrob(p_iso_var),
  ncol = 2,
  bottom = nota_iso
)
grDevices::dev.off()

ggplot2::ggsave(
  filename = ruta_png_iso,
  plot = grafico_iso,
  width = 12,
  height = 7,
  dpi = 200,
  bg = "white"
)

# --- 8. Guardado ---------------------------------------------------

fronteras_tbl <- fronteras |>
  dplyr::select(
    muestra,
    rendimiento_objetivo,
    rendimiento_esperado,
    varianza,
    desviacion_estandar
  )

readr::write_csv(pesos_tbl, ruta_pesos, na = "NA")
readr::write_csv(fronteras_tbl, ruta_frontera, na = "NA")

rutas_salida <- c(
  ruta_pesos,
  ruta_frontera,
  ruta_png_frontera,
  ruta_png_iso
)

for (ruta in rutas_salida) {
  if (!file.exists(ruta) || file.info(ruta)$size == 0) {
    stop("No se generó la salida: ", ruta, call. = FALSE)
  }
}

# --- 9. Resumen en consola -----------------------------------------

top3 <- pesos_tbl[seq_len(3), ]
resto <- sum(pesos_tbl$peso[-seq_len(3)])

cat("\n============================================================\n")
cat(" TABLA 4 — PORTAFOLIO DE MÍNIMA VARIANZA (12 REGIONES)\n")
cat("============================================================\n")
cat(sprintf(
  "Muestra completa: %d retornos mensuales (%s a %s)\n",
  nrow(matriz_ret),
  format(min(ret_ancho$mes), "%Y-%m"),
  format(max(ret_ancho$mes), "%Y-%m")
))

cat("\nLas 3 regiones con mayor peso:\n")
for (i in seq_len(nrow(top3))) {
  cat(sprintf(
    "  %-32s %7.2f %%\n",
    top3$region_ema[i],
    100 * top3$peso[i]
  ))
}
cat(sprintf(
  "  %-32s %7.2f %%\n",
  "Resto de regiones",
  100 * resto
))

cat(sprintf(
  "\nRendimiento esperado del MVP: %+.4f %% mensual\n",
  100 * ret_mvp
))
cat(sprintf(
  "Desviación estándar del MVP:  %.4f %% mensual\n",
  100 * sd_mvp
))
cat(sprintf(
  "Varianza del MVP:             %.8f\n",
  var_mvp
))
cat(sprintf(
  "Pesos: suma %.6f | mínimo %.6f | máximo %.4f\n",
  sum(pesos_tbl$peso), min(pesos_tbl$peso), max(pesos_tbl$peso)
))

cat("\n--- Fronteras eficientes ---\n")
cat(sprintf(
  "Puntos resueltos: normal %d/%d | epidémica %d/%d\n",
  n_validos_normal, n_niveles, n_validos_epi, n_niveles
))

rango_finito <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(c(NA_real_, NA_real_))
  }
  range(x)
}

sd_normal <- rango_finito(frontera_normal$desviacion_estandar)
sd_epi <- rango_finito(frontera_epi$desviacion_estandar)

cat(sprintf(
  "Riesgo (sd mensual): normal %.2f %% a %.2f %% | epidémica %.2f %% a %.2f %%\n",
  100 * sd_normal[1], 100 * sd_normal[2],
  100 * sd_epi[1], 100 * sd_epi[2]
))
cat(sprintf(
  "MVP normal: ret %+.3f %% | sd %.3f %% | MVP epidémica: ret %+.3f %% | sd %.3f %%\n",
  100 * ret_mvp_normal, 100 * sd_mvp_normal,
  100 * ret_mvp_epi, 100 * sd_mvp_epi
))

cat("\n--- Isomedia e isovarianza (3 regiones) ---\n")
cat(sprintf(
  "Regiones: %s\n",
  paste(regiones_iso, collapse = ", ")
))
cat(sprintf(
  "Portafolios en la malla: %s | paso %.4f\n",
  format(nrow(malla), big.mark = ".", decimal.mark = ","),
  paso_malla
))
cat(sprintf(
  "Rendimiento esperado: %.3f %% a %.3f %%\n",
  100 * min(malla$media), 100 * max(malla$media)
))
cat(sprintf(
  "Varianza: %.6f a %.6f\n",
  min(malla$varianza), max(malla$varianza)
))

cat("\nSalidas:\n")
for (ruta in rutas_salida) {
  cat("  ", ruta, "\n", sep = "")
}

cat("\nNotas:\n")
cat("  - 47 rendimientos mensuales es una muestra corta: la frontera es\n")
cat("    in-sample y no incorpora incertidumbre de estimación.\n")
cat("  - No hay activo libre de riesgo: no se calcula CML ni portafolio\n")
cat("    tangente.\n")
cat("  - El tramo de cada frontera con rendimiento menor que el del MVP\n")
cat("    es la rama ineficiente; el rombo marca el MVP de cada muestra.\n")
cat("  - La Sigma epidémica se regularizó con +1e-8 en la diagonal porque\n")
cat("    12 retornos no dan rango completo a una matriz de 12x12.\n")

message("\n✅ Script 20 finalizado.")
