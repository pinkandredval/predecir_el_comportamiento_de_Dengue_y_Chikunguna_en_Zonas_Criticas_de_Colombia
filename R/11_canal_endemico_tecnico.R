# ============================================================
# 11_canal_endemico_tecnico.R
# Canal endémico (método de Bortman) y análisis técnico
# Dengue | regiones EMA | 2020-2025
# ============================================================
#
# INSUMOS
#   data/series/serie_dengue_tasas.csv   (generado por 10_tasas_incidencia.R)
#     columnas: cod_mun_o, anio, semana, fecha_sem, casos, region_ema,
#               poblacion_region, tasa
#   data/series/ranking_tasas_dengue.csv (generado por 10_tasas_incidencia.R)
#     se usa para ordenar las regiones y elegir el Top 6 del análisis técnico.
#
# SALIDAS
#   graficas/11_canal_endemico.png
#   graficas/11_analisis_tecnico.png  (requiere el paquete tidyquant)
#
# MÉTODOS
#   Canal endémico de Bortman:
#     Para cada región y semana epidemiológica se calcula la media geométrica
#     de los casos y su intervalo de confianza al 95% mediante transformación
#     logarítmica. Se usa log(1 + casos) para que las semanas con cero no
#     indefinan la media; al back-transformar con expm1 el canal queda en la
#     escala original de casos.
#     Zonas: éxito (≤ IC inferior), seguridad (IC inferior-media), alerta
#     (media-IC superior) y brote (> IC superior).
#     La semana 53 se estima junto con la 52 porque solo 2020 y 2025 tuvieron
#     53 semanas; el canal resultante se replica en la semana 53.
#
#   Análisis técnico:
#     OHLC mensual a partir de los totales semanales (Open = primer valor del
#     mes, High = máximo, Low = mínimo, Close = último valor).
#     Indicadores semanales: SMA(4), EMA(12), Bandas de Bollinger (20 semanas,
#     ±2 desviaciones estándar) y MACD(12, 26, 9). Las bandas se representan
#     al cierre de cada mes para coincidir con las velas mensuales.
#     Se grafican las 6 regiones más afectadas del ranking del Script 10; las
#     otras 6 se omiten por espacio y se informan en consola.
#
# NOTA SOBRE tidyquant
#   La regla del proyecto exige tidyquant para habilitar el análisis técnico.
#   Si no está instalado, el script guarda primero el canal endémico, imprime
#   la instrucción de instalación y se detiene. Las velas se dibujan con
#   tidyquant::geom_candlestick(); los indicadores se calculan con TTR, que es
#   el motor de indicadores de quantmod.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(scales)
  library(TTR)
})

# --- 0. Rutas -----------------------------------------------------

ruta_graficas <- here("graficas")

dir.create(ruta_graficas, showWarnings = FALSE, recursive = TRUE)

ruta_dengue <- here("data", "series", "serie_dengue_tasas.csv")
ruta_ranking <- here("data", "series", "ranking_tasas_dengue.csv")

if (!file.exists(ruta_dengue)) {
  stop(
    "No existe '", ruta_dengue, "'. ",
    "Ejecute primero R/10_tasas_incidencia.R.",
    call. = FALSE
  )
}

if (!file.exists(ruta_ranking)) {
  stop(
    "No existe '", ruta_ranking, "'. ",
    "Ejecute primero R/10_tasas_incidencia.R.",
    call. = FALSE
  )
}

# --- 1. Lectura y validación de insumos ---------------------------

message("Leyendo series de tasas y ranking del Script 10...")

dengue <- read_csv(
  ruta_dengue,
  show_col_types = FALSE,
  progress = FALSE
)

ranking <- read_csv(
  ruta_ranking,
  show_col_types = FALSE,
  progress = FALSE
)

requeridas_dengue <- c(
  "region_ema", "anio", "semana", "fecha_sem", "casos"
)
faltantes_dengue <- setdiff(requeridas_dengue, names(dengue))

if (length(faltantes_dengue) > 0) {
  stop(
    "Faltan columnas en '", basename(ruta_dengue), "': ",
    paste(faltantes_dengue, collapse = ", "), ".",
    call. = FALSE
  )
}

requeridas_ranking <- c("posicion", "region_ema")
faltantes_ranking <- setdiff(requeridas_ranking, names(ranking))

if (length(faltantes_ranking) > 0) {
  stop(
    "Faltan columnas en '", basename(ruta_ranking), "': ",
    paste(faltantes_ranking, collapse = ", "), ".",
    call. = FALSE
  )
}

dengue <- dengue |>
  mutate(
    region_ema = as.character(region_ema),
    anio = as.integer(anio),
    semana = as.integer(semana),
    fecha_sem = as.Date(fecha_sem),
    casos = as.numeric(casos)
  ) |>
  filter(
    !is.na(region_ema),
    !is.na(anio),
    !is.na(semana),
    !is.na(fecha_sem),
    !is.na(casos),
    anio >= 2020,
    anio <= 2025,
    semana >= 1,
    semana <= 53
  )

if (nrow(dengue) == 0) {
  stop(
    "La serie de dengue quedó vacía después de filtrar 2020-2025 y ",
    "semanas 1-53. Revise '", basename(ruta_dengue), "'.",
    call. = FALSE
  )
}

if (any(dengue$casos < 0)) {
  stop("Hay casos negativos en la serie; revise el Script 10.", call. = FALSE)
}

regiones_ema <- sort(unique(dengue$region_ema))

if (length(regiones_ema) != 12) {
  stop(
    "Se esperaban 12 regiones EMA y se encontraron ",
    length(regiones_ema), ": ",
    paste(regiones_ema, collapse = ", "), ".",
    call. = FALSE
  )
}

ranking <- ranking |>
  mutate(
    posicion = as.integer(posicion),
    region_ema = as.character(region_ema)
  ) |>
  filter(!is.na(posicion), !is.na(region_ema)) |>
  arrange(posicion)

if (nrow(ranking) < 6) {
  stop(
    "El ranking debe tener al menos 6 regiones para elegir el Top 6. ",
    "Revise '", basename(ruta_ranking), "'.",
    call. = FALSE
  )
}

if (anyDuplicated(ranking$region_ema)) {
  warning(
    "Hay regiones repetidas en el ranking; se conserva la primera posición.",
    call. = FALSE
  )
  ranking <- ranking |> distinct(region_ema, .keep_all = TRUE)
}

# Orden de las regiones: el del ranking (más afectadas primero). Si alguna
# región de la serie no está en el ranking, se agrega al final para no perderla.
niveles_region <- ranking$region_ema
faltan_en_ranking <- setdiff(regiones_ema, niveles_region)

if (length(faltan_en_ranking) > 0) {
  warning(
    "Hay regiones de la serie que no están en el ranking: ",
    paste(faltan_en_ranking, collapse = ", "), ". Se agregan al final.",
    call. = FALSE
  )
  niveles_region <- c(niveles_region, faltan_en_ranking)
}

# --- 2. Totales semanales por región ------------------------------

message("Agrupando casos por región y semana epidemiológica...")

serie_semanal <- dengue |>
  group_by(region_ema, anio, semana) |>
  summarise(
    fecha_sem = min(fecha_sem),
    casos = sum(casos, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(region_ema, anio, semana)

if (anyDuplicated(serie_semanal[, c("region_ema", "anio", "semana")])) {
  stop(
    "La agregación semanal produjo filas duplicadas por región-año-semana.",
    call. = FALSE
  )
}

n_semanas <- serie_semanal |>
  count(region_ema, name = "n_semanas")

if (any(n_semanas$n_semanas != 314)) {
  warning(
    "No todas las regiones tienen 314 semanas (2020-2025). ",
    "Semanas por región: ",
    paste(
      sprintf("%s (%d)", n_semanas$region_ema, n_semanas$n_semanas),
      collapse = "; "
    ),
    ".",
    call. = FALSE
  )
}

serie_semanal <- serie_semanal |>
  mutate(region_ema = factor(region_ema, levels = niveles_region))

# --- 3. Canal endémico de Bortman ---------------------------------

# Calcula, para cada región y semana, la media geométrica y el IC 95% con
# transformación logarítmica. La semana 53 se agrupa con la 52 para estimar
# (solo 2020 y 2025 tienen 53 semanas) y luego se replica en la semana 53.
calcular_canal_bortman <- function(serie_semanal, z = 1.96) {
  canal_base <- serie_semanal |>
    mutate(semana_canal = pmin(as.integer(semana), 52L)) |>
    group_by(region_ema, semana_canal) |>
    summarise(
      n_anios = n(),
      media_log = mean(log1p(casos)),
      sd_log = sd(log1p(casos)),
      .groups = "drop"
    ) |>
    mutate(
      # Se usa z = 1,96 (aproximación normal) para evitar intervalos
      # inestables en semanas con pocos años de historia.
      # Si solo hay un año, sd() es NA; se usa un IC de ancho cero.
      se_log = if_else(is.na(sd_log), 0, sd_log / sqrt(n_anios)),
      media_geometrica = expm1(media_log),
      ic_inferior = pmax(expm1(media_log - z * se_log), 0),
      ic_superior = expm1(media_log + z * se_log),
      # El orden de las bandas nunca debe invertirse.
      ic_inferior = pmin(ic_inferior, media_geometrica),
      ic_superior = pmax(ic_superior, media_geometrica)
    ) |>
    select(-media_log, -sd_log, -se_log)

  bind_rows(
    canal_base |> mutate(semana = semana_canal),
    canal_base |>
      filter(semana_canal == 52L) |>
      mutate(semana = 53L)
  ) |>
    arrange(region_ema, semana)
}

message("Calculando el canal endémico de Bortman...")

canal <- calcular_canal_bortman(serie_semanal)

techo_region <- serie_semanal |>
  group_by(region_ema) |>
  summarise(max_casos = max(casos, na.rm = TRUE), .groups = "drop")

canal <- canal |>
  left_join(techo_region, by = "region_ema") |>
  mutate(
    # Techo de cada panel: máximo observado o límite superior, con holgura.
    techo = pmax(max_casos, ic_superior, 1) * 1.05,
    # No se dibuja zona de brote en regiones sin casos (el canal es 0).
    brote_ymax = if_else(max_casos > 0, techo, NA_real_)
  )

anio_actual <- max(serie_semanal$anio, na.rm = TRUE)
anios_historicos <- sort(unique(serie_semanal$anio))
anios_historicos <- anios_historicos[anios_historicos < anio_actual]

observado <- serie_semanal |>
  mutate(
    grupo = if_else(anio == anio_actual, "actual", "historico")
  )

# Regiones sin casos en todo el período (Bogotá D.C. por diseño del análisis).
sin_casos <- canal |>
  group_by(region_ema) |>
  summarise(
    max_casos = max(max_casos, na.rm = TRUE),
    techo = max(techo, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(max_casos == 0) |>
  mutate(
    x = 27,
    y = techo / 2,
    label = "Sin casos notificados"
  )

caption_canal <- c(
  paste0(
    "Fuente: SIVIGILA/INS, serie municipal agregada por región EMA (Script 10). ",
    "Canal estimado con ", min(serie_semanal$anio), "-",
    max(serie_semanal$anio), "."
  ),
  paste(
    "Zonas: éxito ≤ IC inferior; seguridad hasta la media;",
    "alerta hasta el IC superior; brote > IC superior."
  ),
  paste(
    "La semana 53 se estima junto con la 52 porque solo 2020 y 2025",
    "tuvieron 53 semanas."
  )
)

if ("Bogotá D.C." %in% as.character(sin_casos$region_ema)) {
  caption_canal <- c(
    caption_canal,
    paste(
      "Bogotá D.C. no registra casos de dengue por diseño del análisis",
      "(placebo de enfermedad sin miedo)."
    )
  )
}

p_canal <- ggplot(canal, aes(x = semana)) +
  # Fondo neutro para las regiones sin casos.
  geom_rect(
    data = sin_casos,
    aes(xmin = 1, xmax = 53, ymin = 0, ymax = techo),
    fill = "grey93",
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    aes(ymin = 0, ymax = ic_inferior, fill = "Éxito"),
    alpha = 0.45
  ) +
  geom_ribbon(
    aes(ymin = ic_inferior, ymax = media_geometrica, fill = "Seguridad"),
    alpha = 0.45
  ) +
  geom_ribbon(
    aes(ymin = media_geometrica, ymax = ic_superior, fill = "Alerta"),
    alpha = 0.45
  ) +
  geom_ribbon(
    data = filter(canal, !is.na(brote_ymax)),
    aes(ymin = ic_superior, ymax = brote_ymax, fill = "Brote"),
    alpha = 0.45
  ) +
  geom_line(
    data = filter(observado, grupo == "historico"),
    aes(x = semana, y = casos, group = anio),
    color = "grey55",
    linewidth = 0.25,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = filter(observado, grupo == "actual"),
    aes(x = semana, y = casos, group = 1),
    color = "black",
    linewidth = 0.55,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = sin_casos,
    aes(x = x, y = y, label = label),
    size = 3,
    color = "grey35",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ region_ema, ncol = 4, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Éxito" = "#66BB6A",
      "Seguridad" = "#FFD54F",
      "Alerta" = "#FF9800",
      "Brote" = "#E53935"
    ),
    breaks = c("Éxito", "Seguridad", "Alerta", "Brote"),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = seq(1, 53, by = 8),
    expand = c(0, 0)
  ) +
  labs(
    title = "Canal endémico de dengue por regiones EMA (método de Bortman)",
    subtitle = paste0(
      "Media geométrica semanal (log1p) e IC 95% (",
      min(serie_semanal$anio), "-", max(serie_semanal$anio),
      "). Línea negra: ", anio_actual,
      "; líneas grises: ", min(anios_historicos), "-", max(anios_historicos), "."
    ),
    x = "Semana epidemiológica",
    y = "Casos semanales",
    caption = paste(caption_canal, collapse = "\n")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.1, "cm"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.9, "lines"),
    strip.text = element_text(face = "bold", size = 8, color = "grey15"),
    axis.text = element_text(size = 7, color = "grey30"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35", size = 9),
    plot.caption = element_text(hjust = 0, color = "grey45", size = 7.5),
    plot.margin = margin(10, 14, 8, 10)
  )

ruta_canal <- here("graficas", "11_canal_endemico.png")

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = ruta_canal,
    plot = p_canal,
    width = 16,
    height = 11,
    dpi = 300,
    bg = "white",
    device = ragg::agg_png
  )
} else {
  ggsave(
    filename = ruta_canal,
    plot = p_canal,
    width = 16,
    height = 11,
    dpi = 300,
    bg = "white"
  )
}

message(
  "Canal endémico calculado para ", n_distinct(canal$region_ema),
  " regiones y ", n_distinct(canal$semana), " semanas."
)
message("Gráfico guardado en: ", ruta_canal)

# --- 4. Verificación de tidyquant ---------------------------------

if (!requireNamespace("tidyquant", quietly = TRUE)) {
  message("")
  message("============================================================")
  message("El paquete 'tidyquant' no está instalado.")
  message("Para generar el análisis técnico, instálelo en R con:")
  message('  install.packages("tidyquant")')
  message("Después vuelva a ejecutar este script.")
  message("El canal endémico ya quedó guardado en: ", ruta_canal)
  message("============================================================")
  message("")
  stop(
    "Script detenido: falta el paquete 'tidyquant'.",
    call. = FALSE
  )
}

# --- 5. Top 6 de regiones para el análisis técnico ----------------

top6 <- ranking |>
  arrange(posicion) |>
  slice_head(n = 6) |>
  pull(region_ema)

top6 <- intersect(top6, as.character(serie_semanal$region_ema))

if (length(top6) < 6) {
  stop(
    "No se pudieron identificar 6 regiones del ranking en la serie. ",
    "Revise '", basename(ruta_ranking), "'.",
    call. = FALSE
  )
}

regiones_omitidas <- niveles_region[!niveles_region %in% top6]

message("")
message("Regiones incluidas en el análisis técnico (Top 6): ",
        paste(top6, collapse = ", "))
message("Regiones omitidas por espacio: ",
        paste(regiones_omitidas, collapse = ", "))

# --- 6. Serie semanal e indicadores técnicos ----------------------

# Calcula los indicadores sobre la serie semanal de cada región. Se usa TTR
# (el motor de quantmod) y se mantiene el orden temporal de fecha_sem.
calcular_indicadores <- function(datos) {
  macd_mat <- TTR::MACD(
    datos$casos,
    nFast = 12,
    nSlow = 26,
    nSig = 9
  )

  datos |>
    mutate(
      sma_4 = TTR::SMA(casos, n = 4),
      ema_12 = TTR::EMA(casos, n = 12),
      bb_media = TTR::SMA(casos, n = 20),
      bb_sd = TTR::runSD(casos, n = 20),
      bb_sup = bb_media + 2 * bb_sd,
      bb_inf = pmax(bb_media - 2 * bb_sd, 0),
      macd = as.numeric(macd_mat[, "macd"]),
      macd_senal = as.numeric(macd_mat[, "signal"])
    )
}

message("Calculando indicadores técnicos semanales...")

serie_tecnica <- serie_semanal |>
  filter(region_ema %in% top6) |>
  arrange(region_ema, fecha_sem) |>
  group_by(region_ema) |>
  group_modify(~ calcular_indicadores(.x)) |>
  ungroup()

# --- 7. OHLC mensual y bandas al cierre de mes --------------------

message("Construyendo velas mensuales (OHLC)...")

ohlc_mensual <- serie_semanal |>
  filter(region_ema %in% top6) |>
  arrange(region_ema, fecha_sem) |>
  mutate(fecha_mes = as.Date(format(fecha_sem, "%Y-%m-01"))) |>
  group_by(region_ema, fecha_mes) |>
  summarise(
    open = first(casos),
    high = max(casos),
    low = min(casos),
    close = last(casos),
    .groups = "drop"
  ) |>
  ungroup() |>
  mutate(region_ema = factor(region_ema, levels = top6))

if (anyNA(ohlc_mensual$open) || anyNA(ohlc_mensual$close)) {
  stop("Hay meses sin datos para construir las velas OHLC.", call. = FALSE)
}

# tidyquant::geom_candlestick() usa un eje x numérico con ancho de 0.9
# unidades. Se crea un índice de mes global (1..72) para que las velas
# mensuales tengan el ancho correcto y compartan el eje entre facetas.
meses <- sort(unique(ohlc_mensual$fecha_mes))
ohlc_mensual <- ohlc_mensual |>
  mutate(x_mes = match(fecha_mes, meses))

bandas_mes <- serie_tecnica |>
  mutate(fecha_mes = as.Date(format(fecha_sem, "%Y-%m-01"))) |>
  group_by(region_ema, fecha_mes) |>
  slice_max(fecha_sem, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    region_ema = factor(region_ema, levels = top6),
    x_mes = match(fecha_mes, meses)
  ) |>
  filter(!is.na(bb_sup), !is.na(bb_inf))

meses_enero <- which(format(meses, "%m") == "01")
etiquetas_anio <- format(meses[meses_enero], "%Y")

# --- 8. Gráfico de velas y bandas de Bollinger --------------------

message("Generando gráfico de análisis técnico...")

p_tecnico <- ggplot(
  ohlc_mensual,
  aes(x = x_mes, y = close)
) +
  geom_ribbon(
    data = bandas_mes,
    aes(x = x_mes, ymin = bb_inf, ymax = bb_sup),
    inherit.aes = FALSE,
    fill = "#1565C0",
    alpha = 0.12
  ) +
  tidyquant::geom_candlestick(
    aes(open = open, high = high, low = low, close = close),
    colour_up = "#2E7D32",
    colour_down = "#C62828",
    fill_up = "#2E7D32",
    fill_down = "#C62828",
    linewidth = 0.3
  ) +
  geom_line(
    data = bandas_mes,
    aes(x = x_mes, y = bb_sup),
    inherit.aes = FALSE,
    color = "#1565C0",
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_line(
    data = bandas_mes,
    aes(x = x_mes, y = bb_inf),
    inherit.aes = FALSE,
    color = "#1565C0",
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_line(
    data = bandas_mes,
    aes(x = x_mes, y = bb_media),
    inherit.aes = FALSE,
    color = "#1565C0",
    linewidth = 0.35
  ) +
  facet_wrap(~ region_ema, ncol = 3, scales = "free_y") +
  scale_x_continuous(
    breaks = meses_enero,
    labels = etiquetas_anio,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = label_number(
      big.mark = ".",
      decimal.mark = ",",
      accuracy = 1
    ),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Análisis técnico de casos de dengue: velas mensuales y bandas de Bollinger",
    subtitle = paste0(
      "Top 6 regiones EMA por tasa acumulada (Script 10). ",
      "Bandas calculadas sobre la serie semanal (20 semanas)."
    ),
    x = NULL,
    y = "Casos",
    caption = paste(
      "Velas: Open/High/Low/Close mensual de los casos semanales.",
      "Verde: cierre ≥ apertura; rojo: cierre < apertura.",
      "Bandas de Bollinger: 20 semanas, ±2 desviaciones estándar, al cierre de cada mes.",
      "SMA(4), EMA(12) y MACD(12,26,9) se calculan sobre la serie semanal y se reportan en consola.",
      "Se omiten por espacio las 6 regiones restantes.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.text = element_text(size = 7, color = "grey30"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35", size = 9),
    plot.caption = element_text(hjust = 0, color = "grey45", size = 7.2),
    plot.margin = margin(10, 14, 8, 10)
  )

ruta_tecnico <- here("graficas", "11_analisis_tecnico.png")

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = ruta_tecnico,
    plot = p_tecnico,
    width = 15,
    height = 8.5,
    dpi = 300,
    bg = "white",
    device = ragg::agg_png
  )
} else {
  ggsave(
    filename = ruta_tecnico,
    plot = p_tecnico,
    width = 15,
    height = 8.5,
    dpi = 300,
    bg = "white"
  )
}

message("Gráfico guardado en: ", ruta_tecnico)

# --- 9. Resumen de indicadores en consola -------------------------

resumen_tecnico <- serie_tecnica |>
  group_by(region_ema) |>
  slice_max(fecha_sem, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(
    region_ema, fecha_sem, casos,
    sma_4, ema_12,
    bb_inf, bb_media, bb_sup,
    macd, macd_senal
  ) |>
  mutate(
    across(
      c(sma_4, ema_12, bb_inf, bb_media, bb_sup, macd, macd_senal),
      ~ round(.x, 2)
    )
  )

cat("\n--- Últimos indicadores técnicos por región ---\n")
print(
  as.data.frame(resumen_tecnico),
  row.names = FALSE
)

cat("\nRegiones omitidas por espacio (fuera del Top 6 del ranking):\n")
cat(paste0("  - ", regiones_omitidas, collapse = "\n"), "\n")

message("✅ Script 11 finalizado.")
