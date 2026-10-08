# ============================================================
# 04_descriptivo.R
# Análisis descriptivo y exploratorio
# ============================================================

library(tidyverse)
library(here)
library(janitor)
library(scales)

# --- Evitar que MASS::select / MASS::filter enmascaren dplyr ---
select <- dplyr::select
filter <- dplyr::filter

# --- Cargar datos preparados ---
dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

# ============================================================
# 0. UTILIDADES (nombres de departamentos DIVIPOLA)
# ============================================================

nombrar_dpto <- function(codigo) {
  case_when(
    codigo == "05" ~ "Antioquia",
    codigo == "08" ~ "Atlántico",
    codigo == "11" ~ "Bogotá",
    codigo == "13" ~ "Bolívar",
    codigo == "15" ~ "Boyacá",
    codigo == "17" ~ "Caldas",
    codigo == "18" ~ "Caquetá",
    codigo == "19" ~ "Cauca",
    codigo == "20" ~ "Cesar",
    codigo == "23" ~ "Córdoba",
    codigo == "25" ~ "Cundinamarca",
    codigo == "27" ~ "Chocó",
    codigo == "41" ~ "Huila",
    codigo == "44" ~ "La Guajira",
    codigo == "47" ~ "Magdalena",
    codigo == "50" ~ "Meta",
    codigo == "52" ~ "Nariño",
    codigo == "54" ~ "Norte de Santander",
    codigo == "63" ~ "Quindío",
    codigo == "66" ~ "Risaralda",
    codigo == "68" ~ "Santander",
    codigo == "70" ~ "Sucre",
    codigo == "73" ~ "Tolima",
    codigo == "76" ~ "Valle del Cauca",
    codigo == "81" ~ "Arauca",
    codigo == "85" ~ "Casanare",
    codigo == "86" ~ "Putumayo",
    codigo == "88" ~ "San Andrés",
    codigo == "91" ~ "Amazonas",
    codigo == "94" ~ "Guainía",
    codigo == "95" ~ "Guaviare",
    codigo == "97" ~ "Vaupés",
    codigo == "99" ~ "Vichada",
    TRUE ~ paste0("Cod. ", codigo)
  )
}

# ============================================================
# 1. DIMENSIONES
# ============================================================

cat("\n===== DIMENSIONES =====\n")
cat("Dengue:", nrow(dengue), "filas x", ncol(dengue), "columnas\n")
cat("Chikungunya:", nrow(chikungunya), "filas x", ncol(chikungunya), "columnas\n")

# ============================================================
# 2. RESUMEN NUMÉRICO
# ============================================================

cat("\n===== SUMMARY DENGUE (numéricas) =====\n")
print(summary(dengue[, c("edad", "semana", "anio")]))

cat("\n===== SUMMARY CHIKUNGUNYA (numéricas) =====\n")
print(summary(chikungunya[, c("edad", "semana", "anio")]))

resumen_var <- function(df, var, nombre_evento) {
  x <- df[[var]]
  data.frame(
    evento    = nombre_evento,
    variable  = var,
    n         = sum(!is.na(x)),
    media     = round(mean(x, na.rm = TRUE), 2),
    mediana   = median(x, na.rm = TRUE),
    sd        = round(sd(x, na.rm = TRUE), 2),
    min       = min(x, na.rm = TRUE),
    q1        = quantile(x, 0.25, na.rm = TRUE),
    q3        = quantile(x, 0.75, na.rm = TRUE),
    max       = max(x, na.rm = TRUE),
    row.names = NULL
  )
}

resumen_num <- bind_rows(
  resumen_var(dengue,      "edad",   "Dengue"),
  resumen_var(dengue,      "semana", "Dengue"),
  resumen_var(chikungunya, "edad",   "Chikungunya"),
  resumen_var(chikungunya, "semana", "Chikungunya")
)

cat("\n===== RESUMEN NUMÉRICO DETALLADO =====\n")
print(resumen_num)

# ============================================================
# 3. FRECUENCIAS DE CATEGÓRICAS
# ============================================================

vars_cat <- c("sexo", "tip_ss", "area", "tip_cas", "pac_hos",
              "con_fin", "ajuste", "per_etn", "estrato", "uni_med")

frecuencias <- function(df, vars, nombre) {
  cat("\n========================================\n")
  cat("FRECUENCIAS:", nombre, "\n")
  cat("========================================\n")
  for (v in vars) {
    if (v %in% names(df)) {
      cat("\n---", v, "---\n")
      print(df %>% count(.data[[v]], sort = TRUE))
    }
  }
}

frecuencias(dengue, vars_cat, "DENGUE")
frecuencias(chikungunya, vars_cat, "CHIKUNGUNYA")

# ============================================================
# 4. RESUMEN POR EVENTO (comparativo)
# ============================================================

columnas_comunes <- intersect(names(dengue), names(chikungunya))

datos_comp <- bind_rows(
  dengue[, columnas_comunes],
  chikungunya[, columnas_comunes]
)

cat("\n===== RESUMEN COMPARATIVO POR EVENTO =====\n")
print(
  datos_comp %>%
    group_by(evento) %>%
    summarise(
      n            = n(),
      edad_media   = round(mean(edad, na.rm = TRUE), 1),
      edad_mediana = median(edad, na.rm = TRUE),
      prop_fem     = round(mean(sexo == "F", na.rm = TRUE) * 100, 1),
      prop_hosp    = round(mean(pac_hos == "1", na.rm = TRUE) * 100, 1),
      .groups = "drop"
    ) %>%
    as.data.frame()
)

# ============================================================
# 5. TABLAS CRUZADAS
# ============================================================

cat("\n===== DENGUE: SEXO x TIPO DE CASO =====\n")
print(dengue %>% tabyl(sexo, tip_cas))

cat("\n===== DENGUE: SEXO x HOSPITALIZACIÓN =====\n")
print(dengue %>% tabyl(sexo, pac_hos))

cat("\n===== CHIKUNGUNYA: SEXO x TIPO DE CASO =====\n")
print(chikungunya %>% tabyl(sexo, tip_cas))

cat("\n===== CHIKUNGUNYA: SEXO x HOSPITALIZACIÓN =====\n")
print(chikungunya %>% tabyl(sexo, pac_hos))

# ============================================================
# 6. GRÁFICOS
# ============================================================

# --- 6.1 Distribución de edad ---
datos_comp <- datos_comp %>%
  mutate(grupo_edad = cut(edad,
                          breaks = seq(0, 125, by = 5),
                          right = FALSE,
                          labels = paste0(seq(0, 120, by = 5), "-",
                                          seq(4, 124, by = 5))))

edad_barras <- datos_comp %>%
  count(evento, grupo_edad) %>%
  group_by(evento) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup()

ggplot(edad_barras, aes(x = grupo_edad, y = pct, fill = evento)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = ifelse(pct >= 2, paste0(pct, "%"), "")),
            vjust = -0.3, size = 2.8) +
  facet_wrap(~ evento, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Dengue" = "#08519C",
                               "Chikungunya" = "#A50F15")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Distribución porcentual por grupos de edad",
       x = "Grupo de edad (años)", y = "% de casos") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# --- 6.2 Boxplot edad × sexo ---
ggplot(datos_comp, aes(x = sexo, y = edad, fill = sexo)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.5, outlier.alpha = 0.4) +
  facet_wrap(~ evento, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#74A9CF", "M" = "#0570B0")) +
  labs(title = "Edad por sexo y evento",
       x = "Sexo", y = "Edad (años)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 6.3 Tipo de caso (corregido: unificar 2 y 4) ---
datos_comp %>%
  mutate(
    tip_cas_label = case_when(
      tip_cas == "1" ~ "1. Sospechoso",
      tip_cas %in% c("2", "4") ~ "2/4. Conf. clínica",
      tip_cas == "3" ~ "3. Conf. laboratorio",
      tip_cas == "5" ~ "5. Otro",
      TRUE ~ "Sin dato"
    )
  ) %>%
  count(evento, tip_cas_label) %>%
  group_by(evento) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  ggplot(aes(x = tip_cas_label, y = pct, fill = evento)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.3, size = 3.5) +
  facet_wrap(~ evento, scales = "free_y") +
  scale_fill_manual(values = c("Dengue" = "#3182BD",
                               "Chikungunya" = "#E6550D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Tipo de caso por evento",
       x = "Tipo de caso", y = "% de casos") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

# --- 6.4 Sexo ---
datos_comp %>%
  count(evento, sexo) %>%
  group_by(evento) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  ggplot(aes(x = sexo, y = pct, fill = evento)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.3, size = 4) +
  facet_wrap(~ evento, scales = "free_y") +
  scale_fill_manual(values = c("Dengue" = "#6BAED6",
                               "Chikungunya" = "#FD8D3C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Distribución porcentual por sexo",
       x = "Sexo", y = "% de casos") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 6.5 Hospitalización ---
datos_comp %>%
  count(evento, pac_hos) %>%
  group_by(evento) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  ggplot(aes(x = pac_hos, y = pct, fill = evento)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.3, size = 4) +
  facet_wrap(~ evento, scales = "free_y") +
  scale_fill_manual(values = c("Dengue" = "#4292C6",
                               "Chikungunya" = "#E31A1C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Hospitalización por evento",
       x = "¿Hospitalizado? (1=Sí, 2=No)", y = "% de casos") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# --- 6.6 Serie temporal Dengue (log1p) ---
serie_dengue <- datos_comp %>%
  filter(evento == "Dengue") %>%
  count(anio, semana)

ggplot(serie_dengue, aes(x = semana, y = n, color = factor(anio))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("2020" = "#08519C", "2021" = "#2171B5",
               "2022" = "#4292C6", "2023" = "#6BAED6",
               "2024" = "#9ECAE1", "2025" = "#C6DBEF"),
    name = "Año"
  ) +
  scale_x_continuous(breaks = seq(0, 52, by = 4)) +
  scale_y_continuous(trans = "log1p", labels = comma) +
  labs(title = "Dengue: casos por semana epidemiológica (2020-2025)",
       subtitle = "Cada línea representa un año (escala log)",
       x = "Semana epidemiológica", y = "Casos (log)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# --- 6.6 Serie temporal Chikungunya (log1p) ---
serie_chik <- datos_comp %>%
  filter(evento == "Chikungunya") %>%
  count(anio, semana)

ggplot(serie_chik, aes(x = semana, y = n, color = factor(anio))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("2020" = "#A50F15", "2021" = "#CB181D",
               "2022" = "#EF3B2C", "2023" = "#FB6A4A",
               "2024" = "#FC9272", "2025" = "#FCBBA1"),
    name = "Año"
  ) +
  scale_x_continuous(breaks = seq(0, 52, by = 4)) +
  scale_y_continuous(trans = "log1p", labels = comma) +
  labs(title = "Chikungunya: casos por semana epidemiológica (2020-2025)",
       subtitle = "Cada línea representa un año (escala log)",
       x = "Semana epidemiológica", y = "Casos (log)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# --- 6.7 Serie temporal continua ---
colores_evento <- c("Dengue" = "#08519C", "Chikungunya" = "#A50F15")
rango_fechas <- range(datos_comp$fecha_epi, na.rm = TRUE)

datos_comp %>%
  count(fecha_epi, evento) %>%
  ggplot(aes(x = fecha_epi, y = n, color = evento)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ evento, scales = "free", ncol = 1) +
  scale_color_manual(values = colores_evento) +
  scale_y_continuous(labels = comma) +
  scale_x_date(limits = rango_fechas, date_labels = "%Y",
               date_breaks = "1 year") +
  labs(title = "Serie temporal por evento (escala libre)",
       x = "Fecha", y = "Casos") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# 6.8 DISTRIBUCIÓN POR DEPARTAMENTO (TODOS)
# ============================================================

# Detectar columna de departamento
col_dpto <- intersect(
  c("cod_dpto_o", "cod_dpto_r", "cod_dpto_n",
    "departamento_ocurrencia", "departamento_residencia",
    "departamento_notificacion"),
  names(datos_comp)
)[1]

if (!is.na(col_dpto)) {
  
  cat("\n")
  cat("=========================================================\n")
  cat("DISTRIBUCIÓN POR DEPARTAMENTO (TODOS)\n")
  cat("=========================================================\n")
  cat("Columna usada:", col_dpto, "\n\n")
  
  # --- Resumen por departamento y evento (TODOS) ---
  dptos_todos <- datos_comp %>%
    filter(!is.na(.data[[col_dpto]])) %>%
    mutate(dpto = nombrar_dpto(as.character(.data[[col_dpto]]))) %>%
    count(evento, dpto, name = "casos") %>%
    group_by(evento) %>%
    mutate(pct = round(casos / sum(casos) * 100, 2)) %>%
    ungroup()
  
  # --- Gráfico: TODOS los departamentos, escala libre por evento ---
  print(
    ggplot(dptos_todos,
           aes(x = reorder(dpto, casos), y = casos, fill = evento)) +
      geom_col(alpha = 0.9) +
      geom_text(aes(label = comma(casos)),
                hjust = -0.1, size = 2.3) +
      coord_flip() +
      facet_wrap(~ evento, scales = "free", ncol = 2) +
      scale_fill_manual(values = c("Dengue" = "#08519C",
                                   "Chikungunya" = "#A50F15")) +
      scale_y_continuous(
        labels = comma,
        expand = expansion(mult = c(0, 0.15))
      ) +
      labs(
        title    = "Casos por departamento y evento",
        subtitle = "Todos los departamentos. Escala independiente por panel.",
        x        = NULL,
        y        = "Casos"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        axis.text.y     = element_text(size = 7),
        panel.grid.major.y = element_blank()
      )
  )
  
  # --- Tabla completa (todos los departamentos por evento) ---
  cat("\n--- Tabla completa por departamento (ordenada por casos) ---\n")
  tabla_dptos_full <- dptos_todos %>%
    arrange(evento, desc(casos)) %>%
    select(evento, dpto, casos, pct)
  
  print(as.data.frame(tabla_dptos_full), row.names = FALSE)
  
  # --- Guardar tabla ---
  write.csv(tabla_dptos_full,
            here("data", "processed", "casos_por_departamento.csv"),
            row.names = FALSE)
  
} else {
  cat("\n⚠️  No se encontró columna de departamento en datos_comp\n")
  cat("Columnas con 'dpto' o 'departamento':\n")
  print(grep("dpto|departamento", names(datos_comp), value = TRUE))
}

# ============================================================
# 7. ATÍPICOS EN EDAD
# ============================================================

outliers_edad <- function(df, nombre) {
  df %>%
    mutate(
      q1 = quantile(edad, 0.25, na.rm = TRUE),
      q3 = quantile(edad, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      atipico = edad < (q1 - 1.5 * iqr) | edad > (q3 + 1.5 * iqr)
    ) %>%
    filter(atipico) %>%
    count(edad, sort = TRUE) %>%
    mutate(evento = nombre)
}

cat("\n===== ATÍPICOS EN EDAD - DENGUE =====\n")
print(as.data.frame(outliers_edad(dengue, "Dengue")))

cat("\n===== ATÍPICOS EN EDAD - CHIKUNGUNYA =====\n")
print(as.data.frame(outliers_edad(chikungunya, "Chikungunya")))

cat("\n✅ Análisis descriptivo y exploratorio completado.\n")

