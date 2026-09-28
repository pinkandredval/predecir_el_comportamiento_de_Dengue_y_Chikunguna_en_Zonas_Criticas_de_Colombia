# ============================================================
# 04_descriptivo.R
# Análisis descriptivo y exploratorio (solo consola)
# ============================================================

library(tidyverse)
library(here)
library(janitor)
library(scales)   # para comma() y label_number()

# --- Cargar datos preparados ---
dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

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
print(summary(dengue %>% select(edad, semana, anio)))

cat("\n===== SUMMARY CHIKUNGUNYA (numéricas) =====\n")
print(summary(chikungunya %>% select(edad, semana, anio)))

resumen_var <- function(df, var, nombre_evento) {
  x <- df[[var]]
  data.frame(
    evento   = nombre_evento,
    variable = var,
    n        = sum(!is.na(x)),
    media    = round(mean(x, na.rm = TRUE), 2),
    mediana  = median(x, na.rm = TRUE),
    sd       = round(sd(x, na.rm = TRUE), 2),
    min      = min(x, na.rm = TRUE),
    q1       = quantile(x, 0.25, na.rm = TRUE),
    q3       = quantile(x, 0.75, na.rm = TRUE),
    max      = max(x, na.rm = TRUE),
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
  dengue      %>% select(all_of(columnas_comunes)),
  chikungunya %>% select(all_of(columnas_comunes))
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

# --- 6.1 Histograma de edad (facetado por evento) ---
ggplot(datos_comp, aes(x = edad, fill = evento)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  facet_wrap(~ evento, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Distribución de edad por evento",
       x = "Edad", y = "Frecuencia") +
  theme_minimal() +
  theme(legend.position = "none")

# --- 6.2 Boxplot edad × sexo (por evento) ---
ggplot(datos_comp, aes(x = sexo, y = edad, fill = sexo)) +
  geom_boxplot() +
  facet_wrap(~ evento, scales = "free_y") +
  labs(title = "Edad por sexo y evento",
       x = "Sexo", y = "Edad") +
  theme_minimal() +
  theme(legend.position = "none")

# --- 6.3 Barras: tipo de caso (por evento, escala libre) ---
datos_comp %>%
  count(evento, tip_cas) %>%
  ggplot(aes(x = tip_cas, y = n, fill = evento)) +
  geom_col() +
  facet_wrap(~ evento, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Tipo de caso por evento",
       x = "Código", y = "Casos") +
  theme_minimal() +
  theme(legend.position = "none")

# --- 6.4 Barras: sexo (por evento, escala libre) ---
datos_comp %>%
  count(evento, sexo) %>%
  ggplot(aes(x = sexo, y = n, fill = evento)) +
  geom_col() +
  facet_wrap(~ evento, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Distribución por sexo",
       x = "Sexo", y = "Casos") +
  theme_minimal() +
  theme(legend.position = "none")

# --- 6.5 Barras: hospitalización (por evento, escala libre) ---
datos_comp %>%
  count(evento, pac_hos) %>%
  ggplot(aes(x = pac_hos, y = n, fill = evento)) +
  geom_col() +
  facet_wrap(~ evento, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Hospitalización por evento",
       x = "1 = Sí, 2 = No", y = "Casos") +
  theme_minimal() +
  theme(legend.position = "none")

# --- 6.6 Serie temporal por semana (facetado por evento y año) ---
datos_comp %>%
  count(anio, semana, evento) %>%
  ggplot(aes(x = semana, y = n, color = evento)) +
  geom_line(linewidth = 0.7) +
  facet_grid(evento ~ anio, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Casos por semana epidemiológica, año y evento",
       x = "Semana", y = "Casos") +
  theme_minimal() +
  theme(legend.position = "none")

# --- 6.7 Serie temporal con eje X y Y libres por evento ---
datos_comp %>%
  count(fecha_epi, evento) %>%
  ggplot(aes(x = fecha_epi, y = n, color = evento)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ evento, scales = "free", ncol = 1) +
  scale_y_continuous(labels = comma) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(title = "Serie temporal por evento (escala libre)",
       x = "Fecha", y = "Casos") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside"
  )


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