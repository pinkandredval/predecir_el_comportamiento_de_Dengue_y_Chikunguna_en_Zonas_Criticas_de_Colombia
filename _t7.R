options(width = 160)
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(lubridate); library(scales)
})

cat("== select resuelto ==\n")
print(find("select"))
cat("\n== conflictos ==\n")
print(conflicts())

semanas_por_anio <- tibble(
  anio = 2020:2025,
  n_semanas = c(53L, 52L, 52L, 52L, 52L, 53L)
)

calendario_semanal <- expand_grid(anio = 2020:2025, semana = 1:53) %>%
  left_join(semanas_por_anio, by = "anio") %>%
  filter(semana <= n_semanas) %>%
  mutate(fecha_sem = as.Date(paste0(anio, "-01-01")) + (semana - 1) * 7) %>%
  select(anio, semana, fecha_sem) %>%
  arrange(anio, semana)

cat("\n== filas:", nrow(calendario_semanal), "\n")
print(head(calendario_semanal, 3))
print(tail(calendario_semanal, 3))

cat("\n== con MASS cargado despues ==\n")
library(MASS)
print(find("select"))
r <- tryCatch(select(calendario_semanal, anio), error = function(e) conditionMessage(e))
print(r)
