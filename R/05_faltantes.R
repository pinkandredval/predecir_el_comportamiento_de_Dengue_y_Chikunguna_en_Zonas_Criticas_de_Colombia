# ============================================================
# 05_faltantes.R (Versión completa + diagnóstico pre-imputación)
# ============================================================

library(tidyverse)
library(here)
library(janitor)
library(naniar)
library(mice)
library(scales)

# --- Cargar datos ---
dengue      <- readRDS(here("data", "processed", "dengue_prep.rds"))
chikungunya <- readRDS(here("data", "processed", "chikungunya_prep.rds"))

# ============================================================
# 0. VERIFICACIÓN DE NA ESTRUCTURALES
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("0. VERIFICACIÓN DE NA ESTRUCTURALES\n")
cat("=========================================================\n")

verificar_estructurales <- function(df, nombre) {
  cat("\n---", nombre, "---\n")
  cols_est <- names(df)[grepl("^na_estructural_", names(df))]
  
  if (length(cols_est) == 0) {
    cat("No hay columnas na_estructural_\n")
    return(invisible(NULL))
  }
  
  df %>%
    summarise(across(all_of(cols_est), ~ sum(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "variable_estructural", values_to = "n_estructural") %>%
    mutate(
      variable_original = sub("^na_estructural_", "", variable_estructural),
      pct = round(n_estructural / nrow(df) * 100, 2)
    ) %>%
    select(-variable_estructural) %>%
    arrange(desc(pct)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

verificar_estructurales(dengue, "DENGUE")
verificar_estructurales(chikungunya, "CHIKUNGUNYA")

# ============================================================
# 1. FUNCIONES AUXILIARES Y CATEGORÍAS
# ============================================================

clasificar_celdas <- function(df) {
  df_analisis <- df %>% select(-starts_with("na_estructural_"))
  vars_analisis <- names(df_analisis)
  
  estados <- df_analisis %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "faltante", "presente")))
  
  for (v in vars_analisis) {
    col_est <- paste0("na_estructural_", v)
    if (col_est %in% names(df)) {
      idx <- which(!is.na(df[[col_est]]) & df[[col_est]])
      if (length(idx) > 0) estados[[v]][idx] <- "no_aplica"
    }
  }
  estados
}

categorias_vars <- list(
  "Identificación"        = c("consecutive", "cod_eve", "nombre_evento", "evento"),
  "Temporales"            = c("fec_not", "semana", "ano", "anio", "fecha_epi", "ini_sin", "fec_con"),
  "Sociodemográficas"     = c("edad", "uni_med", "sexo", "per_etn", "estrato", "ocupacion"),
  "Geográficas"           = c("cod_pais_o", "cod_dpto_o", "cod_mun_o", "cod_pais_r", "cod_dpto_r", "cod_mun_r", "cod_dpto_n", "cod_mun_n", "area"),
  "Seguridad social"      = c("tip_ss", "cod_ase"),
  "Grupos poblacionales"  = c("gru_pob", "nom_grupo", "gp_discapa", "gp_desplaz", "gp_migrant", "gp_carcela", "gp_gestan", "sem_ges", "gp_indigen", "gp_pobicfb", "gp_mad_com", "gp_desmovi", "gp_psiquia", "gp_vic_vio", "gp_otros"),
  "Clínicas"              = c("tip_cas", "pac_hos", "fec_hos", "con_fin", "fec_def", "ajuste", "fecha_nto", "fec_aju"),
  "Militares"             = c("fm_fuerza", "fm_unidad", "fm_grado"),
  "Administrativas"       = c("fuente", "confirmados", "va_sispro", "estado_final_de_caso", "nom_est_f_caso", "nombre_upgd", "pais_ocurrencia", "departamento_ocurrencia", "municipio_ocurrencia", "pais_residencia", "departamento_residencia", "municipio_residencia", "departamento_notificacion", "municipio_notificacion")
)

asignar_categoria <- function(variable, categorias = categorias_vars) {
  for (cat_nombre in names(categorias)) {
    if (variable %in% categorias[[cat_nombre]]) return(cat_nombre)
  }
  return("Otras")
}

# ============================================================
# 2. RESUMEN DE FALTANTES POR VARIABLE
# ============================================================

resumen_faltantes <- function(df, nombre) {
  df_analisis <- df %>% select(-starts_with("na_estructural_"))
  
  na_total <- df_analisis %>%
    summarise(across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_na_total")
  
  na_estructural <- sapply(names(df_analisis), function(v) {
    col_est <- paste0("na_estructural_", v)
    if (col_est %in% names(df)) sum(df[[col_est]], na.rm = TRUE) else 0
  })
  
  na_total %>%
    mutate(
      n_na_estructural = na_estructural,
      n_na_real        = n_na_total - n_na_estructural,
      pct_total        = round(n_na_total / nrow(df) * 100, 2),
      pct_estructural  = round(n_na_estructural / nrow(df) * 100, 2),
      pct_real         = round(n_na_real / nrow(df) * 100, 2),
      categoria        = sapply(variable, asignar_categoria),
      evento           = nombre
    ) %>%
    arrange(desc(pct_real), desc(pct_estructural))
}

falt_dengue <- resumen_faltantes(dengue, "Dengue")
falt_chik   <- resumen_faltantes(chikungunya, "Chikungunya")

# ============================================================
# 3. FUNCIONES DE VISUALIZACIÓN (Patrones y Mapas corregidos)
# ============================================================

patron_faltantes <- function(df, nombre_evento, categoria_filtro) {
  estados <- clasificar_celdas(df)
  vars_cat <- intersect(categorias_vars[[categoria_filtro]], names(estados))
  
  if (length(vars_cat) < 2) return(invisible(NULL))
  
  estados <- estados[, vars_cat, drop = FALSE]
  n_total_filas <- nrow(estados)
  
  estados_str <- estados %>%
    mutate(
      patron = do.call(paste0, lapply(., function(x) {
        ifelse(x == "presente", "P", ifelse(x == "faltante", "F", "N"))
      }))
    )
  
  tab <- sort(table(estados_str$patron), decreasing = TRUE)
  patrones <- data.frame(patron = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
  
  mat_pat <- do.call(rbind, strsplit(patrones$patron, ""))
  colnames(mat_pat) <- vars_cat
  mat_pat <- as.data.frame(mat_pat, stringsAsFactors = FALSE)
  
  mat_faltantes_real <- as.data.frame(lapply(mat_pat, function(x) as.integer(x == "F")))
  patrones$n_faltantes_fila <- rowSums(mat_faltantes_real)
  faltantes_por_var <- colSums(mat_faltantes_real * patrones$n)
  
  df_long <- cbind(mat_pat, n = patrones$n, fila = seq_len(nrow(patrones))) %>%
    pivot_longer(cols = all_of(vars_cat), names_to = "variable", values_to = "estado_letra") %>%
    mutate(
      variable = factor(variable, levels = vars_cat),
      fila = factor(fila, levels = rev(seq_len(nrow(patrones)))),
      estado = case_when(
        estado_letra == "P" ~ "Presente",
        estado_letra == "N" ~ "No aplica",
        estado_letra == "F" ~ "Faltante real"
      ),
      estado = factor(estado, levels = c("Presente", "No aplica", "Faltante real"))
    )
  
  etiq_fila_izq <- data.frame(fila = factor(seq_len(nrow(patrones)), levels = rev(seq_len(nrow(patrones)))), n = patrones$n)
  etiq_fila_der <- data.frame(fila = factor(seq_len(nrow(patrones)), levels = rev(seq_len(nrow(patrones)))), nf = patrones$n_faltantes_fila)
  etiq_col_abajo <- data.frame(variable = factor(vars_cat, levels = vars_cat), n_falt = faltantes_por_var)
  n_vars <- length(vars_cat)
  
  ggplot(df_long, aes(x = variable, y = fila, fill = estado)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(data = etiq_fila_izq, aes(x = 0.4, y = fila, label = format(n, big.mark = ",")), inherit.aes = FALSE, hjust = 1, size = 3, color = "black") +
    geom_text(data = etiq_fila_der, aes(x = n_vars + 0.6, y = fila, label = nf), inherit.aes = FALSE, hjust = 0, size = 3, color = "black", fontface = "bold") +
    geom_text(data = etiq_col_abajo, aes(x = variable, y = 0.3, label = n_falt), inherit.aes = FALSE, vjust = 1, size = 3, color = "black", fontface = "bold") +
    scale_fill_manual(
      values = c("Presente" = "#2166AC", "No aplica" = "#CCCCCC", "Faltante real" = "#B2182B"),
      name = NULL
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(df_long$fila))) +
    labs(title = paste0("Patrón de faltantes (3 estados) — ", nombre_evento),
         subtitle = paste0("Categoría: ", categoria_filtro, " | Total filas: ", format(n_total_filas, big.mark = ",")),
         x = NULL, y = NULL) +
    coord_cartesian(xlim = c(0.2, n_vars + 1), ylim = c(0.5, nrow(patrones) + 0.5), clip = "off") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 0),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      plot.margin = margin(10, 50, 20, 50)
    )
}

mapa_faltantes_completo <- function(df, nombre_evento, categoria_filtro) {
  estados <- clasificar_celdas(df)
  vars_cat <- intersect(categorias_vars[[categoria_filtro]], names(estados))
  
  if (length(vars_cat) == 0) return(invisible(NULL))
  
  estados <- estados[, vars_cat, drop = FALSE]
  n_total <- nrow(estados)
  
  resumen <- purrr::map_dfr(vars_cat, function(v) {
    x <- estados[[v]]
    tibble(
      variable  = v,
      presente  = sum(x == "presente",  na.rm = TRUE),
      faltante  = sum(x == "faltante",  na.rm = TRUE),
      no_aplica = sum(x == "no_aplica", na.rm = TRUE)
    )
  }) %>%
    mutate(
      pct_faltante  = faltante  / n_total * 100,
      pct_no_aplica = no_aplica / n_total * 100,
      pct_presente  = presente  / n_total * 100,
      etiqueta_txt = case_when(
        pct_faltante > 0     ~ paste0(round(pct_faltante, 1), "%"),
        pct_no_aplica == 100 ~ "100%",
        TRUE                 ~ "0%"
      ),
      color_txt = case_when(
        pct_faltante > 0     ~ "#B2182B",
        pct_no_aplica == 100 ~ "#666666",
        TRUE                 ~ "#2166AC"
      )
    ) %>%
    arrange(desc(pct_faltante), desc(pct_no_aplica))
  
  df_long <- resumen %>%
    select(variable, pct_presente, pct_faltante, pct_no_aplica) %>%
    pivot_longer(cols = starts_with("pct_"), names_to = "estado", values_to = "pct") %>%
    mutate(
      estado   = factor(sub("^pct_", "", estado), levels = c("presente", "faltante", "no_aplica")),
      variable = factor(variable, levels = rev(resumen$variable))
    )
  
  ggplot(df_long, aes(x = variable, y = pct, fill = estado)) +
    geom_col(width = 0.75) +
    geom_text(
      data = resumen,
      aes(x = variable, y = 102, label = etiqueta_txt, color = I(color_txt)),
      inherit.aes = FALSE, size = 2.8, fontface = "bold"
    ) +
    scale_fill_manual(
      values = c("presente" = "#2166AC", "faltante" = "#B2182B", "no_aplica" = "#CCCCCC"),
      labels = c("presente" = "Presente", "faltante" = "Faltante real", "no_aplica" = "No aplica"),
      name = NULL
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.05))) +
    coord_flip(clip = "off") +
    labs(
      title    = paste0("Mapa de faltantes — ", nombre_evento),
      subtitle = paste0("Categoría: ", categoria_filtro, "  |  n = ", format(n_total, big.mark = ","), " filas"),
      x = NULL, y = "% de registros"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 8),
      legend.position    = "top",
      plot.margin        = margin(5, 40, 5, 5)
    )
}

# ============================================================
# 4. EJECUCIÓN DEL BUCLE PRINCIPAL
# ============================================================

dir.create("graficas/faltantes", recursive = TRUE, showWarnings = FALSE)

for (evento in c("dengue", "chikungunya")) {
  df_evt <- get(evento)
  nom_evt <- tools::toTitleCase(evento)
  
  cat("\nProcesando evento:", nom_evt, "\n")
  
  for (cat_actual in names(categorias_vars)) {
    vars_cat <- intersect(categorias_vars[[cat_actual]], names(clasificar_celdas(df_evt)))
    n_vars_cat <- length(vars_cat)
    
    if (n_vars_cat < 2) next
    
    p1 <- patron_faltantes(df_evt, nom_evt, cat_actual)
    if (!is.null(p1)) {
      print(p1)
      ggsave(
        file.path("graficas/faltantes", paste0("patron_md_", evento, "_", cat_actual, ".png")),
        plot = p1,
        width = max(7, 0.8 * n_vars_cat + 3),
        height = max(4, 0.25 * nrow(unique(clasificar_celdas(df_evt)[, vars_cat, drop = FALSE])) + 3),
        dpi = 150, bg = "white"
      )
    }
    
    p2 <- mapa_faltantes_completo(df_evt, nom_evt, cat_actual)
    if (!is.null(p2)) {
      print(p2)
      ggsave(
        file.path("graficas/faltantes", paste0("mapa_", evento, "_", cat_actual, ".png")),
        plot = p2,
        width = 9,
        height = max(4, 0.3 * n_vars_cat + 2),
        dpi = 150, bg = "white"
      )
    }
  }
}

cat("\n✅ ¡Proceso finalizado correctamente! Gráficas sincronizadas y guardadas.\n")

# ============================================================
# 5. RESUMEN DE FALTANTES POR CATEGORÍA (tablas)
# ============================================================

patron_por_categoria <- function(df, nombre) {
  estados <- clasificar_celdas(df)
  resultados <- list()
  
  for (cat_nombre in names(categorias_vars)) {
    vars_cat <- intersect(categorias_vars[[cat_nombre]], names(estados))
    if (length(vars_cat) == 0) next
    
    sub_estados <- estados %>% select(all_of(vars_cat))
    
    resultados[[cat_nombre]] <- tibble(
      categoria = cat_nombre,
      variable  = names(sub_estados),
      n_presente  = colSums(sub_estados == "presente", na.rm = TRUE),
      n_faltante  = colSums(sub_estados == "faltante", na.rm = TRUE),
      n_no_aplica = colSums(sub_estados == "no_aplica", na.rm = TRUE),
      n_total     = nrow(sub_estados)
    ) %>%
      mutate(
        pct_presente  = round(n_presente  / n_total * 100, 2),
        pct_faltante  = round(n_faltante  / n_total * 100, 2),
        pct_no_aplica = round(n_no_aplica / n_total * 100, 2),
        evento = nombre
      )
  }
  
  bind_rows(resultados)
}

patron_dengue <- patron_por_categoria(dengue, "Dengue")
patron_chik   <- patron_por_categoria(chikungunya, "Chikungunya")

# ============================================================
# 6. TOP 15 PATRONES DE FALTANTES REALES (registros)
# ============================================================

patrones_reales <- function(df, nombre) {
  estados <- clasificar_celdas(df) %>%
    select(-any_of(c("consecutive", "cod_eve", "nombre_evento",
                     "evento", "va_sispro", "ano")))
  
  mat <- as.data.frame(lapply(estados, function(x) as.integer(x == "faltante")))
  names(mat) <- names(estados)
  
  patron_str <- do.call(paste, c(mat, sep = ""))
  tab <- sort(table(patron_str), decreasing = TRUE)
  top <- head(tab, 15)
  
  mat_pat <- do.call(rbind, strsplit(names(top), ""))
  colnames(mat_pat) <- names(mat)
  mat_pat <- as.data.frame(mat_pat)
  mat_pat[] <- lapply(mat_pat, as.integer)
  
  data.frame(
    n = as.integer(top),
    n_faltantes = rowSums(mat_pat),
    variables_faltantes = apply(mat_pat, 1, function(r) {
      cols <- names(mat_pat)[r == 1]
      if (length(cols) == 0) "(ninguna)" else paste(cols, collapse = ", ")
    }),
    evento = nombre,
    stringsAsFactors = FALSE
  )
}

cat("\n--- TOP 15 PATRONES DE FALTANTES REALES (Dengue) ---\n")
print(as.data.frame(patrones_reales(dengue, "Dengue")), row.names = FALSE)

cat("\n--- TOP 15 PATRONES DE FALTANTES REALES (Chikungunya) ---\n")
print(as.data.frame(patrones_reales(chikungunya, "Chikungunya")), row.names = FALSE)

# ============================================================
# 7. DIAGNÓSTICO FINAL
# ============================================================

diagnostico <- function(falt_resumen, nombre) {
  cat("\n---", nombre, "---\n")
  cat("NA totales:         ", comma(sum(falt_resumen$n_na_total)), "\n")
  cat("NA estructurales:   ", comma(sum(falt_resumen$n_na_estructural)),
      "(", round(sum(falt_resumen$n_na_estructural) / sum(falt_resumen$n_na_total) * 100, 1), "%)\n")
  cat("NA reales:          ", comma(sum(falt_resumen$n_na_real)), "\n")
  cat("\nVariables con faltantes reales:\n")
  print(
    falt_resumen %>%
      filter(pct_real > 0) %>%
      arrange(desc(pct_real)) %>%
      select(variable, pct_real, categoria) %>%
      as.data.frame(),
    row.names = FALSE
  )
}

diagnostico(falt_dengue, "DENGUE")
diagnostico(falt_chik, "CHIKUNGUNYA")

# ============================================================
# 8. DIAGNÓSTICO PREVIO A LA IMPUTACIÓN (solo consola)
# ============================================================

cat("\n")
cat("=========================================================\n")
cat("8. DIAGNÓSTICO PREVIO A LA IMPUTACIÓN\n")
cat("=========================================================\n")

# ------------------------------------------------------------
# 8.1 Diagnóstico básico
# ------------------------------------------------------------
diagnostico_basico <- function(df, nombre) {
  
  estados <- clasificar_celdas(df)
  
  resumen_var <- tibble(
    variable      = names(estados),
    pct_faltante  = sapply(estados, function(x) mean(x == "faltante") * 100),
    pct_no_aplica = sapply(estados, function(x) mean(x == "no_aplica") * 100),
    pct_presente  = sapply(estados, function(x) mean(x == "presente") * 100),
    n_faltante    = sapply(estados, function(x) sum(x == "faltante"))
  ) %>%
    arrange(desc(pct_faltante))
  
  na_por_fila <- rowSums(estados == "faltante")
  umbral      <- 0.5 * ncol(estados)
  
  list(
    resumen_var       = resumen_var,
    n_filas           = nrow(df),
    n_vars            = ncol(estados),
    max_pct           = max(resumen_var$pct_faltante),
    med_pct           = median(resumen_var$pct_faltante[resumen_var$pct_faltante > 0]),
    n_criticas        = sum(resumen_var$pct_faltante > 40),
    n_filas_malas     = sum(na_por_fila > umbral),
    pct_filas_malas   = round(sum(na_por_fila > umbral) / nrow(df) * 100, 2),
    mediana_na_fila   = median(na_por_fila),
    max_na_fila       = max(na_por_fila)
  )
}

# ------------------------------------------------------------
# 8.2 Test MCAR de Little
# ------------------------------------------------------------
test_mcar_seguro <- function(df, nombre, n_max = 5000) {
  
  df_analisis <- df %>%
    select(-starts_with("na_estructural_")) %>%
    mutate(across(where(is.character), as.factor))
  
  if (nrow(df_analisis) > n_max) {
    set.seed(123)
    df_analisis <- df_analisis %>% slice_sample(n = n_max)
  }
  
  df_analisis <- df_analisis %>% select(where(~ any(is.na(.))))
  
  if (ncol(df_analisis) < 2) return(NULL)
  
  tryCatch(
    naniar::mcar_test(df_analisis),
    error = function(e) NULL
  )
}

# ------------------------------------------------------------
# 8.3 Coocurrencia de faltantes
# ------------------------------------------------------------
coocurrencia_na <- function(df) {
  estados <- clasificar_celdas(df)
  vars    <- names(estados)[sapply(estados, function(x) any(x == "faltante"))]
  if (length(vars) < 2) return(NULL)
  
  mat <- as.data.frame(lapply(estados[vars], function(x) as.integer(x == "faltante")))
  cor(mat, use = "pairwise.complete.obs")
}

# ------------------------------------------------------------
# 8.4 Función maestra: diagnóstico + recomendación
# ------------------------------------------------------------
diagnostico_imputacion <- function(df, nombre, mcar_result, cor_na) {
  
  diag <- diagnostico_basico(df, nombre)
  
  # --- Decisión MCAR ---
  if (is.null(mcar_result)) {
    p_mcar    <- NA
    mech_txt  <- "No calculable"
    mech_code <- "?"
  } else {
    p_mcar    <- as.numeric(mcar_result$p.value)
    if (p_mcar < 0.05) {
      mech_txt  <- "MAR / MNAR (p < 0.05)"
      mech_code <- "MAR"
    } else {
      mech_txt  <- "MCAR (p >= 0.05)"
      mech_code <- "MCAR"
    }
  }
  
  # --- Pares de NA correlacionados ---
  n_pares_corr <- 0
  top_pares    <- character(0)
  if (!is.null(cor_na)) {
    idx <- which(abs(cor_na) > 0.3 & upper.tri(cor_na), arr.ind = TRUE)
    n_pares_corr <- nrow(idx)
    if (n_pares_corr > 0) {
      idx <- idx[order(-abs(cor_na[idx])), , drop = FALSE]
      top_n <- head(idx, 5)
      top_pares <- apply(top_n, 1, function(i) {
        sprintf("%s<->%s (r=%.2f)",
                rownames(cor_na)[as.integer(i[1])],
                colnames(cor_na)[as.integer(i[2])],
                cor_na[as.integer(i[1]), as.integer(i[2])])
      })
    }
  }
  
  # --- DECISIÓN DEL MÉTODO (lógica corregida) ---
  metodo <- case_when(
    # MCAR + pocas filas -> pmm (Chikungunya cae aquí)
    mech_code == "MCAR" & diag$n_filas < 1e4
    ~ "mice + method='pmm' (m=5, maxit=5) — MCAR, datos chicos",
    
    # MCAR + muchas filas -> pmm también sirve
    mech_code == "MCAR" & diag$n_filas >= 1e4
    ~ "mice + method='pmm' (m=5)",
    
    # MAR + pocas filas -> pmm
    mech_code == "MAR" & diag$n_filas < 1e4
    ~ "mice + method='pmm' (m=5) — MAR, datos chicos",
    
    # MAR + muchas filas + poco faltante -> rf
    mech_code == "MAR" & diag$n_filas >= 1e4 & diag$max_pct < 30
    ~ "mice + method='rf' (m=5, maxit=3) — MAR, base grande",
    
    # MAR + muchas filas + mucho faltante -> rf o kNN
    mech_code == "MAR" & diag$n_filas >= 1e4 & diag$max_pct >= 30
    ~ "VIM::kNN o missRanger (base grande + alta falta)",
    
    # Fallback
    TRUE ~ "Revisar manualmente"
  )
  
  # --- Variables sugeridas a eliminar ---
  vars_eliminar <- diag$resumen_var %>%
    filter(pct_faltante > 40) %>%
    pull(variable)
  
  # --- IMPRIMIR REPORTE ---
  cat("\n")
  cat("=====================================================\n")
  cat("  DIAGNÓSTICO DE IMPUTACIÓN:", nombre, "\n")
  cat("=====================================================\n")
  
  cat("\n── DIMENSIONES ──\n")
  cat("  Filas totales:                ", comma(diag$n_filas), "\n")
  cat("  Variables:                    ", diag$n_vars, "\n")
  
  cat("\n── FALTANTES POR VARIABLE ──\n")
  cat("  Máx % faltante real:          ", round(diag$max_pct, 2), "%\n")
  cat("  Mediana % faltante real:      ", round(diag$med_pct, 2), "%\n")
  cat("  Variables con >40% faltante:  ", diag$n_criticas, "\n")
  
  cat("\n── FALTANTES POR FILA ──\n")
  cat("  Mediana NA por fila:          ", diag$mediana_na_fila, "\n")
  cat("  Máx NA en una fila:           ", diag$max_na_fila, "\n")
  cat("  Filas con >50% NA:            ",
      comma(diag$n_filas_malas),
      "(", diag$pct_filas_malas, "%)\n")
  
  cat("\n── MECANISMO DE FALTANTES (Little MCAR test) ──\n")
  if (is.na(p_mcar)) {
    cat("  p-value:                       No calculable\n")
  } else {
    cat("  p-value:                       ", signif(p_mcar, 4), "\n")
  }
  cat("  Mecanismo:                     ", mech_txt, "\n")
  
  cat("\n── COOCURRENCIA DE FALTANTES ──\n")
  cat("  Pares con |r| > 0.3:           ", n_pares_corr, "\n")
  if (n_pares_corr > 0) {
    cat("  Top pares correlacionados:\n")
    for (p in top_pares) cat("    - ", p, "\n", sep = "")
  }
  
  cat("\n── RECOMENDACIÓN ──\n")
  cat("  Método sugerido:               ", metodo, "\n")
  if (length(vars_eliminar) > 0) {
    cat("  Variables a considerar ELIMINAR (>40% NA):\n")
    for (v in vars_eliminar) cat("    - ", v, "\n", sep = "")
  } else {
    cat("  Variables a eliminar:          ninguna\n")
  }
  
  tibble(
    evento                = nombre,
    n_filas               = diag$n_filas,
    n_vars                = diag$n_vars,
    max_pct_faltante      = round(diag$max_pct, 2),
    mediana_pct_faltante  = round(diag$med_pct, 2),
    n_vars_criticas_40    = diag$n_criticas,
    mediana_na_por_fila   = diag$mediana_na_fila,
    pct_filas_malas       = diag$pct_filas_malas,
    mcar_p_value          = ifelse(is.na(p_mcar), NA, signif(p_mcar, 4)),
    mecanismo             = mech_code,
    n_pares_corr_na       = n_pares_corr,
    metodo_recomendado    = metodo,
    vars_eliminar_sugeridas = paste(vars_eliminar, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 8.5 Ejecutar por evento
# ------------------------------------------------------------

mcar_d <- test_mcar_seguro(dengue,      "Dengue",      n_max = 5000)
mcar_c <- test_mcar_seguro(chikungunya, "Chikungunya", n_max = 5000)

cor_d  <- coocurrencia_na(dengue)
cor_c  <- coocurrencia_na(chikungunya)

tabla_d <- diagnostico_imputacion(dengue,      "Dengue",      mcar_d, cor_d)
tabla_c <- diagnostico_imputacion(chikungunya, "Chikungunya", mcar_c, cor_c)

tabla_resumen <- bind_rows(tabla_d, tabla_c)

# ------------------------------------------------------------
# 8.6 Tabla resumen
# ------------------------------------------------------------

cat("\n")
cat("=========================================================\n")
cat("RESUMEN COMPARATIVO FINAL\n")
cat("=========================================================\n")
print(as.data.frame(tabla_resumen), row.names = FALSE)

write.csv(tabla_resumen,
          here("data", "processed", "diagnostico_imputacion.csv"),
          row.names = FALSE)

cat("\n✅ Diagnóstico pre-imputación completado.\n")
