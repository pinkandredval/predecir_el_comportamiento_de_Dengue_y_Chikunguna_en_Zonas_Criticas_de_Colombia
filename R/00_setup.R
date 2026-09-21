# ============================================================
# 00_setup.R
# Configuración inicial y descompresión de datos crudos
# Proyecto: Dengue y Chikunguña en Colombia (2020-2025)
# ============================================================

# --- Descomprimir datos de SIVIGILA ---
# Solo se ejecuta la primera vez. Si los archivos ya están
# descomprimidos, este bloque no hace nada (sobrescribe).

unzip("data/raw/dengue.zip",
      exdir = "data/raw/dengue/")

unzip("data/raw/chikunguya.zip",
      exdir = "data/raw/chikungunya/")

# --- Verificar que se descomprimió correctamente ---
list.files("data/raw/dengue/")
list.files("data/raw/chikungunya/")


# Mover archivos de dengue un nivel arriba
archivos_dengue <- list.files("data/raw/dengue/dengue/",
                              full.names = TRUE)
file.rename(archivos_dengue,
            file.path("data/raw/dengue", basename(archivos_dengue)))

# Mover archivos de chikungunya un nivel arriba
archivos_chik <- list.files("data/raw/chikungunya/chikunguya/",
                            full.names = TRUE)
file.rename(archivos_chik,
            file.path("data/raw/chikungunya", basename(archivos_chik)))

# Borrar las carpetas vacías
unlink("data/raw/dengue/dengue", recursive = TRUE)
unlink("data/raw/chikungunya/chikunguya", recursive = TRUE)

list.files("data/raw/dengue/")
list.files("data/raw/chikungunya/")
