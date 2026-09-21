# Proyecto Dengue y Chikunguña en Colombia (2020-2025)

Análisis epidemiológico y sistema de apoyo a la decisión para autoridades regionales.

## Fuente de datos
- SIVIGILA (Instituto Nacional de Salud)
- Códigos: 210 (Dengue), 217 (Chikunguña)
- Período: 2020-2025

## Datos
Los archivos crudos de SIVIGILA no se incluyen en el repositorio.
Descárgalos desde: [https://portalsivigila.ins.gov.co/Paginas/Buscador.aspx]
- Dengue (código 210): colocar en `data/raw/dengue/`
- Chikunguña (código 217): colocar en `data/raw/chikungunya/`
Fecha de descarga: [la fecha que corresponda]

## Requisitos
- R versión 4.5.x
- RStudio
- Git

## Cómo replicar este proyecto

1. clonar el proyecto

- Abre Rstuudio
- Ve a `File → New Project…`
- Elige Version Control → Git
- En Repository URL, pega: https://github.com/pinkandredval/predecir_el_comportamiento_de_Dengue_y_Chikunguna_en_Zonas_Criticas_de_Colombia
- En Project directory name elige un nombre para el proyecto
- Elige la carpeta donde quieres guardar el proyecto
- Clic en Create Project

RStudio descargará el proyecto y lo abrirá automáticamente.

2. Restaurar los paquetes de R

Al abrir el proyecto, verás en la consola de R algo como: `- Project '...' loaded. [renv 1.2.4]`

- Ahora ejecuta en la consola: 
```bash
renv::restore()
```
   - Te va a preguntar si quieres continuar. Responde Y y presiona Enter.
```bash
renv::status()
```
Debe responder: `No issues found -- the project is in a consistent state.` o podria salir una respuesta con la lista de paquetes descargados pero que no estan en uso aun

3. Descarga los datos crudos

Los archivos de SIVIGILA no vienen incluidos en el repositorio

- Entra al portal: https://portalsivigila.ins.gov.co/Paginas/Buscador.aspx
- Descarga los 6 archivos de Dengue (código 210) y los 6 de Chikunguña (código 217), uno por cada año entre 2020 y 2025.



## Cómo replicar este proyecto

1. Clona el repositorio:
   ```bash
   git clone <URL-del-repo>
   ```
   bhfbh
