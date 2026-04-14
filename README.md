# Webmining - Optimizador de recorridos turisticos

Proyecto en Python para planificar recorridos por categorias de lugares (por ejemplo: museo, restaurante, heladeria), saliendo y regresando a un depot.

La app ofrece:

- interfaz grafica con Flet,
- optimizacion por criterio (distancia, tiempo o puntaje),
- ventanas horarias y tiempos de estadia,
- visualizacion de ruta y exportes (CSV/PDF).

## Estado actual

- La interfaz principal esta en `app.py`.
- La logica de carga de datos esta en `src/data_loader.py`.
- El motor de optimizacion esta en `src/optimizer.py`.
- La salida por consola heredada esta en `main.py` + `src/display.py`.
- Los diagramas estan en `docs/diagrams/`.

## Interfaz grafica (Flet)

La aplicacion incluye una UI moderna para:

- seleccionar depot de salida,
- definir criterio de optimizacion,
- configurar orden de actividades y tiempos de estadia,
- visualizar metricas, ruta optima y tablas de detalle.

### Ejecucion en escritorio (desktop)

```bash
uv run flet run app.py
```

Alternativa equivalente:

```bash
uv run python app.py
```

### Ejecucion en web

```bash
uv run flet run --web --host 0.0.0.0 --port 8550 app.py
```

Luego abrir en navegador:

```text
http://localhost:8550
```

### Ejecucion en movil (modo desarrollo)

1. Verificar dispositivos/emuladores disponibles:

```bash
uv run flet devices
uv run flet emulators
```

2. Ejecutar en Android:

```bash
uv run flet run --android app.py
```

3. Ejecutar en iOS:

```bash
uv run flet run --ios app.py
```

Nota: tambien podes usar el modo web desde el celular abriendo la URL del host en la misma red local.

### Build para distribucion (opcional)

```bash
uv run flet build apk .
uv run flet build ipa .
uv run flet build web .
```

Entrada de consola heredada:

```bash
uv run python main.py
```

## Objetivo

Encontrar el mejor recorrido que:

- comienza y termina en el depot,
- visita una parada por cada tipo seleccionado,
- respeta ventanas horarias (apertura/cierre) si aplica,
- optimiza el criterio elegido (distancia, tiempo o puntaje).

## Arquitectura del proyecto

```text
webmining/
├─ app.py                         # UI Flet
├─ main.py
├─ src/
│  ├─ data_loader.py             # Carga de nodos y matrices
│  ├─ optimizer.py               # Optimizacion y validacion de horarios
│  └─ display.py                 # Impresion por consola
├─ docs/
│  └─ diagrams/
│     ├─ tsp_flow.py             # Genera diagrama con Python Diagrams (usa Graphviz)
│     └─ tsp_flow.dot            # Fuente DOT de Graphviz
├─ data/
│  ├─ depots.txt
│  ├─ depot0/..depot10/          # datasets por depot
│  ├─ g_nodos.txt                # dataset base de nodos
│  ├─ g_distancias.csv           # dataset base de distancias
│  └─ g_tiempos.csv              # dataset base de tiempos
├─ matriz_distancias.csv         # matriz historica de referencia
├─ Conexión_API_GoogleMaps.ipynb
├─ Conexión_API_GoogleMaps+Hoteles.ipynb
├─ requirements.txt
├─ pyproject.toml
└─ Makefile
```

## Diagramas (Diagrams + Graphviz)

El proyecto usa **Graphviz** como motor y ofrece dos formas de mantener diagramas:

- `docs/diagrams/tsp_flow.dot`: definición manual en formato DOT.
- `docs/diagrams/tsp_flow.py`: generación programática con la librería `diagrams` (que internamente usa Graphviz).

### Generar con Graphviz (DOT)

```bash
dot -Tpng docs/diagrams/tsp_flow.dot -o docs/diagrams/tsp_flow.png
dot -Tsvg docs/diagrams/tsp_flow.dot -o docs/diagrams/tsp_flow.svg
```

También podés usar:

```bash
uv run make diagram-dot
```

### Generar con Python Diagrams

Instalación de dependencias (flujo recomendado del proyecto):

```bash
uv sync
```

Ejecutar generador:

```bash
uv run python docs/diagrams/tsp_flow.py
```

También podés usar:

```bash
uv run make diagram-diagrams
```

Salida esperada:

- `docs/diagrams/tsp_flow_diagrams.png`

### Generar todo junto

```bash
uv run make diagrams
```

Este comando ejecuta automáticamente `check-env` antes de generar los diagramas.

### Verificar entorno

```bash
make check-env
```

## Enfoque de optimizacion actual

El optimizador actual implementa una busqueda por combinaciones:

- construye candidatos eligiendo un nodo por cada tipo en el orden definido,
- evalua cada ruta como `[depot, ...paradas..., depot]`,
- filtra rutas infactibles por ventanas horarias,
- selecciona la mejor segun el criterio:
	- `distancia`: minimiza distancia total,
	- `tiempo`: minimiza tiempo total de traslado,
	- `puntaje`: maximiza puntaje total y desempata por distancia.

Nota: no se usa actualmente un modelo MIP/MTZ en runtime.

## Requisitos

- Python 3.12+
- Dependencias:
	- `flet`
	- `staticmap`
	- `pandas`
	- `matplotlib`
	- `scipy`
	- `ortools`
	- `graphviz`
	- `diagrams`
	- `googlemaps`

## Instalación

### Opción recomendada (uv)

```bash
uv sync
```

### Opción alternativa: pip

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Opción 2: entorno con `pyproject.toml` (uv/pip)

```bash
pip install -e .
```

## Uso rápido (API actual)

Ejemplo minimo para ejecutar el optimizador desde script o notebook:

```python
from src.data_loader import load_nodes, load_matrix, get_nodes_by_type
from src.optimizer import optimize

nodes = load_nodes("data/depot0/g_nodos.txt")
dist = load_matrix("data/depot0/g_distancias.csv")
time_mat = load_matrix("data/depot0/g_tiempos.csv")
by_type = get_nodes_by_type(nodes)

type_order = ["Museo", "Restaurante", "Heladeria"]
stay_times = [45, 60, 30]

best_route, best_value, schedule = optimize(
	nodes_df=nodes,
	nodes_by_type=by_type,
	cost_matrix=dist,
	type_order=type_order,
	criterion="distancia",
	dist_matrix=dist,
	stay_times=stay_times,
	time_matrix=time_mat,
	start_minutes=10 * 60,
)

print("Ruta:", best_route)
print("Valor:", best_value)
print("Schedule:", schedule)
```

## Formato de datos

Cada dataset de depot debe incluir:

- `g_nodos.txt` con columnas como `ID`, `nombre`, `tipo`, `puntaje`, `lat`, `lon`, `Apertura`, `Cierre`.
- `g_distancias.csv` con matriz cuadrada de distancias.
- `g_tiempos.csv` con matriz cuadrada de tiempos (segundos).

Los IDs deben coincidir entre nodos y matrices.

## Limitaciones y próximos pasos

- El enfoque por combinaciones escala peor cuando crece mucho la cantidad de nodos por tipo.
- Falta una suite de tests automatizados.
- Puede mejorarse:
	- cacheo de evaluaciones,
	- poda heuristica,
	- benchmark por tamaño de instancia,
	- exporte de reportes con metadatos de configuracion.