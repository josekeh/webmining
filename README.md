# Webmining - Optimización de rutas (TSP) con OR-Tools

Proyecto en Python para resolver un problema de **Traveling Salesman Problem (TSP)** a partir de una matriz de distancias entre puntos (lugares), usando **programación lineal entera mixta** con OR-Tools.

## Estado actual

- **Núcleo de optimización**: implementado en `src/or_tools.py` con módulos refactorizados (`src/optimizer.py`, `src/data_loader.py`, `src/display.py`).
- **Interfaz web**: `app.py` proporciona una aplicación **Streamlit** interactiva.
- **CLI ejecutable**: `main.py` orquesta lectura de datos, solve y reporte desde terminal.
- **Notebooks**: `Conexión_API_GoogleMaps.ipynb` y `test.ipynb` para exploración.
- **Matrices de ejemplo**: `matriz_distancias.csv` (con unidades), `data/distancias_test.csv`, `data/distancias.csv`; también `data/nodos.txt` y `data/tiempos.csv`.

## Objetivo

Encontrar un ciclo Hamiltoniano de costo mínimo que:

- visite todos los nodos exactamente una vez,
- regrese al nodo de origen,
- minimice la suma de distancias.

## Arquitectura del proyecto

```text
webmining/
├─ main.py                       # CLI: orquesta lectura, solve y reporte
├─ app.py                        # Streamlit: interfaz web interactiva para el optimizer
├─ src/
│  ├─ or_tools.py               # Modelo TSP original (variables, restricciones, objetivo)
│  ├─ optimizer.py              # Lógica de optimización refactorizada
│  ├─ data_loader.py            # Carga y preprocesamiento de matrices de distancias
│  └─ display.py                # Componentes de visualización y reporte
├─ docs/
│  └─ diagrams/
│     ├─ tsp_flow.py            # Genera diagrama con Python Diagrams (usa Graphviz)
│     └─ tsp_flow.dot           # Fuente DOT de Graphviz
├─ data/
│  ├─ distancias.csv            # Matriz de distancias (formato principal)
│  ├─ distancias_test.csv       # Matriz numérica de ejemplo
│  ├─ nodos.txt                 # Definición de nodos
│  └─ tiempos.csv               # Matriz de costos/tiempos
├─ matriz_distancias.csv        # Matriz alternativa con unidades de texto ("km", "m")
├─ Conexión_API_GoogleMaps.ipynb # Obtención/experimentación con distancias
├─ test.ipynb                   # Pruebas exploratorias
├─ Makefile                     # Targets para generación de diagramas
├─ requirements.txt
├─ pyproject.toml
├─ uv.lock                      # Dependencias bloqueadas para reproducibilidad
└─ .python-version              # Python 3.12+
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

## Modelo matemático (resumen)

### Variables de decisión

- $x_{ij} \in \{0,1\}$ para $i \neq j$: vale 1 si el tour va de $i$ a $j$.
- Variables auxiliares $u_i$ (MTZ) para eliminar subciclos.

### Función objetivo

Minimizar:

$$
\min \sum_{i}\sum_{j \neq i} c_{ij}x_{ij}
$$

### Restricciones

- **Una entrada por nodo**: $\sum_{i \neq j} x_{ij} = 1$.
- **Una salida por nodo**: $\sum_{j \neq i} x_{ij} = 1$.
- **Eliminación de subciclos (MTZ)** para nodos distintos del depósito/base.

## Requisitos

- Python 3.12+
- Dependencias:
	- `googlemaps`
	- `graphviz`
	- `diagrams`
	- `pandas`
	- `matplotlib`
	- `scipy`
	- `ortools`

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

Ejemplo mínimo para ejecutar el solver desde script o notebook:

```python
import pandas as pd
from src.or_tools import calcular_ruta_optima

# Para el solver conviene usar la versión numérica (sin sufijos "km"/"m")
dist = pd.read_csv("data/distancias_test.csv", index_col=0)

status, modelo, x_ij, nodos, costos = calcular_ruta_optima(dist, dist.values)

print("Status:", status)

arcos = []
for i in nodos:
		for j in x_ij[i]:
				if x_ij[i][j].solution_value() > 0.5:
						arcos.append((i, j))

print("Arcos elegidos:", arcos)
```

## Formato de datos

El solver espera una matriz cuadrada de costos consistente con los índices:

- filas y columnas representan los mismos nodos,
- diagonal puede ser 0 o 1 según preprocesamiento,
- valores recomendados: numéricos (`float`/`int`).

### Archivos incluidos

- `data/distancias_test.csv`: matriz numérica lista para probar.
- `matriz_distancias.csv`: matriz con unidades de texto (`"4.1 km"`, `"1 m"`), útil como fuente cruda para preprocesar.

## Ejecución

### Interfaz web (Streamlit)

```bash
uv run streamlit run app.py
```

La aplicación abrirá en `http://localhost:8501` con UI interactiva para carga de datos, ejecución del optimizer y visualización.

### CLI (Terminal)

```bash
uv run python main.py
```

Orquesta automáticamente: lectura de datos, ejecución del solver y generación de reporte.

## Limitaciones y próximos pasos

- `case="vrp"` en `load_constraints` está declarado pero no implementado.
- Sería útil agregar:
	- validación y normalización automática de matrices con unidades,
	- reconstrucción y visualización explícita del tour ordenado,
	- tests unitarios completos para restricciones y extracción de solución,
	- documentación de API pública para integración en third-party apps.