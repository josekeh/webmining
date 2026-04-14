from diagrams import Diagram
from diagrams.generic.blank import Blank


def build_tsp_flow_diagram() -> None:
    graph_attr = {
        "rankdir": "LR",
        "splines": "spline",
        "fontname": "Helvetica",
        "fontsize": "12",
    }

    node_attr = {
        "shape": "box",
        "style": "rounded,filled",
        "fillcolor": "#F8FAFC",
        "color": "#334155",
        "fontname": "Helvetica",
        "fontsize": "11",
    }

    edge_attr = {
        "color": "#475569",
        "fontname": "Helvetica",
        "fontsize": "10",
    }

    with Diagram(
        name="tsp_flow_diagrams",
        filename="docs/diagrams/tsp_flow_diagrams",
        outformat="png",
        show=False,
        graph_attr=graph_attr,
        node_attr=node_attr,
        edge_attr=edge_attr,
    ):
        input_data = Blank("Dataset depot (nodos + matrices)")
        load = Blank("Carga con data_loader.py")
        choose = Blank("Elegir criterio y orden de tipos")
        combos = Blank("Generar combinaciones de nodos")
        schedule = Blank("Validar horarios y estadias")
        score = Blank("Evaluar costo/puntaje")
        select = Blank("Seleccionar mejor ruta")
        status = Blank("¿Ruta factible encontrada?")
        extract = Blank("Mostrar mapa, tablas y exportes")
        review = Blank("Ajustar parametros / reintentar")

        input_data >> load >> choose >> combos >> schedule >> score >> select >> status
        status >> extract
        status >> review


if __name__ == "__main__":
    build_tsp_flow_diagram()
