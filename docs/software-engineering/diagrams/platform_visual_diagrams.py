from diagrams import Cluster, Diagram
from diagrams.gcp.database import SQL
from diagrams.gcp.network import LoadBalancing
from diagrams.gcp.storage import GCS
from diagrams.gis.geocoding import Pelias
from diagrams.gis.routing import Graphhopper
from diagrams.onprem.client import User
from diagrams.programming.language import Python


def build_platform_visual() -> None:
    graph_attr = {
        "rankdir": "LR",
        "splines": "spline",
        "bgcolor": "#F4F8FF",
        "pad": "0.50",
        "nodesep": "0.75",
        "ranksep": "0.95",
        "fontname": "Helvetica",
        "fontsize": "13",
        "labelloc": "t",
        "label": "Plataforma Webmining (Vista Visual)",
    }

    node_attr = {
        "fontname": "Helvetica",
        "fontsize": "12",
        "fontcolor": "#17324C",
    }

    edge_attr = {
        "color": "#4E6785",
        "penwidth": "1.6",
        "fontname": "Helvetica",
        "fontsize": "10",
    }

    with Diagram(
        name="platform_visual_diagrams",
        filename="docs/software-engineering/diagrams/platform_visual_diagrams",
        outformat="png",
        show=False,
        graph_attr=graph_attr,
        node_attr=node_attr,
        edge_attr=edge_attr,
    ):
        user = User("Usuario")

        with Cluster("Capa de Aplicacion"):
            ui = Python("Flet UI\n(app.py)")
            orchestrator = LoadBalancing("Controlador de flujo")
            exports = GCS("Exportes CSV/PDF")

        with Cluster("Capa de Dominio"):
            optimizer = Python("Motor de optimizacion")
            scheduler = Python("Validador de horarios")

        with Cluster("Capa de Datos"):
            loader = Python("data_loader.py")
            datasets = SQL("Depots + matrices")
            places = Pelias("Google Places API")
            matrix = Graphhopper("Distance Matrix API")

        user >> ui >> orchestrator >> optimizer >> scheduler
        optimizer >> exports
        orchestrator >> loader >> datasets
        places >> datasets
        matrix >> datasets


if __name__ == "__main__":
    build_platform_visual()
