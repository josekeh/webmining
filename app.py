import base64
import asyncio
import io
import os
from datetime import datetime

import flet as ft
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from staticmap import CircleMarker, Line, StaticMap

from src.data_loader import get_nodes_by_type, list_depots, load_matrix, load_nodes
from src.optimizer import optimize


TYPE_COLORS = {
    "depot": "#111827",
    "Cervecería": "#E63946",
    "Cerveceria": "#E63946",
    "Heladería": "#1D4ED8",
    "Heladeria": "#1D4ED8",
    "Museo": "#047857",
    "Restaurante": "#D97706",
}

PLACEHOLDER = "-- Selecciona un punto de partida --"


def load_data(depot_folder: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, list[int]]]:
    nodes = load_nodes(f"{depot_folder}/g_nodos.txt")
    dist = load_matrix(f"{depot_folder}/g_distancias.csv")
    time_mat = load_matrix(f"{depot_folder}/g_tiempos.csv")
    by_type = get_nodes_by_type(nodes)
    return nodes, dist, time_mat, by_type


def parse_hhmm(value: str) -> int | None:
    try:
        hh, mm = value.strip().split(":")
        h = int(hh)
        m = int(mm)
        if h < 0 or h > 23 or m < 0 or m > 59:
            return None
        return h * 60 + m
    except (AttributeError, ValueError):
        return None


def theme_tokens(is_light: bool) -> dict[str, str]:
    if is_light:
        return {
            "row_bg": "#FFFFFF",
            "row_border": "#D7E2F0",
            "input_bg": "#F8FBFF",
            "input_text": "#0F172A",
            "input_label": "#475569",
            "input_border": "#9EB4D0",
            "input_focus": "#2563EB",
            "tile_bg": "#FFFFFF",
            "tile_border": "#D7E2F0",
            "tile_text": "#0F172A",
            "tile_icon": "#1D4ED8",
            "table_head_bg": "#EEF4FC",
            "table_text": "#0F172A",
        }

    return {
        "row_bg": "#0F182D",
        "row_border": "#243350",
        "input_bg": "#0E172D",
        "input_text": "#F8FAFC",
        "input_label": "#BFDBFE",
        "input_border": "#2E3A5B",
        "input_focus": "#60A5FA",
        "tile_bg": "#0F182D",
        "tile_border": "#243350",
        "tile_text": "#E2E8F0",
        "tile_icon": "#BFDBFE",
        "table_head_bg": "#13203A",
        "table_text": "#E2E8F0",
    }


def metric_card(title: str, value: str, accent: str) -> ft.Container:
    return ft.Container(
        width=220,
        padding=18,
        border_radius=16,
        gradient=ft.LinearGradient(
            begin=ft.Alignment.TOP_LEFT,
            end=ft.Alignment.BOTTOM_RIGHT,
            colors=["#101A30", "#0E1628"],
        ),
        border=ft.border.all(1, "#2D3D5F"),
        shadow=ft.BoxShadow(
            blur_radius=18,
            spread_radius=0,
            color="#22000000",
            offset=ft.Offset(0, 8),
        ),
        content=ft.Column(
            spacing=4,
            controls=[
                ft.Text(title, size=12, color="#95A3C7", weight=ft.FontWeight.W_500),
                ft.Text(value, size=26, color=accent, weight=ft.FontWeight.W_700),
            ],
        ),
    )


def build_route_plot(nodes: pd.DataFrame, route: list[int]) -> str:
    # Preferred rendering: OpenStreetMap tiles with route overlay.
    try:
        route_points = nodes.loc[route]
        static_map = StaticMap(1200, 700)

        line_points = [(float(row["lon"]), float(row["lat"])) for _, row in route_points.iterrows()]
        if len(line_points) >= 2:
            static_map.add_line(Line(line_points, "#1D4ED8", 5))

        for step, (_, row) in enumerate(route_points.iterrows()):
            lon = float(row["lon"])
            lat = float(row["lat"])
            static_map.add_marker(CircleMarker((lon, lat), "#FFFFFF", 9))
            if step == 0:
                static_map.add_marker(CircleMarker((lon, lat), "#16A34A", 11))

        image = static_map.render()
        buffer = io.BytesIO()
        image.save(buffer, format="PNG")
        return base64.b64encode(buffer.getvalue()).decode("utf-8")
    except Exception:
        # Fallback when tile download fails (e.g., no internet).
        pass

    fig, ax = plt.subplots(figsize=(8, 5), dpi=130)
    fig.patch.set_facecolor("#0B1223")
    ax.set_facecolor("#0F1A33")

    for tipo, group in nodes.groupby("tipo"):
        color = TYPE_COLORS.get(tipo, "#9CA3AF")
        ax.scatter(group["lon"], group["lat"], c=color, s=45, alpha=0.35, label=tipo)

    route_points = nodes.loc[route]
    ax.plot(route_points["lon"], route_points["lat"], color="#60A5FA", linewidth=2.8, alpha=0.95)
    ax.scatter(route_points["lon"], route_points["lat"], c="#F8FAFC", s=30, zorder=3)

    for step, node_id in enumerate(route):
        x = nodes.loc[node_id, "lon"]
        y = nodes.loc[node_id, "lat"]
        ax.text(x, y, str(step), fontsize=8, color="#93C5FD", ha="center", va="bottom", weight="bold")

    ax.set_title("Recorrido optimizado", color="#E5E7EB", fontsize=13, pad=10)
    ax.set_xlabel("Longitud", color="#CBD5E1")
    ax.set_ylabel("Latitud", color="#CBD5E1")
    ax.tick_params(colors="#94A3B8")
    for spine in ax.spines.values():
        spine.set_color("#334155")

    legend = ax.legend(facecolor="#0B1223", edgecolor="#334155", fontsize=8)
    for text in legend.get_texts():
        text.set_color("#CBD5E1")

    buffer = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buffer, format="png", bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("utf-8")


def table_scroller(table: ft.DataTable) -> ft.Row:
    return ft.Row(
        scroll=ft.ScrollMode.AUTO,
        controls=[table],
    )


def export_pdf_report(
    output_path: str,
    summary: dict[str, str],
    route_df: pd.DataFrame,
    legs_df: pd.DataFrame,
) -> None:
    with PdfPages(output_path) as pdf:
        fig_summary = plt.figure(figsize=(8.27, 11.69))
        fig_summary.patch.set_facecolor("white")
        ax_summary = fig_summary.add_subplot(111)
        ax_summary.axis("off")

        lines = [
            "Optimizador de Recorrido Turistico - Reporte",
            "",
            f"Generado: {summary['generated_at']}",
            f"Depot: {summary['depot']}",
            f"Criterio: {summary['criterion']}",
            f"Distancia total: {summary['distance']}",
            f"Tiempo traslado: {summary['travel']}",
            f"Tiempo estadia: {summary['stay']}",
            f"Tiempo total: {summary['total_time']}",
            f"Puntaje total: {summary['score']}",
        ]
        ax_summary.text(0.05, 0.95, "\n".join(lines), fontsize=12, va="top", family="sans-serif")
        pdf.savefig(fig_summary, bbox_inches="tight")
        plt.close(fig_summary)

        fig_route = plt.figure(figsize=(11.69, 8.27))
        ax_route = fig_route.add_subplot(111)
        ax_route.axis("off")
        ax_route.set_title("Recorrido optimizado", fontsize=14, pad=16)
        route_table = ax_route.table(
            cellText=route_df.astype(str).values,
            colLabels=route_df.columns,
            loc="center",
            cellLoc="center",
        )
        route_table.auto_set_font_size(False)
        route_table.set_fontsize(8)
        route_table.scale(1, 1.35)
        pdf.savefig(fig_route, bbox_inches="tight")
        plt.close(fig_route)

        fig_legs = plt.figure(figsize=(11.69, 8.27))
        ax_legs = fig_legs.add_subplot(111)
        ax_legs.axis("off")
        ax_legs.set_title("Detalle tramo a tramo", fontsize=14, pad=16)
        legs_table = ax_legs.table(
            cellText=legs_df.astype(str).values,
            colLabels=legs_df.columns,
            loc="center",
            cellLoc="center",
        )
        legs_table.auto_set_font_size(False)
        legs_table.set_fontsize(9)
        legs_table.scale(1, 1.4)
        pdf.savefig(fig_legs, bbox_inches="tight")
        plt.close(fig_legs)


def main(page: ft.Page) -> None:
    page.title = "Optimizador de Recorrido Turistico"
    page.theme_mode = ft.ThemeMode.DARK
    page.padding = 0
    page.bgcolor = "#090E1A"
    page.scroll = ft.ScrollMode.AUTO
    page.window_min_width = 360
    page.window_min_height = 640
    page.fonts = {
        "Sora": "https://fonts.googleapis.com/css2?family=Sora:wght@400;500;600;700&display=swap"
    }
    page.theme = ft.Theme(font_family="Sora")

    depots = list_depots("data")
    if not depots:
        page.add(
            ft.Container(
                expand=True,
                    alignment=ft.Alignment.CENTER,
                content=ft.Text(
                    "No se encontraron depots en data/.",
                    color="#FCA5A5",
                    size=20,
                    weight=ft.FontWeight.W_600,
                ),
            )
        )
        return

    depot_paths = {name: path for path, name in depots}
    depot_options = [PLACEHOLDER] + sorted(depot_paths.keys())

    cache: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, list[int]]]] = {}

    def get_depot_data(depot_name: str):
        if depot_name not in cache:
            cache[depot_name] = load_data(depot_paths[depot_name])
        return cache[depot_name]

    initial_nodes, _, _, initial_by_type = get_depot_data(depot_options[1])
    available_types = list(initial_by_type.keys())

    depot_dropdown = ft.Dropdown(
        label="Punto de partida",
        width=360,
        value=PLACEHOLDER,
        options=[ft.dropdown.Option(name) for name in depot_options],
        border_color="#3B4B71",
        focused_border_color="#60A5FA",
    )

    criterion_picker = ft.SegmentedButton(
        selected=["distancia"],
        show_selected_icon=False,
        style=ft.ButtonStyle(
            bgcolor={ft.ControlState.SELECTED: "#1D4ED8", ft.ControlState.DEFAULT: "#0F182D"},
            color={ft.ControlState.SELECTED: "#F8FAFC", ft.ControlState.DEFAULT: "#BFDBFE"},
            side={ft.ControlState.DEFAULT: ft.BorderSide(1, "#2E3A5B")},
            shape={ft.ControlState.DEFAULT: ft.RoundedRectangleBorder(radius=10)},
        ),
        segments=[
            ft.Segment(
                value="distancia",
                label=ft.Text("Distancia", size=13),
                icon=ft.Icon(ft.Icons.ROUTE),
            ),
            ft.Segment(
                value="tiempo",
                label=ft.Text("Tiempo", size=13),
                icon=ft.Icon(ft.Icons.TIMER),
            ),
            ft.Segment(
                value="puntaje",
                label=ft.Text("Puntaje", size=13),
                icon=ft.Icon(ft.Icons.STAR),
            ),
        ],
    )

    theme_switch = ft.Switch(
        label="Modo claro",
        value=False,
        active_color="#0EA5E9",
    )

    stops_input = ft.TextField(
        label="Cantidad de paradas",
        value=str(max(1, len(available_types))),
        width=190,
        keyboard_type=ft.KeyboardType.NUMBER,
        border_color="#3B4B71",
        focused_border_color="#60A5FA",
    )

    start_time_input = ft.TextField(
        label="Hora de salida (HH:MM)",
        value="10:00",
        width=200,
        border_color="#3B4B71",
        focused_border_color="#60A5FA",
    )

    status_text = ft.Text(color="#FCA5A5", size=14, visible=False)
    stops_column = ft.Column(spacing=10)
    result_column = ft.Column(
        spacing=14,
        visible=False,
        opacity=0,
        offset=ft.Offset(0, 0.03),
        animate_opacity=ft.Animation(300, ft.AnimationCurve.EASE_OUT),
        animate_offset=ft.Animation(350, ft.AnimationCurve.DECELERATE),
    )

    export_status = ft.Text(color="#93C5FD", size=13, visible=False)
    latest_result: dict[str, pd.DataFrame | dict[str, str]] = {}

    def selected_types() -> list[str]:
        return [c.data["type"].value for c in stops_column.controls if c.data and c.data["type"].value]

    def selected_stays() -> list[int]:
        values: list[int] = []
        for c in stops_column.controls:
            raw = c.data["stay"].value if c.data else "30"
            try:
                values.append(max(0, int(raw)))
            except ValueError:
                values.append(30)
        return values

    def style_stop_rows() -> None:
        tokens = theme_tokens(theme_switch.value)
        for row in stops_column.controls:
            if not row.data:
                continue

            type_dropdown = row.data["type"]
            stay_field = row.data["stay"]

            row.bgcolor = tokens["row_bg"]
            row.border = ft.border.all(1, tokens["row_border"])

            type_dropdown.bgcolor = tokens["input_bg"]
            type_dropdown.color = tokens["input_text"]
            type_dropdown.text_style = ft.TextStyle(color=tokens["input_text"])
            type_dropdown.label_style = ft.TextStyle(color=tokens["input_label"])
            type_dropdown.border_color = tokens["input_border"]
            type_dropdown.focused_border_color = tokens["input_focus"]
            type_dropdown.filled = True

            stay_field.bgcolor = tokens["input_bg"]
            stay_field.fill_color = tokens["input_bg"]
            stay_field.color = tokens["input_text"]
            stay_field.text_style = ft.TextStyle(color=tokens["input_text"])
            stay_field.label_style = ft.TextStyle(color=tokens["input_label"])
            stay_field.border_color = tokens["input_border"]
            stay_field.focused_border_color = tokens["input_focus"]
            stay_field.filled = True

    def export_csv_click(_: ft.ControlEvent) -> None:
        if "route_df" not in latest_result or "legs_df" not in latest_result:
            export_status.value = "Primero optimiza un recorrido para poder exportar."
            export_status.color = "#FCA5A5"
            export_status.visible = True
            page.update()
            return

        os.makedirs("exports", exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        route_path = os.path.join("exports", f"recorrido_{stamp}.csv")
        legs_path = os.path.join("exports", f"tramos_{stamp}.csv")

        latest_result["route_df"].to_csv(route_path, index=False)
        latest_result["legs_df"].to_csv(legs_path, index=False)

        export_status.value = f"CSV exportado en {route_path} y {legs_path}"
        export_status.color = "#86EFAC"
        export_status.visible = True
        page.update()

    def export_pdf_click(_: ft.ControlEvent) -> None:
        if "route_df" not in latest_result or "legs_df" not in latest_result or "summary" not in latest_result:
            export_status.value = "Primero optimiza un recorrido para poder exportar."
            export_status.color = "#FCA5A5"
            export_status.visible = True
            page.update()
            return

        os.makedirs("exports", exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        pdf_path = os.path.join("exports", f"reporte_recorrido_{stamp}.pdf")

        export_pdf_report(
            output_path=pdf_path,
            summary=latest_result["summary"],
            route_df=latest_result["route_df"],
            legs_df=latest_result["legs_df"],
        )

        export_status.value = f"PDF exportado en {pdf_path}"
        export_status.color = "#86EFAC"
        export_status.visible = True
        page.update()

    def refresh_stop_rows(count: int) -> None:
        stops_column.controls.clear()
        safe_count = max(1, min(20, count))
        tokens = theme_tokens(theme_switch.value)

        active_depot = depot_dropdown.value
        types = available_types
        if active_depot and active_depot != PLACEHOLDER:
            _, _, _, by_type = get_depot_data(active_depot)
            types = list(by_type.keys())

        for idx in range(safe_count):
            type_dropdown = ft.Dropdown(
                label=f"Parada {idx + 1}",
                expand=True,
                value=types[idx % len(types)] if types else None,
                options=[ft.dropdown.Option(t) for t in types],
                bgcolor=tokens["input_bg"],
                color=tokens["input_text"],
                text_style=ft.TextStyle(color=tokens["input_text"]),
                label_style=ft.TextStyle(color=tokens["input_label"]),
                border_color=tokens["input_border"],
                focused_border_color=tokens["input_focus"],
                filled=True,
            )
            stay_field = ft.TextField(
                label="Estadia (min)",
                width=150,
                value="30",
                keyboard_type=ft.KeyboardType.NUMBER,
                bgcolor=tokens["input_bg"],
                fill_color=tokens["input_bg"],
                color=tokens["input_text"],
                text_style=ft.TextStyle(color=tokens["input_text"]),
                label_style=ft.TextStyle(color=tokens["input_label"]),
                border_color=tokens["input_border"],
                focused_border_color=tokens["input_focus"],
                filled=True,
            )
            row = ft.Container(
                padding=12,
                border_radius=12,
                bgcolor=tokens["row_bg"],
                border=ft.border.all(1, tokens["row_border"]),
                content=ft.Row(
                    spacing=12,
                    controls=[type_dropdown, stay_field],
                ),
                data={"type": type_dropdown, "stay": stay_field},
            )
            stops_column.controls.append(row)
        style_stop_rows()
        page.update()

    def on_depot_change(_: ft.ControlEvent) -> None:
        if depot_dropdown.value and depot_dropdown.value != PLACEHOLDER:
            _, _, _, by_type = get_depot_data(depot_dropdown.value)
            count = int(stops_input.value) if stops_input.value and stops_input.value.isdigit() else len(by_type)
            refresh_stop_rows(count)
        result_column.visible = False
        result_column.opacity = 0
        result_column.offset = ft.Offset(0, 0.03)
        export_status.visible = False
        status_text.visible = False
        page.update()

    def on_stop_count_change(_: ft.ControlEvent) -> None:
        value = stops_input.value.strip() if stops_input.value else ""
        count = int(value) if value.isdigit() else len(stops_column.controls) or 1
        refresh_stop_rows(count)
        result_column.visible = False
        result_column.opacity = 0
        result_column.offset = ft.Offset(0, 0.03)
        export_status.visible = False

    def optimize_click(_: ft.ControlEvent) -> None:
        status_text.visible = False
        export_status.visible = False
        result_column.visible = False
        result_column.opacity = 0
        result_column.offset = ft.Offset(0, 0.03)

        selected_depot = depot_dropdown.value
        if not selected_depot or selected_depot == PLACEHOLDER:
            status_text.value = "Selecciona un punto de partida antes de optimizar."
            status_text.visible = True
            page.update()
            return

        start_minutes = parse_hhmm(start_time_input.value or "")
        if start_minutes is None:
            status_text.value = "Hora invalida. Usa formato HH:MM."
            status_text.visible = True
            page.update()
            return

        type_order = selected_types()
        stay_times = selected_stays()
        if not type_order:
            status_text.value = "Configura al menos una parada."
            status_text.visible = True
            page.update()
            return

        nodes, dist, time_mat, by_type = get_depot_data(selected_depot)
        criterion = next(iter(criterion_picker.selected), "distancia")
        cost_matrix = dist if criterion == "distancia" else time_mat if criterion == "tiempo" else dist

        best_route, _, schedule = optimize(
            nodes,
            by_type,
            cost_matrix,
            type_order,
            criterion,
            dist_matrix=dist,
            stay_times=stay_times,
            time_matrix=time_mat,
            start_minutes=start_minutes,
        )

        if best_route is None:
            status_text.value = "No existe ruta factible con los horarios y tiempos indicados."
            status_text.visible = True
            page.update()
            return

        total_dist = sum(dist.loc[best_route[i], best_route[i + 1]] for i in range(len(best_route) - 1))
        total_travel = sum(time_mat.loc[best_route[i], best_route[i + 1]] for i in range(len(best_route) - 1)) / 60
        total_stay = sum(stay_times)
        total_time = total_travel + total_stay
        total_score = sum(nodes.loc[n, "puntaje"] for n in best_route[1:-1])
        route_plot = build_route_plot(nodes, best_route)
        is_light = theme_switch.value
        result_panel_bg = "#FFFFFF" if is_light else "#0D162A"
        result_panel_border = "#D7E2F0" if is_light else "#253658"
        result_title_color = "#0F172A" if is_light else "#E2E8F0"
        result_head_bg = "#EEF4FC" if is_light else "#13203A"
        result_head_text = "#1E3A8A" if is_light else "#BFDBFE"
        result_cell_text = "#0F172A" if is_light else "#E2E8F0"

        route_rows: list[ft.DataRow] = []
        for step, node_id in enumerate(best_route):
            is_depot = nodes.loc[node_id, "tipo"] == "depot"
            stay_idx = step - 1
            stay_val = 0 if is_depot else stay_times[stay_idx] if 0 <= stay_idx < len(stay_times) else 0
            arrival = schedule[step]["llegada"] if schedule else "-"
            departure = schedule[step]["salida"] if schedule else "-"
            wait_flag = "si" if schedule and schedule[step]["espera"] else ""
            route_rows.append(
                ft.DataRow(
                    cells=[
                        ft.DataCell(ft.Text(str(step), color=result_cell_text)),
                        ft.DataCell(ft.Text(str(node_id), color=result_cell_text)),
                        ft.DataCell(ft.Text(str(nodes.loc[node_id, "nombre"]), color=result_cell_text)),
                        ft.DataCell(ft.Text(str(nodes.loc[node_id, "tipo"]), color=result_cell_text)),
                        ft.DataCell(ft.Text(f"{nodes.loc[node_id, 'puntaje']:.1f}", color=result_cell_text)),
                        ft.DataCell(ft.Text(str(stay_val), color=result_cell_text)),
                        ft.DataCell(ft.Text(arrival or "-", color=result_cell_text)),
                        ft.DataCell(ft.Text(departure or "-", color=result_cell_text)),
                        ft.DataCell(ft.Text(wait_flag, color=result_cell_text)),
                    ]
                )
            )

        leg_rows: list[ft.DataRow] = []
        legs_data: list[dict[str, str | float]] = []
        for i in range(len(best_route) - 1):
            orig = best_route[i]
            dest = best_route[i + 1]
            dist_km = round(float(dist.loc[orig, dest] / 1000), 2)
            time_min = round(float(time_mat.loc[orig, dest] / 60), 1)
            leg_rows.append(
                ft.DataRow(
                    cells=[
                        ft.DataCell(ft.Text(str(nodes.loc[orig, "nombre"]), color=result_cell_text)),
                        ft.DataCell(ft.Text(str(nodes.loc[dest, "nombre"]), color=result_cell_text)),
                        ft.DataCell(ft.Text(f"{dist_km:.2f}", color=result_cell_text)),
                        ft.DataCell(ft.Text(f"{time_min:.0f}", color=result_cell_text)),
                    ]
                )
            )
            legs_data.append(
                {
                    "Desde": str(nodes.loc[orig, "nombre"]),
                    "Hasta": str(nodes.loc[dest, "nombre"]),
                    "Distancia (km)": dist_km,
                    "Tiempo (min)": time_min,
                }
            )

        route_data: list[dict[str, str | float | int]] = []
        for step, node_id in enumerate(best_route):
            is_depot = nodes.loc[node_id, "tipo"] == "depot"
            stay_idx = step - 1
            stay_val = 0 if is_depot else stay_times[stay_idx] if 0 <= stay_idx < len(stay_times) else 0
            arrival = schedule[step]["llegada"] if schedule else "-"
            departure = schedule[step]["salida"] if schedule else "-"
            wait_flag = "si" if schedule and schedule[step]["espera"] else ""
            route_data.append(
                {
                    "Paso": step,
                    "ID": int(node_id),
                    "Nombre": str(nodes.loc[node_id, "nombre"]),
                    "Tipo": str(nodes.loc[node_id, "tipo"]),
                    "Puntaje": float(nodes.loc[node_id, "puntaje"]),
                    "Estadia (min)": int(stay_val),
                    "Llegada": arrival or "-",
                    "Salida": departure or "-",
                    "Espera": wait_flag,
                }
            )

        route_df = pd.DataFrame(route_data)
        legs_df = pd.DataFrame(legs_data)
        latest_result["route_df"] = route_df
        latest_result["legs_df"] = legs_df
        latest_result["summary"] = {
            "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "depot": selected_depot,
            "criterion": criterion,
            "distance": f"{(total_dist / 1000):.1f} km",
            "travel": f"{total_travel:.0f} min",
            "stay": f"{total_stay:.0f} min",
            "total_time": f"{total_time:.0f} min",
            "score": f"{total_score:.1f} pts",
        }

        by_type_controls: list[ft.Control] = []
        tile_tokens = theme_tokens(theme_switch.value)
        for tipo in type_order:
            node_ids = by_type[tipo]
            selected = [n for n in best_route[1:-1] if nodes.loc[n, "tipo"] == tipo]
            selected_names = ", ".join(str(nodes.loc[n, "nombre"]) for n in selected) or "-"
            type_rows = [
                ft.DataRow(
                    cells=[
                        ft.DataCell(ft.Text(str(nodes.loc[nid, "nombre"]), color=tile_tokens["table_text"])),
                        ft.DataCell(ft.Text(f"{nodes.loc[nid, 'puntaje']:.1f}", color=tile_tokens["table_text"])),
                    ]
                )
                for nid in sorted(node_ids, key=lambda nid: nodes.loc[nid, "puntaje"], reverse=True)
            ]
            by_type_controls.append(
                ft.ExpansionTile(
                    title=ft.Text(f"{tipo} | seleccionados: {selected_names}", color=tile_tokens["tile_text"], size=13),
                    bgcolor=tile_tokens["tile_bg"],
                    collapsed_bgcolor=tile_tokens["tile_bg"],
                    text_color=tile_tokens["tile_text"],
                    collapsed_text_color=tile_tokens["tile_text"],
                    icon_color=tile_tokens["tile_icon"],
                    collapsed_icon_color=tile_tokens["tile_icon"],
                    controls=[
                        ft.DataTable(
                            heading_row_color=tile_tokens["table_head_bg"],
                            columns=[
                                ft.DataColumn(ft.Text("Nombre", color=tile_tokens["tile_text"])),
                                ft.DataColumn(ft.Text("Puntaje", color=tile_tokens["tile_text"])),
                            ],
                            rows=type_rows,
                        )
                    ],
                )
            )

        result_column.controls = [
            ft.Row(
                wrap=True,
                spacing=10,
                controls=[
                    metric_card("Distancia total", f"{(total_dist / 1000):.1f} km", "#A5F3FC"),
                    metric_card("Tiempo de traslado", f"{total_travel:.0f} min", "#C4B5FD"),
                    metric_card("Tiempo en lugares", f"{total_stay:.0f} min", "#FBCFE8"),
                    metric_card("Tiempo total", f"{total_time:.0f} min", "#FDE68A"),
                    metric_card("Puntaje total", f"{total_score:.1f} pts", "#86EFAC"),
                ],
            ),
            ft.Row(
                wrap=True,
                spacing=10,
                controls=[
                    ft.OutlinedButton(
                        content=ft.Text("Exportar CSV"),
                        icon=ft.Icons.DOWNLOAD,
                        on_click=export_csv_click,
                    ),
                    ft.OutlinedButton(
                        content=ft.Text("Exportar PDF"),
                        icon=ft.Icons.PICTURE_AS_PDF,
                        on_click=export_pdf_click,
                    ),
                    export_status,
                ],
            ),
            ft.Container(
                margin=ft.margin.only(top=8),
                padding=14,
                border_radius=16,
                bgcolor=result_panel_bg,
                border=ft.border.all(1, result_panel_border),
                content=ft.Column(
                    spacing=12,
                    controls=[
                        ft.Text("Visual del recorrido", size=18, weight=ft.FontWeight.W_700, color=result_title_color),
                        ft.Container(
                            height=340,
                            border_radius=12,
                            clip_behavior=ft.ClipBehavior.HARD_EDGE,
                            content=ft.Image(
                                src=base64.b64decode(route_plot),
                                fit=ft.BoxFit.FIT_WIDTH,
                                width=1200,
                                height=340,
                                anti_alias=True,
                            ),
                        ),
                    ],
                ),
            ),
            ft.Container(
                padding=14,
                border_radius=16,
                bgcolor=result_panel_bg,
                border=ft.border.all(1, result_panel_border),
                content=ft.Column(
                    spacing=8,
                    controls=[
                        ft.Text("Recorrido optimizado", size=18, weight=ft.FontWeight.W_700, color=result_title_color),
                        table_scroller(
                            ft.DataTable(
                                data_row_min_height=42,
                                data_row_max_height=48,
                                heading_row_color=result_head_bg,
                                columns=[
                                    ft.DataColumn(ft.Text("Paso", color=result_head_text)),
                                    ft.DataColumn(ft.Text("ID", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Nombre", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Tipo", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Puntaje", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Estadia", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Llegada", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Salida", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Espera", color=result_head_text)),
                                ],
                                rows=route_rows,
                            )
                        ),
                    ],
                ),
            ),
            ft.Container(
                padding=14,
                border_radius=16,
                bgcolor=result_panel_bg,
                border=ft.border.all(1, result_panel_border),
                content=ft.Column(
                    spacing=8,
                    controls=[
                        ft.Text("Detalle tramo a tramo", size=18, weight=ft.FontWeight.W_700, color=result_title_color),
                        table_scroller(
                            ft.DataTable(
                                heading_row_color=result_head_bg,
                                columns=[
                                    ft.DataColumn(ft.Text("Desde", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Hasta", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Distancia (km)", color=result_head_text)),
                                    ft.DataColumn(ft.Text("Tiempo (min)", color=result_head_text)),
                                ],
                                rows=leg_rows,
                            )
                        ),
                    ],
                ),
            ),
            ft.Container(
                padding=14,
                border_radius=16,
                bgcolor=result_panel_bg,
                border=ft.border.all(1, result_panel_border),
                content=ft.Column(
                    spacing=8,
                    controls=[
                        ft.Text("Nodos disponibles por tipo", size=18, weight=ft.FontWeight.W_700, color=result_title_color),
                        ft.Column(spacing=8, controls=by_type_controls),
                    ],
                ),
            ),
        ]

        result_column.visible = True
        result_column.opacity = 1
        result_column.offset = ft.Offset(0, 0)
        page.update()

    depot_dropdown.on_change = on_depot_change
    stops_input.on_blur = on_stop_count_change

    refresh_stop_rows(max(1, len(available_types)))

    hero = ft.Container(
        padding=28,
        border_radius=22,
        gradient=ft.LinearGradient(
            begin=ft.Alignment.TOP_LEFT,
            end=ft.Alignment.BOTTOM_RIGHT,
            colors=["#0891B2", "#1E40AF", "#111827"],
        ),
        animate_opacity=ft.Animation(400, ft.AnimationCurve.EASE_OUT),
        animate_offset=ft.Animation(500, ft.AnimationCurve.DECELERATE),
        opacity=0,
        offset=ft.Offset(0, 0.03),
        content=ft.Column(
            spacing=8,
            controls=[
                ft.Text(
                    "Optimizador de Recorrido Turistico",
                    size=34,
                    weight=ft.FontWeight.W_700,
                    color="#F8FAFC",
                ),
                ft.Text(
                    "Planifica una ruta eficiente por Rosario con una interfaz renovada en Flet.",
                    size=15,
                    color="#DCEBFF",
                ),
            ],
        ),
    )

    config_card = ft.Container(
        margin=ft.margin.only(top=16),
        padding=18,
        border_radius=18,
        gradient=ft.LinearGradient(
            begin=ft.Alignment.TOP_LEFT,
            end=ft.Alignment.BOTTOM_RIGHT,
            colors=["#0E172D", "#0A1324"],
        ),
        border=ft.border.all(1, "#2B3C60"),
        shadow=ft.BoxShadow(
            blur_radius=18,
            spread_radius=0,
            color="#18000000",
            offset=ft.Offset(0, 8),
        ),
        animate_opacity=ft.Animation(500, ft.AnimationCurve.EASE_OUT),
        animate_offset=ft.Animation(550, ft.AnimationCurve.DECELERATE),
        opacity=0,
        offset=ft.Offset(0, 0.04),
        content=ft.Column(
            spacing=14,
            controls=[
                ft.Text("Configuracion", size=22, weight=ft.FontWeight.W_700, color="#E2E8F0"),
                ft.Row(wrap=True, spacing=14, controls=[depot_dropdown, stops_input, start_time_input]),
                ft.Text("Criterio de optimizacion", size=15, weight=ft.FontWeight.W_600, color="#C7D2FE"),
                ft.Row(wrap=True, spacing=16, controls=[criterion_picker, theme_switch]),
                ft.Divider(color="#1F2A44"),
                ft.Text("Orden de actividades y tiempo de estadia", size=16, weight=ft.FontWeight.W_600, color="#E2E8F0"),
                stops_column,
                ft.ElevatedButton(
                    "Optimizar recorrido",
                    height=46,
                    style=ft.ButtonStyle(
                        bgcolor="#2563EB",
                        color="#F8FAFC",
                        shape=ft.RoundedRectangleBorder(radius=10),
                    ),
                    on_click=optimize_click,
                ),
                status_text,
            ],
        ),
    )

    content_column = ft.Column(
        spacing=14,
        controls=[hero, config_card, result_column, ft.Container(height=16)],
    )

    content_container = ft.Container(
        padding=ft.padding.symmetric(horizontal=24, vertical=16),
        content=content_column,
    )

    page.add(content_container)

    def apply_theme(is_light: bool) -> None:
        if is_light:
            page.theme_mode = ft.ThemeMode.LIGHT
            page.bgcolor = "#F3F6FB"
            hero.gradient = ft.LinearGradient(
                begin=ft.Alignment.TOP_LEFT,
                end=ft.Alignment.BOTTOM_RIGHT,
                colors=["#38BDF8", "#2563EB", "#DCE9FF"],
            )
            hero.content.controls[0].color = "#0F172A"
            hero.content.controls[1].color = "#1E293B"
            config_card.gradient = ft.LinearGradient(
                begin=ft.Alignment.TOP_LEFT,
                end=ft.Alignment.BOTTOM_RIGHT,
                colors=["#FFFFFF", "#ECF3FF"],
            )
            config_card.border = ft.border.all(1, "#BFD2F2")
            config_card.content.controls[0].color = "#0F172A"
            config_card.content.controls[2].color = "#1E3A8A"
            config_card.content.controls[4].color = "#C9D9F0"
            config_card.content.controls[5].color = "#0F172A"
            status_text.color = "#B91C1C"
        else:
            page.theme_mode = ft.ThemeMode.DARK
            page.bgcolor = "#090E1A"
            hero.gradient = ft.LinearGradient(
                begin=ft.Alignment.TOP_LEFT,
                end=ft.Alignment.BOTTOM_RIGHT,
                colors=["#0891B2", "#1E40AF", "#111827"],
            )
            hero.content.controls[0].color = "#F8FAFC"
            hero.content.controls[1].color = "#DCEBFF"
            config_card.gradient = ft.LinearGradient(
                begin=ft.Alignment.TOP_LEFT,
                end=ft.Alignment.BOTTOM_RIGHT,
                colors=["#0E172D", "#0A1324"],
            )
            config_card.border = ft.border.all(1, "#2B3C60")
            config_card.content.controls[0].color = "#E2E8F0"
            config_card.content.controls[2].color = "#C7D2FE"
            config_card.content.controls[4].color = "#1F2A44"
            config_card.content.controls[5].color = "#E2E8F0"
            status_text.color = "#FCA5A5"

        style_stop_rows()

        page.update()

    def on_theme_change(_: ft.ControlEvent) -> None:
        apply_theme(theme_switch.value)

    def on_resize(_: ft.ControlEvent | None = None) -> None:
        width = page.width or 1100
        is_mobile = width < 920
        content_container.padding = ft.padding.symmetric(
            horizontal=12 if is_mobile else 24,
            vertical=12 if is_mobile else 16,
        )
        hero.padding = 20 if is_mobile else 28
        hero.content.controls[0].size = 24 if is_mobile else 34
        hero.content.controls[1].size = 13 if is_mobile else 15
        depot_dropdown.width = width - 60 if is_mobile else 360
        stops_input.width = width - 60 if is_mobile else 190
        start_time_input.width = width - 60 if is_mobile else 200
        page.update()

    async def intro_animation() -> None:
        await asyncio.sleep(0.05)
        hero.opacity = 1
        hero.offset = ft.Offset(0, 0)
        page.update()
        await asyncio.sleep(0.06)
        config_card.opacity = 1
        config_card.offset = ft.Offset(0, 0)
        page.update()

    page.on_resize = on_resize
    theme_switch.on_change = on_theme_change
    apply_theme(theme_switch.value)
    on_resize()
    page.run_task(intro_animation)


if __name__ == "__main__":
    ft.app(target=main)
