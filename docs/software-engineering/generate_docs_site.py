from pathlib import Path
from shutil import copy2

import markdown


ROOT = Path(__file__).resolve().parent
SITE_DIR = ROOT / "site"
DIAGRAMS_DIR = ROOT / "diagrams"
SITE_DIAGRAMS_DIR = SITE_DIR / "diagrams"
PRESENTATIONS = [
    "Webmining_Arquitectura_y_Datos_Google_20min.pptx",
    "Webmining_Presentacion_Docente_20min.pptx",
    "Webmining_Presentacion_Tribunal_Tecnico_20min.pptx",
]

DOC_ORDER = [
    "README.md",
    "01-Vision-y-Alcance.md",
    "02-Requerimientos-y-Casos-de-Uso.md",
    "03-Arquitectura.md",
    "04-DFD-Nivel-0-a-2.md",
    "05-Modelo-de-Datos-y-DER.md",
    "06-Flujos-y-Procesos.md",
    "07-Calidad-Seguridad-y-Operacion.md",
    "08-Plan-de-Evolucion.md",
    "09-Presentacion-20min.md",
  "10-UML-y-Casos-de-Uso.md",
]


def html_template(title: str, nav: str, body: str) -> str:
    return f"""<!doctype html>
<html lang=\"es\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>{title}</title>
  <style>
    :root {{
      --bg: #0f172a;
      --bg-soft: #111c34;
      --card: rgba(12, 30, 64, 0.72);
      --card-border: rgba(127, 178, 255, 0.25);
      --title: #e0ecff;
      --text: #d4deef;
      --muted: #9ab0d2;
      --link: #7dd3fc;
      --link-hover: #a5f3fc;
      --accent: #fb7185;
      --shadow: 0 20px 45px rgba(0, 0, 0, 0.35);
    }}

    * {{ box-sizing: border-box; }}

    body {{
      margin: 0;
      font-family: "Avenir Next", "Trebuchet MS", "Segoe UI", sans-serif;
      color: var(--text);
      background:
        radial-gradient(circle at 15% 12%, #24427a 0%, transparent 32%),
        radial-gradient(circle at 84% 18%, #5a2f5b 0%, transparent 26%),
        radial-gradient(circle at 55% 82%, #113749 0%, transparent 32%),
        linear-gradient(130deg, var(--bg) 0%, var(--bg-soft) 100%);
      min-height: 100vh;
    }}

    .layout {{
      display: grid;
      grid-template-columns: 320px 1fr;
      min-height: 100vh;
      backdrop-filter: blur(2px);
    }}

    aside {{
      border-right: 1px solid var(--card-border);
      background: rgba(2, 7, 20, 0.52);
      padding: 24px 18px;
      position: sticky;
      top: 0;
      height: 100vh;
      overflow: auto;
    }}

    .brand {{
      font-weight: 800;
      font-size: 1.2rem;
      letter-spacing: 0.2px;
      color: var(--title);
      margin-bottom: 6px;
    }}

    .subtitle {{
      color: var(--muted);
      font-size: 0.9rem;
      margin-bottom: 20px;
    }}

    nav a {{
      display: block;
      color: var(--text);
      text-decoration: none;
      padding: 10px 12px;
      border: 1px solid transparent;
      border-radius: 10px;
      margin: 3px 0;
      transition: transform 180ms ease, border-color 180ms ease, background 180ms ease;
    }}

    nav a:hover {{
      border-color: var(--card-border);
      background: rgba(25, 51, 95, 0.5);
      transform: translateX(3px);
      color: var(--link-hover);
    }}

    .downloads {{
      margin-top: 22px;
      padding-top: 16px;
      border-top: 1px solid var(--card-border);
      font-size: 0.92rem;
    }}

    .downloads h3 {{
      margin: 0 0 8px 0;
      color: var(--title);
      font-size: 0.96rem;
    }}

    .downloads a {{
      color: var(--link);
      text-decoration: none;
      display: block;
      margin: 4px 0;
    }}

    main {{
      padding: 28px;
    }}

    article {{
      max-width: 980px;
      margin: 0 auto;
      background: var(--card);
      border: 1px solid var(--card-border);
      border-radius: 20px;
      box-shadow: var(--shadow);
      padding: 34px 38px;
      animation: reveal 500ms ease;
    }}

    @keyframes reveal {{
      from {{ opacity: 0; transform: translateY(10px); }}
      to {{ opacity: 1; transform: translateY(0); }}
    }}

    h1, h2, h3 {{ color: var(--title); }}
    h1 {{ font-size: 2rem; margin-top: 0; }}
    h2 {{ margin-top: 2rem; border-left: 4px solid var(--accent); padding-left: 10px; }}
    p, li {{ line-height: 1.62; }}

    a {{ color: var(--link); }}
    a:hover {{ color: var(--link-hover); }}

    img {{
      max-width: 100%;
      border-radius: 12px;
      border: 1px solid var(--card-border);
      margin: 10px 0;
      background: rgba(255, 255, 255, 0.03);
    }}

    code {{
      background: rgba(17, 36, 76, 0.82);
      border: 1px solid rgba(151, 197, 255, 0.3);
      color: #eaf4ff;
      border-radius: 6px;
      padding: 2px 6px;
    }}

    pre code {{
      display: block;
      padding: 12px;
      overflow-x: auto;
    }}

    @media (max-width: 980px) {{
      .layout {{ grid-template-columns: 1fr; }}
      aside {{ position: static; height: auto; }}
      main {{ padding: 16px; }}
      article {{ padding: 22px; }}
    }}
  </style>
</head>
<body>
  <div class=\"layout\">
    <aside>
      <div class=\"brand\">Webmining Docs Hub</div>
      <div class=\"subtitle\">Ingenieria de software, arquitectura y datos Google</div>
      <nav>{nav}</nav>
      <section class=\"downloads\">
        <h3>Presentaciones</h3>
        <a href=\"Webmining_Arquitectura_y_Datos_Google_20min.pptx\">Version ejecutiva</a>
        <a href=\"Webmining_Presentacion_Docente_20min.pptx\">Version docente</a>
        <a href=\"Webmining_Presentacion_Tribunal_Tecnico_20min.pptx\">Version tribunal tecnico</a>
      </section>
    </aside>
    <main>
      <article>{body}</article>
    </main>
  </div>
</body>
</html>
"""


def slug(md_name: str) -> str:
    if md_name == "README.md":
        return "index.html"
    return md_name.replace(".md", ".html")


def build_nav() -> str:
    links = []
    for name in DOC_ORDER:
        target = slug(name)
        label = "Inicio" if name == "README.md" else name.replace(".md", "")
        links.append(f'<a href="{target}">{label}</a>')
    return "\n".join(links)


def copy_assets() -> None:
    SITE_DIAGRAMS_DIR.mkdir(parents=True, exist_ok=True)
    for path in DIAGRAMS_DIR.glob("*.png"):
        copy2(path, SITE_DIAGRAMS_DIR / path.name)
    for path in DIAGRAMS_DIR.glob("*.svg"):
        copy2(path, SITE_DIAGRAMS_DIR / path.name)

    for filename in PRESENTATIONS:
        src = ROOT / filename
        if src.exists():
            copy2(src, SITE_DIR / filename)


def render_document(md_name: str, nav: str) -> None:
    md_path = ROOT / md_name
    html_path = SITE_DIR / slug(md_name)

    md_text = md_path.read_text(encoding="utf-8")
    body = markdown.markdown(md_text, extensions=["fenced_code", "tables", "toc"])

    for source_name in DOC_ORDER:
      body = body.replace(f'href="{source_name}"', f'href="{slug(source_name)}"')

    body = body.replace('href="site/index.html"', 'href="index.html"')
    body = body.replace('src="diagrams/', 'src="diagrams/')
    body = body.replace('href="diagrams/', 'href="diagrams/')

    html = html_template(f"Webmining | {md_name}", nav, body)
    html_path.write_text(html, encoding="utf-8")


def build_site() -> None:
    SITE_DIR.mkdir(parents=True, exist_ok=True)
    nav = build_nav()
    copy_assets()
    for md_name in DOC_ORDER:
        render_document(md_name, nav)
    print(f"Sitio generado en: {SITE_DIR}")


if __name__ == "__main__":
    build_site()
