# Webmining - Documentacion Integral de Ingenieria de Software

Este paquete contiene la documentacion completa del proyecto con foco en arquitectura, analisis funcional, datos y operacion.

## Alcance

- DFD nivel 0, nivel 1 y nivel 2.
- Arquitectura del sistema (contexto, contenedores, componentes y despliegue).
- Diagramas de flujo funcional y tecnico.
- UML: casos de uso, clases de dominio y secuencia principal.
- DER (modelo entidad-relacion) del dominio de datos.
- Pipeline de obtencion de datos con Google Maps.
- Presentacion ejecutiva de 20 minutos.

## Estructura

- [01-Vision-y-Alcance.md](01-Vision-y-Alcance.md)
- [02-Requerimientos-y-Casos-de-Uso.md](02-Requerimientos-y-Casos-de-Uso.md)
- [03-Arquitectura.md](03-Arquitectura.md)
- [04-DFD-Nivel-0-a-2.md](04-DFD-Nivel-0-a-2.md)
- [05-Modelo-de-Datos-y-DER.md](05-Modelo-de-Datos-y-DER.md)
- [06-Flujos-y-Procesos.md](06-Flujos-y-Procesos.md)
- [07-Calidad-Seguridad-y-Operacion.md](07-Calidad-Seguridad-y-Operacion.md)
- [08-Plan-de-Evolucion.md](08-Plan-de-Evolucion.md)
- [09-Presentacion-20min.md](09-Presentacion-20min.md)
- [10-UML-y-Casos-de-Uso.md](10-UML-y-Casos-de-Uso.md)

## Diagramas incluidos

Carpeta: [diagrams](diagrams)

- DFD: nivel 0, nivel 1, nivel 2.
- Arquitectura: contexto, contenedores, componentes, despliegue.
- Datos: DER.
- Procesos: flujo de optimizacion, flujo de obtencion de datos Google.
- UML: casos de uso, clases de dominio, secuencia de optimizacion.
- Diagrama visual de plataforma con Python Diagrams.

## Generacion de artefactos

Desde la raiz del proyecto:

1. Generar todos los diagramas modernos y a color:

   make software-docs-diagrams

2. Generar presentaciones PowerPoint (ejecutiva, docente y tribunal tecnico):

   uv run python docs/software-engineering/generate_presentation.py

3. Generar sitio HTML navegable:

   make software-docs-site

4. Generar todo junto:

   make software-docs

## Artefactos finales esperados

- Imagenes SVG/PNG en [docs/software-engineering/diagrams](diagrams).
- Presentaciones en:
   - [docs/software-engineering/Webmining_Arquitectura_y_Datos_Google_20min.pptx](Webmining_Arquitectura_y_Datos_Google_20min.pptx)
   - [docs/software-engineering/Webmining_Presentacion_Docente_20min.pptx](Webmining_Presentacion_Docente_20min.pptx)
   - [docs/software-engineering/Webmining_Presentacion_Tribunal_Tecnico_20min.pptx](Webmining_Presentacion_Tribunal_Tecnico_20min.pptx)
- Sitio HTML en [docs/software-engineering/site/index.html](site/index.html).

## Audiencia objetivo

- Equipo tecnico (desarrollo, datos, arquitectura).
- Stakeholders de negocio.
- Docentes/evaluadores del proceso de ingenieria de software.
