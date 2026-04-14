# 07 - Calidad, Seguridad y Operacion

## Calidad de software

- Convenciones de codigo y modularidad por responsabilidad.
- Validacion de entradas de usuario.
- Validacion de consistencia de matrices.
- Regeneracion automatica de diagramas.

## Estrategia de testing sugerida

- Unit tests para:
  - parseo de horarios,
  - carga de matrices,
  - reglas de factibilidad,
  - desempates por criterio.
- Integration tests para flujo app -> optimizador -> render.
- Golden tests para dataset de referencia.

## Seguridad

- Gestion de API keys via variables de entorno.
- No hardcodear credenciales en notebooks o codigo.
- Sanitizacion basica de entradas de usuario.
- Registro de auditoria para jobs de carga de datos.

## Observabilidad

- Log estructurado por etapa del pipeline.
- Metricas de performance de optimizacion.
- Metricas de consumo de Google APIs.

## Operacion

- Modos de ejecucion: desktop, web, movil.
- Exportes operativos: CSV y PDF.
- Procedimiento de refresh de dataset por depot.
