# 02 - Requerimientos y Casos de Uso

## Requerimientos funcionales

- RF-01: Cargar nodos y matrices de un depot.
- RF-02: Permitir seleccionar criterio de optimizacion.
- RF-03: Permitir definir orden de tipos a visitar.
- RF-04: Permitir definir tiempo de estadia por parada.
- RF-05: Validar hora de salida y formato.
- RF-06: Calcular ruta optimizada.
- RF-07: Verificar factibilidad con ventanas horarias.
- RF-08: Mostrar mapa y tablas de detalle.
- RF-09: Exportar resultado en CSV.
- RF-10: Exportar resultado en PDF.

## Requerimientos no funcionales

- RNF-01: Interfaz clara en escritorio, web y movil.
- RNF-02: Tiempo de respuesta aceptable para datasets academicos.
- RNF-03: Trazabilidad de insumos y salidas.
- RNF-04: Codigo mantenible y modular.
- RNF-05: Diagramas y documentacion actualizados.

## Casos de uso principales

- CU-01: Configurar recorrido.
- CU-02: Optimizar recorrido.
- CU-03: Analizar detalle por tramo.
- CU-04: Exportar reportes.
- CU-05: Ajustar parametros y reoptimizar.

Diagrama UML asociado: [10-UML-y-Casos-de-Uso.md](10-UML-y-Casos-de-Uso.md)

## Actores

- Actor primario: Usuario operador.
- Actor secundario: Sistema de datos Google (fuente de enriquecimiento).

## Reglas de negocio clave

- RB-01: La ruta inicia y termina en depot.
- RB-02: Se selecciona una parada por tipo configurado.
- RB-03: En criterio puntaje, desempata por distancia.
- RB-04: Si una parada viola horario, la ruta se descarta.
