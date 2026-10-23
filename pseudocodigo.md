# Pseudocódigo del ajuste (fit.py)

Cada paso del programa corresponde a un paso del desarrollo analítico.

```
ALGORITMO AjustePendulo
  ENTRADA: archivo datos.csv con L y tres tiempos de 10 oscilaciones por fila

  // 1. Cargar datos y construir x, y
  PARA cada fila:
      T  ← promedio(t1, t2, t3) / 10          // periodo de una oscilación
      x_i ← L                                  // variable independiente (m)
      y_i ← T²                                 // variable dependiente (s²)

  // 2. Sumas que aparecen en ∇S = 0
  n ← número de datos
  Sx ← Σ x_i ;  Sy ← Σ y_i ;  Sxx ← Σ x_i² ;  Sxy ← Σ x_i y_i

  // 3. Ecuaciones normales (de ∂S/∂a = 0 y ∂S/∂b = 0):
  //      a·Sxx + b·Sx = Sxy
  //      a·Sx  + b·n  = Sy
  D  ← n·Sxx − Sx²                             // determinante del sistema
  SI D = 0 ENTONCES error "todos los x iguales: no hay recta única"
  a* ← (n·Sxy − Sx·Sy) / D                     // regla de Cramer
  b* ← (Sxx·Sy − Sx·Sxy) / D

  // 4. Error mínimo y verificación del punto crítico
  S_min ← Σ (y_i − a*·x_i − b*)²
  grad  ← ( −2 Σ x_i(y_i − a*x_i − b*),  −2 Σ (y_i − a*x_i − b*) )   // debe ser ≈ (0, 0)

  // 5. Hessiano y criterio de la segunda derivada
  H ← [[2·Sxx, 2·Sx], [2·Sx, 2·n]]
  SI H[0][0] > 0 Y det(H) > 0 ENTONCES (a*, b*) es MÍNIMO
  autovalores, autovectores ← eig(H)           // ejes y curvaturas de las elipses de nivel

  // 6. Calidad del ajuste y física
  e_i  ← y_i − (a*·x_i + b*)                   // residuos
  RMSE ← raíz( Σ e_i² / n ) ;  R² ← 1 − Σ e_i² / Σ (y_i − ȳ)²
  g    ← 4π² / a*                              // comparar con 9.78 m/s²
  d    ← b* / a*                               // corrimiento de L implícito (teoría: b = 0)

  // 7. Verificación con bibliotecas
  (a_lib, b_lib) ← numpy.polyfit(x, y, 1)  y  scipy.stats.linregress(x, y)
  tabla de diferencias máximas entre los tres métodos

  // 8. Análisis adicional
  Opción B: PARA cada i: reajustar sin el punto i y registrar a, b, g
  Opción C: ajustar y = c2·x² + c1·x + c0 y comparar su S con S_min

  // 9. Gráficas
  dispersión + recta óptima + 2 rectas peores, superficie S(a,b), trazas,
  curvas de nivel con (a*, b*) y autovectores, secciones en direcciones
  principales, residuos vs x e histograma de residuos
FIN
```
