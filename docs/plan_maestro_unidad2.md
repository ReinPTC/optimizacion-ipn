# Plan Maestro Curricular: Unidad 2 — Programación Lineal

**Unidad de Aprendizaje:** Optimización  
**Institución:** Instituto Politécnico Nacional (IPN) — Escuela Superior de Física y Matemáticas (ESFM)  
**Licenciatura:** Matemática Algorítmica  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming*, 4ta Edición, Springer (2016).  
**Unidad de Competencia:** Formula, analiza y resuelve problemas de optimización lineal empleando el método simplex, la teoría de dualidad y el análisis de sensibilidad con rigor algebraico y computacional.

---

## 1. Mapeo Oficial: Temario del IPN $\longleftrightarrow$ Capítulos de Luenberger & Ye (4ta Ed.)

```
Unidad 2: Programación Lineal

2.1 Propiedades Básicas de los Programas Lineales (Luenberger & Ye, Capítulo 2)
├── 2.1.1 Introducción y Formas Canónica y Estándar (Sec. 2.1, pp. 11–14)
├── 2.1.2 Ejemplos Canónicos de Programación Lineal (Sec. 2.2, pp. 14–19)
├── 2.1.3 Soluciones Básicas y Soluciones Básicas Factibles (Sec. 2.3, pp. 19–21)
├── 2.1.4 El Teorema Fundamental de la Programación Lineal (Sec. 2.4, pp. 21–24)
└── 2.1.5 Relaciones con la Convexidad y Puntos Extremos (Sec. 2.5, pp. 24–28)

2.2 El Método Simplex (Luenberger & Ye, Capítulo 3)
├── 2.2.1 Pivotes y Transformaciones Elementales (Sec. 3.1, pp. 33–37)
├── 2.2.2 Puntos Extremos Adyacentes y Vértices Vecinos (Sec. 3.2, pp. 38–41)
├── 2.2.3 Determinación de una Solución Factible Mínima y Costos Reducidos (Sec. 3.3, pp. 41–45)
├── 2.2.4 Procedimiento Computacional: Algoritmo Simplex Tabular (Sec. 3.4, pp. 45–52)
├── 2.2.5 Determinación de una SBF Inicial: Dos Fases y Gran M (Sec. 3.5, pp. 52–57)
├── 2.2.6 Forma Matricial del Método Simplex (Sec. 3.6, pp. 57–60)
└── 2.2.7 *El Método Simplex Revisado y Descomposición LU (Sec. 3.7, pp. 60–66)

2.3 Dualidad y Complementariedad (Luenberger & Ye, Capítulo 4)
├── 2.3.1 Programas Lineales Duales y Tablas de Transformación (Sec. 4.1, pp. 75–81)
├── 2.3.2 El Teorema de Dualidad: Débil, Fuerte y Lema de Farkas (Sec. 4.2, pp. 81–86)
├── 2.3.3 Relaciones con el Procedimiento Simplex y Multiplicadores (Sec. 4.3, pp. 86–90)
├── 2.3.4 Sensibilidad, Precios Sombra y Holgura Complementaria (Sec. 4.4, pp. 90–97)
├── 2.3.5 El Método Simplex Dual (Sec. 4.6, pp. 102–108)
└── 2.3.6 *El Algoritmo Primal-Dual (Sec. 4.7, pp. 108–114)

2.4 Pendiente
└── Contenido reservado para asignación y desarrollo específico del curso
```

---

## 2. Desglose Teórico, Teoremas y Bibliografía Oficial por Subtema

| Subtema Oficial | Contenido Teórico y Teoremas Clave | Referencia Exacta Luenberger (4ta Ed.) | Textos Alternativos Oficiales (IPN) |
| :--- | :--- | :--- | :--- |
| **2.1.1 Introducción y Formas Canónica y Estándar** | Formulación general, paso de desigualdades a igualdades ($\mathbf{A}\mathbf{x}=\mathbf{b}$ con $\mathbf{x}\ge\mathbf{0}$) mediante holguras $s_i$ y excesos $e_i$, descomposición de variables irrestrictas $x_j = x_j^+ - x_j^-$, inversión de objetivos $\max \mathbf{c}^T\mathbf{x} = -\min(-\mathbf{c})^T\mathbf{x}$. | Cap. 2, Sec. 2.1 (pp. 11–14) | Bazaraa, Jarvis & Sherali (2011) Cap. 1; Cottle & Thapa (2017) Cap. 3 |
| **2.1.2 Ejemplos Canónicos de Programación Lineal** | Modelo de la dieta de Stigler (1945), problema de transporte de Hitchcock-Koopmans, mezcla de hidrocarburos, planeación de manufactura multi-período con inventarios dinámicos. | Cap. 2, Sec. 2.2 (pp. 14–19) | Bazaraa, Jarvis & Sherali (2011) Cap. 1; Dantzig (1963) Cap. 3 |
| **2.1.3 Soluciones Básicas y SBF** | Definición formal de base $\mathbf{B}$ de $\mathbf{A} \in \mathbb{R}^{m \times n}$ con $\operatorname{rango}(\mathbf{A})=m$, solución básica $\mathbf{x}_B = \mathbf{B}^{-1}\mathbf{b}, \mathbf{x}_N = \mathbf{0}$; Solución Básica Factible (SBF) si $\mathbf{x}_B \ge \mathbf{0}$; soluciones degeneradas ($x_{Bi} = 0$). | Cap. 2, Sec. 2.3 (pp. 19–21) | Cottle & Thapa (2017) Cap. 3; Bazaraa et al. (2011) Cap. 2 |
| **2.1.4 Teorema Fundamental de la PL** | **Teorema Fundamental de la Programación Lineal:** (1) Si existe una solución factible, existe al menos una SBF; (2) Si existe una solución factible óptima, existe una SBF que es óptima. Demostración analítica por inducción y reducción de variables activas. | Cap. 2, Sec. 2.4 (pp. 21–24) | Bazaraa, Jarvis & Sherali (2011) Cap. 2; Luenberger (2016) pp. 21–24 |
| **2.1.5 Relaciones con la Convexidad** | **Teorema de Equivalencia:** Un vector $\mathbf{x}$ es una SBF del poliedro $P = \{\mathbf{x} \in \mathbb{R}^n : \mathbf{A}\mathbf{x}=\mathbf{b}, \mathbf{x}\ge\mathbf{0}\}$ si y solo si $\mathbf{x}$ es un punto extremo de $P$. Corolarios sobre existencia de óptimos en politopos compactos. | Cap. 2, Sec. 2.5 (pp. 24–28) | Cottle & Thapa (2017) Cap. 4; Bazaraa et al. (2011) Cap. 2 |
| **2.2.1 Pivotes y Transformaciones** | Álgebra de operaciones elementales por fila, actualización de la matriz aumentada $[\mathbf{B}^{-1}\mathbf{A} \mid \mathbf{B}^{-1}\mathbf{b}]$, intercambio de columnas básicas, propiedades de no singularidad de la nueva base. | Cap. 3, Sec. 3.1 (pp. 33–37) | Bazaraa, Jarvis & Sherali (2011) Cap. 3; Cottle & Thapa (2017) Cap. 3 |
| **2.2.2 Puntos Extremos Adyacentes** | Movimiento a lo largo de aristas del poliedro $P$, vector de dirección simplex $\mathbf{d}$, adyacencia de bases que comparten $m-1$ columnas, caracterización geométrica de las transiciones simplex. | Cap. 3, Sec. 3.2 (pp. 38–41) | Bazaraa, Jarvis & Sherali (2011) Cap. 3; Nocedal & Wright (2006) Cap. 13 |
| **2.2.3 Determinación de una SBF Mínima** | Costos reducidos relativos $r_j = c_j - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{a}_j = c_j - z_j$. Condición de parada y optimalidad ($r_j \ge 0 \ \forall j$). Detección analítica de problema no acotado ($r_j < 0$ con $\mathbf{y}_j \le \mathbf{0}$). | Cap. 3, Sec. 3.3 (pp. 41–45) | Cottle & Thapa (2017) Cap. 5; Bazaraa et al. (2011) Cap. 3 |
| **2.2.4 Algoritmo Simplex Tabular** | Construcción del tablero canónico, regla del costo reducido más negativo de Dantzig, regla de la razón mínima $\theta^* = \min_{i: y_{ij}>0} \{x_{Bi} / y_{ij}\}$, resolución de empates, prevención de ciclaje y regla del subíndice mínimo de Bland. | Cap. 3, Sec. 3.4 (pp. 45–52) | Bazaraa, Jarvis & Sherali (2011) Cap. 3; Cottle & Thapa (2017) Cap. 5 |
| **2.2.5 Métodos de Inicialización (Dos Fases y Gran M)** | Introducción de variables artificiales $\mathbf{x}_a$. Fase I: Minimizar $\sum x_{ai}$; detección de infactibilidad ($W^* > 0$) y transición a Fase II. Método de la Gran M con penalización $M \gg 0$. | Cap. 3, Sec. 3.5 (pp. 52–57) | Bazaraa, Jarvis & Sherali (2011) Cap. 4; Cottle & Thapa (2017) Cap. 5 |
| **2.2.6 Forma Matricial del Método Simplex** | Representación compacta del tablero mediante $\mathbf{B}^{-1}$, cálculo explícito de multiplicadores $\mathbf{y}^T = \mathbf{c}_B^T\mathbf{B}^{-1}$ y vector de recursos actualizado $\bar{\mathbf{b}} = \mathbf{B}^{-1}\mathbf{b}$. | Cap. 3, Sec. 3.6 (pp. 57–60) | Nocedal & Wright (2006) Cap. 13; Cottle & Thapa (2017) Cap. 5 |
| **2.2.7 Simplex Revisado y Factorización LU** | Ahorro en complejidad temporal y espacial, almacenamiento de factores $\mathbf{B} = \mathbf{L}\mathbf{U}$, matrices eta de actualización elemental, estabilidad numérica frente al Simplex estándar. | Cap. 3, Sec. 3.7 (pp. 60–66) | Nocedal & Wright (2006) Cap. 13; Bazaraa et al. (2011) Cap. 5 |
| **2.3.1 Programas Lineales Duales** | Relaciones simétricas y asimétricas primal-dual, formulación canónica: Primal $\min \mathbf{c}^T\mathbf{x}$ s.a. $\mathbf{A}\mathbf{x} \ge \mathbf{b}, \mathbf{x}\ge\mathbf{0} \Longleftrightarrow$ Dual $\max \mathbf{b}^T\mathbf{y}$ s.a. $\mathbf{A}^T\mathbf{y} \le \mathbf{c}, \mathbf{y}\ge\mathbf{0}$. Tabla formal de transformaciones de restricciones y variables. | Cap. 4, Sec. 4.1 (pp. 75–81) | Bazaraa, Jarvis & Sherali (2011) Cap. 6; Cottle & Thapa (2017) Cap. 6 |
| **2.3.2 El Teorema de Dualidad** | **Lema de Dualidad Débil:** $\mathbf{b}^T\mathbf{y} \le \mathbf{c}^T\mathbf{x}$ para cualquier par factible $(\mathbf{x}, \mathbf{y})$. **Teorema de Dualidad Fuerte:** Si uno tiene solución óptima, ambos tienen solución óptima y $\mathbf{c}^T\mathbf{x}^* = \mathbf{b}^T\mathbf{y}^*$. Corolarios sobre infactibilidad y no acotamiento. Lema de Farkas. | Cap. 4, Sec. 4.2 (pp. 81–86) | Cottle & Thapa (2017) Cap. 6; Bazaraa et al. (2011) Cap. 6 |
| **2.3.3 Relaciones con el Procedimiento Simplex** | Multiplicadores simplex $\mathbf{y}^T = \mathbf{c}_B^T\mathbf{B}^{-1}$ como solución dual óptima. Los costos reducidos $r_j = c_j - \mathbf{y}^T\mathbf{a}_j$ corresponden a las variables de holgura duales. Lectura directa del óptimo dual desde el tablero final primal. | Cap. 4, Sec. 4.3 (pp. 86–90) | Bazaraa, Jarvis & Sherali (2011) Cap. 6; Luenberger (2016) pp. 86–90 |
| **2.3.4 Sensibilidad y Holgura Complementaria** | **Teorema de Holgura Complementaria:** Condición necesaria y suficiente de optimalidad: $x_j (c_j - \mathbf{y}^T\mathbf{a}_j) = 0$ y $y_i (\mathbf{a}_i^T\mathbf{x} - b_i) = 0$. Precios sombra (shadow prices): $\nabla_{\mathbf{b}} z^*(\mathbf{b}) = \mathbf{y}^*$. Intervalos de variación permitida para $\mathbf{b}$ y $\mathbf{c}$. | Cap. 4, Sec. 4.4 (pp. 90–97) | Bazaraa, Jarvis & Sherali (2011) Cap. 6; Cottle & Thapa (2017) Cap. 7 |
| **2.3.5 El Método Simplex Dual** | Fundamentación algorítmica: mantiene dual-factibilidad ($r_j \ge 0$) mientras busca primal-factibilidad ($\mathbf{x}_B \ge \mathbf{0}$). Criterio de salida (fila más negativa $x_{Br} < 0$), criterio de entrada (razón mínima dual $\min_{j: y_{rj}<0} \{r_j / |y_{rj}|\}$). Adición de nuevas restricciones sin recomputar desde cero. | Cap. 4, Sec. 4.6 (pp. 102–108) | Bazaraa, Jarvis & Sherali (2011) Cap. 6; Cottle & Thapa (2017) Cap. 6 |
| **2.3.6 *El Algoritmo Primal-Dual** | Algoritmo clásico de Dantzig-Ford-Fulkerson: solución dual factible inicial, construcción del primal restringido con variables artificiales para satisfacer complementariedad, actualización de multiplicadores duales. | Cap. 4, Sec. 4.7 (pp. 108–114) | Bazaraa, Jarvis & Sherali (2011) Cap. 6; Papadimitriou & Steiglitz (1998) Cap. 5 |
| **2.4 Pendiente** | Tema reservado según requerimiento del docente/alumno. | — | — |

---

## 3. Bibliografía Oficial del IPN (Libros de Consulta y Texto Base)

1. **Bazaraa, M. S., Jarvis, J. J. & Sherali, H. D. (2011)**. *Linear Programming and Network Flows* (4th ed.), John Wiley & Sons.
2. **Cottle, R. W. & Thapa, M. N. (2017)**. *Linear and Nonlinear Optimization*, Springer.
3. **Dantzig, G. B. (1963)**. *Linear Programming and Extensions*, Princeton University Press.
4. **Luenberger, D. G. & Ye, Y. (2016)**. *Linear and Nonlinear Programming* (4th ed.), Springer. *(Texto Base Obligatorio)*
5. **Nocedal, J. & Wright, S. (2006)**. *Numerical Optimization* (2nd ed.), Springer.
6. **Papadimitriou, C. H. & Steiglitz, K. (1998)**. *Combinatorial Optimization: Algorithms and Complexity*, Dover Publications.
