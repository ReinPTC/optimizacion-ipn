# Plan de Trabajo y Dosificación de Notas para Alumnos — Unidad 2

**Asignatura:** Optimización  
**Unidad Temática II:** Programación Lineal  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming* (4ta Ed., Springer 2016).  
**Institución:** Instituto Politécnico Nacional (IPN) — Escuela Superior de Física y Matemáticas (ESFM)

---

## 1. Lineamientos Generales para la Elaboración de Notas

### 🎯 Propósito Pedagógico
El objetivo de las notas individuales para la Unidad 2 es consolidar en el estudiante de Matemática Algorítmica el dominio teórico-algebraico del método simplex, la interpretación geométrica de los poliedros convexos, la fundamentación rigurosa de la dualidad y la competencia para programar algoritmos de optimización lineal sin depender de librerías de "caja negra" (e.g. Scipy optimize linprog).

### 📋 Reglas Estrictas de Entrega
1. **Consulta de Textos Alternativos Obligatoria:**  
   Cada reporte debe redactarse consultando **al menos un texto diferente a Luenberger** de la Bibliografía Oficial del IPN (e.g. Bazaraa, Jarvis & Sherali; Cottle & Thapa; Dantzig; Nocedal & Wright). Se debe contrastar la notación y el enfoque pedagógico.
2. **Estructura Requerida de Cada Reporte (Formato Markdown / PDF / Notebook):**
   * **Encabezado Académico:** Nombre del alumno, matrícula, número de reporte, fecha y temas cubiertos.
   * **Marco Teórico y Definiciones Formales:** Definiciones precisas con fórmulas renderizadas en LaTeX.
   * **Demostración Analítica Rigurosa:** Desarrollo paso a paso del teorema asignado sin omitir justificaciones algebraicas.
   * **Ejemplo Práctico Resuelto Paso a Paso:** Resolución numérica completa a mano (tableros simplex, pivoteo explícito).
   * **Implementación Computacional en Python:** Código limpio, vectorizado con NumPy, documentado con *docstrings*, pruebas y visualizaciones.
   * **Pregunta de Reflexión Teórica:** Respuesta argumentada a la pregunta detonadora del reporte.
   * **Bibliografía:** Citas en formato IEEE o APA especificando capítulos y páginas consultadas.

---

## 2. Dosificación y Calendario de Entregas (Reportes R07 a R12)

```
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                             CALENDARIO Y DOSIFICACIÓN DE REPORTES (UNIDAD 2)                                     │
├────────┬──────────────────────────────────────────┬────────────────────────────────────────┬─────────────────────┤
│Reporte │ Temas de la Sesión (Luenberger 4ta Ed.)  │ Textos Oficiales Alternativos          │ Entregable Requerido│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R07. (N)│ 2.1.1 Forma Estándar & Canónica          │ Bazaraa, Jarvis & Sherali (2011) Cap. 1│ Nota teórica +      │
│        │ 2.1.2 Ejemplos Clásicos de Modelado      │ Cottle & Thapa (2017) Cap. 3           │ Formulación de 3    │
│        │ 2.1.3 Soluciones Básicas y Degeneración  │                                        │ modelos canónicos   │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R08. (N)│ 2.1.4 Teorema Fundamental de la PL       │ Bazaraa, Jarvis & Sherali (2011) Cap. 2│ Demostración TFPL + │
│        │ 2.1.5 Equivalencia SBF ↔ Puntos Extremos │ Cottle & Thapa (2017) Cap. 4           │ Script Python de    │
│        │                                          │                                        │ cálculo de todas SBF│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R09. (N)│ 2.2.1-2.2.3 Pivotes, Adyacencia y Costos │ Bazaraa, Jarvis & Sherali (2011) Cap. 3│ Nota teórica +      │
│        │       Reducidos                          │ Cottle & Thapa (2017) Cap. 5           │ Ejercicio de        │
│        │ 2.2.4 Algoritmo Simplex Tabular & Bland  │                                        │ pivoteo a mano      │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R10. (N)│ 2.2.5 Método de Dos Fases y Gran M       │ Bazaraa, Jarvis & Sherali (2011) Cap. 4│ Script Simplex      │
│        │ 2.2.6-2.2.7 Simplex Matricial & Revisado │ Nocedal & Wright (2006) Cap. 13        │ Tabular en Python + │
│        │       con Factorización LU               │                                        │ Dos Fases           │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R11. (N)│ 2.3.1 Formulación de Programas Duales    │ Bazaraa, Jarvis & Sherali (2011) Cap. 6│ Demostración de     │
│        │ 2.3.2 Teoremas de Dualidad Débil y Fuerte│ Cottle & Thapa (2017) Cap. 6           │ Dualidad Fuerte +   │
│        │ 2.3.3 Relaciones con el Simplex          │                                        │ Lema de Farkas      │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R12. (N)│ 2.3.4 Sensibilidad y Holgura Complem.    │ Bazaraa, Jarvis & Sherali (2011) Cap. 6│ Script Python de    │
│        │ 2.3.5 Método Simplex Dual                │ Cottle & Thapa (2017) Cap. 6 & 7       │ Simplex Dual y      │
│        │ 2.3.6 Algoritmo Primal-Dual              │                                        │ Análisis Sensibilid.│
└────────┴──────────────────────────────────────────┴────────────────────────────────────────┴─────────────────────┘
```

---

## 3. Guía Específica por Reporte

### Reporte R07: Formulación Estándar, Modelado y Soluciones Básicas
* **Subtemas:** 2.1.1, 2.1.2 y 2.1.3.
* **Teorema/Proposición a Desarrollar:** Proposición sobre la equivalencia entre cualquier programa lineal general y su forma estándar $(\min \mathbf{c}^T\mathbf{x} \text{ s.a. } \mathbf{A}\mathbf{x}=\mathbf{b}, \mathbf{x}\ge\mathbf{0})$. Demostración de que el rango completo por filas de $\mathbf{A}$ ($m \le n$) es una condición sin pérdida de generalidad (eliminación de redundancias o detección de inconsistencias).
* **Ejercicio Numérico a Mano:** Dado un sistema con $m=2, n=4$:
  $$\mathbf{A} = \begin{pmatrix} 1 & 1 & 1 & 0 \\ 1 & 2 & 0 & 1 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 4 \\ 6 \end{pmatrix}$$
  Calcular exhaustivamente todas las $\binom{4}{2} = 6$ posibles soluciones básicas, clasificándolas en factibles (SBF), infactibles y degeneradas.
* **Pregunta Detonadora:** *¿Por qué una solución óptima en un problema no degenerado tiene exactamente $m$ componentes estrictamente positivas, y qué anomalías geométricas provoca la degeneración durante la optimización?*

---

### Reporte R08: El Teorema Fundamental de la PL y Geometría de Puntos Extremos
* **Subtemas:** 2.1.4 y 2.1.5.
* **Teorema a Demostrar:** Teorema Fundamental de la Programación Lineal (Luenberger pp. 21–23):
  1. Si existe una solución factible, existe una solución básica factible.
  2. Si existe una solución óptima factible, existe una solución básica factible óptima.
* **Demostración Geométrica:** Teorema de equivalencia entre vértices del poliedro $P = \{\mathbf{x}: \mathbf{A}\mathbf{x}=\mathbf{b}, \mathbf{x}\ge\mathbf{0}\}$ y Soluciones Básicas Factibles.
* **Implementación en Python:** Script `06_soluciones_basicas_factibles.py` que reciba $\mathbf{A}, \mathbf{b}, \mathbf{c}$, genere todas las combinaciones de columnas básicas $\mathbf{B}$, calcule $\mathbf{x}_B = \mathbf{B}^{-1}\mathbf{b}$, evalúe $\mathbf{c}^T\mathbf{x}$ y determine por fuerza bruta el óptimo global comprobando el Teorema Fundamental.
* **Pregunta Detonadora:** *Si un poliedro no tiene puntos extremos, ¿puede formularse como un problema de programación lineal en forma estándar $\mathbf{A}\mathbf{x}=\mathbf{b}, \mathbf{x}\ge\mathbf{0}$? Justifique con base en los conos de recesión.*

---

### Reporte R09: Álgebra del Simplex Tabular y Regla de Bland
* **Subtemas:** 2.2.1, 2.2.2, 2.2.3 y 2.2.4.
* **Teorema a Demostrar:** Teorema de Reducción Estricta de Costo: Si para la columna que entra $j$ el costo reducido es $r_j < 0$ y la solución no es degenerada ($\theta^* > 0$), entonces el nuevo costo satisface $z_{nuevo} = z_{ant} + \theta^* r_j < z_{ant}$.
* **Demostración de Parada:** Teorema de Optimalidad del Simplex y Teorema de Bland para prevención de ciclos infinitos.
* **Ejercicio Numérico a Mano:** Resolución completa paso a paso de un problema de minimización con 3 restricciones y 5 variables mediante el tablero Simplex, registrando explícitamente variables básicas, vector $\mathbf{c}_B$, columna pivote, fila pivote y cálculo de $\theta^*$.
* **Pregunta Detonadora:** *¿Bajo qué circunstancias específicas el simplex puede permanecer en el mismo vértice durante una o más iteraciones y cómo garantiza la regla de Bland la convergencia finita?*

---

### Reporte R10: Inicialización (Dos Fases y Gran M) y Simplex Revisado LU
* **Subtemas:** 2.2.5, 2.2.6 y 2.2.7.
* **Teorema a Demostrar:** Teorema de Factibilidad de Fase I: El problema original tiene una solución factible si y solo si el valor óptimo de la Fase I es $W^* = 0$. Procedimiento de expulsión de variables artificiales básicas con valor cero al finalizar la Fase I.
* **Implementación en Python:** Script `07_simplex_primal_tabular.py` y `08_dos_fases_y_gran_m.py` que resuelvan automáticamente problemas lineales desde la entrada estándar o archivos matriciales, imprimiendo el tablero en cada iteración.
* **Pregunta Detonadora:** *¿Cuáles son las desventajas de condicionamiento numérico del método de la Gran M frente al de Dos Fases en aritmética de punto flotante computacional?*

---

### Reporte R11: Teoría de Dualidad y Teorema Fundamental de Dualidad
* **Subtemas:** 2.3.1, 2.3.2 y 2.3.3.
* **Teoremas a Demostrar:**
  1. Lema de Dualidad Débil: $\mathbf{b}^T\mathbf{y} \le \mathbf{c}^T\mathbf{x}$.
  2. Teorema de Dualidad Fuerte: $\mathbf{b}^T\mathbf{y}^* = \mathbf{c}^T\mathbf{x}^*$ utilizando la base óptima primal $\mathbf{B}$ y definiendo $\mathbf{y}^{*T} = \mathbf{c}_B^T \mathbf{B}^{-1}$.
  3. Lema de Farkas como consecuencia de la dualidad.
* **Ejercicio Numérico a Mano:** Construir el dual exacto de un problema con restricciones mixtas ($\le, \ge, =$) y variables con signo restringido e irrestricto. Resolver el primal por el Simplex y verificar que los multiplicadores simplex del último tablero coinciden con la solución óptima del dual.
* **Pregunta Detonadora:** *¿Qué relación geométrica existe entre el Teorema del Hiperplano Separador (estudiado en la Unidad 1) y el Teorema de Dualidad Fuerte de la Programación Lineal?*

---

### Reporte R12: Holgura Complementaria, Sensibilidad y Simplex Dual
* **Subtemas:** 2.3.4, 2.3.5 y 2.3.6.
* **Teoremas a Demostrar:**
  1. Teorema de Holgura Complementaria: $(\mathbf{c} - \mathbf{A}^T\mathbf{y}^*)^T \mathbf{x}^* = 0$.
  2. Teorema de Sensibilidad y Precios Sombra: $\frac{\partial z^*}{\partial b_i} = y_i^*$.
* **Implementación en Python:** Script `10_dualidad_holgura_complementaria.py` y `11_simplex_dual.py` para reoptimizar un problema cuando se agrega una nueva restricción al problema óptimo original sin recalcular desde el inicio.
* **Pregunta Detonadora:** *En un problema de producción económica, si un recurso no se agota por completo en el óptimo ($s_i > 0$), ¿por qué su precio sombra dual correspondiente debe ser estrictamente cero? Explique desde la perspectiva de la teoría de precios de mercado.*
