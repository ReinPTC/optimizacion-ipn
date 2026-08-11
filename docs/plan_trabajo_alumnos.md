# Plan de Trabajo y Dosificación de Notas para Alumnos — Unidad 1

**Asignatura:** Optimización  
**Unidad Temática I:** Herramientas para la optimización  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming* (4ta Ed., Springer 2016).  
**Institución:** Instituto Politécnico Nacional (IPN)

---

## 1. Lineamientos Generales para la Elaboración de Notas

### 🎯 Propósito Pedagógico
El objetivo de las notas individuales es fomentar en el estudiante el rigor matemático, la capacidad de síntesis, la contrastación bibliográfica crítica y la competencia computacional para implementar los conceptos de optimización en código ejecutable de Python.

### 📋 Reglas Estrictas de Entrega
1. **Consulta de Textos Alternativos Obligatoria:**  
   La nota debe redactarse consultando **al menos un texto diferente a Luenberger** seleccionado exclusivamente de la **Bibliografía Oficial del Programa del IPN**. Se debe contrastar la notación, definiciones y enfoque pedagógico del autor alternativo frente a Luenberger.
2. **Estructura Requerida de Cada Reporte (Formato Markdown / PDF / Notebook):**
   * **Encabezado Académico:** Nombre del alumno, número de reporte, fecha y temas cubiertos.
   * **Marco Teórico y Definiciones Formales:** Definiciones precisas con fórmulas renderizadas en LaTeX.
   * **Demostración Analítica Rigurosa:** Desarrollo paso a paso de al menos una proposición/teorema clave (sin omitir justificaciones algebraicas).
   * **Ejemplo Práctico Resuelto Paso a Paso:** Resolución numérica completa a mano.
   * **Implementación Computacional en Python:** Código limpio, documentado con *docstrings*, gráficos y visualizaciones pertinentes.
   * **Pregunta de Reflexión Teórica:** Respuesta argumentada a la pregunta detonadora del reporte.
   * **Bibliografía:** Citas en formato IEEE o APA especificando capítulos y páginas consultadas tanto de Luenberger como del texto alternativo oficial.

---

## 2. Dosificación y Calendario de Entregas (Reportes R01 a R06)

```
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                             CALENDARIO Y DOSIFICACIÓN DE REPORTES (UNIDAD 1)                                     │
├────────┬──────────────────────────────────────────┬────────────────────────────────────────┬─────────────────────┤
│Reporte │ Temas de la Sesión (Luenberger 4ta Ed.)  │ Textos Oficiales Alternativos          │ Entregable Requerido│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R01. (N)│ 1.1.1 Notación Matricial & Particiones   │ Bazaraa, Jarvis & Sherali (2011) Cap. 2│ Nota teórica +      │
│        │ 1.1.2 Espacios En, Producto Escalar,     │ Cottle & Thapa (2017) Cap. 1           │ 3 ejercicios de     │
│        │       Cauchy-Schwarz y Proy. Ortogonal   │                                        │ proyección ortogonal│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R02. (N)│ 1.1.3 Autovalores & Formas Cuadráticas   │ Nocedal & Wright (2006) Apéndice A     │ Nota teórica +      │
│        │ 1.1.5 Eliminación Gaussiana, Fact. LU,   │ Cottle & Thapa (2017) Cap. 2           │ Script Python de    │
│        │       Matrices Elementales y Pivoteo     │                                        │ Factorización LU    │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R03. (N)│ 1.1.4 Topología en En (Abiertos/Cerrados,│ Bazaraa, Sherali & Shetty (2013) Cap. 1│ Nota teórica +      │
│        │       Weierstrass), Gradiente, Hessiano, │ Bonnans et al. (2006) Cap. 1           │ Demostración TFI y  │
│        │       Taylor y Teorema Función Implícita │                                        │ Desarrollo Taylor   │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R04. (N)│ 1.2.1 Variedades Lineales & Hiperplanos  │ Bazaraa, Sherali & Shetty (2013) Cap. 2│ Nota teórica +      │
│        │ 1.2.2 Hiperplanos Separadores y de Apoyo │ Bazaraa, Jarvis & Sherali (2011) Cap. 2│ Visualización de    │
│        │       (Teoremas 1, 2 y 3 de Separación)  │ Cottle & Thapa (2017) Cap. 4           │ Separación Convexa  │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R05. (N)│ 1.2.3 Puntos Extremos & Teorema de       │ Bazaraa, Jarvis & Sherali (2011) Cap. 2│ Nota teórica +      │
│        │       Representación de Politopos        │ Bazaraa, Sherali & Shetty (2013) Cap. 2│ Tabla y cálculo     │
│        │ 1.3.5 Soluciones Básicas Factibles (SBF) │ Cottle & Thapa (2017) Cap. 3           │ de todas las SBF    │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R06. (N)│ 1.3.1-1.3.3 Filosofía, Escala, Dificultad│ Nocedal & Wright (2006) Cap. 1 & 2     │ Nota teórica +      │
│        │       y Convergencia de Algoritmos       │ Öchsner & Makvandi (2020) Cap. 1 & 2   │ Formulación a Forma │
│        │ 1.3.4 Los 7 Modelos Clásicos de PL       │ Bazaraa, Jarvis & Sherali (2011) Cap. 1│ Estándar en Python  │
└────────┴──────────────────────────────────────────┴────────────────────────────────────────┴─────────────────────┘
```

---

## 3. Guía Detallada por Reporte

### 📄 Reporte R01. (N): Álgebra Lineal, Espacios $E^n$ y Proyecciones Ortogonales
* **Contenido Temático:**
  * Notación matricial de Luenberger, partición por bloques $\mathbf{A} = [\mathbf{B}, \mathbf{C}]$.
  * Espacio euclidiano $E^n$, producto interior $\mathbf{x}^T\mathbf{y}$, norma euclidiana $\|\mathbf{x}\| = (\mathbf{x}^T\mathbf{x})^{1/2}$, desigualdad de Cauchy-Schwarz.
  * Subespacios $M$, complemento ortogonal $M^\perp$ y Teorema de Descomposición Ortogonal única $\mathbf{x} = \mathbf{a} + \mathbf{b}$ con $\mathbf{a} \in M$, $\mathbf{b} \in M^\perp$.
* **Textos de Consulta:**
  * *Base:* Luenberger & Ye (2016), Apéndice A (Sec. A.1–A.3, pp. 495–498).
  * *Alternativo Oficial:* Bazaraa, Jarvis & Sherali (2011), Cap. 2; Cottle & Thapa (2017), Cap. 1.
* **Entregable Específico:**
  1. Resumen teórico riguroso en prosa académica.
  2. Demostración analítica de la desigualdad de Cauchy-Schwarz a partir de $\|\mathbf{x} - \alpha \mathbf{y}\|^2 \ge 0$.
  3. Cálculo manual de la proyección ortogonal del vector $\mathbf{x} = (3, 1, 4)^T$ sobre el subespacio $M = \operatorname{span}\{(1, 1, 0)^T, (0, 1, 1)^T\}$.
  4. Script en Python que calcule y grafique la proyección en 3D con Matplotlib.
* **Pregunta de Reflexión:** ¿Por qué la descomposición ortogonal $\mathbf{x} = \mathbf{a} + \mathbf{b}$ garantiza que $\mathbf{a}$ es el punto de $M$ más cercano a $\mathbf{x}$?

---

### 📄 Reporte R02. (N): Autovalores, Formas Cuadráticas y Factorización LU
* **Contenido Temático:**
  * Ecuación característica $\det(\mathbf{A}-\lambda\mathbf{I})=0$, ortogonalidad de autovectores de matrices simétricas, diagonalización ortogonal $\mathbf{Q}^T\mathbf{A}\mathbf{Q} = \mathbf{\Lambda}$.
  * Formas cuadráticas $\mathbf{x}^T\mathbf{A}\mathbf{x}$ y clasificación por signatura de autovalores (definida positiva, semidefinida, indefinida).
  * Raíz cuadrada simétrica $\mathbf{A}^{1/2} = \mathbf{Q}\mathbf{\Lambda}^{1/2}\mathbf{Q}^T$.
  * Eliminación gaussiana, matrices elementales $\mathbf{M}_k$, factorización $\mathbf{A} = \mathbf{L}\mathbf{U}$, sustitución hacia adelante/atrás y estrategia de pivoteo parcial.
* **Textos de Consulta:**
  * *Base:* Luenberger & Ye (2016), Apéndice A (Sec. A.4, pp. 498–499) y Apéndice C (pp. 513–515).
  * *Alternativo Oficial:* Nocedal & Wright (2006), Apéndice A; Cottle & Thapa (2017), Cap. 2.
* **Entregable Específico:**
  1. Desarrollo analítico de la signatura de la matriz $\mathbf{A} = \begin{pmatrix} 4 & 2 \\ 2 & 3 \end{pmatrix}$ y cálculo de $\mathbf{A}^{1/2}$.
  2. Factorización $\mathbf{A} = \mathbf{L}\mathbf{U}$ paso a paso para un sistema $3 \times 3$ mostrando explícitamente las matrices elementales $\mathbf{M}_1, \mathbf{M}_2$ y $\mathbf{L} = \mathbf{M}_1^{-1}\mathbf{M}_2^{-1}$.
  3. Script en Python con función propia `lu_factorization(A)` que implemente el algoritmo de Luenberger y compare el tiempo de cómputo frente a la inversión matricial explícita $\mathbf{A}^{-1}\mathbf{b}$.
* **Pregunta de Reflexión:** ¿Por qué en optimización numérica se prefiere resolver dos sistemas triangulares $\mathbf{L}\mathbf{y}=\mathbf{b}$ y $\mathbf{U}\mathbf{x}=\mathbf{y}$ en lugar de calcular explícitamente la inversa $\mathbf{A}^{-1}$?

---

### 📄 Reporte R03. (N): Topología en $E^n$, Weierstrass, Derivadas, Taylor y Función Implícita
* **Contenido Temático:**
  * Sucesiones en $E^n$, bolas abiertas, conjuntos abiertos/cerrados, interior $\mathring{S}$, frontera y compacidad.
  * Teorema de Weierstrass para funciones continuas en conjuntos compactos.
  * Gradiente $\nabla f(\mathbf{x})$, matriz Hessiana $\nabla^2 f(\mathbf{x}) = \mathbf{F}(\mathbf{x})$ y simetría de derivadas cruzadas.
  * Teorema de Taylor de orden 1 y 2 con término de residuo matricial.
  * Teorema de la Función Implícita (TFI) y su papel como generalización no lineal de sistemas lineales $\mathbf{A}\mathbf{x}=\mathbf{b}$.
  * Notación de Landau $O(x)$ y $o(x)$.
* **Textos de Consulta:**
  * *Base:* Luenberger & Ye (2016), Apéndice A (Sec. A.5–A.6, pp. 499–503).
  * *Alternativo Oficial:* Bazaraa, Sherali & Shetty (2013), Cap. 1; Bonnans et al. (2006), Cap. 1.
* **Entregable Específico:**
  1. Explicación de la importancia del Teorema de Weierstrass en la existencia de soluciones óptimas en optimización no lineal.
  2. Expansión de Taylor cuadrática de la función $f(x_1, x_2) = x_1^3 + 2x_1 x_2^2 - 3x_2^2$ alrededor del punto $(1, 2)$.
  3. Aplicación del Teorema de la Función Implícita para determinar si el sistema $x_1^2 + x_2^2 + x_3^2 - 1 = 0$, $x_1 + x_2 + x_3 = 0$ puede resolverse localmente para $(x_1, x_2)$ en términos de $x_3$ cerca de $(1/\sqrt{2}, -1/\sqrt{2}, 0)$.
  4. Script en Python con SymPy para verificar analíticamente el gradiente, Hessiano y aproximación de Taylor.
* **Pregunta de Reflexión:** ¿Qué ocurre con la aproximación local de un problema de optimización cuando la matriz Hessiana evaluada en un punto estacionario es indefinida o singular?

---

### 📄 Reporte R04. (N): Conjuntos Convexos, Hiperplanos y Teoremas de Separación
* **Contenido Temático:**
  * Definición de conjunto convexo, preservación por combinaciones, suma de Minkowski e intersecciones (Proposición 1).
  * Conos, conos convexos y cápsula convexa $\operatorname{co}(S)$.
  * Variedades lineales, equivalencia algebraica y geométrica de hiperplanos $H = \{\mathbf{x}: \mathbf{a}^T\mathbf{x} = c\}$ (Proposiciones 2 y 3).
  * Semiespacios $H^+, H^-$, politopos y poliedros.
  * Demostración completa del **Teorema 1 de Luenberger** (Separación estricta de punto exterior usando $\delta = \inf \|\mathbf{x}-\mathbf{y}\|$ y $\mathbf{a} = \mathbf{x}_0 - \mathbf{y}$).
  * Teorema 2 (Hiperplano de apoyo en frontera) y Teorema 3 (Separación de conjuntos convexos disjuntos).
* **Textos de Consulta:**
  * *Base:* Luenberger & Ye (2016), Apéndice B (Sec. B.1–B.3, pp. 505–511).
  * *Alternativo Oficial:* Bazaraa, Sherali & Shetty (2013), Cap. 2 (Sec. 2.4); Bazaraa, Jarvis & Sherali (2011), Cap. 2; Cottle & Thapa (2017), Cap. 4.
* **Entregable Específico:**
  1. Transcripción y desglose paso a paso de la demostración del Teorema 1 de Luenberger, justificando por qué el límite $\alpha \to 0^+$ produce $(\mathbf{x}_0 - \mathbf{y})^T(\mathbf{x} - \mathbf{x}_0) \ge 0$.
  2. Construcción analítica del hiperplano de separación entre el elipsoide $C = \{(x_1, x_2) : 2x_1^2 + x_2^2 \le 4\}$ y el punto exterior $\mathbf{y} = (3, 2)$.
  3. Script interactivo en Python con Matplotlib/Plotly que dibuje el conjunto convexo $C$, el punto $\mathbf{y}$, la proyección $\mathbf{x}_0$, el vector normal $\mathbf{a}$ y el hiperplano separador.
* **Pregunta de Reflexión:** ¿Por qué la convexidad del conjunto $C$ es indispensable para garantizar que el hiperplano separa a todo el conjunto y no solo localmente cerca de $\mathbf{x}_0$?

---

### 📄 Reporte R05. (N): Puntos Extremos, Teorema de Representación y Soluciones Básicas
* **Contenido Temático:**
  * Definición geométrica de punto extremo.
  * Lema 1 de Luenberger (puntos extremos de la intersección con hiperplano de apoyo).
  * Demostración por inducción en la dimensión $n$ del **Teorema 4** ($C = \operatorname{co}(\operatorname{ext}(C))$ para conjuntos convexos cerrados y acotados).
  * Teorema 5: Equivalencia entre la $\mathcal{H}$-representación (intersección finita de semiespacios) y la $\mathcal{V}$-representación (cápsula convexa de vértices) de poliedros.
  * Sistemas $\mathbf{A}\mathbf{x} = \mathbf{b}$, submatrices base $\mathbf{B}_{m \times m}$, variables básicas y no básicas.
  * Definición algebraica de Solución Básica y Solución Básica Factible (SBF), y su correspondencia 1-a-1 con los puntos extremos del poliedro factible.
* **Textos de Consulta:**
  * *Base:* Luenberger & Ye (2016), Apéndice B (Sec. B.4, pp. 511–512) y Cap. 2 (Sec. 2.3, pp. 19–20).
  * *Alternativo Oficial:* Bazaraa, Jarvis & Sherali (2011), Caps. 2 y 3; Bazaraa, Sherali & Shetty (2013), Cap. 2; Cottle & Thapa (2017), Cap. 3.
* **Entregable Específico:**
  1. Demostración analítica del Lema 1 de Luenberger.
  2. Para el politopo definido por $x_1 + 2x_2 \le 6$, $2x_1 + x_2 \le 6$, $x_1, x_2 \ge 0$:
     - Convertir a forma estándar con variables de holgura $\mathbf{A}\mathbf{x} = \mathbf{b}$.
     - Enumerar todas las $\binom{4}{2} = 6$ posibles bases $\mathbf{B}$, calcular sus soluciones básicas $\mathbf{x}_B = \mathbf{B}^{-1}\mathbf{b}$ y determinar cuáles son factibles (SBF).
     - Identificar las coordenadas geométricas $(x_1, x_2)$ de cada SBF y comprobar que corresponden exactamente a los vértices del poliedro.
  3. Script en Python que automatice el cálculo de todas las combinaciones de columnas de $\mathbf{A}$, determine la invertibilidad de $\mathbf{B}$, obtenga las SBF y grafique el poliedro resaltando los puntos extremos.
* **Pregunta de Reflexión:** ¿Qué ventaja computacional proporciona saber que el óptimo de un programa lineal siempre se encuentra en una Solución Básica Factible (punto extremo)?

---

### 📄 Reporte R06. (N): Modelado en PL, Escala, Convergencia y los 7 Problemas Canónicos
* **Contenido Temático:**
  * Filosofía de optimización de Luenberger, balance precisión vs. tratabilidad computacional.
  * Clasificación por escala (pequeña, mediana, gran escala) y la importancia crítica de la esparcidad.
  * Algoritmos iterativos: convergencia global vs. local, tasa canónica de convergencia (*"One good theory is worth a thousand computer runs"*).
  * Transformación a forma estándar (holguras, exceso, variables libres por descomposición $x = u - v$ o sustitución).
  * Los **7 modelos canónicos de Luenberger**:
    1. Problema de la Dieta (Stigler).
    2. Problema de Manufactura.
    3. Problema de Transporte (matriz de coeficientes rala).
    4. Problema de Flujo Máximo en redes capacitadas.
    5. Problema de Almacenamiento multitemporal (Warehousing).
    6. Clasificador Lineal y Support Vector Machine (SVM).
    7. Subasta Combinatoria Parimutuel (formulación minimax linealizada).
* **Textos de Consulta:**
  * *Base:* Luenberger & Ye (2016), Cap. 1 (Sec. 1.1–1.4, pp. 1–8) y Cap. 2 (Sec. 2.1–2.2, pp. 11–19).
  * *Alternativo Oficial:* Nocedal & Wright (2006), Caps. 1 y 2; Bazaraa, Jarvis & Sherali (2011), Cap. 1; Öchsner & Makvandi (2020), Caps. 1 y 2.
* **Entregable Específico:**
  1. Ensayo crítico (1–2 páginas) sobre la clasificación de escala y esparcidad de Luenberger y cómo los métodos de búsqueda iterativa superan la resolución analítica de condiciones de Lagrange.
  2. Formulación matemática formal en forma estándar matricial $(\mathbf{c}, \mathbf{A}, \mathbf{b})$ de al menos 3 de los 7 modelos canónicos de Luenberger (incluyendo obligatoriamente el Clasificador SVM y el Problema de Almacenamiento).
  3. Script en Python utilizando `scipy.optimize.linprog` que resuelva una instancia numérica del Problema de Almacenamiento ($n=4$ periodos) y del Clasificador SVM con puntos 2D separables, graficando el hiperplano resultante.
* **Pregunta de Reflexión:** ¿Cómo se linealiza la función objetivo minimax $\max [\boldsymbol{\pi}^T\mathbf{x} - \max_i (\mathbf{A}\mathbf{x})_i]$ en el problema de la subasta combinatoria para transformarlo en un programa lineal estándar?

---

## 4. Rúbrica de Evaluación Oficial (Escala de 100 Puntos)

| Criterio | Descripción | Puntaje Máximo |
| :--- | :--- | :---: |
| **Rigor Matemático y Demostraciones** | Precisión en las definiciones formales, notación vectorial/matricial correcta y desarrollo paso a paso sin saltos lógicos en las demostraciones analíticas. | **25 pts** |
| **Fidelidad y Consulta Bibliográfica** | Apego conceptual a Luenberger & Ye (2016) y contrastación explícita y documentada con al menos un **texto alternativo oficial del IPN**, incluyendo citas formales (IEEE/APA). | **20 pts** |
| **Ejercicios Numéricos y Procedimientos** | Resolución detallada, paso a paso, con operaciones matriciales y vectoriales completas y verificación de resultados. | **20 pts** |
| **Implementación Computacional en Python** | Código funcional, modular, con *docstrings*, buenas prácticas PEP 8, reproducible y con visualizaciones gráficas claras (Matplotlib/Plotly). | **20 pts** |
| **Análisis Crítico y Pregunta Detonadora** | Profundidad conceptual, claridad argumentativa y justificación sólida de la pregunta de reflexión teórica. | **15 pts** |
| **Total** | | **100 pts** |

---

## 5. Bibliografía Oficial del IPN para Consulta de Alumnos

1. **Bazaraa, M. S., Jarvis, J. J. & Sherali, H. D. (2011)**. *Linear Programming and Network Flows* (4th ed.). Hoboken, NJ: John Wiley & Sons. ISBN: 978-1-118-21132-8. [Tipo B]
2. **Bazaraa, M. S., Sherali, H. D. & Shetty, C. M. (2013)**. *Nonlinear Programming: Theory and Algorithms* (3rd ed.). Hoboken, NJ: John Wiley & Sons. ISBN: 978-1-118-85756-4. [Tipo B]
3. **Bonnans, J. F., Gilbert, J. C., Lemarechal, C. & Sagastizábal, C. (2006)**. *Numerical Optimization: Theoretical and Practical Aspects* (2nd ed.). Berlin: Springer-Verlag. ISBN: 978-3-540-35445-1. [Tipo C]
4. **Cottle, R. W. & Thapa, M. N. (2017)**. *Linear and Nonlinear Optimization*. New York: Springer. ISBN: 978-1-4939-7053-7. [Tipo B]
5. **Luenberger, D. G. & Ye, Y. (2016)**. *Linear and Nonlinear Programming* (4th ed.). Cham, Switzerland: Springer International Publishing. ISBN: 978-3-319-18841-6. [Tipo B — **Texto Base Obligatorio**]
6. **Nocedal, J. & Wright, S. (2006)**. *Numerical Optimization* (2nd ed.). New York: Springer. ISBN: 978-0-387-30303-1. [Tipo B]
7. **Öchsner, A. & Makvandi, R. (2020)**. *Numerical Engineering Optimization: Continuous and Discrete Problems*. Cham, Switzerland: Springer. ISBN: 978-3-030-43387-1. [Tipo C]
