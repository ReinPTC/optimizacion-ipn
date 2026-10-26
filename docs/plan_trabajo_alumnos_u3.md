# Plan de Trabajo y Dosificación de Notas para Alumnos — Unidad 3

**Asignatura:** Optimización  
**Unidad Temática III:** Optimización Sin Restricciones  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming* (4ta Ed., Springer 2016).  
**Institución:** Instituto Politécnico Nacional (IPN) — Escuela Superior de Física y Matemáticas (ESFM)  
**Licenciatura:** Matemática Algorítmica  

---

## 1. Lineamientos Pedagógicos y Reglas de Elaboración de Notas

### 🎯 Propósito Formativo
El objetivo de las notas individuales para la Unidad 3 es desarrollar en el estudiante de Matemática Algorítmica el dominio analítico del cálculo diferencial multivariable aplicado, la habilidad de demostrar rigurosamente tasas y órdenes de convergencia asintótica, la comprensión geométrica del acondicionamiento espectral de matrices de curvatura y la competencia para implementar esquemas algorítmicos no lineales vectorizados desde cero en Python (NumPy) sin recurrir a funciones prefabricadas de caja negra (e.g. `scipy.optimize.minimize`).

### 📋 Reglas Estrictas de Elaboración y Entrega
1. **Consulta Obligatoria de Textos Alternativos:**  
   Cada reporte debe contrastar las demostraciones de Luenberger & Ye con **al menos un texto alternativo** de la bibliografía oficial del IPN:
   - Bazaraa, Sherali & Shetty (2013), *Nonlinear Programming: Theory and Algorithms* (3ra Ed.).
   - Nocedal & Wright (2006), *Numerical Optimization* (2da Ed.).
   - Bonnans, Gilbert, Lemarechal & Sagastizábal (2006), *Numerical Optimization*.
   - Cottle & Thapa (2017), *Linear and Nonlinear Optimization*.
2. **Estructura Requerida de Cada Reporte (Formato Markdown / PDF / Notebook):**
   * **Encabezado Académico:** Nombre del alumno, número de boleta, número de reporte (R13 a R18), fecha y temas cubiertos.
   * **Marco Teórico y Definiciones Rigurosas:** Notación matemática formal con expresiones en LaTeX.
   * **Demostración Analítica Asignada:** Desarrollo paso a paso sin omitir pasos intermedios ni justificaciones de hipótesis (e.g. compacidad, diferenciabilidad $C^2/C^3$, definición positiva $\mathbf{Q} \succ 0$).
   * **Ejercicio Numérico a Mano:** Resolución completa con cálculos explícitos de gradientes, matrices Hessianas, descomposiciones espectrales o pasos de búsqueda de línea.
   * **Implementación Computacional en Python:** Código limpio, vectorizado en NumPy, modular, con *docstrings*, pruebas unitarias y gráficos comparativos de convergencia.
   * **Pregunta de Reflexión Teórica:** Análisis crítico argumentado sobre la pregunta detonadora del reporte.
   * **Referencias Bibliográficas:** Citas formales en formato APA/IEEE especificando capítulos, secciones y páginas.

---

## 2. Dosificación y Calendario de Entregas (Reportes R13 a R18)

```
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                             CALENDARIO Y DOSIFICACIÓN DE REPORTES (UNIDAD 3)                                     │
├────────┬──────────────────────────────────────────┬────────────────────────────────────────┬─────────────────────┤
│Reporte │ Temas de la Sesión (Luenberger 4ta Ed.)  │ Textos Oficiales Alternativos          │ Entregable Requerido│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R13. (N)│ 3.2 Condiciones de Optimalidad (Cap. 7)  │ Bazaraa, Sherali & Shetty (2013) Cap. 3│ Demostración FONC / │
│        │ • FONC en abiertos y conos factibles     │ Nocedal & Wright (2006) Cap. 2         │ SOSC + Análisis de  │
│        │ • SONC, SOSC y Criterio de Sylvester     │                                        │ signatura espectral │
│        │ • Convexidad y Puntos de Ensilladura     │                                        │ en Python (Mód. 12) │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R14. (N)│ 3.1 Búsqueda sobre una Línea (Cap. 8.1)  │ Bazaraa, Sherali & Shetty (2013) Cap. 8│ Comparativa de      │
│        │ • Unimodalidad, Fibonacci y Sección Áurea│ Cottle & Thapa (2017) Cap. 7           │ contracción Áurea vs│
│        │ • Ajustes Cuadrático y Cúbico (Davidon)  │                                        │ Armijo + Módulo 13  │
│        │ • Búsqueda Inexacta de Armijo y Wolfe    │                                        │ en Python           │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R15. (N)│ 3.3 Descenso Más Pronunciado (Cap. 8.2)  │ Nocedal & Wright (2006) Cap. 3         │ Demostración de     │
│        │ • Deducción paso de Rayleigh y ortogonal.│ Bonnans et al. (2006) Cap. 3           │ Kantórovich + Cota  │
│        │ • Fenómeno de zig-zag en elipsoides      │                                        │ ((κ-1)/(κ+1))² +    │
│        │ • Desigualdad de Kantórovich y tasa      │                                        │ Módulo 14 en Python │
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R16. (N)│ 3.3 Método de Newton Multidim. (Cap. 8.5)│ Nocedal & Wright (2006) Cap. 3 & 4     │ Demostración orden  │
│        │ • Aproximación de Taylor y paso Newton   │ Bazaraa, Sherali & Shetty (2013) Cap. 8│ cuadrático (p = 2) +│
│        │ • Teorema de convergencia cuadrática local│                                        │ Invarianza afín +   │
│        │ • Invarianza afín y patologías           │                                        │ Benchmark Rosenbrock│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R17. (N)│ 3.4 Métodos de Direcciones Conjugadas    │ Nocedal & Wright (2006) Cap. 5         │ Demostración del    │
│        │ • Q-ortogonalidad y Gram-Schmidt         │ Cottle & Thapa (2017) Cap. 8           │ Subespacio en       │
│        │ • Teorema del Subespacio en Expansión    │                                        │ Expansión + C-G     │
│        │ • Algoritmo C-G, Krylov y Chebyshev      │                                        │ en n pasos (Mód. 15)│
├────────┼──────────────────────────────────────────┼────────────────────────────────────────┼─────────────────────┤
│R18. (N)│ 3.5 Métodos Cuasi-Newton (Capítulo 10)   │ Nocedal & Wright (2006) Cap. 6 & 7     │ Demostración H_{k+1}│
│        │ • Ecuación de la secante y Rango 1/2     │ Bazaraa, Sherali & Shetty (2013) Cap.10│ def. positiva en    │
│        │ • Algoritmos DFP y BFGS                  │ Bonnans et al. (2006) Cap. 4           │ BFGS + Algoritmo    │
│        │ • L-BFGS sin memoria y Región Confianza  │                                        │ L-BFGS (Módulo 16)  │
└────────┴──────────────────────────────────────────┴────────────────────────────────────────┴─────────────────────┘
```

---

## 3. Catálogo de Problemas Analíticos Obligatorios (Problem Set U3)

### Bloque A: Optimalidad, Curvatura y Búsqueda 1D (Reportes R13 y R14)
1. **Problema 3.1 (Punto de Silla y Signatura):**  
   Dada la función $f(x_1, x_2) = x_1^3 - x_1^2 x_2 + 2x_2^2$ (Luenberger Ejemplo 7.3.2, p. 186):
   - Determine todos los puntos estacionarios en $\mathbb{R}^2$.
   - Calcule la matriz Hessiana $\mathbf{F}(\mathbf{x})$ en cada punto.
   - Demuestre que el punto $(6, 9)$ es un punto de ensilladura evaluando su determinante ($\det = -72$) y determinando explícitamente las direcciones de ascenso y descenso local.
2. **Problema 3.2 (Cerradura del Mapeo de Búsqueda Lineal):**  
   Sea $f(x) = (x - 1)^2$ en $\mathbb{R}$ y consideremos la sucesión de direcciones $d_k = 1/k \to 0$ con punto inicial $x_k = 0$ (Luenberger p. 228).
   - Demuestre que el mapeo de búsqueda exacta $\mathbf{S}(x_k, d_k)$ genera $y_k = 1$ para todo $k$, mientras que $\mathbf{S}(0, 0) = \{0\}$.
   - Concluya formalmente por qué la hipótesis $\mathbf{d} \neq \mathbf{0}$ es indispensable en el Teorema de Cerradura de Zangwill.
3. **Problema 3.3 (Minimax de Fibonacci vs. Sección Áurea):**  
   Deseamos reducir el intervalo de incertidumbre inicial $[0, 1]$ a una longitud no mayor a $\varepsilon = 10^{-3}$ en una función unimodal continua:
   - Calcule el número mínimo de evaluaciones de función requeridas mediante Búsqueda de Fibonacci ($1/F_{N+1} \le 10^{-3}$).
   - Compare dicho valor con el número de evaluaciones requeridas mediante Sección Áurea ($(0.618034)^{N-1} \le 10^{-3}$).

---

### Bloque B: Descenso Más Pronunciado y Desigualdad de Kantórovich (Reporte R15)
4. **Problema 3.4 (Ortogonalidad Sucesiva de Gradientes):**  
   Sea $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T \mathbf{Q} \mathbf{x} - \mathbf{b}^T \mathbf{x}$ con $\mathbf{Q} \succ 0$.
   - Demuestre analíticamente que bajo el paso de búsqueda exacta $\alpha_k = \frac{\mathbf{g}_k^T\mathbf{g}_k}{\mathbf{g}_k^T\mathbf{Q}\mathbf{g}_k}$, el nuevo gradiente $\mathbf{g}_{k+1} = \mathbf{g}_k - \alpha_k \mathbf{Q}\mathbf{g}_k$ satisface $\mathbf{g}_{k+1}^T \mathbf{g}_k = 0$.
   - Deduzca la relación para la norma del nuevo gradiente:
     $$\|\mathbf{g}_{k+1}\|^2 = \|\mathbf{g}_k\|^2 \left[ \frac{(\mathbf{g}_k^T\mathbf{Q}^2\mathbf{g}_k)(\mathbf{g}_k^T\mathbf{g}_k) - (\mathbf{g}_k^T\mathbf{Q}\mathbf{g}_k)^2}{(\mathbf{g}_k^T\mathbf{Q}\mathbf{g}_k)^2} \right]$$
5. **Problema 3.5 (Vector Extremal de Kantórovich):**  
   Dada una matriz $\mathbf{Q} \succ 0$ con autovalores extremos $a = \lambda_{\min}$ y $A = \lambda_{\max}$ asociados a los autovectores ortonormales $\mathbf{v}_1$ y $\mathbf{v}_n$:
   - Demuestre que el vector $\mathbf{x}_{\text{ext}} = \frac{1}{\sqrt{2}}(\mathbf{v}_1 + \mathbf{v}_n)$ alcanza la cota inferior de Kantórovich con igualdad exacta:
     $$\frac{(\mathbf{x}_{\text{ext}}^T \mathbf{x}_{\text{ext}})^2}{(\mathbf{x}_{\text{ext}}^T \mathbf{Q} \mathbf{x}_{\text{ext}})(\mathbf{x}_{\text{ext}}^T \mathbf{Q}^{-1} \mathbf{x}_{\text{ext}})} = \frac{4 a A}{(a + A)^2}$$
   - Grafique la tasa de convergencia asintótica $\beta(\kappa) = \left(\frac{\kappa-1}{\kappa+1}\right)^2$ en función del número de condición $\kappa \in [1, 1000]$.

---

### Bloque C: Método de Newton, Direcciones Conjugadas y Cuasi-Newton (Reportes R16, R17 y R18)
6. **Problema 3.6 (Invarianza Afín de Newton):**  
   Sea $f: \mathbb{R}^n \to \mathbb{R}$ diferenciable y sea $\mathbf{T} \in \mathbb{R}^{n \times n}$ una matriz no singular. Definamos $\mathbf{y} = \mathbf{T}\mathbf{x}$ y $\hat{f}(\mathbf{y}) = f(\mathbf{T}^{-1}\mathbf{y})$.
   - Demuestre que el paso de Newton en la variable $\mathbf{y}$, $\Delta \mathbf{y} = -[\nabla^2 \hat{f}(\mathbf{y})]^{-1}\nabla \hat{f}(\mathbf{y})$, satisface idénticamente $\Delta \mathbf{y} = \mathbf{T} \Delta \mathbf{x}$, demostrando que el método de Newton es invariante ante transformaciones afines del espacio.
7. **Problema 3.7 (Teorema del Subespacio en Expansión de C-G):**  
   En el método de Direcciones Conjugadas aplicado a $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T\mathbf{Q}\mathbf{x} - \mathbf{b}^T\mathbf{x}$:
   - Demuestre que para todo $k < n$, el vector $\mathbf{x}_k$ minimiza la función cuadrática sobre la variedad afín $\mathbf{x}_0 + \operatorname{span}\{\mathbf{d}_0, \mathbf{d}_1, \dots, \mathbf{d}_{k-1}\}$.
   - Concluya por qué el algoritmo encuentra la solución exacta $\mathbf{x}^*$ en a lo más $n$ iteraciones.
8. **Problema 3.8 (Preservación de Definición Positiva en BFGS):**  
   Sea $\mathbf{H}_k \succ 0$ y sea la actualización BFGS:
   $$\mathbf{H}_{k+1} = (\mathbf{I} - \rho_k \mathbf{s}_k \mathbf{y}_k^T) \mathbf{H}_k (\mathbf{I} - \rho_k \mathbf{y}_k \mathbf{s}_k^T) + \rho_k \mathbf{s}_k \mathbf{s}_k^T, \quad \text{donde } \rho_k = \frac{1}{\mathbf{y}_k^T \mathbf{s}_k}$$
   - Demuestre que si el paso cumple la condición de curvatura estricta $\mathbf{y}_k^T \mathbf{s}_k > 0$, entonces $\mathbf{H}_{k+1} \succ 0$ (es estrictamente definida positiva).

---

## 4. Preguntas Detonadoras de Examen Oral / Defensa de Notas

1. *¿Por qué el método de Descenso Más Pronunciado avanza en ángulos rectos exactos ($90^\circ$) bajo búsqueda de línea exacta, y de qué forma la Desigualdad de Kantórovich explica analíticamente el colapso de su velocidad en valles alargados?*
2. *¿Cuál es la diferencia fundamental entre el orden de convergencia ($p$) y la tasa asintótica de convergencia ($\beta$)? Proporcione ejemplos con Descenso Más Pronunciado, Secante y Newton.*
3. *¿Por qué el método de Newton puro puede divergir o converger hacia un máximo local si la matriz Hessiana no es definida positiva, y cómo resuelven esta patología la modificación de Levenberg-Marquardt y el método de Newton Amortiguado?*
4. *¿Por qué el método de Gradiente Conjugado es considerado un 'puente óptimo' entre primer y segundo orden? ¿Cuál es su costo computacional por iteración y su requerimiento de memoria respecto a Newton?*
5. *¿Qué ventajas computacionales ofrece la formulación L-BFGS (Memoryless / Two-Loop Recursion) sobre el método BFGS estándar para problemas con $n \ge 10^5$ variables?*
