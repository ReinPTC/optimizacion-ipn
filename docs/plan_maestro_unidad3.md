# Plan Maestro Curricular: Unidad 3 — Optimización Sin Restricciones

**Unidad de Aprendizaje:** Optimización  
**Institución:** Instituto Politécnico Nacional (IPN) — Escuela Superior de Física y Matemáticas (ESFM)  
**Licenciatura:** Matemática Algorítmica  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming*, 4ta Edición, Springer (2016).  
**Unidad de Competencia 3:** Resuelve problemas de optimización no lineal sin restricciones aplicando caracterizaciones analíticas de primer y segundo orden (FONC, SONC/SOSC), análisis espectral de curvatura, algoritmos de búsqueda lineal unidimensional, métodos de gradiente con análisis asintótico de convergencia (Kantórovich), métodos de direcciones conjugadas y esquemas cuasi-Newton con rigor matemático y computacional.

---

## 1. Arquitectura Pedagógica e Inversión Curricular Didáctica

En el programa oficial de la asignatura, el temario de la Unidad 3 se despliega tradicionalmente iniciando por *Búsqueda sobre una Línea* (Tema 3.1) y continuando con *Condiciones de Optimalidad* (Tema 3.2). Sin embargo, conforme a la arquitectura conceptual canónica de David G. Luenberger & Yinyu Ye:

$$\boxed{\text{Teoría Analítica de Soluciones (Cap. 7)} \longrightarrow \text{Búsqueda 1D (Cap. 8.1)} \longrightarrow \text{Descenso y Newton (Cap. 8.2, 8.5)} \longrightarrow \text{Direcciones Conjugadas (Cap. 9)} \longrightarrow \text{Cuasi-Newton (Cap. 10)}}$$

1. **Prioridad Teórica del Tema 3.2 (Cap. 7):** No es pedagógicamente viable exigir al estudiante comprender por qué una dirección es de descenso ($\nabla f(\mathbf{x})^T \mathbf{d} < 0$) o analizar la convergencia de una búsqueda lineal $\min_\alpha \phi(\alpha)$ sin haber formalizado previamente el gradiente como aproximación afín, las condiciones necesarias de primer orden (FONC), la signatura espectral Hessiana (SONC/SOSC) y la convexidad.
2. **Fundamentación Algorítmica del Tema 3.1 (Cap. 8.1):** Una vez dominadas las condiciones de optimalidad, la búsqueda sobre una línea se comprende con absoluta claridad como el operador de reducción unidimensional nuclear $\mathbf{S}(\mathbf{x}, \mathbf{d})$ necesario para implementar cualquier método iterativo multidimensional.
3. **Escalamiento Progresivo a Métodos Multidimensionales (Temas 3.3, 3.4 y 3.5):** Los métodos multidimensionales avanzan en orden creciente de sofisticación matemática:
   - **Tema 3.3 (Cap. 8.2 & 8.5):** Métodos clásicos de primer orden (*Steepest Descent* y cota de Kantórovich) frente al método exacto de segundo orden (*Newton Multidimensional* con orden cuadrático $p=2$).
   - **Tema 3.4 (Cap. 9):** Métodos de Direcciones Conjugadas y Gradiente Conjugado (C-G), puente óptimo que alcanza convergencia cuadrática en $n$ pasos utilizando únicamente información de primer orden ($\mathcal{O}(n)$ memoria).
   - **Tema 3.5 (Cap. 10):** Métodos Cuasi-Newton (DFP, BFGS y L-BFGS), que reconstruyen inductivamente la inversa Hessiana resolviendo las patologías de Newton.

---

## 2. Mapeo Oficial: Temario del IPN $\longleftrightarrow$ Capítulos de Luenberger & Ye (4ta Ed.)

```
Unidad 3: Optimización Sin Restricciones

3.2 Condiciones de Optimalidad (Luenberger & Ye, Capítulo 7)
├── 3.2.1 Condiciones Necesarias de Primer Orden (FONC) (Sec. 7.1, pp. 180–183)
├── 3.2.2 Ejemplos de Problemas sin Restricciones y Modelos Canónicos (Sec. 7.2, pp. 183–185)
├── 3.2.3 Condiciones de Segundo Orden (SONC / SOSC) y Signatura Hessiana (Sec. 7.3, pp. 185–189)
├── 3.2.4 Funciones Convexas y Cóncavas (Sec. 7.4, pp. 189–192)
├── 3.2.5 Minimización y Maximización de Funciones Convexas (Sec. 7.5, pp. 192–195)
└── 3.2.6 Condiciones de Orden Cero, Epígrafes e Hiperplanos de Soporte (Sec. 7.6, pp. 195–198)

3.1 Búsqueda sobre una Línea (Luenberger & Ye, Capítulo 8, Sección 8.1)
├── 3.1.1 Unimodalidad y Métodos de Orden Cero (Fibonacci y Sección Áurea) (Sec. 8.1, pp. 214–218)
├── 3.1.2 Métodos de Ajuste de Curvas (Lagrange Cuadrático y Davidon Cúbico) (Sec. 8.1, pp. 219–223)
├── 3.1.3 Métodos de Orden 1 y 2 (Secante, Falsa Posición y Newton 1D) (Sec. 8.1, pp. 223–226)
├── 3.1.4 Convergencia Global y Cerradura del Mapeo de Búsqueda S(x,d) (Sec. 8.1, pp. 227–229)
└── 3.1.5 Búsqueda Inexacta de Línea: Regla de Armijo y Condiciones de Wolfe (Sec. 8.1, pp. 230–232)

3.3 Métodos Básicos de Descenso: Descenso Más Pronunciado y Método de Newton (Luenberger & Ye, Cap. 8.2 y 8.5)
├── 3.3.1 El Método del Descenso Más Pronunciado (Steepest Descent) (Sec. 8.2, pp. 233–235)
├── 3.3.2 Ortogonalidad Sucesiva de Gradientes (g_{k+1}^T g_k = 0) y Zig-Zag (Sec. 8.2, pp. 234–235)
├── 3.3.3 Desigualdad de Kantórovich y Tasa de Convergencia Asintótica ((A-a)/(A+a))² (Sec. 8.2, pp. 235–239)
├── 3.3.4 El Método de Newton Multidimensional y Dirección de Newton (Sec. 8.5, pp. 247–249)
└── 3.3.5 Orden Cuadrático (p = 2), Invarianza Afín, Ventajas y Patologías (Sec. 8.5, pp. 248–251)

3.4 Métodos de Direcciones Conjugadas (Luenberger & Ye, Capítulo 9)
├── 3.4.1 Definición de Direcciones Q-Conjugadas y Gram-Schmidt Conjugado (Sec. 9.1, pp. 261–264)
├── 3.4.2 Propiedades de Descenso y Teorema del Subespacio en Expansión (Sec. 9.2, pp. 264–268)
├── 3.4.3 El Algoritmo de Gradiente Conjugado (C-G) para Formas Cuadráticas (Sec. 9.3, pp. 268–272)
├── 3.4.4 C-G como Proceso Óptimo, Subespacios de Krylov y Polinomios de Chebyshev (Sec. 9.4, pp. 272–276)
├── 3.4.5 El Método de Gradiente Conjugado Parcial y Reinicios Cíclicos (Sec. 9.5, pp. 276–279)
└── 3.4.6 Extensión a Problemas No Cuadráticos (Fletcher-Reeves y Polak-Ribière) (Sec. 9.6, pp. 279–284)

3.5 Métodos Cuasi-Newton (Luenberger & Ye, Capítulo 10)
├── 3.5.1 Métodos de Newton Modificados y Regularización Levenberg-Marquardt (Sec. 10.1, pp. 291–293)
├── 3.5.2 Construcción de la Inversa, Ecuación de la Secante y Corrección SR1 (Sec. 10.2, pp. 293–296)
├── 3.5.3 Algoritmos Davidon-Fletcher-Powell (DFP) y BFGS de Rango 2 (Sec. 10.3, pp. 296–301)
├── 3.5.4 Métodos Cuasi-Newton sin Memoria (L-BFGS y Two-Loop Recursion) (Sec. 10.7, pp. 312–315)
└── 3.5.5 Combinación de Descenso Más Pronunciado y Newton: Región de Confianza y Dogleg (Sec. 10.8, pp. 315–319)
```

---

## 3. Matriz de Desglose Teórico, Teoremas Clave y Bibliografía Oficial Comparativa

| Subtema | Contenido Teórico y Teoremas Medulares | Referencia Canónica Luenberger & Ye (4ta Ed.) | Textos Oficiales Alternativos (IPN) |
| :--- | :--- | :--- | :--- |
| **3.2 Condiciones de Optimalidad** | FONC en abiertos $\nabla f(\mathbf{x}^*)=\mathbf{0}$; FONC con restricciones de cono factible $\nabla f(\mathbf{x}^*)^T \mathbf{d} \ge 0$; SONC $\mathbf{F}(\mathbf{x}^*) \succeq 0$; SOSC $\mathbf{F}(\mathbf{x}^*) \succ 0$; Criterio de Sylvester; Desigualdad de Jensen; Teorema de equivalencia de mínimo local $\implies$ global en convexidad; Caracterizaciones de orden cero con epígrafes e hiperplanos de soporte. | Capítulo 7 (pp. 180–198) | Bazaraa, Sherali & Shetty (2013) Cap. 3 & 4; Nocedal & Wright (2006) Cap. 2; Bonnans et al. (2006) Cap. 1 |
| **3.1 Búsqueda sobre una Línea** | Unimodalidad; Búsqueda de Fibonacci (minimax óptimo $\frac{1}{F_{N+1}}$); Sección Áurea con contracción constante $\tau \approx 0.618034$; Ajuste cuadrático de 3 puntos (orden $p \approx 1.325$); Ajuste cúbico de Hermite (Davidon, $p=2$); Método de la Secante ($p \approx 1.618$); Newton 1D ($p=2$); Teorema de convergencia global por patrón de 3 puntos; Cerradura del mapeo de búsqueda $\mathbf{S}(\mathbf{x}, \mathbf{d})$ con $\mathbf{d} \neq \mathbf{0}$; Regla de Armijo con Backtracking y Condiciones de Wolfe. | Capítulo 8, Sec. 8.1 (pp. 214–232) | Bazaraa, Sherali & Shetty (2013) Cap. 8; Nocedal & Wright (2006) Cap. 3; Cottle & Thapa (2017) Cap. 7 |
| **3.3 Métodos Básicos de Descenso** | Derivación de la dirección de máximo descenso por Cauchy-Schwarz; Paso óptimo de Rayleigh $\alpha_k = \frac{\mathbf{g}_k^T\mathbf{g}_k}{\mathbf{g}_k^T\mathbf{Q}\mathbf{g}_k}$; Teorema de ortogonalidad sucesiva $\mathbf{g}_{k+1}^T\mathbf{g}_k = 0$; Fenómeno de zig-zag; Error de energía $E(\mathbf{x})$; **Teorema de la Desigualdad de Kantórovich** $\frac{(\mathbf{x}^T\mathbf{x})^2}{(\mathbf{x}^T\mathbf{Q}\mathbf{x})(\mathbf{x}^T\mathbf{Q}^{-1}\mathbf{x})} \ge \frac{4aA}{(a+A)^2}$; Cota de convergencia lineal $\left(\frac{\kappa-1}{\kappa+1}\right)^2$; Método de Newton Multidimensional $\mathbf{x}_{k+1} = \mathbf{x}_k - [\mathbf{F}(\mathbf{x}_k)]^{-1}\nabla f(\mathbf{x}_k)$; Teorema de convergencia cuadrática local ($p=2$); Invarianza afín de Newton; Patologías: pérdida de definición positiva y costo $\mathcal{O}(n^3)$. | Capítulo 8, Sec. 8.2 & 8.5 (pp. 233–251) | Nocedal & Wright (2006) Cap. 3 & 4; Bazaraa, Sherali & Shetty (2013) Cap. 8; Bonnans et al. (2006) Cap. 3 |
| **3.4 Métodos de Direcciones Conjugadas** | Definición geométrica de $Q$-ortogonalidad $\mathbf{d}_i^T\mathbf{Q}\mathbf{d}_j = 0$; Independencia lineal de direcciones conjugadas; Generación por Gram-Schmidt; **Teorema del Subespacio en Expansión** (minimización sobre la variedad afín $\mathbf{x}_0 + \mathcal{B}_k$); Convergencia en a lo más $n$ iteraciones; Algoritmo de Gradiente Conjugado (C-G) para formas cuadráticas con actualización recursiva $\beta_k = \frac{\|\mathbf{g}_{k+1}\|^2}{\|\mathbf{g}_k\|^2}$; Ortogonalidad de residuos $\mathbf{g}_i^T\mathbf{g}_j = 0$; Relación con subespacios de Krylov $\mathcal{K}_k(\mathbf{Q}, \mathbf{g}_0)$; Polinomios de Chebyshev y cota de convergencia $\|\mathbf{x}_k - \mathbf{x}^*\|_{\mathbf{Q}} \le 2\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^k \|\mathbf{x}_0 - \mathbf{x}^*\|_{\mathbf{Q}}$; C-G parcial y reinicios ante clusters de eigenvalores; Extensión no lineal de Fletcher-Reeves, Polak-Ribière y Hestenes-Stiefel con reinicios de Powell. | Capítulo 9, Sec. 9.1–9.6 (pp. 261–284) | Nocedal & Wright (2006) Cap. 5; Bazaraa, Sherali & Shetty (2013) Cap. 9; Cottle & Thapa (2017) Cap. 8 |
| **3.5 Métodos Cuasi-Newton** | Modificación de Cholesky / Levenberg-Marquardt $\mathbf{F}(\mathbf{x}_k) + \mu_k \mathbf{I}$; Ecuación de la secante $\mathbf{H}_{k+1}\mathbf{y}_k = \mathbf{s}_k$ (con $\mathbf{s}_k = \mathbf{x}_{k+1}-\mathbf{x}_k, \mathbf{y}_k = \mathbf{g}_{k+1}-\mathbf{g}_k$); Deducción y patologías de la actualización simétrica de rango 1 (SR1); Algoritmo Davidon-Fletcher-Powell (DFP) de rango 2; Demostración inductiva de preservación de simetría y definición positiva ($\mathbf{H}_{k+1} \succ 0$); Fórmula BFGS y su superioridad numérica dual; Métodos sin memoria (Memoryless Quasi-Newton / L-BFGS) con almacenamiento acotado de $m$ vectores y algoritmo *two-loop recursion*; Métodos de Región de Confianza, trayectoria continua de Levenberg-Marquardt y método Dogleg de Powell. | Capítulo 10, Sec. 10.1–10.8 (pp. 291–319) | Nocedal & Wright (2006) Cap. 6, 7 & 8; Bazaraa, Sherali & Shetty (2013) Cap. 10; Bonnans et al. (2006) Cap. 4 |

---

## 4. Suite Computacional en Python Asociada a la Unidad 3

Cada bloque teórico está respaldado por un módulo ejecutable en `python_src/` con tipado estricto, cálculos analíticos y benchmarks reproducibles:

1. [`12_condiciones_primer_segundo_orden.py`](file:///c:/Users/TCER/OneDrive/antigravity/RafaelIPN/Optimización/python_src/12_condiciones_primer_segundo_orden.py): Validación de FONC, Hessiana $\mathbf{F}(\mathbf{x}^*)$, signatura espectral (SOSC), análisis de puntos de ensilladura (Ejemplo Luenberger $(6, 9)$ con $\det = -72$) y mínimos cuadrados polinomiales.
2. [`13_busqueda_lineal_unidimensional.py`](file:///c:/Users/TCER/OneDrive/antigravity/RafaelIPN/Optimización/python_src/13_busqueda_lineal_unidimensional.py): Algoritmos de búsqueda 1D: Fibonacci, Sección Áurea, Ajuste Cuadrático de Lagrange, Secante, Newton 1D y regla de Armijo con Backtracking.
3. [`14_metodo_gradiente_descenso_kantorovich.py`](file:///c:/Users/TCER/OneDrive/antigravity/RafaelIPN/Optimización/python_src/14_metodo_gradiente_descenso_kantorovich.py): Descenso Más Pronunciado con ortogonalidad de pasos $\mathbf{g}_{k+1}^T\mathbf{g}_k=0$, verificación analítica de la Desigualdad de Kantórovich $\frac{4aA}{(a+A)^2}$, tasa asintótica $\left(\frac{\kappa-1}{\kappa+1}\right)^2$, Método de Newton Multidimensional con orden cuadrático ($p=2$) y benchmark comparativo en Rosenbrock.
4. [`15_gradiente_conjugado_fletcher_reeves.py`](file:///c:/Users/TCER/OneDrive/antigravity/RafaelIPN/Optimización/python_src/15_gradiente_conjugado_fletcher_reeves.py): Algoritmo C-G cuadrático canónico (terminación en $n$ pasos), evaluación de la cota de Chebyshev, subespacios de Krylov, y C-G no lineal (Fletcher-Reeves / Polak-Ribière) con reinicio de Powell.
5. [`16_metodos_cuasi_newton_bfgs.py`](file:///c:/Users/TCER/OneDrive/antigravity/RafaelIPN/Optimización/python_src/16_metodos_cuasi_newton_bfgs.py): Algoritmos DFP y BFGS con preservación de $\mathbf{H}_k \succ 0$, implementación L-BFGS (*two-loop recursion*), Newton modificado Levenberg-Marquardt y aproximación Dogleg de Powell.

---

## 5. Rúbrica de Evaluación de Competencias (Unidad 3)

| Dimensión de Aprendizaje | Nivel Inicial / Receptivo (6.0 – 7.0) | Nivel Autónomo / Competente (8.0 – 9.0) | Nivel Experto / Matemática Algorítmica (10.0) |
| :--- | :--- | :--- | :--- |
| **Demostraciones Analíticas** | Aplica fórmulas de memoria sin justificar hipótesis (e.g. omite comprobar si $\mathbf{Q} \succ 0$ o si el gradiente es no nulo). | Demuestra rigurosamente teoremas fundamentales (FONC, SOSC, Kantórovich, Ortogonalidad $\mathbf{g}_{k+1}^T\mathbf{g}_k=0$) justificando cada paso algebraico. | Demuestra teoremas complejos con soltura (Subespacio en expansión de C-G, cota de Chebyshev, preservación de definición positiva en BFGS e invarianza afín de Newton). |
| **Análisis de Convergencia** | Confunde orden de convergencia $p$ con tasa asintótica $\beta$. No relaciona $\kappa$ con la velocidad del algoritmo. | Calcula empírica y teóricamente tasas de convergencia lineales y cuadráticas; acota el error mediante la cota de Kantórovich. | Realiza análisis espectral fino; identifica el vector extremal que alcanza la cota inferior exacta de Kantórovich y explica la aceleración por clusters de autovalores. |
| **Diseño e Implementación de Código** | Depende de cajas negras (`scipy.optimize`); código monolítico y poco documentado. | Implementa algoritmos desde cero vectorizados en NumPy con control de paso (Armijo/Wolfe) y pruebas de convergencia. | Desarrolla arquitecturas modulares y eficientes: *two-loop recursion* en L-BFGS, reinicios inteligentes y benchmarks comparativos exhaustivos con profiling temporal. |
