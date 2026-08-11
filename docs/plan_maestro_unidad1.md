# Plan Maestro Curricular: Unidad 1 — Herramientas para la Optimización

**Unidad de Aprendizaje:** Optimización  
**Institución:** Instituto Politécnico Nacional (IPN)  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming*, 4th Edition, Springer (2016).  
**Unidad de Competencia:** Identifica los elementos necesarios para el estudio de los algoritmos de optimización a partir del álgebra lineal y los conjuntos convexos.

---

## 1. Mapeo Oficial: Temario del IPN $\longleftrightarrow$ Apéndices y Capítulos de Luenberger

```
1.1 Revisión de álgebra lineal
├── 1.1.1 Notación matricial (Luenberger Apéndice A, Secciones A.1 y A.2)
├── 1.1.2 Espacios (Luenberger Apéndice A, Sección A.3)
├── 1.1.3 Valores propios y formas cuadráticas (Luenberger Apéndice A, Sección A.4)
├── 1.1.4 Conceptos topológicos y funciones (Luenberger Apéndice A, Secciones A.5 y A.6)
└── 1.1.5 Eliminación gaussiana (Luenberger Apéndice C)

1.2 Conjuntos convexos
├── 1.2.1 Hiperplanos y politopos (Luenberger Apéndice B, Secciones B.1 y B.2)
├── 1.2.2 Hiperplanos separadores y de apoyo (Luenberger Apéndice B, Sección B.3)
└── 1.2.3 Puntos extremos (Luenberger Apéndice B, Sección B.4)

1.3 Elementos para la optimización
├── 1.3.1 Tipos de problemas y espacios de búsqueda (Luenberger Cap. 1, Secciones 1.1 y 1.2)
├── 1.3.2 Tamaño y dificultad de los problemas (Luenberger Cap. 1, Sección 1.3)
└── 1.3.3 Algoritmos iterativos y su convergencia (Luenberger Cap. 1, Sección 1.4)
```

---

## 2. Desglose Teórico y Bibliográfico por Subtema Oficial

| Subtema Oficial | Contenido Teórico y Teoremas | Referencia Exacta Luenberger (4ta Ed.) | Texto Alternativo Oficial (IPN) |
| :--- | :--- | :--- | :--- |
| **1.1.1 Notación matricial** | Conjuntos, $\min$, $\operatorname{argmin}$, supremo, ínfimo; matrices, operaciones, determinantes, matrices similares y particiones por bloques $[\mathbf{B}, \mathbf{C}]$. | Apéndice A, Sec. A.1–A.2 (pp. 495–497) | Bazaraa, Jarvis & Sherali (2011) Cap. 2; Cottle & Thapa (2017) Cap. 1 |
| **1.1.2 Espacios** | Espacios $E^n$, producto escalar $\mathbf{x}^T\mathbf{y}$, norma euclidiana, desigualdad de Cauchy-Schwarz, subespacios $M$, complemento ortogonal $M^\perp$, y Teorema de Descomposición Ortogonal Única ($\mathbf{x} = \mathbf{a} + \mathbf{b}$). | Apéndice A, Sec. A.3 (pp. 497–498) | Bazaraa, Jarvis & Sherali (2011) Cap. 2; Cottle & Thapa (2017) Cap. 1 |
| **1.1.3 Valores propios y formas cuadráticas** | Ecuación característica, propiedades espectrales de matrices simétricas, diagonalización $\mathbf{Q}^T\mathbf{A}\mathbf{Q}=\mathbf{\Lambda}$, signaturas (PD, PSD, ND, Indefinida) y raíz cuadrada $\mathbf{A}^{1/2}$. | Apéndice A, Sec. A.4 (pp. 498–499) | Nocedal & Wright (2006) Apéndice A; Cottle & Thapa (2017) Cap. 2 |
| **1.1.4 Conceptos topológicos y funciones** | Sucesiones $\{\mathbf{x}_k\} \to \mathbf{x}$, bolas abiertas, compacidad, Teorema de Weierstrass, gradiente $\nabla f(\mathbf{x})$, Hessiano $\nabla^2 f(\mathbf{x})$, Teorema de Taylor de orden 1 y 2, Teorema de la Función Implícita (TFI) y notación $o, O$. | Apéndice A, Sec. A.5–A.6 (pp. 499–503) | Bazaraa, Sherali & Shetty (2013) Cap. 1; Bonnans et al. (2006) Cap. 1 |
| **1.1.5 Eliminación gaussiana** | Sistemas triangulares, matrices elementales $\mathbf{M}_k$, factorización $\mathbf{A}=\mathbf{L}\mathbf{U}$, sustitución hacia adelante ($\mathbf{L}\mathbf{y}=\mathbf{b}$) y atrás ($\mathbf{U}\mathbf{x}=\mathbf{y}$), y pivoteo parcial ($|m_{ij}| \le 1$). | Apéndice C (pp. 513–515) | Cottle & Thapa (2017) Cap. 2; Nocedal & Wright (2006) Apéndice A |
| **1.2.1 Hiperplanos y politopos** | Conjuntos convexos, preservación de convexidad (Proposición 1), cápsula convexa $\operatorname{co}(S)$, conos, variedades lineales, equivalencia algebraica y geométrica de hiperplanos (Proposiciones 2 y 3: $\mathbf{a}^T\mathbf{x}=c$), semiespacios, politopos y poliedros. | Apéndice B, Sec. B.1–B.2 (pp. 505–509) | Bazaraa, Sherali & Shetty (2013) Cap. 2; Bazaraa, Jarvis & Sherali (2011) Cap. 2 |
| **1.2.2 Hiperplanos separadores y de apoyo** | **Teorema 1** (Separación estricta de punto exterior $\mathbf{y} \notin \bar{C}$ con $\delta = \inf \|\mathbf{x}-\mathbf{y}\|$ y $\mathbf{a} = \mathbf{x}_0 - \mathbf{y}$); **Teorema 2** (Hiperplano de apoyo en la frontera); **Teorema 3** (Separación de convexos disjuntos). | Apéndice B, Sec. B.3 (pp. 509–511) | Bazaraa, Sherali & Shetty (2013) Cap. 2; Cottle & Thapa (2017) Cap. 4 |
| **1.2.3 Puntos extremos** | Definición de punto extremo, **Lema 1** (puntos extremos de la intersección con hiperplano de apoyo), **Teorema 4** ($C = \operatorname{co}(\operatorname{ext}(C))$ para compactos convexos), y **Teorema 5** (Equivalencia de $\mathcal{H}$-rep y $\mathcal{V}$-rep de poliedros). | Apéndice B, Sec. B.4 (pp. 511–512) | Bazaraa, Jarvis & Sherali (2011) Cap. 2; Bazaraa, Sherali & Shetty (2013) Cap. 2 |
| **1.3.1 Tipos de problemas y espacios de búsqueda** | Filosofía de optimización, compensación precisión vs. tratabilidad computacional, formulación general de programación matemática con restricciones de igualdad y desigualdad, y espacios factibles continuos. | Cap. 1, Sec. 1.1–1.2 (pp. 1–5) | Nocedal & Wright (2006) Cap. 1; Bazaraa, Jarvis & Sherali (2011) Cap. 1 |
| **1.3.2 Tamaño y dificultad de los problemas** | Clasificación de escala (pequeña escala $\le 5$, mediana escala $5-1\,000$, gran escala $>1\,000$ a millones), esparcidad matricial ($\rho(\mathbf{A}) \ll 1\%$) y almacenamiento ralo. | Cap. 1, Sec. 1.3 (pp. 5–6) | Nocedal & Wright (2006) Cap. 1; Öchsner & Makvandi (2020) Cap. 1 |
| **1.3.3 Algoritmos iterativos y su convergencia** | Los 3 pilares teóricos (diseño, convergencia global y velocidad local), clasificación formal de órdenes de convergencia (lineal, superlineal, cuadrática), y Teorema de Tasa Canónica Lineal de Luenberger $\beta = \frac{\kappa(\mathbf{Q})-1}{\kappa(\mathbf{Q})+1}$. | Cap. 1, Sec. 1.4 (pp. 6–8) | Nocedal & Wright (2006) Cap. 2; Bazaraa, Jarvis & Sherali (2011) Cap. 3 |

---

## 3. Bibliografía Oficial del IPN

1. **Bazaraa, M. S., Jarvis, J. J. & Sherali, H. D. (2011)**. *Linear Programming and Network Flows* (4th ed.), Wiley.
2. **Bazaraa, M. S., Sherali, H. D. & Shetty, C. M. (2013)**. *Nonlinear Programming: Theory and Algorithms* (3rd ed.), Wiley.
3. **Bonnans, J. F., Gilbert, J. C., Lemarechal, C. & Sagastizábal, C. (2006)**. *Numerical Optimization*, Springer.
4. **Cottle, R. W. & Thapa, M. N. (2017)**. *Linear and Nonlinear Optimization*, Springer.
5. **Luenberger, D. G. & Ye, Y. (2016)**. *Linear and Nonlinear Programming* (4th ed.), Springer. *(Texto Base Obligatorio)*
6. **Nocedal, J. & Wright, S. (2006)**. *Numerical Optimization* (2nd ed.), Springer.
7. **Öchsner, A. & Makvandi, R. (2020)**. *Numerical Engineering Optimization*, Springer.
