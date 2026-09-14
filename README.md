# Libro Digital Interactivo: Optimización — Unidades 1 y 2

**Institución:** Instituto Politécnico Nacional (IPN)  
**Unidad de Aprendizaje:** Optimización  
**Texto Base Obligatorio:** David G. Luenberger & Yinyu Ye, *Linear and Nonlinear Programming* (4ta Edición, Springer 2016).

---

## 🎯 Descripción del Proyecto
Este repositorio contiene la plataforma digital interactiva (Progressive Web App - PWA) y la suite computacional de notas para:
- **Unidad 1: Herramientas para la optimización**
- **Unidad 2: Programación Lineal (Algoritmo Simplex, Dos Fases, Dualidad)**

Desarrollada bajo rigurosos estándares académicos y pedagógicos del IPN.

El contenido presenta:
- Explicaciones teóricas rigurosas y demostraciones exclusivas de Luenberger & Ye.
- Renderizado matemático vectorial en LaTeX (KaTeX).
- Herramientas y simuladores matemáticos interactivos ejecutables directamente en el navegador (Formas cuadráticas, Factorización LU, Simplex interactivo, Visualizador de politopos 2D).
- Módulos en Python 3 y Jupyter Notebooks para la experimentación computacional.
- Planes maestros y planes de trabajo oficiales para estudiantes con dosificación de reportes.

---

## 📁 Estructura del Repositorio

```
Optimización/
│
├── index.html                   # Aplicación Web Interactiva Principal (PWA)
├── manifest.json                # Web App Manifest (PWA instalable)
├── sw.js                        # Service Worker para funcionamiento 100% offline
├── requirements.txt             # Dependencias de Python (NumPy, SciPy, Matplotlib, Jupyter)
│
├── css/
│   └── main.css                 # Sistema de diseño, temas dark/light, tablas simplex, cajas KaTeX
│
├── js/
│   ├── app.js                   # Lógica de navegación, búsqueda en tiempo real y temas
│   └── interactive_tools.js     # Simuladores Canvas (Formas cuadráticas, LU Solver, Simplex Tabular)
│
├── docs/
│   ├── plan_maestro_unidad1.md  # Plan curricular U1 y matriz de dependencias
│   ├── plan_trabajo_alumnos.md  # Guía de entrega de reportes U1
│   ├── plan_maestro_unidad2.md  # Plan curricular U2 (Programación Lineal)
│   └── plan_trabajo_alumnos_u2.md # Guía de entrega de reportes U2
│
├── python_src/                  # Módulos computacionales en Python
│   ├── 01_algebra_lineal_formas_cuadraticas.py
│   ├── 02_factorizacion_lu_gauss.py
│   ├── 03_geometria_convexa_separacion.py
│   ├── 04_modelado_7_problemas_luenberger.py
│   ├── 05_convergencia_tasas_canonicas.py
│   ├── 06_soluciones_basicas_factibles.py
│   ├── 07_simplex_primal_tabular.py
│   ├── 08_dos_fases_y_gran_m.py
│   ├── 09_simplex_revisado_lu.py
│   ├── 10_dualidad_holgura_complementaria.py
│   └── 11_simplex_dual.py
│
├── notebooks/                   # Cuaderno de Jupyter interactivo
│   └── Unidad1_Herramientas_Optimizacion_Luenberger.ipynb
│
├── Material U1/                 # Material bibliográfico y notas de apoyo Unidad 1
├── Material U2/                 # Material gráfico y láminas de apoyo Unidad 2
│
└── .github/
    └── workflows/
        └── deploy.yml           # Workflow para despliegue automático en GitHub Pages
```

---

## 🚀 Despliegue y Ejecución

### 1. Visualización Web (PWA)
Para visualizar el libro interactivo en el navegador de manera local:
- Abre directamente el archivo `index.html` en tu navegador, o inicia un servidor local simple:
```bash
python -m http.server 8000
```
Luego ingresa a `http://localhost:8000`.

### 2. Ejecución de Módulos de Python
Instala las dependencias necesarias:
```bash
pip install -r requirements.txt
```
Ejecuta cualquiera de los módulos computacionales (ejemplos):
```bash
# Unidad 1
python python_src/01_algebra_lineal_formas_cuadraticas.py
python python_src/02_factorizacion_lu_gauss.py

# Unidad 2
python python_src/06_soluciones_basicas_factibles.py
python python_src/07_simplex_primal_tabular.py
python python_src/10_dualidad_holgura_complementaria.py
```

### 3. Cuaderno Interactivo en Jupyter / Google Colab
Inicia el entorno de Jupyter:
```bash
jupyter lab notebooks/Unidad1_Herramientas_Optimizacion_Luenberger.ipynb
```

---

## 📚 Bibliografía Oficial del IPN

1. **Bazaraa, M. S., Jarvis, J. J. & Sherali, H. D. (2011)**. *Linear Programming and Network Flows* (4th ed.), Wiley.
2. **Bazaraa, M. S., Sherali, H. D. & Shetty, C. M. (2013)**. *Nonlinear Programming: Theory and Algorithms* (3rd ed.), Wiley.
3. **Bonnans, J. F., Gilbert, J. C., Lemarechal, C. & Sagastizábal, C. (2006)**. *Numerical Optimization*, Springer.
4. **Cottle, R. W. & Thapa, M. N. (2017)**. *Linear and Nonlinear Optimization*, Springer.
5. **Luenberger, D. G. & Ye, Y. (2016)**. *Linear and Nonlinear Programming* (4th ed.), Springer. *(Texto Base Obligatorio)*
6. **Nocedal, J. & Wright, S. (2006)**. *Numerical Optimization* (2nd ed.), Springer.
7. **Öchsner, A. & Makvandi, R. (2020)**. *Numerical Engineering Optimization*, Springer.
