"""
=============================================================================
ESFM - Instituto Politécnico Nacional
Licenciatura en Matemática Algorítmica
Materia: Optimización

MÓDULO 12: CONDICIONES DE PRIMER Y SEGUNDO ORDEN EN OPTIMIZACIÓN NO LINEAL
Referencia: David G. Luenberger & Yinyu Ye, "Linear and Nonlinear Programming"
            (4ta Edición), Capítulo 7: Basic Properties of Solutions & Algorithms
=============================================================================
Este módulo computacional implementa y valida:
  1. Verificación analítica y numérica de Condiciones Necesarias de Primer Orden
     (FONC): Gradiente nulo ∇f(x*) = 0.
  2. Condiciones de Segundo Orden (SONC y SOSC):
     - Cálculo y simetrización de la matriz Hessiana F(x*) = ∇²f(x*).
     - Descomposición espectral y signatura de autovalores (Criterio de Sylvester).
     - Clasificación rigurosa: Mínimo local estricto, Máximo local estricto,
       Punto de Ensilladura (Saddle Point) o Punto Degenerado / Semidefinido.
  3. Ejemplos canónicos del texto de Luenberger & Ye:
     - Ejemplo 7.1.1 (p. 182): f(x1, x2) = x1² - x1*x2 + x2² - 3*x2
     - Ejemplo 7.3.2 (p. 186): f(x1, x2) = x1³ - x1²*x2 + 2*x2² (Punto de silla en (6,9))
     - Ejemplo 7.2.3 (p. 184): Mínimos Cuadrados Polinomiales y Ecuaciones Normales.
=============================================================================
"""

import sys
import numpy as np

# Asegurar codificación UTF-8 en consola de Windows
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def gradiente_numerico(f, x, h=1e-6):
    """
    Aproxima el gradiente ∇f(x) mediante diferencias finitas centrales.
    Error de aproximación: O(h²).
    """
    n = len(x)
    grad = np.zeros(n)
    for i in range(n):
        e_i = np.zeros(n)
        e_i[i] = h
        grad[i] = (f(x + e_i) - f(x - e_i)) / (2.0 * h)
    return grad


def hessiano_numerico(f, x, h=1e-5):
    """
    Aproxima la matriz Hessiana F(x) = ∇²f(x) mediante diferencias finitas.
    Garantiza simetría numérica F = (F + F^T) / 2.
    """
    n = len(x)
    H = np.zeros((n, n))
    f_x = f(x)
    for i in range(n):
        e_i = np.zeros(n)
        e_i[i] = h
        for j in range(i, n):
            e_j = np.zeros(n)
            e_j[j] = h
            if i == j:
                # Segunda derivada pura: (f(x+h) - 2f(x) + f(x-h)) / h²
                H[i, i] = (f(x + e_i) - 2.0 * f_x + f(x - e_i)) / (h ** 2)
            else:
                # Derivada cruzada: (f(x+h_i+h_j) - f(x+h_i-h_j) - f(x-h_i+h_j) + f(x-h_i-h_j)) / (4h²)
                f_pp = f(x + e_i + e_j)
                f_pm = f(x + e_i - e_j)
                f_mp = f(x - e_i + e_j)
                f_mm = f(x - e_i - e_j)
                val = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h ** 2)
                H[i, j] = val
                H[j, i] = val
    return H


def clasificar_punto_estacionario(f, x_star, tol_grad=1e-4, tol_eig=1e-5):
    """
    Evalúa FONC, SONC y SOSC para clasificar formalmente un punto candidato x*.
    
    Retorna:
      dict con: gradiente, hessiano, autovalores, clasificación y detalles matemáticos.
    """
    grad = gradiente_numerico(f, x_star)
    norm_grad = np.linalg.norm(grad)
    es_estacionario = norm_grad < tol_grad

    H = hessiano_numerico(f, x_star)
    autovalores, autovectores = np.linalg.eigh(H)
    det_H = np.linalg.det(H)
    traza_H = np.trace(H)

    pos_eigs = np.sum(autovalores > tol_eig)
    neg_eigs = np.sum(autovalores < -tol_eig)
    cero_eigs = len(autovalores) - pos_eigs - neg_eigs

    if not es_estacionario:
        clasificacion = "NO ESTACIONARIO (∇f ≠ 0; punto ordinario)"
        descripcion = f"La norma del gradiente ||∇f|| = {norm_grad:.4e} supera la tolerancia {tol_grad}."
    elif pos_eigs == len(autovalores):
        clasificacion = "MÍNIMO LOCAL ESTRICTO (SOSC Satisfecha: F ≻ 0)"
        descripcion = (
            f"Todos los autovalores son estrictamente positivos ({autovalores}). "
            f"La forma cuadrática d^T F d >= λ_min ||d||² > 0 garantiza mínimo local aislado."
        )
    elif neg_eigs == len(autovalores):
        clasificacion = "MÁXIMO LOCAL ESTRICTO (SOSC Satisfecha: F ≺ 0)"
        descripcion = f"Todos los autovalores son estrictamente negativos ({autovalores})."
    elif pos_eigs > 0 and neg_eigs > 0:
        clasificacion = "PUNTO DE ENSILLADURA (Saddle Point: F Indefinida)"
        descripcion = (
            f"Existen autovalores de signos opuestos (positivos: {pos_eigs}, negativos: {neg_eigs}). "
            f"det(F) = {det_H:.4e}. Hay direcciones de descenso y de ascenso."
        )
    else:
        clasificacion = "PUNTO DEGENERADO / SEMIDEFINIDO (Criterio Inconcluso: det(F) ≈ 0)"
        descripcion = (
            f"Al menos un autovalor es prácticamente nulo ({cero_eigs} autovalores ≈ 0). "
            f"Se requiere análisis de derivadas de orden superior o curvas direccionales."
        )

    return {
        "x": x_star,
        "f_val": f(x_star),
        "gradiente": grad,
        "norma_grad": norm_grad,
        "es_estacionario": es_estacionario,
        "hessiano": H,
        "autovalores": autovalores,
        "det_hessiano": det_H,
        "traza_hessiano": traza_H,
        "clasificacion": clasificacion,
        "descripcion": descripcion,
    }


def resolver_minimos_cuadrados_polinomial(t_puntos, y_puntos, grado=2):
    """
    Ejemplo 7.2.3 de Luenberger: Ajuste de curvas polinomiales mediante
    resolución de Ecuaciones Normales Q a = b, donde Q = V^T V es definida positiva.
    """
    m = len(t_puntos)
    # Matriz de Vandermonde: V_ij = t_i^j
    V = np.vander(t_puntos, N=grado + 1, increasing=True)
    Q = V.T @ V
    b = V.T @ y_puntos
    
    # Resolver ecuaciones normales
    a_coef = np.linalg.solve(Q, b)
    
    # Calcular residuo cuadrático mínimo
    y_pred = V @ a_coef
    residuo = np.sum((y_pred - y_puntos) ** 2)
    
    # Autovalores del Hessiano Q (la función de costo es f(a) = 1/2 a^T (2Q) a - 2b^T a + c)
    autovals_Q = np.linalg.eigvalsh(2 * Q)
    
    return {
        "grado": grado,
        "coeficientes_a": a_coef,
        "matriz_Q": Q,
        "vector_b": b,
        "autovalores_Hessiano": autovals_Q,
        "residuo_minimo": residuo
    }


# =============================================================================
# EJECUCIÓN Y DEMOSTRACIÓN PEDAGÓGICA (LUENBERGER CAP. 7)
# =============================================================================
if __name__ == "__main__":
    print("=" * 80)
    print("DEMOSTRACIÓN DE CONDICIONES DE PRIMER Y SEGUNDO ORDEN (LUENBERGER CAP. 7)")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # CASO 1: Ejemplo 7.1.1 (p. 182)
    # f(x1, x2) = x1^2 - x1*x2 + x2^2 - 3*x2
    # ∇f = [2*x1 - x2, -x1 + 2*x2 - 3] = [0, 0]  ==> x* = [1, 2]
    # F = [[2, -1], [-1, 2]]  ==> det(F) = 3 > 0, traza = 4 > 0 ==> Mínimo estricto
    # -------------------------------------------------------------------------
    print("\n--- CASO 1: Ejemplo 7.1.1 (Luenberger p. 182) ---")
    f1 = lambda x: x[0]**2 - x[0]*x[1] + x[1]**2 - 3.0*x[1]
    x_candidato_1 = np.array([1.0, 2.0])
    res1 = clasificar_punto_estacionario(f1, x_candidato_1)
    
    print(f"Punto analizado: x* = {res1['x']}")
    print(f"Valor f(x*): {res1['f_val']:.4f}")
    print(f"Gradiente ∇f(x*): {res1['gradiente']} (||∇f|| = {res1['norma_grad']:.2e})")
    print(f"Matriz Hessiana F(x*):\n{res1['hessiano']}")
    print(f"Autovalores λ(F): {res1['autovalores']}")
    print(f"det(F): {res1['det_hessiano']:.4f}, traza(F): {res1['traza_hessiano']:.4f}")
    print(f"Diagnóstico: {res1['clasificacion']}")
    print(f"Fundamento: {res1['descripcion']}")

    # -------------------------------------------------------------------------
    # CASO 2: Ejemplo 7.3.2 (p. 186)
    # f(x1, x2) = x1^3 - x1^2*x2 + 2*x2^2
    # Puntos estacionarios: (6, 9) y (0, 0)
    # En (6, 9): F = [[18, -12], [-12, 4]] ==> det(F) = 72 - 144 = -72 < 0 ==> Silla!
    # En (0, 0): F = [[0, 0], [0, 4]] ==> det(F) = 0 ==> Degenerado
    # -------------------------------------------------------------------------
    print("\n--- CASO 2: Ejemplo 7.3.2 (Luenberger p. 186) ---")
    f2 = lambda x: x[0]**3 - (x[0]**2)*x[1] + 2.0*(x[1]**2)
    
    # Análisis en (6, 9)
    x_silla = np.array([6.0, 9.0])
    res2_silla = clasificar_punto_estacionario(f2, x_silla)
    print(f"\n[Subcaso A: Punto (6, 9)]")
    print(f"Punto: {res2_silla['x']}, f(x) = {res2_silla['f_val']}")
    print(f"Gradiente: {res2_silla['gradiente']}, ||∇f|| = {res2_silla['norma_grad']:.2e}")
    print(f"Hessiana F(6,9):\n{res2_silla['hessiano']}")
    print(f"Autovalores: {res2_silla['autovalores']}")
    print(f"det(F): {res2_silla['det_hessiano']:.4f}")
    print(f"Diagnóstico: {res2_silla['clasificacion']}")
    print(f"Fundamento: {res2_silla['descripcion']}")

    # Análisis en (0, 0)
    x_origen = np.array([0.0, 0.0])
    res2_origen = clasificar_punto_estacionario(f2, x_origen)
    print(f"\n[Subcaso B: Origen (0, 0)]")
    print(f"Punto: {res2_origen['x']}, f(x) = {res2_origen['f_val']}")
    print(f"Gradiente: {res2_origen['gradiente']}, ||∇f|| = {res2_origen['norma_grad']:.2e}")
    print(f"Hessiana F(0,0):\n{res2_origen['hessiano']}")
    print(f"Autovalores: {res2_origen['autovalores']}")
    print(f"det(F): {res2_origen['det_hessiano']:.4f}")
    print(f"Diagnóstico: {res2_origen['clasificacion']}")
    print(f"Fundamento: {res2_origen['descripcion']}")

    # -------------------------------------------------------------------------
    # CASO 3: Ejemplo 7.2.3 (p. 184)
    # Mínimos Cuadrados y Ecuaciones Normales
    # -------------------------------------------------------------------------
    print("\n--- CASO 3: Ejemplo 7.2.3 - Mínimos Cuadrados Polinomiales ---")
    t_data = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    # y = 1.5 - 0.8*t + 0.5*t^2 + ruido
    y_data = np.array([1.52, 1.18, 1.95, 3.71, 6.48])
    res_lsq = resolver_minimos_cuadrados_polinomial(t_data, y_data, grado=2)
    print(f"Puntos observados (t): {t_data}")
    print(f"Valores observados (y): {y_data}")
    print(f"Coeficientes óptimos a*: {res_lsq['coeficientes_a']} (a0 + a1*t + a2*t²)")
    print(f"Matriz de Ecuaciones Normales Q = V^T V:\n{res_lsq['matriz_Q']}")
    print(f"Autovalores del Hessiano 2Q: {res_lsq['autovalores_Hessiano']}")
    print(f"¿Hessiano estrictamente convexo?: {np.all(res_lsq['autovalores_Hessiano'] > 0)}")
    print(f"Residuo mínimo suma de cuadrados: {res_lsq['residuo_minimo']:.6f}")
    print("\n" + "=" * 80)
    print("FIN DE LA VALIDACIÓN COMPUTACIONAL DEL MÓDULO 12")
    print("=" * 80)
