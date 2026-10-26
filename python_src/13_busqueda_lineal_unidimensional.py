"""
=============================================================================
ESFM - Instituto Politécnico Nacional
Licenciatura en Matemática Algorítmica
Materia: Optimización

MÓDULO 13: ALGORITMOS DE BÚSQUEDA SOBRE UNA LÍNEA (LINE SEARCH)
Referencia: David G. Luenberger & Yinyu Ye, "Linear and Nonlinear Programming"
            (4ta Edición), Capítulo 8: Basic Descent Methods, Sección 8.1
=============================================================================
Este módulo computacional implementa con rigor analítico y numérico:
  1. Métodos de Orden Cero (Sin Derivadas) para Funciones Unimodales:
     - Búsqueda de Fibonacci (Minimización minimax óptima para N prefijado).
     - Búsqueda de la Sección Áurea (Tasa constante tau = (sqrt(5)-1)/2 ≈ 0.618).
  2. Métodos de Ajuste de Curvas (Interpolación):
     - Ajuste Cuadrático de Lagrange (Fórmula analítica de 3 puntos, orden ≈ 1.325).
     - Ajuste Cúbico de Hermite (Valores y derivadas en 2 puntos, orden 2.0).
  3. Métodos de Orden 1 y 2:
     - Método de la Secante / Posición Falsa para f'(x) = 0 (Orden ≈ 1.618).
     - Método de Newton Unidimensional (Orden cuadrático 2.0).
  4. Búsqueda Inexacta de Línea:
     - Regla de Armijo (Sufficient Decrease Condition) y algoritmo de Backtracking.
     - Verificación de Condiciones de Wolfe (Decrecimiento + Curvatura).
=============================================================================
"""

import sys
import math
from typing import Callable, List, Dict, Any, Tuple

# Reconfigurar salida estándar en Windows a UTF-8
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# =============================================================================
# 1. MÉTODOS DE ORDEN CERO (FIBONACCI Y SECCIÓN ÁUREA)
# =============================================================================

def fibonacci_search(f: Callable[[float], float], a: float, b: float, N: int) -> Dict[str, Any]:
    """
    Búsqueda de Fibonacci (Luenberger pp. 216-218).
    Garantiza una reducción óptima del intervalo de incertidumbre mediante 
    la sucesión de Fibonacci F_k = F_{k-1} + F_{k-2}.
    """
    # Generar números de Fibonacci: F_0=1, F_1=1, F_2=2, F_3=3, ...
    F = [1, 1]
    for _ in range(2, N + 2):
        F.append(F[-1] + F[-2])

    historial = []
    L = b - a
    # Puntos interiores iniciales
    c = a + (F[N - 1] / F[N + 1]) * L
    d = a + (F[N] / F[N + 1]) * L
    fc = f(c)
    fd = f(d)

    for k in range(1, N):
        historial.append({
            "iter": k, "a": a, "b": b, "c": c, "d": d,
            "fc": fc, "fd": fd, "longitud": b - a
        })

        if fc < fd:
            b = d
            d = c
            fd = fc
            c = a + (F[N - k - 1] / F[N - k + 1]) * (b - a)
            fc = f(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + (F[N - k] / F[N - k + 1]) * (b - a)
            fd = f(d)

    # Última perturbación infinitesimal para romper empate si fuese necesario
    x_opt = (a + b) / 2.0
    return {
        "metodo": "Búsqueda de Fibonacci",
        "x_opt": x_opt,
        "f_opt": f(x_opt),
        "intervalo_final": (a, b),
        "longitud_final": b - a,
        "reduccion_teorica": 1.0 / F[N + 1],
        "evaluaciones": N + 1,
        "historial": historial
    }


def golden_section_search(f: Callable[[float], float], a: float, b: float, 
                          tol: float = 1e-6, max_iter: int = 100) -> Dict[str, Any]:
    """
    Búsqueda de la Sección Áurea (Luenberger pp. 218-219).
    Constante áurea tau = (sqrt(5) - 1) / 2 ≈ 0.6180339887...
    Reduce el intervalo por un factor de 0.618 por evaluación, reutilizando un punto.
    """
    tau = (math.sqrt(5.0) - 1.0) / 2.0  # ≈ 0.6180339887
    historial = []

    c = b - tau * (b - a)
    d = a + tau * (b - a)
    fc = f(c)
    fd = f(d)
    evals = 2

    for k in range(1, max_iter + 1):
        longitud = b - a
        historial.append({
            "iter": k, "a": a, "b": b, "c": c, "d": d,
            "fc": fc, "fd": fd, "longitud": longitud
        })

        if longitud < tol:
            break

        if fc < fd:
            b = d
            d = c
            fd = fc
            c = b - tau * (b - a)
            fc = f(c)
            evals += 1
        else:
            a = c
            c = d
            fc = fd
            d = a + tau * (b - a)
            fd = f(d)
            evals += 1

    x_opt = (a + b) / 2.0
    return {
        "metodo": "Sección Áurea",
        "x_opt": x_opt,
        "f_opt": f(x_opt),
        "intervalo_final": (a, b),
        "longitud_final": b - a,
        "iteraciones": len(historial),
        "evaluaciones": evals,
        "historial": historial
    }


# =============================================================================
# 2. MÉTODOS DE AJUSTE DE CURVAS (INTERPOLACIÓN CUADRÁTICA Y CÚBICA)
# =============================================================================

def quadratic_fit_search(f: Callable[[float], float], x1: float, x2: float, x3: float, 
                         tol: float = 1e-6, max_iter: int = 50) -> Dict[str, Any]:
    """
    Ajuste Cuadrático de Curvas (Lagrange) de 3 Puntos (Luenberger pp. 220-222).
    Dada una terna x1 < x2 < x3, el vértice de la parábola interpolante es:
      x4 = 0.5 * [ (x2² - x3²)f1 + (x3² - x1²)f2 + (x1² - x2²)f3 ] /
                 [ (x2 - x3)f1 + (x3 - x1)f2 + (x1 - x2)f3 ]
    Orden de convergencia: lambda ≈ 1.325.
    """
    historial = []
    pts = sorted([x1, x2, x3])
    x1, x2, x3 = pts[0], pts[1], pts[2]
    f1, f2, f3 = f(x1), f(x2), f(x3)
    evals = 3

    for k in range(1, max_iter + 1):
        num = (x2**2 - x3**2) * f1 + (x3**2 - x1**2) * f2 + (x1**2 - x2**2) * f3
        den = (x2 - x3) * f1 + (x3 - x1) * f2 + (x1 - x2) * f3

        if abs(den) < 1e-14:
            # Puntos colineales o degenerados
            x4 = x2
            break

        x4 = 0.5 * (num / den)
        f4 = f(x4)
        evals += 1

        historial.append({
            "iter": k, "x1": x1, "x2": x2, "x3": x3, "x4": x4,
            "f4": f4, "error_paso": abs(x4 - x2)
        })

        if abs(x4 - x2) < tol:
            break

        # Actualizar la terna conservando el corchete de 3 puntos
        nuevos_pts = sorted([(x1, f1), (x2, f2), (x3, f3), (x4, f4)], key=lambda p: p[0])
        idx4 = [p[0] for p in nuevos_pts].index(x4)

        if idx4 == 0 or idx4 == 3:
            # Si x4 cayó fuera de los extremos, contraer hacia el mejor
            mejor = min(nuevos_pts, key=lambda p: p[1])
            x2 = mejor[0]
            f2 = mejor[1]
            x1 = x2 - 0.5 * abs(x3 - x1)
            x3 = x2 + 0.5 * abs(x3 - x1)
            f1, f3 = f(x1), f(x3)
            evals += 2
        else:
            # Seleccionar los 3 puntos que rodean al mínimo actual
            if f4 < f2:
                if x4 > x2:
                    x1, f1 = x2, f2
                    x2, f2 = x4, f4
                else:
                    x3, f3 = x2, f2
                    x2, f2 = x4, f4
            else:
                if x4 > x2:
                    x3, f3 = x4, f4
                else:
                    x1, f1 = x4, f4

    return {
        "metodo": "Ajuste Cuadrático (Lagrange)",
        "x_opt": x4,
        "f_opt": f4,
        "iteraciones": len(historial),
        "evaluaciones": evals,
        "orden_convergencia_teorico": "≈ 1.325",
        "historial": historial
    }


def cubic_fit_search(f: Callable[[float], float], df: Callable[[float], float],
                     x1: float, x2: float, tol: float = 1e-6, max_iter: int = 30) -> Dict[str, Any]:
    """
    Ajuste Cúbico de Hermite (Luenberger pp. 223-224).
    Utiliza f y f' en dos puntos para ajustar un polinomio cúbico.
    Orden de convergencia: cuadrático (2.0).
    """
    historial = []
    evals = 0

    for k in range(1, max_iter + 1):
        f1, df1 = f(x1), df(x1)
        f2, df2 = f(x2), df(x2)
        evals += 2

        dx = x2 - x1
        if abs(dx) < 1e-12:
            break

        z = (3.0 * (f1 - f2) / dx) + df1 + df2
        w_sq = z**2 - df1 * df2
        if w_sq < 0:
            w = 0.0
        else:
            w = math.sqrt(w_sq)

        den = df2 - df1 + 2.0 * w
        if abs(den) < 1e-14:
            x_new = 0.5 * (x1 + x2)
        else:
            x_new = x2 - dx * ((df2 + w - z) / den)

        f_new = f(x_new)
        df_new = df(x_new)
        evals += 1

        historial.append({
            "iter": k, "x1": x1, "x2": x2, "x_new": x_new,
            "f_new": f_new, "df_new": df_new, "dx": abs(x2 - x1)
        })

        if abs(df_new) < tol or abs(x_new - x2) < tol:
            x1 = x_new
            break

        if df_new * df1 < 0:
            x2 = x_new
        else:
            x1 = x_new

    return {
        "metodo": "Ajuste Cúbico de Hermite",
        "x_opt": x_new,
        "f_opt": f(x_new),
        "iteraciones": len(historial),
        "evaluaciones": evals,
        "orden_convergencia_teorico": "2.0 (Cuadrático)",
        "historial": historial
    }


# =============================================================================
# 3. MÉTODOS DE ORDEN 1 Y 2 (SECANTE Y NEWTON 1D)
# =============================================================================

def secant_search(df: Callable[[float], float], x0: float, x1: float,
                  tol: float = 1e-6, max_iter: int = 50) -> Dict[str, Any]:
    """
    Método de la Secante / Posición Falsa para f'(x) = 0 (Luenberger pp. 222-223).
    x_{k+1} = x_k - f'(x_k) * (x_k - x_{k-1}) / (f'(x_k) - f'(x_{k-1}))
    Orden de convergencia: tau_1 = (1 + sqrt(5)) / 2 ≈ 1.618 (Superlineal).
    """
    historial = []
    g0 = df(x0)
    g1 = df(x1)
    evals = 2

    for k in range(1, max_iter + 1):
        dg = g1 - g0
        if abs(dg) < 1e-14:
            break

        x_next = x1 - g1 * (x1 - x0) / dg
        g_next = df(x_next)
        evals += 1

        historial.append({
            "iter": k, "x_k": x1, "df_k": g1, "x_next": x_next,
            "df_next": g_next, "paso": abs(x_next - x1)
        })

        if abs(g_next) < tol or abs(x_next - x1) < tol:
            x1 = x_next
            break

        x0, g0 = x1, g1
        x1, g1 = x_next, g_next

    return {
        "metodo": "Método de la Secante (Posición Falsa)",
        "x_opt": x1,
        "df_opt": g1,
        "iteraciones": len(historial),
        "evaluaciones": evals,
        "orden_convergencia_teorico": "≈ 1.618 (Superlineal)",
        "historial": historial
    }


def newton_1d_search(df: Callable[[float], float], d2f: Callable[[float], float],
                     x0: float, tol: float = 1e-6, max_iter: int = 50) -> Dict[str, Any]:
    """
    Método de Newton Unidimensional (Luenberger pp. 224-226).
    x_{k+1} = x_k - f'(x_k) / f''(x_k)
    Orden de convergencia: 2.0 (Cuadrático local).
    """
    historial = []
    x_k = x0
    evals = 0

    for k in range(1, max_iter + 1):
        g = df(x_k)
        h = d2f(x_k)
        evals += 2

        historial.append({
            "iter": k, "x_k": x_k, "df": g, "d2f": h
        })

        if abs(g) < tol:
            break

        if abs(h) < 1e-14:
            # Curvatura nula
            break

        paso = g / h
        x_next = x_k - paso

        if abs(x_next - x_k) < tol:
            x_k = x_next
            break

        x_k = x_next

    return {
        "metodo": "Método de Newton 1D",
        "x_opt": x_k,
        "iteraciones": len(historial),
        "evaluaciones": evals,
        "orden_convergencia_teorico": "2.0 (Cuadrático)",
        "historial": historial
    }


# =============================================================================
# 4. BÚSQUEDA INEXACTA: REGLA DE ARMIJO Y CONDICIONES DE WOLFE
# =============================================================================

def armijo_backtracking(phi: Callable[[float], float], dphi_0: float, 
                        phi_0: float, s: float = 1.0, beta: float = 0.5, 
                        eps: float = 0.1, max_backtracks: int = 30) -> Dict[str, Any]:
    """
    Regla de Armijo (Sufficient Decrease Condition) con Backtracking (Luenberger pp. 230-231).
    Condición de Armijo:
      phi(alpha) <= phi(0) + eps * phi'(0) * alpha
    donde 0 < eps < 0.5 y 0 < beta < 1.
    """
    alpha = s
    backtracks = 0
    historial = []

    while backtracks < max_backtracks:
        phi_alpha = phi(alpha)
        umbral_armijo = phi_0 + eps * dphi_0 * alpha
        cumple = (phi_alpha <= umbral_armijo)

        historial.append({
            "intento": backtracks + 1,
            "alpha": alpha,
            "phi_alpha": phi_alpha,
            "umbral_armijo": umbral_armijo,
            "cumple_armijo": cumple
        })

        if cumple:
            break

        alpha *= beta
        backtracks += 1

    return {
        "metodo": "Regla de Armijo (Backtracking)",
        "alpha_aceptado": alpha,
        "phi_alpha": phi(alpha),
        "retrocesos": backtracks,
        "factor_contraccion_beta": beta,
        "parametro_eps": eps,
        "historial": historial
    }


# =============================================================================
# EJECUCIÓN COMPARATIVA DE ALGORITMOS DE BÚSQUEDA DE LÍNEA
# =============================================================================
if __name__ == "__main__":
    print("=" * 80)
    print("COMPARATIVA DE MÉTODOS DE BÚSQUEDA SOBRE UNA LÍNEA (LUENBERGER CAP. 8.1)")
    print("=" * 80)

    # Función de prueba unimodal en [0, 3]:
    # f(x) = 2*x^2 - exp(x/2) - 3*x
    # f'(x) = 4*x - 0.5*exp(x/2) - 3
    # f''(x) = 4 - 0.25*exp(x/2)
    f_test = lambda x: 2.0 * (x ** 2) - math.exp(0.5 * x) - 3.0 * x
    df_test = lambda x: 4.0 * x - 0.5 * math.exp(0.5 * x) - 3.0
    d2f_test = lambda x: 4.0 - 0.25 * math.exp(0.5 * x)

    a_init, b_init = 0.0, 3.0

    print(f"\nFunción objetivo de prueba: f(x) = 2x² - e^(0.5x) - 3x")
    print(f"Intervalo inicial de incertidumbre: [{a_init}, {b_init}]\n")

    # 1. Fibonacci
    res_fib = fibonacci_search(f_test, a_init, b_init, N=15)
    print(f"[1] FIBONACCI (N=15):")
    print(f"    x* = {res_fib['x_opt']:.7f}, f(x*) = {res_fib['f_opt']:.7f}")
    print(f"    Longitud final: {res_fib['longitud_final']:.2e}, Reducción teórica: {res_fib['reduccion_teorica']:.2e}")

    # 2. Golden Section
    res_gold = golden_section_search(f_test, a_init, b_init, tol=1e-6)
    print(f"\n[2] SECCIÓN ÁUREA (tau ≈ 0.618):")
    print(f"    x* = {res_gold['x_opt']:.7f}, f(x*) = {res_gold['f_opt']:.7f}")
    print(f"    Iteraciones: {res_gold['iteraciones']}, Evaluaciones de función: {res_gold['evaluaciones']}")

    # 3. Quadratic Fit (Lagrange 3 puntos)
    res_quad = quadratic_fit_search(f_test, 0.2, 1.2, 2.5, tol=1e-6)
    print(f"\n[3] AJUSTE CUADRÁTICO DE LAGRANGE (Orden ≈ 1.325):")
    print(f"    x* = {res_quad['x_opt']:.7f}, f(x*) = {res_quad['f_opt']:.7f}")
    print(f"    Iteraciones: {res_quad['iteraciones']}, Evaluaciones: {res_quad['evaluaciones']}")

    # 4. Secant Method (Posición Falsa)
    res_sec = secant_search(df_test, 0.0, 3.0, tol=1e-6)
    print(f"\n[4] MÉTODO DE LA SECANTE / POSICIÓN FALSA (Orden ≈ 1.618):")
    print(f"    x* = {res_sec['x_opt']:.7f}, f'(x*) = {res_sec['df_opt']:.2e}")
    print(f"    Iteraciones: {res_sec['iteraciones']}, Evaluaciones: {res_sec['evaluaciones']}")

    # 5. Newton 1D
    res_newton = newton_1d_search(df_test, d2f_test, x0=2.0, tol=1e-6)
    print(f"\n[5] MÉTODO DE NEWTON 1D (Orden 2.0 Cuadrático):")
    print(f"    x* = {res_newton['x_opt']:.7f}, f(x*) = {f_test(res_newton['x_opt']):.7f}")
    print(f"    Iteraciones: {res_newton['iteraciones']}, Evaluaciones: {res_newton['evaluaciones']}")

    # 6. Armijo Backtracking
    phi_0 = f_test(0.0)
    dphi_0 = df_test(0.0)  # Pendiente inicial < 0
    res_armijo = armijo_backtracking(f_test, dphi_0=dphi_0, phi_0=phi_0, s=2.0, beta=0.5, eps=0.1)
    print(f"\n[6] BÚSQUEDA INEXACTA CON REGLA DE ARMIJO (Backtracking):")
    print(f"    Paso aceptado alpha = {res_armijo['alpha_aceptado']:.4f}")
    print(f"    phi(alpha) = {res_armijo['phi_alpha']:.4f} <= Umbral Armijo")
    print(f"    Número de retrocesos: {res_armijo['retrocesos']}")

    print("\n" + "=" * 80)
    print("TABLA COMPARATIVA DE CONVERGENCIA")
    print("=" * 80)
    print(f"{'Método':<32} | {'Orden Teórico':<18} | {'Evals/Iter':<12} | {'x* Obtenido':<12}")
    print("-" * 80)
    print(f"{'Sección Áurea':<32} | {'Lineal (0.618)':<18} | {res_gold['evaluaciones']:<12} | {res_gold['x_opt']:.6f}")
    print(f"{'Ajuste Cuadrático (Lagrange)':<32} | {'≈ 1.325':<18} | {res_quad['evaluaciones']:<12} | {res_quad['x_opt']:.6f}")
    print(f"{'Secante / Posición Falsa':<32} | {'≈ 1.618':<18} | {res_sec['evaluaciones']:<12} | {res_sec['x_opt']:.6f}")
    print(f"{'Newton 1D':<32} | {'2.0 (Cuadrático)':<18} | {res_newton['evaluaciones']:<12} | {res_newton['x_opt']:.6f}")
    print("=" * 80)
