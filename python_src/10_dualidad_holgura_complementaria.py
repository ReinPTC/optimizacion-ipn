"""
================================================================================
Unidad 2 — Programación Lineal
Módulo 10: Teoría de Dualidad, Holgura Complementaria y Análisis de Sensibilidad
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 4: Dualidad y Complementariedad
           Secciones 4.1 - 4.4 (pp. 75-97)
================================================================================
"""

import numpy as np
from typing import Dict, Any, Tuple


def verify_duality_and_slackness(
    c: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    tol: float = 1e-7
) -> Dict[str, Any]:
    """
    Verifica las propiedades matemáticas fundamentales del par Primal-Dual:
    
    Primal Canónico:
        min  z = c^T x
        s.a. A x >= b,   x >= 0
    
    Dual Canónico:
        max  w = b^T y
        s.a. A^T y <= c,  y >= 0
    
    Propiedades a verificar:
    1. Factibilidad Primal: A x >= b y x >= 0.
    2. Factibilidad Dual: A^T y <= c y y >= 0.
    3. Lema de Dualidad Débil: b^T y <= c^T x.
    4. Teorema de Dualidad Fuerte: Si ambos son óptimos, c^T x = b^T y.
    5. Teorema de Holgura Complementaria:
       - x_j * (c_j - (A^T y)_j) = 0  para todo j.
       - y_i * ((A x)_i - b_i) = 0    para todo i.
    """
    c = np.asarray(c, dtype=float).flatten()
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).flatten()
    x = np.asarray(x, dtype=float).flatten()
    y = np.asarray(y, dtype=float).flatten()

    m, n = A.shape

    # 1. Factibilidad Primal
    primal_slacks = A @ x - b  # debe ser >= 0
    primal_feasible = np.all(primal_slacks >= -tol) and np.all(x >= -tol)

    # 2. Factibilidad Dual
    dual_slacks = c - (A.T @ y)  # debe ser >= 0
    dual_feasible = np.all(dual_slacks >= -tol) and np.all(y >= -tol)

    # 3. Valores objetivos
    primal_obj = float(np.dot(c, x))
    dual_obj = float(np.dot(b, y))
    duality_gap = primal_obj - dual_obj

    weak_duality_holds = bool(dual_obj <= primal_obj + tol)
    strong_duality_holds = bool(abs(duality_gap) <= tol and primal_feasible and dual_feasible)

    # 4. Condiciones de Holgura Complementaria
    # x_j * dual_slack_j = 0
    comp_slack_x = np.abs(x * dual_slacks)
    # y_i * primal_slack_i = 0
    comp_slack_y = np.abs(y * primal_slacks)

    comp_slack_satisfied = bool(np.all(comp_slack_x <= tol) and np.all(comp_slack_y <= tol))

    return {
        "primal_feasible": bool(primal_feasible),
        "dual_feasible": bool(dual_feasible),
        "primal_objective": primal_obj,
        "dual_objective": dual_obj,
        "duality_gap": float(duality_gap),
        "weak_duality_holds": weak_duality_holds,
        "strong_duality_holds": strong_duality_holds,
        "complementary_slackness_holds": comp_slack_satisfied,
        "primal_slacks": primal_slacks,
        "dual_slacks": dual_slacks,
        "comp_products_x": comp_slack_x,
        "comp_products_y": comp_slack_y
    }


def shadow_price_sensitivity(
    c: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    optimal_basis: list,
    delta_b: np.ndarray
) -> Dict[str, Any]:
    """
    Calcula el impacto en el costo óptimo de un cambio en los recursos b -> b + delta_b
    usando los Precios Sombra (Shadow Prices) y = c_B^T B^(-1):
        Delta z* = y^T delta_b
    y verifica el rango de validez donde la base B sigue siendo factible.
    """
    c = np.asarray(c, dtype=float).flatten()
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).flatten()
    delta_b = np.asarray(delta_b, dtype=float).flatten()

    B = A[:, optimal_basis]
    inv_B = np.linalg.inv(B)
    c_B = c[optimal_basis]

    # Precios sombra (vector dual óptimo)
    y_star = c_B @ inv_B
    predicted_delta_z = float(np.dot(y_star, delta_b))

    # Nueva solución básica: x_B_new = B^(-1) (b + delta_b)
    x_B_orig = inv_B @ b
    x_B_new = inv_B @ (b + delta_b)
    remains_feasible = bool(np.all(x_B_new >= -1e-9))

    return {
        "shadow_prices_y": y_star,
        "predicted_delta_z": predicted_delta_z,
        "x_B_original": x_B_orig,
        "x_B_new": x_B_new,
        "basis_remains_feasible": remains_feasible,
        "interpretation": "La base óptima original sigue siendo factible; el precio sombra predice con exactitud el nuevo óptimo."
        if remains_feasible else "El cambio en recursos expulsó a la base de la región factible (se requiere reoptimizar con Simplex Dual)."
    }


if __name__ == "__main__":
    print("=" * 75)
    print("DEMOSTRACIÓN: TEOREMA DE DUALIDAD Y HOLGURA COMPLEMENTARIA")
    print("=" * 75)

    # Primal:
    # min z = 3 x1 + 2 x2
    # s.a.   2 x1 +   x2 >= 6
    #        x1 + 2 x2 >= 8
    #        x1, x2 >= 0
    c_ej = np.array([3.0, 2.0])
    A_ej = np.array([
        [2.0, 1.0],
        [1.0, 2.0]
    ])
    b_ej = np.array([6.0, 8.0])

    # Solución analítica óptima:
    # Intersección: 2x1 + x2 = 6, x1 + 2x2 = 8 => x1* = 4/3, x2* = 10/3
    x_star = np.array([4.0 / 3.0, 10.0 / 3.0])
    # Solución dual óptima:
    # 2y1 + y2 = 3, y1 + 2y2 = 2 => y1* = 4/3, y2* = 1/3
    y_star = np.array([4.0 / 3.0, 1.0 / 3.0])

    res_dual = verify_duality_and_slackness(c_ej, A_ej, b_ej, x_star, y_star)

    print("\nResultados de la Verificación Primal-Dual:")
    print(f"-> Factibilidad Primal: {res_dual['primal_feasible']}")
    print(f"-> Factibilidad Dual:   {res_dual['dual_feasible']}")
    print(f"-> Valor Primal c^T x*: {res_dual['primal_objective']:.4f}")
    print(f"-> Valor Dual   b^T y*: {res_dual['dual_objective']:.4f}")
    print(f"-> Brecha de Dualidad:  {res_dual['duality_gap']:.6e}")
    print(f"-> ¿Cumple Dualidad Fuerte?:        {res_dual['strong_duality_holds']}")
    print(f"-> ¿Cumple Holgura Complementaria?: {res_dual['complementary_slackness_holds']}")
    print(f"-> Holguras Primales (Ax - b):      {res_dual['primal_slacks']}")
    print(f"-> Holguras Duales   (c - A^T y):   {res_dual['dual_slacks']}")

    print("\n" + "-" * 75)
    print("ANÁLISIS DE SENSIBILIDAD Y PRECIOS SOMBRA (SHADOW PRICES)")
    print("-" * 75)
    # Forma estándar para A con variables de exceso: A_std = [A | -I]
    A_std = np.hstack([A_ej, -np.eye(2)])
    c_std = np.array([3.0, 2.0, 0.0, 0.0])
    base_opt = [0, 1]  # Variables x1 y x2 son básicas
    cambio_recurso = np.array([0.5, 0.0])  # Incrementar b1 en 0.5 unidades

    sens = shadow_price_sensitivity(c_std, A_std, b_ej, base_opt, cambio_recurso)
    print(f"-> Precios Sombra y*: {sens['shadow_prices_y']}")
    print(f"-> Cambio en b: {cambio_recurso}")
    print(f"-> Cambio predicho en costo mínimo Delta z*: {sens['predicted_delta_z']:.4f}")
    print(f"-> Nuevo valor x_B: {sens['x_B_new']}")
    print(f"-> Diagnóstico: {sens['interpretation']}")
    print("=" * 75)
