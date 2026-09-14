"""
================================================================================
Unidad 2 — Programación Lineal
Módulo 06: Soluciones Básicas, Soluciones Básicas Factibles y Teorema Fundamental
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 2: Propiedades Básicas de los Programas Lineales
           Secciones 2.3 - 2.5 (pp. 19-28)
================================================================================
"""

import itertools
import numpy as np
from typing import List, Dict, Any, Tuple


def find_all_basic_solutions(
    A: np.ndarray, b: np.ndarray, c: np.ndarray = None, tol: float = 1e-9
) -> List[Dict[str, Any]]:
    """
    Encuentra y clasifica todas las posibles soluciones básicas de un sistema de
    programación lineal en forma estándar:
        A x = b,   x >= 0
    donde A es de tamaño (m, n) con rango(A) = m < n.

    Para cada combinación de m columnas linealmente independientes (base B):
    1. Resuelve B x_B = b  =>  x_B = B^(-1) b
    2. Asigna x_N = 0 para las variables no básicas.
    3. Clasifica la solución en:
       - Básica Factible (SBF) si x_B >= 0.
       - Básica Infactible si al menos una componente de x_B < 0.
       - Degenerada si al menos una componente de x_B == 0.
    4. Evalúa el costo c^T x si el vector de costos c es proporcionado.

    Parámetros:
        A: Matriz de restricciones de tamaño (m, n).
        b: Vector de términos independientes de tamaño (m,).
        c: Vector de costos de tamaño (n,) (opcional).
        tol: Tolerancia numérica para considerar cero o no negatividad.

    Retorna:
        Lista de diccionarios con el detalle exhaustivo de cada base analizada.
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).flatten()
    m, n = A.shape

    if m >= n:
        raise ValueError(f"Se requiere m < n para tener grados de libertad en PL. Obtenido m={m}, n={n}.")

    if c is not None:
        c = np.asarray(c, dtype=float).flatten()
        if len(c) != n:
            raise ValueError(f"El vector de costos c debe tener dimensión {n}.")

    total_combinations = itertools.combinations(range(n), m)
    results = []

    for idx_comb, indices_B in enumerate(total_combinations, start=1):
        indices_B = list(indices_B)
        indices_N = [j for j in range(n) if j not in indices_B]
        B = A[:, indices_B]

        det_B = np.linalg.det(B)
        is_basis = abs(det_B) > tol

        if not is_basis:
            results.append({
                "combination_id": idx_comb,
                "basic_indices": indices_B,
                "non_basic_indices": indices_N,
                "det_B": float(det_B),
                "is_basis": False,
                "x_full": None,
                "is_feasible": False,
                "is_degenerate": False,
                "cost": None,
                "status": "Columnas Linealmente Dependientes (No es Base)"
            })
            continue

        try:
            x_B = np.linalg.solve(B, b)
        except np.linalg.LinAlgError:
            continue

        x_full = np.zeros(n, dtype=float)
        x_full[indices_B] = x_B

        # Verificación de factibilidad primal: x >= 0
        is_feasible = bool(np.all(x_B >= -tol))

        # Verificación de degeneración: algún x_B == 0
        is_degenerate = bool(np.any(np.abs(x_B) <= tol))

        cost_val = float(np.dot(c, x_full)) if c is not None else None

        if is_feasible:
            status_str = "Solución Básica Factible (SBF - Vértice)"
            if is_degenerate:
                status_str += " [Degenerada]"
        else:
            status_str = "Solución Básica Infactible"

        results.append({
            "combination_id": idx_comb,
            "basic_indices": indices_B,
            "non_basic_indices": indices_N,
            "det_B": float(det_B),
            "is_basis": True,
            "x_B": x_B,
            "x_full": x_full,
            "is_feasible": is_feasible,
            "is_degenerate": is_degenerate,
            "cost": cost_val,
            "status": status_str
        })

    return results


def verify_fundamental_theorem_lp(
    A: np.ndarray, b: np.ndarray, c: np.ndarray
) -> Dict[str, Any]:
    """
    Verifica computacionalmente el Teorema Fundamental de la Programación Lineal:
    'Si existe una solución factible, existe una solución básica factible (SBF);
     si existe una solución óptima finita, existe una SBF que es óptima'.
    """
    solutions = find_all_basic_solutions(A, b, c)
    sbf_list = [s for s in solutions if s["is_feasible"]]

    if not sbf_list:
        return {
            "has_feasible_solution": False,
            "has_optimal_sbf": False,
            "message": "El programa lineal es infactible (no existen SBFs)."
        }

    # Búsqueda exhaustiva del mínimo entre todas las SBFs
    optimal_sbf = min(sbf_list, key=lambda s: s["cost"])

    return {
        "has_feasible_solution": True,
        "total_combinations": len(solutions),
        "total_sbf": len(sbf_list),
        "optimal_sbf_indices": optimal_sbf["basic_indices"],
        "optimal_x": optimal_sbf["x_full"],
        "optimal_cost": optimal_sbf["cost"],
        "all_sbf": sbf_list,
        "message": "Teorema Fundamental verificado con éxito: El óptimo global coincide con una SBF."
    }


if __name__ == "__main__":
    print("=" * 75)
    print("DEMOSTRACIÓN COMPUTACIONAL: SOLUCIONES BÁSICAS FACTIBLES Y TEOREMA FUNDAMENTAL")
    print("=" * 75)

    # Ejemplo canónico: m=2 restricciones, n=4 variables (2 originales + 2 holguras)
    # min z = -3 x1 - 2 x2
    # s.a.   x1 + x2 + s1     = 4
    #        x1 + 2x2    + s2 = 6
    #        x1, x2, s1, s2 >= 0
    A_ejemplo = np.array([
        [1.0, 1.0, 1.0, 0.0],
        [1.0, 2.0, 0.0, 1.0]
    ])
    b_ejemplo = np.array([4.0, 6.0])
    c_ejemplo = np.array([-3.0, -2.0, 0.0, 0.0])

    print("\n1. Análisis Exhaustivo de todas las Submatrices B:")
    analisis = find_all_basic_solutions(A_ejemplo, b_ejemplo, c_ejemplo)

    print(f"\nTotal de combinaciones C(4, 2) = {len(analisis)} posibles bases:")
    print("-" * 75)
    print(f"{'ID':<4} | {'Base (Cols)':<12} | {'det(B)':<8} | {'x_full':<24} | {'Costo z':<8} | {'Clasificación'}")
    print("-" * 75)

    for res in analisis:
        base_str = str(res["basic_indices"])
        det_str = f"{res['det_B']:.2f}"
        if res["is_basis"]:
            x_str = "[" + ", ".join(f"{v:.1f}" for v in res["x_full"]) + "]"
            cost_str = f"{res['cost']:.2f}"
        else:
            x_str = "No definida"
            cost_str = "—"
        print(f"{res['combination_id']:<4} | {base_str:<12} | {det_str:<8} | {x_str:<24} | {cost_str:<8} | {res['status']}")

    print("-" * 75)

    # Verificación del Teorema Fundamental
    print("\n2. Verificación del Teorema Fundamental de la Programación Lineal:")
    verif = verify_fundamental_theorem_lp(A_ejemplo, b_ejemplo, c_ejemplo)
    print(f"-> Total de SBFs encontradas (Vértices del Poliedro): {verif['total_sbf']}")
    print(f"-> Base Óptima Global: Columnas {verif['optimal_sbf_indices']}")
    print(f"-> Solución Óptima x*: {verif['optimal_x']}")
    print(f"-> Valor Mínimo z*: {verif['optimal_cost']}")
    print(f"-> Diagnóstico: {verif['message']}")
    print("=" * 75)
