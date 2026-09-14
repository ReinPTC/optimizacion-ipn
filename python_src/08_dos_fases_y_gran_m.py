"""
================================================================================
Unidad 2 — Programación Lineal
Módulo 08: Inicialización del Simplex: Método de las Dos Fases y Método de la Gran M
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 3: El Método Simplex
           Sección 3.5: Búsqueda de una Solución Básica Factible Inicial (pp. 52-57)
================================================================================
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional


def two_phase_simplex(
    c: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    tol: float = 1e-9,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Resuelve un problema de programación lineal general en forma estándar:
        min  z = c^T x
        s.a. A x = b,   x >= 0
    utilizando el Método de las Dos Fases con variables artificiales.

    Fase I:
        Minimiza la suma de variables artificiales w = sum(x_a).
        - Si w* > 0, el conjunto factible original es VACÍO (Infactible).
        - Si w* = 0, se eliminan las columnas artificiales y se procede a la Fase II.

    Fase II:
        Restaura el vector de costos original c y optimiza a partir de la SBF obtenida.
    """
    c = np.asarray(c, dtype=float).flatten()
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).flatten()
    m, n = A.shape

    # Normalizar restricciones para asegurar b >= 0
    for i in range(m):
        if b[i] < -tol:
            b[i] *= -1.0
            A[i, :] *= -1.0

    if verbose:
        print("\n" + "=" * 60)
        print("INICIANDO MÉTODO DE LAS DOS FASES")
        print("=" * 60)

    # ---------------------------------------------------------
    # FASE I: Creación del problema auxiliar con variables artificiales
    # ---------------------------------------------------------
    # Variables originales: x_1 .. x_n
    # Variables artificiales: a_1 .. a_m
    A_phase1 = np.hstack([A, np.eye(m)])
    c_phase1 = np.zeros(n + m)
    c_phase1[n:] = 1.0  # Costo 1 para cada variable artificial
    basis_phase1 = list(range(n, n + m))

    # Tablero de Fase I
    tableau_p1 = np.zeros((m + 1, n + m + 1), dtype=float)
    tableau_p1[1:, :n+m] = A_phase1
    tableau_p1[1:, -1] = b

    # Hacer canónica la fila de costos de Fase I sumando -fila_i
    # w = sum(a_i) = sum(b_i - sum_j a_ij x_j)  =>  w + sum_j (sum_i a_ij) x_j = sum b_i
    sum_rows = np.sum(tableau_p1[1:, :], axis=0)
    tableau_p1[0, :n+m] = c_phase1 - sum_rows[:n+m]
    tableau_p1[0, -1] = -sum_rows[-1]

    # Iteraciones de Fase I
    iter_p1 = 0
    while True:
        red_costs = tableau_p1[0, :n+m]
        candidates = np.where(red_costs < -tol)[0]
        if len(candidates) == 0:
            break

        # Regla de Bland: menor índice
        col_in = candidates[0]
        y_q = tableau_p1[1:, col_in]
        rhs = tableau_p1[1:, -1]

        pos = np.where(y_q > tol)[0]
        if len(pos) == 0:
            break

        ratios = np.full(m, np.inf)
        for i in pos:
            ratios[i] = rhs[i] / y_q[i]

        min_r = np.min(ratios)
        candidates_out = np.where(np.abs(ratios - min_r) < tol)[0]
        row_out = min(candidates_out, key=lambda i: basis_phase1[i])

        # Pivoteo
        pivot = tableau_p1[row_out + 1, col_in]
        tableau_p1[row_out + 1, :] /= pivot
        for i in range(m + 1):
            if i != row_out + 1:
                tableau_p1[i, :] -= tableau_p1[i, col_in] * tableau_p1[row_out + 1, :]

        basis_phase1[row_out] = col_in
        iter_p1 += 1

    opt_w = -tableau_p1[0, -1]
    if verbose:
        print(f"-> Fin de Fase I en {iter_p1} iteraciones.")
        print(f"-> Valor óptimo de Fase I (Suma de artificiales w*): {opt_w:.6f}")

    if opt_w > tol:
        return {
            "status": "Infactible",
            "message": "El problema lineal no tiene solución factible (w* > 0 en Fase I).",
            "optimal_value": None,
            "optimal_solution": None
        }

    # ---------------------------------------------------------
    # TRANSICIÓN A FASE II: Expulsar artificiales y restaurar costos
    # ---------------------------------------------------------
    # Si alguna variable artificial sigue en la base con valor 0 (degeneración),
    # se pivotea con una columna original independiente.
    for i in range(m):
        if basis_phase1[i] >= n:
            # Buscar una variable original no básica con elemento no nulo en esa fila
            found_pivot = False
            for j in range(n):
                if j not in basis_phase1 and abs(tableau_p1[i + 1, j]) > tol:
                    # Pivotear para expulsar la artificial
                    p_val = tableau_p1[i + 1, j]
                    tableau_p1[i + 1, :] /= p_val
                    for k in range(m + 1):
                        if k != i + 1:
                            tableau_p1[k, :] -= tableau_p1[k, j] * tableau_p1[i + 1, :]
                    basis_phase1[i] = j
                    found_pivot = True
                    break

    # Construir Tablero de Fase II eliminando las columnas artificiales
    basis_phase2 = [idx for idx in basis_phase1 if idx < n]
    if len(basis_phase2) != m:
        # Fila redundante en A, recortar base
        basis_phase2 = basis_phase1[:m]

    tableau_p2 = np.zeros((m + 1, n + 1), dtype=float)
    tableau_p2[1:, :n] = tableau_p1[1:, :n]
    tableau_p2[1:, -1] = tableau_p1[1:, -1]

    # Restaurar costos originales en la fila 0: c_j - c_B^T * y_j
    c_B = c[basis_phase2]
    red_costs_p2 = c - (c_B @ tableau_p2[1:, :n])
    curr_z_p2 = float(np.dot(c_B, tableau_p2[1:, -1]))

    tableau_p2[0, :n] = red_costs_p2
    tableau_p2[0, -1] = -curr_z_p2

    if verbose:
        print("\n" + "-" * 60)
        print(f"INICIANDO FASE II (Base Inicial SBF: {basis_phase2})")
        print("-" * 60)

    # Iteraciones de Fase II
    iter_p2 = 0
    while True:
        red_costs = tableau_p2[0, :n]
        candidates = np.where(red_costs < -tol)[0]
        if len(candidates) == 0:
            break

        col_in = candidates[0]
        y_q = tableau_p2[1:, col_in]
        rhs = tableau_p2[1:, -1]

        pos = np.where(y_q > tol)[0]
        if len(pos) == 0:
            return {
                "status": "No Acotado",
                "message": "El problema lineal es no acotado hacia -infinito.",
                "optimal_value": -np.inf,
                "optimal_solution": None
            }

        ratios = np.full(m, np.inf)
        for i in pos:
            ratios[i] = rhs[i] / y_q[i]

        min_r = np.min(ratios)
        candidates_out = np.where(np.abs(ratios - min_r) < tol)[0]
        row_out = min(candidates_out, key=lambda i: basis_phase2[i])

        pivot = tableau_p2[row_out + 1, col_in]
        tableau_p2[row_out + 1, :] /= pivot
        for i in range(m + 1):
            if i != row_out + 1:
                tableau_p2[i, :] -= tableau_p2[i, col_in] * tableau_p2[row_out + 1, :]

        basis_phase2[row_out] = col_in
        iter_p2 += 1

    x_opt = np.zeros(n, dtype=float)
    for i in range(m):
        if basis_phase2[i] < n:
            x_opt[basis_phase2[i]] = tableau_p2[i + 1, -1]

    z_opt = -tableau_p2[0, -1]

    if verbose:
        print(f"-> Fin de Fase II en {iter_p2} iteraciones.")
        print(f"-> Solución Óptima x*: {x_opt}")
        print(f"-> Costo Mínimo Óptimo z*: {z_opt:.4f}")

    return {
        "status": "Óptimo",
        "optimal_solution": x_opt,
        "optimal_value": z_opt,
        "optimal_basis": basis_phase2,
        "iterations_phase1": iter_p1,
        "iterations_phase2": iter_p2,
        "total_iterations": iter_p1 + iter_p2
    }


if __name__ == "__main__":
    print("=" * 75)
    print("DEMOSTRACIÓN: MÉTODO DE LAS DOS FASES PARA ENCONTRAR SBF INICIAL")
    print("=" * 75)

    # Problema con restricciones de tipo >= e = (sin base identidad inmediata):
    # min  z = 2 x1 + 3 x2
    # s.a.   x1 + 2 x2 >= 4
    #        x1 +   x2 == 3
    #        x1, x2 >= 0
    #
    # En forma estándar con exceso e1:
    # min  z = 2 x1 + 3 x2 + 0 e1
    # s.a.   x1 + 2 x2 - e1 = 4
    #        x1 +   x2      = 3
    c_ex = np.array([2.0, 3.0, 0.0])
    A_ex = np.array([
        [1.0, 2.0, -1.0],
        [1.0, 1.0,  0.0]
    ])
    b_ex = np.array([4.0, 3.0])

    res_dos_fases = two_phase_simplex(c_ex, A_ex, b_ex, verbose=True)

    print("\nResumen Final:")
    print(f"-> Estado: {res_dos_fases['status']}")
    print(f"-> Solución x*: {res_dos_fases['optimal_solution']}")
    print(f"-> Valor Óptimo z*: {res_dos_fases['optimal_value']}")
    print(f"-> Iteraciones Totales: {res_dos_fases['total_iterations']}")
    print("=" * 75)
