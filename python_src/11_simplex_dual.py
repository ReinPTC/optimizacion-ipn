"""
================================================================================
Unidad 2 — Programación Lineal
Módulo 11: Algoritmo Simplex Dual y Reoptimización Post-Óptima
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 4: Dualidad y Complementariedad
           Sección 4.6: El Método Simplex Dual (pp. 102-108)
================================================================================
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional


class SimplexDual:
    """
    Implementación del Algoritmo Simplex Dual de Lemke.
    
    Filosofía del Simplex Dual (Luenberger Sec. 4.6):
    A diferencia del Simplex Primal (que mantiene x_B >= 0 y busca r_j >= 0),
    el Simplex Dual inicia con una base DUAL-FACTIBLE (todos los costos reducidos
    r_j >= 0), pero que puede ser PRIMAL-INFACTIBLE (algunos x_Bi < 0).
    
    Reglas de Pivoteo:
    1. Variable que Sale (Fila Pivote p):
       Se elige la fila con la mayor violación de no negatividad primal:
           p = argmin_i { x_Bi : x_Bi < 0 }
    2. Variable que Entra (Columna Pivote q):
       Prueba de la Razón Mínima Dual para preservar la dual-factibilidad (r_j >= 0):
           q = argmin_{j: y_{pj} < 0} { r_j / |y_{pj}| }
    """

    def __init__(
        self,
        c: np.ndarray,
        A: np.ndarray,
        b: np.ndarray,
        initial_basis: List[int],
        tol: float = 1e-9
    ):
        self.c = np.asarray(c, dtype=float).flatten()
        self.A = np.asarray(A, dtype=float)
        self.b = np.asarray(b, dtype=float).flatten()
        self.m, self.n = self.A.shape
        self.basis = list(initial_basis)
        self.tol = tol

        # Construcción del Tablero:
        # Fila 0: [ r_1 .. r_n | -z_0 ]
        # Filas 1..m: [ y_ij | x_Bi ]
        self.tableau = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self._initialize_tableau()

    def _initialize_tableau(self) -> None:
        B = self.A[:, self.basis]
        inv_B = np.linalg.inv(B)
        
        y_A = inv_B @ self.A
        x_B = inv_B @ self.b
        
        self.tableau[1:, :self.n] = y_A
        self.tableau[1:, self.n] = x_B

        c_B = self.c[self.basis]
        red_costs = self.c - (c_B @ y_A)
        current_z = float(np.dot(c_B, x_B))

        self.tableau[0, :self.n] = red_costs
        self.tableau[0, self.n] = -current_z

        # Verificar condición inicial de Dual-Factibilidad: r_j >= 0
        if np.any(red_costs < -self.tol):
            raise ValueError("El tablero inicial NO es dual-factible (existen r_j < 0). Se requiere r >= 0 para el Simplex Dual.")

    def print_tableau(self, iteration: int) -> None:
        print(f"\n--- Tablero Simplex Dual [Iteración {iteration}] ---")
        header = ["Base"] + [f"x{j+1}" for j in range(self.n)] + ["RHS (x_B)"]
        print(f"{header[0]:<8} | " + " | ".join(f"{h:>8}" for h in header[1:]))
        print("-" * (10 + 11 * (self.n + 1)))

        c_row = ["-z / r"] + [f"{self.tableau[0, j]:8.3f}" for j in range(self.n)] + [f"{-self.tableau[0, self.n]:8.3f}"]
        print(f"{c_row[0]:<8} | " + " | ".join(c_row[1:]))
        print("-" * (10 + 11 * (self.n + 1)))

        for i in range(self.m):
            row_label = f"x{self.basis[i]+1}"
            row_vals = [f"{self.tableau[i+1, j]:8.3f}" for j in range(self.n)] + [f"{self.tableau[i+1, self.n]:8.3f}"]
            print(f"{row_label:<8} | " + " | ".join(row_vals))
        print("-" * (10 + 11 * (self.n + 1)))

    def solve(self, max_iter: int = 50, verbose: bool = True) -> Dict[str, Any]:
        iteration = 0

        while iteration < max_iter:
            rhs = self.tableau[1:, self.n]
            red_costs = self.tableau[0, :self.n]

            # 1. Criterio de Optimalidad Primal:
            # Si todos los x_Bi >= 0, la solución es primal-factible y por ende ÓPTIMA.
            neg_rhs_indices = np.where(rhs < -self.tol)[0]
            if len(neg_rhs_indices) == 0:
                x_opt = np.zeros(self.n, dtype=float)
                for i in range(self.m):
                    x_opt[self.basis[i]] = self.tableau[i + 1, self.n]
                z_opt = -self.tableau[0, self.n]

                if verbose:
                    print(f"-> Óptimo global alcanzado en {iteration} iteraciones del Simplex Dual.")

                return {
                    "status": "Óptimo Encontrado",
                    "optimal_solution": x_opt,
                    "optimal_value": z_opt,
                    "optimal_basis": self.basis.copy(),
                    "iterations": iteration
                }

            # 2. Selección de la Variable que Sale (Fila más negativa: p)
            leaving_row = neg_rhs_indices[np.argmin(rhs[neg_rhs_indices])]
            leaving_var = self.basis[leaving_row]

            # 3. Selección de la Variable que Entra (Prueba de la Razón Mínima Dual)
            pivot_row = self.tableau[leaving_row + 1, :self.n]
            neg_elements = np.where(pivot_row < -self.tol)[0]

            if len(neg_elements) == 0:
                # Si todos los elementos de la fila pivote son >= 0, el problema primal es INFACTIBLE
                return {
                    "status": "Infactible",
                    "message": "El programa primal es infactible (la fila pivote no tiene elementos negativos).",
                    "optimal_solution": None,
                    "optimal_value": None,
                    "iterations": iteration
                }

            # Razón mínima: theta = r_j / |y_{pj}|
            ratios = np.full(self.n, np.inf)
            for j in neg_elements:
                ratios[j] = red_costs[j] / np.abs(pivot_row[j])

            entering_col = np.argmin(ratios)
            pivot_val = self.tableau[leaving_row + 1, entering_col]

            if verbose:
                print(f"Iter {iteration}: Sale x{leaving_var+1} (RHS = {rhs[leaving_row]:.3f}) | Entra x{entering_col+1} (Razón = {ratios[entering_col]:.3f})")

            # 4. Operación de Pivoteo de Gauss-Jordan
            self.tableau[leaving_row + 1, :] /= pivot_val
            for i in range(self.m + 1):
                if i != leaving_row + 1:
                    factor = self.tableau[i, entering_col]
                    self.tableau[i, :] -= factor * self.tableau[leaving_row + 1, :]

            self.basis[leaving_row] = entering_col
            iteration += 1

        return {
            "status": "Límite Máximo de Iteraciones Alcanzado",
            "iterations": iteration,
            "optimal_solution": None
        }


if __name__ == "__main__":
    print("=" * 75)
    print("DEMOSTRACIÓN: MÉTODO SIMPLEX DUAL DE LEMKE")
    print("=" * 75)

    # Problema canónico con base dual-factible inmediata:
    # min  z = 2 x1 + 1 x2
    # s.a.   3 x1 +   x2 >= 3   =>  -3 x1 -   x2 + s1 = -3
    #        4 x1 + 3 x2 >= 6   =>  -4 x1 - 3 x2 + s2 = -6
    #          x1 + 2 x2 >= 2   =>  -  x1 - 2 x2 + s3 = -2
    #        x1, x2, s1, s2, s3 >= 0
    #
    # Nótese que con c = [2, 1, 0, 0, 0] y base = [s1, s2, s3],
    # los costos reducidos iniciales son r = [2, 1, 0, 0, 0] >= 0 (¡Dual-Factible!),
    # pero el lado derecho es [-3, -6, -2] < 0 (¡Primal-Infactible!).
    c_dual = np.array([2.0, 1.0, 0.0, 0.0, 0.0])
    A_dual = np.array([
        [-3.0, -1.0, 1.0, 0.0, 0.0],
        [-4.0, -3.0, 0.0, 1.0, 0.0],
        [-1.0, -2.0, 0.0, 0.0, 1.0]
    ])
    b_dual = np.array([-3.0, -6.0, -2.0])
    base_init = [2, 3, 4]  # Columnas s1, s2, s3

    simplex_d = SimplexDual(c_dual, A_dual, b_dual, base_init)
    print("\nTablero Simplex Dual Inicial:")
    simplex_d.print_tableau(0)

    res_dual = simplex_d.solve(verbose=True)

    print("\nTablero Final:")
    simplex_d.print_tableau(res_dual["iterations"])

    print("\nResultados Simplex Dual:")
    print(f"-> Estado: {res_dual['status']}")
    print(f"-> Solución Primal Óptima x*: {res_dual['optimal_solution']}")
    print(f"-> Costo Mínimo Óptimo z*: {res_dual['optimal_value']:.4f}")
    print(f"-> Base Óptima Final: {res_dual['optimal_basis']}")
    print("=" * 75)
