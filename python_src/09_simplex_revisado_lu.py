"""
================================================================================
Unidad 2 — Programación Lineal
Módulo 09: Método Simplex Revisado y Actualización con Descomposición LU
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 3: El Método Simplex
           Sección 3.7: The Revised Simplex Method and LU Decomposition (pp. 60-66)
================================================================================
"""

import numpy as np
import scipy.linalg as la
from typing import List, Dict, Any, Tuple, Optional


class SimplexRevisadoLU:
    """
    Implementación del Método Simplex Revisado empleando factorización LU
    para resolver sistemas triangulares en cada iteración:
        1. B^T y = c_B           =>  L U y = P c_B   (Multiplicadores duales y)
        2. B d = a_q             =>  P L U d = a_q   (Dirección simplex d)
        3. B x_B = b             =>  P L U x_B = b   (Valores de variables básicas)
    
    Ventaja computacional de Luenberger (Sec. 3.7):
    No almacena ni actualiza el tablero completo (m x n), sino únicamente los factores LU
    de la base (m x m), lo cual es ideal para problemas a gran escala y matrices ralas.
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

        if len(self.basis) != self.m:
            raise ValueError(f"La base debe tener longitud m={self.m}.")

    def _solve_B(self, P: np.ndarray, L: np.ndarray, U: np.ndarray, rhs: np.ndarray) -> np.ndarray:
        """Resuelve B x = rhs mediante P L U x = rhs => L y = P^T rhs, U x = y."""
        # P L U x = rhs  => L U x = P^T rhs
        rhs_p = P.T @ rhs
        y = la.solve_triangular(L, rhs_p, lower=True, unit_diagonal=True)
        x = la.solve_triangular(U, y, lower=False)
        return x

    def _solve_BT(self, P: np.ndarray, L: np.ndarray, U: np.ndarray, rhs: np.ndarray) -> np.ndarray:
        """Resuelve B^T y = rhs mediante U^T L^T P y = rhs."""
        # B = P L U  =>  B^T = U^T L^T P^T
        # U^T w = rhs,  L^T v = w,  y = P v
        w = la.solve_triangular(U.T, rhs, lower=True)
        v = la.solve_triangular(L.T, w, lower=False, unit_diagonal=True)
        y = P @ v
        return y

    def solve(self, max_iter: int = 100, verbose: bool = True) -> Dict[str, Any]:
        """Ejecuta el Simplex Revisado."""
        iteration = 0

        while iteration < max_iter:
            # 1. Extraer submatriz base y factorizar B = P L U
            B = self.A[:, self.basis]
            P, L, U = la.lu(B)

            # 2. Calcular valores de las variables básicas: B x_B = b
            x_B = self._solve_B(P, L, U, self.b)

            # 3. Calcular multiplicadores simplex: B^T y = c_B
            c_B = self.c[self.basis]
            y = self._solve_BT(P, L, U, c_B)

            # 4. Calcular costos reducidos para variables no básicas: r_j = c_j - y^T a_j
            non_basis = [j for j in range(self.n) if j not in self.basis]
            reduced_costs = self.c[non_basis] - (y @ self.A[:, non_basis])

            # Criterio de optimalidad
            candidate_entering = [non_basis[k] for k, r in enumerate(reduced_costs) if r < -self.tol]

            if len(candidate_entering) == 0:
                x_opt = np.zeros(self.n, dtype=float)
                x_opt[self.basis] = x_B
                z_opt = float(np.dot(c_B, x_B))

                if verbose:
                    print(f"-> Óptimo alcanzado en {iteration} iteraciones del Simplex Revisado.")

                return {
                    "status": "Óptimo Encontrado",
                    "optimal_solution": x_opt,
                    "optimal_value": z_opt,
                    "dual_solution": y,
                    "optimal_basis": self.basis.copy(),
                    "iterations": iteration
                }

            # Regla de Dantzig: costo reducido más negativo
            idx_min = np.argmin([self.c[j] - np.dot(y, self.A[:, j]) for j in candidate_entering])
            q = candidate_entering[idx_min]

            # 5. Calcular columna transformada: B d = a_q
            a_q = self.A[:, q]
            d = self._solve_B(P, L, U, a_q)

            # 6. Prueba de la razón mínima
            pos_indices = [i for i, d_i in enumerate(d) if d_i > self.tol]
            if len(pos_indices) == 0:
                return {
                    "status": "No Acotado",
                    "optimal_value": -np.inf,
                    "optimal_solution": None,
                    "iterations": iteration
                }

            ratios = [x_B[i] / d[i] for i in pos_indices]
            min_ratio_idx = np.argmin(ratios)
            p = pos_indices[min_ratio_idx]

            if verbose:
                print(f"Iter {iteration:2d}: z = {np.dot(c_B, x_B):8.3f} | Entra x{q+1} | Sale x{self.basis[p]+1} (Paso = {ratios[min_ratio_idx]:.3f})")

            # 7. Actualización de la base
            self.basis[p] = q
            iteration += 1

        return {
            "status": "Máximo de iteraciones alcanzado",
            "iterations": iteration,
            "optimal_solution": None
        }


if __name__ == "__main__":
    print("=" * 75)
    print("DEMOSTRACIÓN: MÉTODO SIMPLEX REVISADO CON FACTORIZACIÓN LU")
    print("=" * 75)

    # Problema canónico:
    # min  z = -3 x1 - 5 x2
    # s.a.    x1          + s1      = 4
    #              2 x2        + s2 = 12
    #        3 x1 + 2 x2            + s3 = 18
    #        x1, x2, s1, s2, s3 >= 0
    c_rev = np.array([-3.0, -5.0, 0.0, 0.0, 0.0])
    A_rev = np.array([
        [1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0, 1.0, 0.0],
        [3.0, 2.0, 0.0, 0.0, 1.0]
    ])
    b_rev = np.array([4.0, 12.0, 18.0])
    base_inicial = [2, 3, 4]  # Columnas s1, s2, s3 (identidad)

    rev_simplex = SimplexRevisadoLU(c_rev, A_rev, b_rev, base_inicial)
    res_rev = rev_simplex.solve(verbose=True)

    print("\nResultados Simplex Revisado LU:")
    print(f"-> Estado: {res_rev['status']}")
    print(f"-> Solución Primal Óptima x*: {res_rev['optimal_solution']}")
    print(f"-> Solución Dual Óptima y*: {res_rev['dual_solution']}")
    print(f"-> Valor Mínimo z*: {res_rev['optimal_value']:.4f}")
    print(f"-> Base Óptima Final: {res_rev['optimal_basis']}")
    print("=" * 75)
