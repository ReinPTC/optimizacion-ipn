"""
================================================================================
Unidad 2 — Programación Lineal
Módulo 07: Algoritmo Simplex Primal Tabular y Regla Anticlaje de Bland
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 3: El Método Simplex
           Secciones 3.1 - 3.4 (pp. 33-52)
================================================================================
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional


class SimplexTabular:
    """
    Implementación rigurosa del Algoritmo Simplex Primal en forma tabular
    para problemas de minimización en forma estándar:
        min  z = c^T x
        s.a. A x = b,   x >= 0
    donde b >= 0 y se dispone de una base factible inicial identidad B = I.
    """

    def __init__(
        self,
        c: np.ndarray,
        A: np.ndarray,
        b: np.ndarray,
        basic_vars: Optional[List[int]] = None,
        use_bland_rule: bool = True,
        tol: float = 1e-9
    ):
        self.c = np.asarray(c, dtype=float).flatten()
        self.A = np.asarray(A, dtype=float)
        self.b = np.asarray(b, dtype=float).flatten()
        self.m, self.n = self.A.shape
        self.use_bland = use_bland_rule
        self.tol = tol

        if np.any(self.b < -self.tol):
            raise ValueError("El vector de recursos b debe ser no negativo (b >= 0) para iniciar el Simplex primal.")

        if basic_vars is None:
            # Asume por defecto que las últimas m columnas forman la base canónica (holguras)
            self.basis = list(range(self.n - self.m, self.n))
        else:
            if len(basic_vars) != self.m:
                raise ValueError(f"La base debe contener exactamente {self.m} índices de columnas.")
            self.basis = list(basic_vars)

        # Construcción del Tablero Simplex Canónico de tamaño (m + 1, n + 1)
        # Fila 0: [ r_1, r_2, ..., r_n | -z_0 ]
        # Filas 1..m: [ y_{i1}, y_{i2}, ..., y_{in} | x_{Bi} ]
        self.tableau = np.zeros((self.m + 1, self.n + 1), dtype=float)
        self._initialize_tableau()

    def _initialize_tableau(self) -> None:
        """Inicializa el tablero canónico aplicando eliminación para costos reducidos."""
        B = self.A[:, self.basis]
        inv_B = np.linalg.inv(B)
        
        # Actualización de restricciones: inv(B) * A y inv(B) * b
        y_A = inv_B @ self.A
        x_B = inv_B @ self.b
        
        self.tableau[1:, :self.n] = y_A
        self.tableau[1:, self.n] = x_B

        # Multiplicadores simplex y costos reducidos: r_j = c_j - c_B^T * y_j
        c_B = self.c[self.basis]
        reduced_costs = self.c - (c_B @ y_A)
        current_z = float(np.dot(c_B, x_B))

        self.tableau[0, :self.n] = reduced_costs
        self.tableau[0, self.n] = -current_z

    def print_tableau(self, iteration: int) -> None:
        """Imprime el tablero Simplex de forma formateada y académica."""
        print(f"\n--- Tablero Simplex [Iteración {iteration}] ---")
        header = ["Base"] + [f"x{j+1}" for j in range(self.n)] + ["RHS (b)"]
        print(f"{header[0]:<8} | " + " | ".join(f"{h:>8}" for h in header[1:]))
        print("-" * (10 + 11 * (self.n + 1)))

        # Fila de costos reducidos (Fila 0)
        c_row = ["-z / r"] + [f"{self.tableau[0, j]:8.3f}" for j in range(self.n)] + [f"{-self.tableau[0, self.n]:8.3f}"]
        print(f"{c_row[0]:<8} | " + " | ".join(c_row[1:]))
        print("-" * (10 + 11 * (self.n + 1)))

        # Filas de restricciones
        for i in range(self.m):
            row_label = f"x{self.basis[i]+1}"
            row_vals = [f"{self.tableau[i+1, j]:8.3f}" for j in range(self.n)] + [f"{self.tableau[i+1, self.n]:8.3f}"]
            print(f"{row_label:<8} | " + " | ".join(row_vals))
        print("-" * (10 + 11 * (self.n + 1)))

    def solve(self, max_iter: int = 100) -> Dict[str, Any]:
        """
        Ejecuta el ciclo iterativo del Algoritmo Simplex hasta alcanzar
        optimalidad, no acotamiento o agotar el límite de iteraciones.
        """
        iteration = 0
        history = []

        while iteration < max_iter:
            reduced_costs = self.tableau[0, :self.n]

            # 1. Criterio de Optimalidad (Minimización):
            # Si todos los costos reducidos r_j >= 0, la solución actual es óptima.
            candidate_entering = np.where(reduced_costs < -self.tol)[0]
            if len(candidate_entering) == 0:
                x_opt = np.zeros(self.n, dtype=float)
                for i in range(self.m):
                    x_opt[self.basis[i]] = self.tableau[i + 1, self.n]
                z_opt = -self.tableau[0, self.n]

                return {
                    "status": "Óptimo Encontrado",
                    "optimal_solution": x_opt,
                    "optimal_value": z_opt,
                    "iterations": iteration,
                    "optimal_basis": self.basis.copy(),
                    "history": history
                }

            # 2. Selección de la Variable que Entra a la Base (Columna Pivote q)
            if self.use_bland:
                # Regla de Bland: menor subíndice con costo reducido negativo
                entering_col = candidate_entering[0]
            else:
                # Regla de Dantzig: costo reducido más negativo
                entering_col = candidate_entering[np.argmin(reduced_costs[candidate_entering])]

            # 3. Detección de Problema No Acotado y Prueba de la Razón Mínima (Fila Pivote p)
            y_q = self.tableau[1:, entering_col]
            rhs = self.tableau[1:, self.n]

            positive_pivot_indices = np.where(y_q > self.tol)[0]
            if len(positive_pivot_indices) == 0:
                return {
                    "status": "Problema No Acotado (Costo Infinito a -inf)",
                    "optimal_solution": None,
                    "optimal_value": -np.inf,
                    "iterations": iteration,
                    "entering_col": entering_col,
                    "direction": -y_q
                }

            # Cálculo del paso theta = rhs_i / y_{iq}
            ratios = np.full(self.m, np.inf)
            for idx in positive_pivot_indices:
                ratios[idx] = rhs[idx] / y_q[idx]

            min_ratio = np.min(ratios)
            candidate_leaving = np.where(np.abs(ratios - min_ratio) < self.tol)[0]

            if self.use_bland:
                # Regla de Bland: menor subíndice de la variable básica saliente
                leaving_row = min(candidate_leaving, key=lambda i: self.basis[i])
            else:
                leaving_row = candidate_leaving[0]

            pivot_val = self.tableau[leaving_row + 1, entering_col]
            leaving_var = self.basis[leaving_row]

            history.append({
                "iteration": iteration,
                "current_z": -self.tableau[0, self.n],
                "entering_var": entering_col,
                "leaving_var": leaving_var,
                "min_ratio": min_ratio,
                "pivot_element": pivot_val
            })

            # 4. Operación de Pivoteo de Gauss-Jordan sobre el Tablero
            # Normalizar fila pivote
            self.tableau[leaving_row + 1, :] /= pivot_val

            # Eliminar en todas las demás filas (incluyendo la fila 0 de costos)
            for i in range(self.m + 1):
                if i != leaving_row + 1:
                    factor = self.tableau[i, entering_col]
                    self.tableau[i, :] -= factor * self.tableau[leaving_row + 1, :]

            # Actualizar índice de la base
            self.basis[leaving_row] = entering_col
            iteration += 1

        return {
            "status": "Límite Máximo de Iteraciones Alcanzado",
            "iterations": iteration,
            "optimal_solution": None,
            "optimal_value": None
        }


if __name__ == "__main__":
    print("=" * 75)
    print("MÉTODO SIMPLEX PRIMAL TABULAR (CON REGLA DE BLAND ANTICLAJE)")
    print("=" * 75)

    # Problema canónico:
    # min z = -4 x1 - 3 x2
    # s.a.   2 x1 + 3 x2 + s1          = 6
    #       -3 x1 + 2 x2      + s2     = 3
    #        2 x2                 + s3 = 5
    #        2 x1 +   x2               + s4 = 4
    #        x1, x2, s1, s2, s3, s4 >= 0
    c_prob = np.array([-4.0, -3.0, 0.0, 0.0, 0.0, 0.0])
    A_prob = np.array([
        [ 2.0, 3.0, 1.0, 0.0, 0.0, 0.0],
        [-3.0, 2.0, 0.0, 1.0, 0.0, 0.0],
        [ 0.0, 2.0, 0.0, 0.0, 1.0, 0.0],
        [ 2.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    ])
    b_prob = np.array([6.0, 3.0, 5.0, 4.0])

    simplex = SimplexTabular(c_prob, A_prob, b_prob, use_bland_rule=True)
    print("\nTablero Inicial:")
    simplex.print_tableau(0)

    resultado = simplex.solve(max_iter=10)

    print("\nTablero Final Óptimo:")
    simplex.print_tableau(resultado["iterations"])

    print("\nResultados del Algoritmo:")
    print(f"-> Estado: {resultado['status']}")
    print(f"-> Iteraciones realizadas: {resultado['iterations']}")
    print(f"-> Solución Óptima x*: {resultado['optimal_solution']}")
    print(f"-> Valor Mínimo z*: {resultado['optimal_value']:.4f}")
    print(f"-> Base Óptima Final: {resultado['optimal_basis']}")
    print("=" * 75)
