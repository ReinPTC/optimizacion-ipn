"""
================================================================================
Unidad 1 — Herramientas para la Optimización
Módulo 02: Eliminación Gaussiana, Factorización LU y Estabilidad Numérica
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Apéndice C (pp. 513-515)
================================================================================
"""

import numpy as np
import time
from typing import Tuple, List, Dict, Any


def lu_decomposition_luenberger(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray]]:
    """
    Implementa la factorización A = L U siguiendo estrictamente la construcción
    teórica de Luenberger (Apéndice C, pp. 513-515) mediante la secuencia de matrices
    elementales de eliminación M_k:
    
        A^(n) = M_{n-1} ... M_1 A = U
        L = M^(-1) = M_1^(-1) ... M_{n-1}^(-1)
        
    Parámetros:
        A: Matriz cuadrada invertible de tamaño (n x n).
        
    Retorna:
        L: Matriz triangular inferior con unos en la diagonal principal.
        U: Matriz triangular superior.
        elementary_matrices: Lista con las matrices M_k generadas en cada paso.
    """
    A_curr = np.array(A, dtype=float)
    n = A_curr.shape[0]
    assert A_curr.shape[0] == A_curr.shape[1], "La matriz A debe ser cuadrada."
    
    L = np.eye(n, dtype=float)
    elementary_matrices = []
    
    for k in range(n - 1):
        pivot = A_curr[k, k]
        if np.abs(pivot) < 1e-14:
            raise ValueError(f"Pivote nulo o casi nulo detectado en fila {k+1}. Se requiere pivoteo.")
            
        Mk = np.eye(n, dtype=float)
        for i in range(k + 1, n):
            multiplier = A_curr[i, k] / pivot
            Mk[i, k] = -multiplier
            L[i, k] = multiplier  # M_k^(-1) tiene signo positivo en los coeficientes fuera de la diagonal
            
        elementary_matrices.append(Mk)
        A_curr = Mk @ A_curr
        
    U = A_curr
    return L, U, elementary_matrices


def forward_substitution(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Resuelve el sistema triangular inferior L y = b mediante sustitución hacia adelante:
        y_1 = b_1 / l_11
        y_i = (b_i - sum_{j=1}^{i-1} l_ij y_j) / l_ii
    """
    n = len(b)
    y = np.zeros(n, dtype=float)
    for i in range(n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
    return y


def back_substitution(U: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Resuelve el sistema triangular superior U x = y mediante sustitución hacia atrás:
        x_n = y_n / u_nn
        x_i = (y_i - sum_{j=i+1}^n u_ij x_j) / u_ii
    """
    n = len(y)
    x = np.zeros(n, dtype=float)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i + 1:], x[i + 1:])) / U[i, i]
    return x


def lu_solve(A: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Resuelve A x = b mediante el proceso bifásico de Luenberger (Ecs. C.1 y C.2):
        1) L y = b  (Forward elimination)
        2) U x = y  (Back substitution)
    """
    L, U, _ = lu_decomposition_luenberger(A)
    y = forward_substitution(L, b)
    x = back_substitution(U, y)
    return x, L, U


def lu_pivoting_scipy_compare(A: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    """
    Compara la solución LU con el solver estándar de NumPy y evalúa el residuo ||Ax - b||.
    """
    t0 = time.perf_counter()
    x_lu, L, U = lu_solve(A, b)
    t_lu = time.perf_counter() - t0
    
    t1 = time.perf_counter()
    x_np = np.linalg.solve(A, b)
    t_np = time.perf_counter() - t1
    
    residual_lu = np.linalg.norm(A @ x_lu - b)
    residual_np = np.linalg.norm(A @ x_np - b)
    
    return {
        "x_lu": x_lu,
        "x_np": x_np,
        "residual_lu": float(residual_lu),
        "residual_np": float(residual_np),
        "time_lu_sec": t_lu,
        "time_np_sec": t_np,
        "L": L,
        "U": U
    }


if __name__ == "__main__":
    print("=" * 75)
    print("FACTORIZACIÓN LU SEGÚN LUENBERGER (APÉNDICE C, pp. 513-515)")
    print("=" * 75)
    
    # Sistema de prueba 3x3
    A = np.array([[2.0, 1.0, 1.0],
                  [4.0, 3.0, 3.0],
                  [8.0, 7.0, 9.0]])
    b = np.array([5.0, 13.0, 37.0])
    
    print("Matriz A:")
    print(A)
    print(f"\nVector b: {b}")
    
    L, U, matrices_M = lu_decomposition_luenberger(A)
    print("\n1. Matrices Elementales de Eliminación M_k:")
    for idx, M_k in enumerate(matrices_M, start=1):
        print(f"M_{idx}:")
        print(M_k)
        
    print("\n2. Matriz Triangular Inferior L:")
    print(L)
    print("\n3. Matriz Triangular Superior U:")
    print(U)
    print("\n4. Comprobación del Producto L * U:")
    print(L @ U)
    print(f"¿L @ U == A?: {np.allclose(L @ U, A)}")
    
    # Solución bifásica
    x, _, _ = lu_solve(A, b)
    print(f"\n5. Solución x calculada mediante L y = b y U x = y:")
    print(f"x = {x}")
    print(f"Comprobación A x = {A @ x} (Esperado: {b})")
    print(f"Residuo ||A x - b|| = {np.linalg.norm(A @ x - b):.2e}")
