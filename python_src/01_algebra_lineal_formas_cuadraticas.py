"""
================================================================================
Unidad 1 — Herramientas para la Optimización
Módulo 01: Álgebra Lineal, Espacios Euclidianos y Formas Cuadráticas
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Apéndice A (Secciones A.1 - A.4, pp. 495-499)
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict, Any


def verify_cauchy_schwarz(x: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
    """
    Verifica la Desigualdad de Cauchy-Schwarz en E^n:
        |x^T y| <= ||x|| * ||y||
    
    Parámetros:
        x, y: Vectores en E^n (arreglos 1D de NumPy).
    
    Retorna:
        Diccionario con el producto escalar, normas y la verificación de la desigualdad.
    """
    x = np.asarray(x, dtype=float).flatten()
    y = np.asarray(y, dtype=float).flatten()
    
    inner_product = np.dot(x, y)
    norm_x = np.linalg.norm(x)
    norm_y = np.linalg.norm(y)
    upper_bound = norm_x * norm_y
    
    satisfied = np.abs(inner_product) <= upper_bound + 1e-12
    cos_theta = inner_product / (norm_x * norm_y) if upper_bound > 0 else 1.0
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    
    return {
        "x": x,
        "y": y,
        "inner_product": float(inner_product),
        "norm_x": float(norm_x),
        "norm_y": float(norm_y),
        "upper_bound": float(upper_bound),
        "satisfied": bool(satisfied),
        "angle_degrees": float(angle_deg)
    }


def orthogonal_projection(x: np.ndarray, basis_M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calcula la descomposición ortogonal única de x en E^n sobre un subespacio M y su
    complemento ortogonal M_perp:
        x = a + b,  donde a in M,  b in M_perp (Luenberger Apéndice A.3).
    
    Parámetros:
        x: Vector a proyectar en E^n (tamaño n).
        basis_M: Matriz de tamaño (n, k) cuyas columnas son una base linealmente independiente de M.
        
    Retorna:
        a: Proyección ortogonal de x sobre M.
        b: Proyección ortogonal de x sobre M_perp (b = x - a).
    """
    x = np.asarray(x, dtype=float).flatten()
    V = np.asarray(basis_M, dtype=float)
    if V.ndim == 1:
        V = V[:, np.newaxis]
        
    # Matriz de proyección ortogonal P = V (V^T V)^(-1) V^T
    VtV = V.T @ V
    proj_matrix = V @ np.linalg.inv(VtV) @ V.T
    
    a = proj_matrix @ x
    b = x - a
    
    # Verificación de ortogonalidad: a^T b = 0
    orthogonality_error = np.abs(np.dot(a, b))
    assert orthogonality_error < 1e-10, f"Error de ortogonalidad: {orthogonality_error}"
    
    return a, b


def spectral_decomposition(A: np.ndarray) -> Dict[str, Any]:
    """
    Realiza la descomposición espectral de una matriz simétrica real A:
        Q^T A Q = Lambda  ==>  A = Q Lambda Q^T
    y clasifica su forma cuadrática q(x) = x^T A x (Luenberger Apéndice A.4).
    
    Parámetros:
        A: Matriz cuadrada simétrica (n x n).
        
    Retorna:
        Diccionario con autovalores, matriz ortogonal Q, signatura y raíz cuadrada A^(1/2).
    """
    A = np.asarray(A, dtype=float)
    assert A.shape[0] == A.shape[1], "La matriz debe ser cuadrada."
    assert np.allclose(A, A.T), "La matriz debe ser simétrica (A = A^T)."
    
    # eigh garantiza autovalores reales y base ortonormal para matrices simétricas
    eigenvalues, Q = np.linalg.eigh(A)
    
    # Clasificación de la forma cuadrática
    tol = 1e-10
    if np.all(eigenvalues > tol):
        signature = "Definida Positiva (Positive Definite)"
    elif np.all(eigenvalues >= -tol):
        signature = "Semidefinida Positiva (Positive Semidefinite)"
    elif np.all(eigenvalues < -tol):
        signature = "Definida Negativa (Negative Definite)"
    elif np.all(eigenvalues <= tol):
        signature = "Semidefinida Negativa (Negative Semidefinite)"
    else:
        signature = "Indefinida (Indefinite)"
        
    # Raíz cuadrada simétrica A^(1/2) = Q Lambda^(1/2) Q^T (si es semidefinida positiva)
    if np.all(eigenvalues >= -tol):
        sqrt_lambda = np.sqrt(np.maximum(eigenvalues, 0.0))
        A_sqrt = Q @ np.diag(sqrt_lambda) @ Q.T
    else:
        A_sqrt = None
        
    return {
        "A": A,
        "eigenvalues": eigenvalues,
        "Q": Q,
        "is_orthogonal": bool(np.allclose(Q.T @ Q, np.eye(len(A)))),
        "signature": signature,
        "A_sqrt": A_sqrt
    }


def plot_quadratic_forms_2d(A: np.ndarray, title: str = "Forma Cuadrática x^T A x", save_path: str = None):
    """
    Visualiza las curvas de nivel y la superficie 3D de la forma cuadrática q(x) = x^T A x.
    """
    A = np.asarray(A, dtype=float)
    x1 = np.linspace(-3, 3, 200)
    x2 = np.linspace(-3, 3, 200)
    X1, X2 = np.meshgrid(x1, x2)
    
    # q(x) = a11 x1^2 + 2 a12 x1 x2 + a22 x2^2
    Z = A[0, 0] * X1**2 + (A[0, 1] + A[1, 0]) * X1 * X2 + A[1, 1] * X2**2
    
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    
    # 1. Curvas de nivel
    cs = axes[0].contour(X1, X2, Z, levels=15, cmap="viridis")
    axes[0].clabel(cs, inline=1, fontsize=9)
    axes[0].axhline(0, color="grey", linestyle="--", alpha=0.6)
    axes[0].axvline(0, color="grey", linestyle="--", alpha=0.6)
    axes[0].set_title(f"Curvas de Nivel: {title}")
    axes[0].set_xlabel("$x_1$")
    axes[0].set_ylabel("$x_2$")
    axes[0].grid(True, alpha=0.3)
    axes[0].axis("equal")
    
    # 2. Superficie 3D
    ax3d = fig.add_subplot(1, 2, 2, projection="3d")
    surf = ax3d.plot_surface(X1, X2, Z, cmap="plasma", alpha=0.85, edgecolor="none")
    ax3d.set_title("Superficie $q(x_1, x_2)$")
    ax3d.set_xlabel("$x_1$")
    ax3d.set_ylabel("$x_2$")
    ax3d.set_zlabel("$q(x)$")
    fig.colorbar(surf, ax=ax3d, shrink=0.5, aspect=10)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Gráfico guardado en: {save_path}")
    plt.show()


if __name__ == "__main__":
    print("=" * 70)
    print("EJEMPLO 1: Desigualdad de Cauchy-Schwarz y Proyección Ortogonal")
    print("=" * 70)
    x = np.array([3.0, 1.0, 4.0])
    y = np.array([1.0, 2.0, 2.0])
    res_cs = verify_cauchy_schwarz(x, y)
    print(f"x = {x}, y = {y}")
    print(f"Producto escalar x^T y = {res_cs['inner_product']:.4f}")
    print(f"||x|| ||y|| = {res_cs['upper_bound']:.4f}")
    print(f"Cumple Cauchy-Schwarz: {res_cs['satisfied']} (Ángulo = {res_cs['angle_degrees']:.2f}°)")
    
    # Subespacio M en E^3 generado por dos vectores linealmente independientes
    V = np.array([[1.0, 0.0],
                  [1.0, 1.0],
                  [0.0, 1.0]])
    a, b = orthogonal_projection(x, V)
    print(f"\nProyección de x sobre M (a)      = {np.round(a, 4)}")
    print(f"Componente ortogonal M_perp (b) = {np.round(b, 4)}")
    print(f"Reconstrucción a + b           = {np.round(a + b, 4)} (Original: {x})")
    print(f"Producto ortogonal a^T b       = {np.dot(a, b):.2e}")
    
    print("\n" + "=" * 70)
    print("EJEMPLO 2: Descomposición Espectral y Signatura de Formas Cuadráticas")
    print("=" * 70)
    # Matriz simétrica definida positiva
    A_dp = np.array([[4.0, 2.0],
                     [2.0, 3.0]])
    res_dp = spectral_decomposition(A_dp)
    print("Matriz A (Definida Positiva):")
    print(A_dp)
    print(f"Autovalores: {res_dp['eigenvalues']}")
    print(f"Clasificación: {res_dp['signature']}")
    print(f"Matriz Ortogonal Q (Q^T Q = I): {res_dp['is_orthogonal']}")
    print("Raíz cuadrada simétrica A^(1/2):")
    print(np.round(res_dp['A_sqrt'], 4))
    print("Comprobación A^(1/2) * A^(1/2) = A:")
    print(np.round(res_dp['A_sqrt'] @ res_dp['A_sqrt'], 4))
