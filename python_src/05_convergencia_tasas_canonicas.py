"""
================================================================================
Unidad 1 — Herramientas para la Optimización
Módulo 05: Algoritmos Iterativos y Análisis de Tasas Canónicas de Convergencia
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 1 (Sección 1.4, pp. 6-8)
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Any, List, Tuple


def gradient_descent_quadratic(
    Q: np.ndarray,
    x0: np.ndarray,
    max_iter: int = 50,
    tol: float = 1e-12
) -> Dict[str, Any]:
    """
    Ejecuta el método de descenso de gradiente con paso óptimo exacto sobre
    la forma cuadrática pura f(x) = (1/2) x^T Q x (mínimo global en x* = 0):
    
        x_{k+1} = x_k - alpha_k g_k,  donde g_k = Q x_k, alpha_k = (g_k^T g_k) / (g_k^T Q g_k)
        
    Demuestra la Tasa Canónica de Luenberger (Cap. 1.4 y Cap. 8):
        ||x_{k+1} - x*||_Q <= ((kappa - 1) / (kappa + 1)) * ||x_k - x*||_Q
        donde kappa(Q) = lambda_max / lambda_min es el número de condición de la matriz Q.
    """
    Q = np.asarray(Q, dtype=float)
    x = np.asarray(x0, dtype=float).flatten()
    
    eigenvalues = np.linalg.eigvalsh(Q)
    lambda_min = np.min(eigenvalues)
    lambda_max = np.max(eigenvalues)
    kappa = lambda_max / lambda_min
    canonical_rate = (kappa - 1.0) / (kappa + 1.0)
    
    trajectory = [x.copy()]
    errors_euclidean = [np.linalg.norm(x)]
    errors_Q = [np.sqrt(x.T @ Q @ x)]
    
    for k in range(max_iter):
        g = Q @ x
        if np.linalg.norm(g) < tol:
            break
            
        alpha = np.dot(g, g) / np.dot(g, Q @ g)
        x = x - alpha * g
        
        trajectory.append(x.copy())
        errors_euclidean.append(np.linalg.norm(x))
        errors_Q.append(np.sqrt(np.maximum(x.T @ Q @ x, 0.0)))
        
    return {
        "trajectory": np.array(trajectory),
        "errors_euclidean": np.array(errors_euclidean),
        "errors_Q": np.array(errors_Q),
        "lambda_min": float(lambda_min),
        "lambda_max": float(lambda_max),
        "condition_number_kappa": float(kappa),
        "canonical_rate_beta": float(canonical_rate),
        "iterations": len(trajectory) - 1
    }


def compare_convergence_orders(n_steps: int = 15) -> Dict[str, np.ndarray]:
    """
    Genera secuencias teóricas representativas para contrastar los órdenes de convergencia:
        1. Lineal (tasa beta = 0.5): e_{k+1} = 0.5 * e_k
        2. Superlineal (convergencia acelerada): e_{k+1} = e_k / (k + 1)
        3. Cuadrática (Método de Newton): e_{k+1} = e_k^2
    """
    e_lin = np.zeros(n_steps)
    e_sup = np.zeros(n_steps)
    e_qua = np.zeros(n_steps)
    
    e_lin[0] = 0.8
    e_sup[0] = 0.8
    e_qua[0] = 0.8
    
    for k in range(n_steps - 1):
        e_lin[k + 1] = 0.5 * e_lin[k]
        e_sup[k + 1] = e_sup[k] / (k + 2.0)
        e_qua[k + 1] = min(e_qua[k] ** 2, 1e-16)
        
    return {
        "k": np.arange(n_steps),
        "lineal": e_lin,
        "superlineal": e_sup,
        "cuadratica": e_qua
    }


def plot_convergence_comparison(save_path: str = None):
    """
    Genera el gráfico comparativo de órdenes de convergencia en escala semilogarítmica.
    """
    data = compare_convergence_orders(12)
    k = data["k"]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(k, data["lineal"], "o-", color="#3b82f6", lw=2, label="Convergencia Lineal (Tasa Canónica $\\beta = 0.5$)")
    ax.semilogy(k, data["superlineal"], "s--", color="#10b981", lw=2, label="Convergencia Superlineal (Quasi-Newton)")
    ax.semilogy(k, data["cuadratica"], "^-.", color="#ef4444", lw=2, label="Convergencia Cuadrática (Newton)")
    
    ax.grid(True, which="both", alpha=0.3)
    ax.set_title("Comparación de Órdenes de Convergencia (Luenberger Cap. 1.4)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Iteración $k$")
    ax.set_ylabel("Error asintótico $\\|\\mathbf{x}_k - \\mathbf{x}^*\\|$ (Escala Log)")
    ax.legend(loc="upper right")
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Gráfico guardado en: {save_path}")
    plt.show()


if __name__ == "__main__":
    print("=" * 75)
    print("ANÁLISIS DE TASAS CANÓNICAS DE CONVERGENCIA (LUENBERGER CAP. 1.4)")
    print("=" * 75)
    
    # Matriz con número de condición moderado kappa = 10 / 1 = 10
    Q = np.array([[10.0, 0.0],
                  [0.0, 1.0]])
    x0 = np.array([1.0, 10.0])
    
    res = gradient_descent_quadratic(Q, x0, max_iter=20)
    print(f"Autovalores de Q: lambda_min = {res['lambda_min']}, lambda_max = {res['lambda_max']}")
    print(f"Número de condición kappa(Q)       = {res['condition_number_kappa']:.2f}")
    print(f"Tasa canónica teórica beta = (k-1)/(k+1) = {res['canonical_rate_beta']:.4f}")
    print(f"Iteraciones requeridas para converger: {res['iterations']}\n")
    
    print("Evolución del error en las primeras 5 iteraciones:")
    for k in range(min(6, len(res["errors_Q"]))):
        print(f"  Iteración {k:2d}: ||x_k||_Q = {res['errors_Q'][k]:.6e}")
