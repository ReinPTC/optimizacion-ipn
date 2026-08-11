"""
================================================================================
Unidad 1 — Herramientas para la Optimización
Módulo 03: Geometría de Conjuntos Convexos, Separación y Puntos Extremos
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Apéndice B (Secciones B.1 - B.4, pp. 505-512)
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from typing import Tuple, Dict, Any, List


def separating_hyperplane_point_convex(
    point_y: np.ndarray,
    convex_set_evaluator: callable,
    initial_guess: np.ndarray
) -> Dict[str, Any]:
    """
    Implementa el Teorema 1 de Luenberger (Apéndice B.3, p. 510):
    "Sea C un conjunto convexo y sea y un punto exterior a la clausura de C.
     Entonces existe un vector a tal que a^T y < inf_{x in C} a^T x."
     
    Demostración constructiva de Luenberger:
        1. delta = inf_{x in C} ||x - y|| > 0 alcanzado en x_0 in frontera(C).
        2. a = x_0 - y.
        3. a^T x >= a^T x_0 = a^T y + delta^2  para todo x in C.
        
    Parámetros:
        point_y: Vector exterior y en E^n.
        convex_set_evaluator: Función que define las restricciones g(x) <= 0 del conjunto C.
        initial_guess: Punto inicial dentro de C.
        
    Retorna:
        Diccionario con x_0, delta, vector normal a, constante c y verificación estricta.
    """
    y = np.asarray(point_y, dtype=float).flatten()
    
    # 1. Encontrar la proyección de mínima distancia x_0 = argmin_{x in C} ||x - y||^2
    def distance_sq(x):
        return np.sum((x - y) ** 2)
        
    constraints = {"type": "ineq", "fun": lambda x: -convex_set_evaluator(x)}
    res = minimize(distance_sq, initial_guess, constraints=constraints, method="SLSQP")
    
    if not res.success:
        raise RuntimeError(f"Error al proyectar sobre el conjunto convexo: {res.message}")
        
    x_0 = res.x
    delta = np.linalg.norm(x_0 - y)
    
    # 2. Vector normal del hiperplano separador (Luenberger: a = x_0 - y)
    a = x_0 - y
    
    # 3. Constante del hiperplano de soporte en x_0: c_0 = a^T x_0
    c_0 = np.dot(a, x_0)
    c_y = np.dot(a, y)
    
    # Hiperplano medio de separación: c_mid = (c_0 + c_y) / 2
    c_mid = (c_0 + c_y) / 2.0
    
    return {
        "y": y,
        "x_0": x_0,
        "delta": float(delta),
        "a": a,
        "c_0": float(c_0),
        "c_y": float(c_y),
        "c_mid": float(c_mid),
        "strictly_separated": bool(c_y < c_0)
    }


def find_polytope_extreme_points_2d(A_ineq: np.ndarray, b_ineq: np.ndarray) -> np.ndarray:
    """
    Calcula los puntos extremos de un politopo 2D definido por A_ineq x <= b_ineq (Teorema 5).
    """
    A_ineq = np.asarray(A_ineq, dtype=float)
    b_ineq = np.asarray(b_ineq, dtype=float).flatten()
    m = len(b_ineq)
    
    intersections = []
    # Intersección de pares de rectas a_i^T x = b_i
    for i in range(m):
        for j in range(i + 1, m):
            A_pair = np.vstack([A_ineq[i], A_ineq[j]])
            b_pair = np.array([b_ineq[i], b_ineq[j]])
            
            if np.abs(np.linalg.det(A_pair)) > 1e-8:
                pt = np.linalg.solve(A_pair, b_pair)
                # Verificar factibilidad: A_ineq pt <= b_ineq
                if np.all(A_ineq @ pt <= b_ineq + 1e-7):
                    intersections.append(pt)
                    
    if not intersections:
        return np.empty((0, 2))
        
    pts = np.unique(np.round(np.array(intersections), 6), axis=0)
    
    # Ordenar vértices en sentido antihorario mediante ConvexHull
    if len(pts) >= 3:
        hull = ConvexHull(pts)
        pts = pts[hull.vertices]
        
    return pts


def plot_separation_luenberger(save_path: str = None):
    """
    Genera el gráfico interactivo del Teorema 1 de Separación de Luenberger.
    """
    # Conjunto convexo C: Elipsoide (x1/2)^2 + x2^2 <= 1
    def ellipsoid_constraint(x):
        return (x[0] / 2.0)**2 + x[1]**2 - 1.0

    y_ext = np.array([2.5, 1.8])
    sep = separating_hyperplane_point_convex(y_ext, ellipsoid_constraint, initial_guess=np.array([0.0, 0.0]))
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Dibujar el conjunto convexo C
    theta = np.linspace(0, 2 * np.pi, 200)
    c_x1 = 2.0 * np.cos(theta)
    c_x2 = 1.0 * np.sin(theta)
    ax.fill(c_x1, c_x2, color="#3b82f6", alpha=0.3, label="Conjunto Convexo $C$")
    ax.plot(c_x1, c_x2, color="#1d4ed8", lw=2)
    
    # Punto exterior y, proyección x_0
    y = sep["y"]
    x0 = sep["x_0"]
    a = sep["a"]
    c_mid = sep["c_mid"]
    
    ax.plot(y[0], y[1], "ro", markersize=8, label="Punto exterior $\\mathbf{y}$")
    ax.plot(x0[0], x0[1], "go", markersize=8, label="Proyección $\\mathbf{x}_0$ (Mín. Distancia)")
    
    # Segmento de distancia mínima y vector a = x0 - y
    ax.annotate("", xy=x0, xytext=y, arrowprops=dict(arrowstyle="<-", color="crimson", lw=2, linestyle="--"))
    ax.text((y[0] + x0[0]) / 2 + 0.1, (y[1] + x0[1]) / 2, f"$\\delta = {sep['delta']:.2f}$", color="crimson", fontsize=11)
    
    # Recta del hiperplano separador: a1 x1 + a2 x2 = c_mid
    t_line = np.linspace(-1, 4, 100)
    if np.abs(a[1]) > 1e-5:
        line_x2 = (c_mid - a[0] * t_line) / a[1]
        ax.plot(t_line, line_x2, "k--", lw=2, label="Hiperplano Separador $H = \\{\\mathbf{x} : \\mathbf{a}^T\\mathbf{x} = c\\}$")
    
    ax.set_xlim(-3, 4)
    ax.set_ylim(-2, 3)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_title("Teorema 1 de Luenberger: Hiperplano Separador de Punto Exterior", fontsize=12, fontweight="bold")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.legend(loc="upper left")
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Gráfico guardado en {save_path}")
    plt.show()


if __name__ == "__main__":
    print("=" * 75)
    print("TEOREMA 1 DE LUENBERGER: SEPARACIÓN DE PUNTO EXTERIOR (APÉNDICE B.3)")
    print("=" * 75)
    
    def ellipsoid(x):
        return (x[0] / 2.0)**2 + x[1]**2 - 1.0
        
    y_ext = np.array([2.5, 1.8])
    res_sep = separating_hyperplane_point_convex(y_ext, ellipsoid, initial_guess=np.array([0.0, 0.0]))
    
    print(f"Punto exterior y                     = {res_sep['y']}")
    print(f"Punto frontera x_0 (Mínima distancia) = {np.round(res_sep['x_0'], 4)}")
    print(f"Distancia mínima delta = ||x_0 - y|| = {res_sep['delta']:.4f}")
    print(f"Vector normal a = x_0 - y            = {np.round(res_sep['a'], 4)}")
    print(f"Valor a^T y                          = {res_sep['c_y']:.4f}")
    print(f"Valor inf_{{x in C}} a^T x = a^T x_0   = {res_sep['c_0']:.4f}")
    print(f"¿Separación estricta a^T y < a^T x_0?: {res_sep['strictly_separated']}")
    
    print("\n" + "=" * 75)
    print("TEOREMAS 4 Y 5: PUNTOS EXTREMOS Y REPRESENTACIÓN DE POLÍTOPOS")
    print("=" * 75)
    # Politopo: x1 + 2*x2 <= 6, 2*x1 + x2 <= 6, x1 >= 0, x2 >= 0
    A_ineq = np.array([[1.0, 2.0],
                       [2.0, 1.0],
                       [-1.0, 0.0],
                       [0.0, -1.0]])
    b_ineq = np.array([6.0, 6.0, 0.0, 0.0])
    
    vertices = find_polytope_extreme_points_2d(A_ineq, b_ineq)
    print("Puntos extremos (vértices) del poliedro calculados:")
    for idx, v in enumerate(vertices, start=1):
        print(f"  Vértice {idx}: {v}")
