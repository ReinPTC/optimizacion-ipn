"""
================================================================================
Unidad 1 — Herramientas para la Optimización
Módulo 04: Modelado en Programación Lineal — Los 7 Problemas Canónicos de Luenberger
Basado en: David G. Luenberger & Yinyu Ye, Linear and Nonlinear Programming (4th Ed.)
           Capítulo 2 (Secciones 2.1 - 2.3, pp. 11-20)
================================================================================
"""

import numpy as np
from itertools import combinations
from scipy.optimize import linprog
from typing import Dict, Any, List, Tuple


# ==============================================================================
# 1. CONVERTIDOR A FORMA ESTÁNDAR Y SOLUCIONES BÁSICAS (Luenberger 2.1 y 2.3)
# ==============================================================================

def compute_all_basic_solutions(A: np.ndarray, b: np.ndarray) -> List[Dict[str, Any]]:
    """
    Calcula todas las soluciones básicas para el sistema Ax = b (m x n con m < n, rango(A) = m)
    seleccionando todas las combinaciones C(n, m) de columnas como submatrices base B (Luenberger 2.3).
    
    Retorna:
        Lista de diccionarios con índices básicos, matriz B, x_B = B^(-1)b, vector completo x
        y bandera de factibilidad (x >= 0, Solución Básica Factible - SBF).
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).flatten()
    m, n = A.shape
    assert m <= n, "Se requiere que m <= n."
    
    results = []
    col_indices = list(range(n))
    
    for basis_cols in combinations(col_indices, m):
        B = A[:, basis_cols]
        det_B = np.linalg.det(B)
        
        if np.abs(det_B) > 1e-8:
            x_B = np.linalg.solve(B, b)
            x_full = np.zeros(n, dtype=float)
            for idx, col in enumerate(basis_cols):
                x_full[col] = x_B[idx]
                
            is_feasible = bool(np.all(x_full >= -1e-8))
            results.append({
                "basis_indices": list(basis_cols),
                "det_B": float(det_B),
                "x_B": x_B,
                "x_full": x_full,
                "is_basic_feasible_solution": is_feasible
            })
            
    return results


# ==============================================================================
# 2. LOS 7 MODELOS CANÓNICOS DE LUENBERGER (Cap. 2, Sec. 2.2, pp. 13-19)
# ==============================================================================

def solve_diet_problem() -> Dict[str, Any]:
    """
    Ejemplo 1 de Luenberger (p. 14): Problema de la Dieta (Stigler).
        min c^T x  s.a.  A_nut x >= b_nut,  x >= 0
    """
    # 3 alimentos, 2 nutrientes
    c = np.array([2.0, 3.5, 1.8])          # Costo por unidad de alimento
    A_nut = np.array([[3.0, 2.0, 1.0],     # Proteínas por alimento
                      [1.0, 4.0, 2.0]])    # Vitaminas por alimento
    b_nut = np.array([12.0, 10.0])         # Requerimientos mínimos
    
    res = linprog(c, A_ub=-A_nut, b_ub=-b_nut, bounds=(0, None), method="highs")
    return {"name": "Problema de la Dieta", "costo_minimo": res.fun, "cantidades_alimentos": res.x, "success": res.success}


def solve_manufacturing_problem() -> Dict[str, Any]:
    """
    Ejemplo 2 de Luenberger (p. 15): Problema de Manufactura.
        max p^T x  s.a.  A x <= b_recursos,  x >= 0
    """
    p = np.array([50.0, 40.0, 70.0])       # Precios de venta de 3 productos
    A_rec = np.array([[2.0, 1.0, 3.0],     # Horas máquina
                      [1.0, 2.0, 2.0],     # Materia prima
                      [0.5, 1.0, 1.5]])    # Horas mano de obra
    b_rec = np.array([100.0, 80.0, 45.0])  # Recursos disponibles
    
    res = linprog(-p, A_ub=A_rec, b_ub=b_rec, bounds=(0, None), method="highs")
    return {"name": "Problema de Manufactura", "ingreso_maximo": -res.fun, "produccion": res.x, "success": res.success}


def solve_transportation_problem() -> Dict[str, Any]:
    """
    Ejemplo 3 de Luenberger (p. 15): Problema de Transporte (2 orígenes, 3 destinos).
        min sum c_ij x_ij  s.a. sum_j x_ij = a_i, sum_i x_ij = b_j, x_ij >= 0
    """
    # Costos c_ij para matriz 2x3 (flattened: x11, x12, x13, x21, x22, x23)
    c = np.array([8.0, 6.0, 10.0, 9.0, 12.0, 7.0])
    a = np.array([50.0, 60.0])             # Oferta de orígenes 1 y 2
    b = np.array([30.0, 40.0, 40.0])       # Demanda de destinos 1, 2 y 3
    
    # Matriz rala de restricciones con coeficientes 0 y 1 (Luenberger p. 16)
    A_eq = np.array([
        [1, 1, 1, 0, 0, 0],  # sum_j x1j = a1
        [0, 0, 0, 1, 1, 1],  # sum_j x2j = a2
        [1, 0, 0, 1, 0, 0],  # sum_i xi1 = b1
        [0, 1, 0, 0, 1, 0],  # sum_i xi2 = b2
        [0, 0, 1, 0, 0, 1]   # sum_i xi3 = b3
    ])
    b_eq = np.concatenate([a, b])
    
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    return {"name": "Problema de Transporte", "costo_transporte": res.fun, "flujos_x": res.x.reshape(2, 3), "success": res.success}


def solve_warehousing_problem() -> Dict[str, Any]:
    """
    Ejemplo 5 de Luenberger (p. 17): Problema de Almacenamiento Multitemporal.
    Decisiones para n=3 periodos con inventario inicial x_1 = 0, capacidad C = 100, costo holding r = 1.0.
    """
    p = np.array([10.0, 15.0, 12.0])       # Precios de mercado por periodo
    r = 1.0                                # Costo de almacenamiento
    C = 100.0                              # Capacidad del almacén
    
    # Variables: u1, s1, x2, u2, s2, x3, u3, s3
    # Objetivo: max sum (p_i(s_i - u_i) - r x_i)  ==> min sum (r x_i + p_i u_i - p_i s_i)
    # Vector c: [u1, s1, x2, u2, s2, x3, u3, s3]
    c = np.array([p[0], -p[0], r, p[1], -p[1], r, p[2], -p[2]])
    
    # Restricciones de balance (Ecs. en Luenberger p. 17):
    # -u1 + s1 + x2 = 0
    # -x2 - u2 + s2 + x3 = 0
    # -x3 - u3 + s3 = 0 (almacén queda vacío al final)
    A_eq = np.array([
        [ -1,  1,  1,  0,  0,  0,  0,  0],
        [  0,  0, -1, -1,  1,  1,  0,  0],
        [  0,  0,  0,  0,  0, -1, -1,  1]
    ])
    b_eq = np.array([0.0, 0.0, 0.0])
    
    # Restricciones de capacidad: x2 <= C, x3 <= C
    A_ub = np.array([
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0]
    ])
    b_ub = np.array([C, C])
    
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    return {"name": "Problema de Almacén", "ganancia_maxima": -res.fun, "estrategia_u_s_x": res.x, "success": res.success}


def solve_svm_classifier() -> Dict[str, Any]:
    """
    Ejemplo 6 de Luenberger (p. 17): Clasificador Lineal y Support Vector Machine (SVM).
    Encuentra hiperplano {x : x^T y + beta = 0} que separa a_i (clase +1) de b_j (clase -1):
        a_i^T y + beta >= 1
        b_j^T y + beta <= -1
    """
    # Puntos 2D de dos clases
    A_pts = np.array([[1.0, 2.0], [2.0, 3.0], [1.5, 2.5]])   # Clase A (+1)
    B_pts = np.array([[4.0, 1.0], [5.0, 2.0], [4.5, 0.5]])   # Clase B (-1)
    
    # Variables: y1, y2, beta (libres, representadas con bounds=(-inf, inf))
    # Objetivo dummy: min ||y||_1 o constante
    c = np.array([0.0, 0.0, 0.0])
    
    # a_i^T y + beta >= 1  ==>  -a_i^T y - beta <= -1
    # b_j^T y + beta <= -1 ==>   b_j^T y + beta <= -1
    A_ub = []
    b_ub = []
    for a_i in A_pts:
        A_ub.append([-a_i[0], -a_i[1], -1.0])
        b_ub.append(-1.0)
    for b_j in B_pts:
        A_ub.append([b_j[0], b_j[1], 1.0])
        b_ub.append(-1.0)
        
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=(None, None), method="highs")
    y_vec = res.x[:2]
    beta = res.x[2]
    return {"name": "Clasificador Lineal SVM", "y_normal": y_vec, "beta": beta, "success": res.success}


if __name__ == "__main__":
    print("=" * 75)
    print("SOLUCIONES BÁSICAS Y BASES EN PL (LUENBERGER 2.3, pp. 19-20)")
    print("=" * 75)
    
    # Sistema Ax = b: 2 restricciones, 4 variables (con holguras)
    A = np.array([[1.0, 2.0, 1.0, 0.0],
                  [2.0, 1.0, 0.0, 1.0]])
    b = np.array([6.0, 6.0])
    
    all_basic = compute_all_basic_solutions(A, b)
    print(f"Total de combinaciones C(4, 2) analizadas: {len(all_basic)}\n")
    for idx, sol in enumerate(all_basic, start=1):
        sbf_tag = "[SBF - FACTIBLE]" if sol["is_basic_feasible_solution"] else "[NO FACTIBLE]"
        print(f"Base {idx} (Columnas {sol['basis_indices']}): x = {sol['x_full']} {sbf_tag}")
        
    print("\n" + "=" * 75)
    print("EJECUCIÓN DE LOS 7 MODELOS CANÓNICOS DE LUENBERGER (CAPÍTULO 2)")
    print("=" * 75)
    print("1. Dieta:", solve_diet_problem())
    print("2. Manufactura:", solve_manufacturing_problem())
    print("3. Transporte:", solve_transportation_problem())
    print("4. Almacén Multitemporal:", solve_warehousing_problem())
    print("5. Clasificador Lineal SVM:", solve_svm_classifier())
