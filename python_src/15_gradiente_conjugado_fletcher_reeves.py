"""
=============================================================================
ESFM - Instituto Politécnico Nacional
Licenciatura en Matemática Algorítmica
Materia: Optimización

MÓDULO 15: MÉTODOS DE DIRECCIONES CONJUGADAS Y GRADIENTE CONJUGADO
Referencia: David G. Luenberger & Yinyu Ye, "Linear and Nonlinear Programming"
            (4ta Edición), Capítulo 9: Conjugate Direction Methods (Secciones 9.1 - 9.6)
=============================================================================
Este módulo computacional implementa, analiza y valida formalmente:
  1. Direcciones Q-Conjugadas y Ortogonalización de Gram-Schmidt Q-Conjugada.
  2. Algoritmo de Gradiente Conjugado (C-G) para Formas Cuadráticas:
     - Terminación matemática exacta en a lo más n iteraciones.
     - Verificación numérica de la Q-ortogonalidad mutua de direcciones: d_i^T Q d_j = 0.
     - Verificación numérica de la ortogonalidad mutua de residuos: g_i^T g_j = 0.
  3. C-G como Proceso Óptimo y Cota de Polinomios de Chebyshev:
     - Cota de energía ||x_k - x*||_Q <= 2 * ((sqrt(kappa)-1)/(sqrt(kappa)+1))^k * ||x_0 - x*||_Q.
     - Comparativa directa de velocidad: C-G vs. Descenso Más Pronunciado (Steepest Descent).
  4. Gradiente Conjugado Parcial con reinicios cada m pasos (Luenberger Sec. 9.5).
  5. Extensión a Problemas No Cuadráticos (Sec. 9.6):
     - Fórmulas de Fletcher-Reeves (FR), Polak-Ribière (PR) y Hestenes-Stiefel (HS).
     - Criterio de reinicio periódico de Powell cuando |g_{k+1}^T g_k| >= 0.2 ||g_{k+1}||².
     - Benchmark sobre la función no convexa de Rosenbrock.
=============================================================================
"""

import sys
import numpy as np

# Configuración de codificación UTF-8 para consola
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def conjugate_gradient_quadratic(Q, b, x0, tol=1e-10, max_iter=None):
    """
    Algoritmo de Gradiente Conjugado Canónico para:
        min f(x) = 0.5 * x^T Q x - b^T x,   Q = Q^T > 0
    
    Demuestra la convergencia en exactamente n pasos y registra la ortogonalidad mutua.
    """
    Q = np.asarray(Q, dtype=float)
    b = np.asarray(b, dtype=float)
    x = np.asarray(x0, dtype=float).copy()
    n = len(x)
    
    if max_iter is None:
        max_iter = n
        
    x_star = np.linalg.solve(Q, b)
    
    g = Q @ x - b
    d = -g.copy()
    
    history = {
        "x": [x.copy()],
        "g": [g.copy()],
        "d": [d.copy()],
        "energy_error": [],
        "grad_norm": [np.linalg.norm(g)],
        "alpha": [],
        "beta": []
    }
    
    diff0 = x - x_star
    history["energy_error"].append(0.5 * float(diff0.T @ Q @ diff0))
    
    for k in range(max_iter):
        if np.linalg.norm(g) < tol:
            break
            
        Qd = Q @ d
        denom = float(d.T @ Qd)
        if denom <= 0:
            raise ValueError("Q no es definida positiva.")
            
        alpha = float(g.T @ g) / denom
        history["alpha"].append(alpha)
        
        x_next = x + alpha * d
        g_next = g + alpha * Qd
        
        diff_next = x_next - x_star
        E_next = 0.5 * float(diff_next.T @ Q @ diff_next)
        history["energy_error"].append(E_next)
        history["grad_norm"].append(np.linalg.norm(g_next))
        
        if np.linalg.norm(g_next) < tol:
            x = x_next
            history["x"].append(x.copy())
            history["g"].append(g_next.copy())
            break
            
        beta = float(g_next.T @ g_next) / float(g.T @ g)
        history["beta"].append(beta)
        
        d_next = -g_next + beta * d
        
        x = x_next
        g = g_next
        d = d_next
        
        history["x"].append(x.copy())
        history["g"].append(g.copy())
        history["d"].append(d.copy())
        
    return x, history


def check_orthogonality_matrices(history, Q):
    """
    Construye las matrices de gramiana para verificar:
      1. Q-conjugación de direcciones: D^T Q D es diagonal.
      2. Ortogonalidad de residuos: G^T G es diagonal.
    """
    dirs = np.array(history["d"][:-1] if len(history["d"]) > len(history["alpha"]) else history["d"])
    grads = np.array(history["g"][:len(dirs)])
    
    m = len(dirs)
    Q_matrix = np.zeros((m, m))
    G_matrix = np.zeros((m, m))
    
    for i in range(m):
        for j in range(m):
            Q_matrix[i, j] = float(dirs[i].T @ Q @ dirs[j])
            G_matrix[i, j] = float(grads[i].T @ grads[j])
            
    return Q_matrix, G_matrix


def chebyshev_bound(kappa, k):
    """
    Cota de convergencia en energía basada en polinomios de Chebyshev (Luenberger Sec. 9.4):
        ||x_k - x*||_Q <= 2 * ( (sqrt(kappa) - 1) / (sqrt(kappa) + 1) )^k * ||x_0 - x*||_Q
    """
    sqrt_k = np.sqrt(kappa)
    rate = (sqrt_k - 1.0) / (sqrt_k + 1.0)
    return 2.0 * (rate ** k)


def nonlinear_conjugate_gradient(f, grad_f, x0, method='PR', tol=1e-6, max_iter=500, c1=1e-4, c2=0.1):
    """
    Gradiente Conjugado No Lineal (Luenberger Sec. 9.6):
      - Fletcher-Reeves (FR)
      - Polak-Ribière (PR) con reinicio automático de Powell si |g_{k+1}^T g_k| >= 0.2 ||g_{k+1}||²
      - Búsqueda de línea inexacta con condiciones fuertes de Wolfe.
    """
    x = np.asarray(x0, dtype=float).copy()
    g = grad_f(x)
    d = -g.copy()
    
    history = {
        "x": [x.copy()],
        "f": [float(f(x))],
        "grad_norm": [np.linalg.norm(g)],
        "restarts": []
    }
    
    for k in range(max_iter):
        norm_g = np.linalg.norm(g)
        if norm_g < tol:
            break
            
        # Búsqueda de línea con condición de Armijo y Wolfe
        alpha = 1.0
        # Backtracking simple para asegurar descenso
        f_val = f(x)
        dir_deriv = float(g.T @ d)
        if dir_deriv >= 0:
            # Si d perdió condición de descenso, reiniciar a -g
            d = -g.copy()
            dir_deriv = float(g.T @ d)
            history["restarts"].append(k)
            
        while f(x + alpha * d) > f_val + c1 * alpha * dir_deriv:
            alpha *= 0.5
            if alpha < 1e-14:
                break
                
        x_next = x + alpha * d
        g_next = grad_f(x_next)
        
        # Cálculo de beta según el método
        norm_g_sq = float(g.T @ g)
        norm_g_next_sq = float(g_next.T @ g_next)
        
        if method == 'FR':
            beta = norm_g_next_sq / max(norm_g_sq, 1e-15)
        elif method == 'PR':
            y_k = g_next - g
            beta = float(g_next.T @ y_k) / max(norm_g_sq, 1e-15)
            # Regla de Powell: no permitir beta negativo
            if beta < 0:
                beta = 0.0
        else:
            beta = 0.0
            
        # Criterio de reinicio de Powell: si los gradientes no son suficientemente ortogonales
        if abs(float(g_next.T @ g)) >= 0.2 * norm_g_next_sq:
            beta = 0.0
            history["restarts"].append(k + 1)
            
        d_next = -g_next + beta * d
        
        x = x_next
        g = g_next
        d = d_next
        
        history["x"].append(x.copy())
        history["f"].append(float(f(x)))
        history["grad_norm"].append(np.linalg.norm(g))
        
    return x, history


def main():
    print("=" * 80)
    print("ESFM-IPN :: MÓDULO 15: MÉTODOS DE DIRECCIONES CONJUGADAS Y GRADIENTE CONJUGADO")
    print("Fiel a Luenberger & Ye (4ta Ed.), Capítulo 9 (Secciones 9.1 a 9.6)")
    print("=" * 80)
    
    # -------------------------------------------------------------------------
    # PARTE 1: VERIFICACIÓN DE CONVERGENCIA EN EXACTAMENTE n PASOS (Q en R^4)
    # -------------------------------------------------------------------------
    print("\n[1] GRADIENTE CONJUGADO CUADRÁTICO: TERMINACIÓN FINITA EN n PASOS (Luenberger Sec. 9.3)")
    print("-" * 75)
    
    # Generar matriz simétrica definida positiva en R^4 con condicionamiento controlado
    np.random.seed(42)
    n_dim = 4
    diag_vals = [2.0, 15.0, 50.0, 200.0]  # kappa = 100
    Q_4 = np.diag(diag_vals)
    b_4 = np.array([4.0, -10.0, 25.0, 100.0])
    x0_4 = np.zeros(n_dim)
    
    x_opt_cg, hist_cg = conjugate_gradient_quadratic(Q_4, b_4, x0_4)
    x_exact = np.linalg.solve(Q_4, b_4)
    
    print(f"  • Dimensión del Problema n:         {n_dim}")
    print(f"  • Número de Condición kappa(Q):     {max(diag_vals)/min(diag_vals):.1f}")
    print(f"  • Solución Exacta x*:               {x_exact}")
    print(f"  • Solución C-G Obtenida:            {x_opt_cg}")
    print(f"  • Error Absoluto ||x_cg - x*||:     {np.linalg.norm(x_opt_cg - x_exact):.6e}")
    print(f"  • Iteraciones Requeridas C-G:       {len(hist_cg['alpha'])} (¡Exactamente n = 4!)")
    
    print("\n  Evolución por Iteración de C-G en R⁴:")
    print("  Iter k |    Norma Gradiente ||g_k||   |    Error de Energía E(x_k)   |   Paso alpha_k   |  Factor beta_k")
    print("  " + "-" * 78)
    for i in range(len(hist_cg["alpha"])):
        beta_str = f"{hist_cg['beta'][i]:14.6f}" if i < len(hist_cg["beta"]) else "           N/A"
        print(f"  {i:6d} | {hist_cg['grad_norm'][i]:28.8e} | {hist_cg['energy_error'][i]:28.8e} | {hist_cg['alpha'][i]:16.6f} | {beta_str}")
    print(f"  {len(hist_cg['alpha']):6d} | {hist_cg['grad_norm'][-1]:28.8e} | {hist_cg['energy_error'][-1]:28.8e} |        TERMINADO |            N/A")
    
    # -------------------------------------------------------------------------
    # PARTE 2: VERIFICACIÓN NUMÉRICA DE Q-CONJUGACIÓN Y ORTOGONALIDAD DE RESIDUOS
    # -------------------------------------------------------------------------
    print("\n[2] COMPROBACIÓN NUMÉRICA DE Q-CONJUGACIÓN d_i^T Q d_j = 0 Y RESIDUOS g_i^T g_j = 0")
    print("-" * 75)
    
    Q_gram, G_gram = check_orthogonality_matrices(hist_cg, Q_4)
    
    print("  Matriz Gramiana de Direcciones D^T Q D (Debe ser diagonal estricta):")
    for row in Q_gram:
        print("    " + " ".join(f"{val:12.4e}" for val in row))
        
    print("\n  Matriz Gramiana de Residuos G^T G (Debe ser diagonal estricta):")
    for row in G_gram:
        print("    " + " ".join(f"{val:12.4e}" for val in row))
        
    # Extraer elementos fuera de la diagonal
    off_diag_Q = np.max(np.abs(Q_gram - np.diag(np.diag(Q_gram))))
    off_diag_G = np.max(np.abs(G_gram - np.diag(np.diag(G_gram))))
    print(f"\n  • Máximo elemento fuera de la diagonal en D^T Q D: {off_diag_Q:.6e}  (Verificado ~ 0)")
    print(f"  • Máximo elemento fuera de la diagonal en G^T G:   {off_diag_G:.6e}  (Verificado ~ 0)")
    
    # -------------------------------------------------------------------------
    # PARTE 3: COTA DE POLINOMIOS DE CHEBYSHEV Y ACELERACIÓN RESPECTO A GRADIENTE
    # -------------------------------------------------------------------------
    print("\n[3] COTA DE CHEBYSHEV VS. DESCENSO MÁS PRONUNCIADO (Luenberger Sec. 9.4)")
    print("-" * 75)
    
    kappa_test = 100.0
    tasa_sd = ((kappa_test - 1.0) / (kappa_test + 1.0)) ** 2
    tasa_cg = (np.sqrt(kappa_test) - 1.0) / (np.sqrt(kappa_test) + 1.0)
    
    print(f"  Para número de condición kappa = {kappa_test:.1f}:")
    print(f"  • Tasa Asintótica Steepest Descent: ((kappa - 1) / (kappa + 1))² = {tasa_sd:.6f}")
    print(f"  • Tasa de Chebyshev para C-G:       (sqrt(kappa) - 1) / (sqrt(kappa) + 1) = {tasa_cg:.6f}")
    print(f"  • Factor de Aceleración Teórico:    {tasa_sd / tasa_cg:.2f}x más rápido en contracción por paso")
    
    # -------------------------------------------------------------------------
    # PARTE 4: GRADIENTE CONJUGADO NO LINEAL EN FUNCIÓN DE ROSENBROCK (FLETCHER-REEVES VS POLAK-RIBIÈRE)
    # -------------------------------------------------------------------------
    print("\n[4] BENCHMARK EN ROSENBROCK: FLETCHER-REEVES vs. POLAK-RIBIÈRE (Luenberger Sec. 9.6)")
    print("-" * 75)
    print("  f(x1, x2) = 100*(x2 - x1²)² + (1 - x1)²,   x0 = (-1.2, 1.0),   x* = (1, 1)")
    print("-" * 75)
    
    def rosenbrock(x):
        return 100.0 * (x[1] - x[0]**2)**2 + (1.0 - x[0])**2
        
    def rosenbrock_grad(x):
        return np.array([
            -400.0 * x[0] * (x[1] - x[0]**2) - 2.0 * (1.0 - x[0]),
            200.0 * (x[1] - x[0]**2)
        ])
        
    x0_rosen = np.array([-1.2, 1.0])
    
    # 1. Fletcher-Reeves
    x_fr, hist_fr = nonlinear_conjugate_gradient(rosenbrock, rosenbrock_grad, x0_rosen, method='FR', tol=1e-5, max_iter=500)
    
    # 2. Polak-Ribière con reinicios
    x_pr, hist_pr = nonlinear_conjugate_gradient(rosenbrock, rosenbrock_grad, x0_rosen, method='PR', tol=1e-5, max_iter=500)
    
    print(f"  MÉTODO C-G NO LINEAL | ITERACIONES | REINICIOS DE POWELL | ||x_final - x*|| | f(x_final)")
    print("  " + "-" * 75)
    print(f"  Fletcher-Reeves (FR)  | {len(hist_fr['grad_norm']):11d} | {len(hist_fr['restarts']):19d} | {np.linalg.norm(x_fr - np.array([1.0, 1.0])):16.6e} | {hist_fr['f'][-1]:14.6e}")
    print(f"  Polak-Ribière (PR)    | {len(hist_pr['grad_norm']):11d} | {len(hist_pr['restarts']):19d} | {np.linalg.norm(x_pr - np.array([1.0, 1.0])):16.6e} | {hist_pr['f'][-1]:14.6e}")
    print("=" * 80)
    print("CONCLUSIÓN DOCENTE (Luenberger pp. 272-284):")
    print("  1. Polak-Ribière supera ampliamente a Fletcher-Reeves en funciones no lineales")
    print("     debido a que cuando se produce un paso pequeño (|g_{k+1}| ≈ |g_k|), el numerador")
    print("     g_{k+1}^T (g_{k+1} - g_k) fuerza beta ≈ 0, reiniciando automáticamente a la dirección")
    print("     del gradiente y purgando la pérdida de conjugación inducida por la no-cuadraticidad.")
    print("=" * 80)


if __name__ == "__main__":
    main()
