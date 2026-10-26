"""
=============================================================================
ESFM - Instituto Politécnico Nacional
Licenciatura en Matemática Algorítmica
Materia: Optimización

MÓDULO 14: MÉTODO DE DESCENSO MÁS PRONUNCIADO Y MÉTODO DE NEWTON MULTIDIMENSIONAL
Referencia: David G. Luenberger & Yinyu Ye, "Linear and Nonlinear Programming"
            (4ta Edición), Capítulo 8: Basic Descent Methods (Secciones 8.2 y 8.5)
=============================================================================
Este módulo computacional implementa, analiza y valida formalmente:
  1. Método de Descenso Más Pronunciado (Steepest Descent) para formas cuadráticas:
     - Cálculo analítico del paso óptimo alpha_k = (g_k^T g_k) / (g_k^T Q g_k).
     - Comprobación matemática de la ortogonalidad sucesiva: g_{k+1}^T g_k = 0.
     - Trayectoria oscilatoria en zig-zag.
  2. Teorema de la Desigualdad de Kantórovich y Tasa de Convergencia:
     - Demostración numérica de la cota inferior 4*a*A / (a + A)².
     - Identificación del vector extremal que alcanza la cota exacta (combinación
       equitativa de autovectores extremos v_min y v_max).
     - Validación empírica de la tasa asintótica: E(x_{k+1}) <= ((kappa-1)/(kappa+1))² E(x_k).
  3. Método de Newton Multidimensional:
     - Paso de Newton puro y amortiguado (Damped Newton con búsqueda de Armijo).
     - Demostración experimental del orden de convergencia cuadrático (p = 2).
  4. Comparativa Experimental Exhaustiva (Steepest Descent vs. Newton):
     - Prueba sobre elipsoide mal condicionado (kappa = 100).
     - Prueba sobre la función no convexa de Rosenbrock.
=============================================================================
"""

import sys
import numpy as np

# Configuración de codificación UTF-8 para consola
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def steepest_descent_quadratic(Q, b, x0, tol=1e-8, max_iter=2000):
    """
    Método de Descenso Más Pronunciado con búsqueda de línea exacta para:
        min f(x) = 0.5 * x^T Q x - b^T x
    
    Parámetros:
        Q: Matriz hessiana simétrica y definida positiva (n x n).
        b: Vector término lineal (n,).
        x0: Vector de inicio (n,).
        tol: Tolerancia en la norma del gradiente ||g_k||.
        max_iter: Límite máximo de iteraciones.
    
    Retorna:
        historial: Diccionario con trayectoria, errores de energía, tasas y ortogonalidad.
    """
    Q = np.asarray(Q, dtype=float)
    b = np.asarray(b, dtype=float)
    x = np.asarray(x0, dtype=float).copy()
    
    # Solución analítica exacta x* = Q^(-1) b
    x_star = np.linalg.solve(Q, b)
    f_star = 0.5 * float(x_star.T @ Q @ x_star) - float(b.T @ x_star)
    
    # Autovalores y número de condición
    eigvals = np.linalg.eigvalsh(Q)
    a, A = float(np.min(eigvals)), float(np.max(eigvals))
    kappa = A / a
    tasa_kantorovich = ((kappa - 1.0) / (kappa + 1.0)) ** 2
    
    history = {
        "x": [x.copy()],
        "grad_norm": [],
        "energy_error": [],
        "step_sizes": [],
        "orthogonalities": [],
        "observed_rates": [],
        "kantorovich_bound": tasa_kantorovich,
        "kappa": kappa
    }
    
    for k in range(max_iter):
        g = Q @ x - b
        norm_g = np.linalg.norm(g)
        
        # Error de energía E(x_k) = 0.5 * (x - x*)^T Q (x - x*)
        diff = x - x_star
        E_k = 0.5 * float(diff.T @ Q @ diff)
        
        history["grad_norm"].append(norm_g)
        history["energy_error"].append(E_k)
        
        if norm_g < tol:
            break
            
        # Paso óptimo analítico en dirección d_k = -g_k
        denom = float(g.T @ Q @ g)
        if denom <= 0:
            raise ValueError("La matriz Q no es definida positiva.")
        alpha = float(g.T @ g) / denom
        history["step_sizes"].append(alpha)
        
        # Siguiente iteración
        x_next = x - alpha * g
        g_next = Q @ x_next - b
        
        # Verificación de ortogonalidad sucesiva g_{k+1}^T g_k
        prod_ortogonal = float(g_next.T @ g)
        history["orthogonalities"].append(prod_ortogonal)
        
        # Tasa observada de error de energía E_{k+1} / E_k
        diff_next = x_next - x_star
        E_next = 0.5 * float(diff_next.T @ Q @ diff_next)
        if E_k > 1e-15:
            history["observed_rates"].append(E_next / E_k)
        else:
            history["observed_rates"].append(0.0)
            
        x = x_next
        history["x"].append(x.copy())
        
    return x, history


def verify_kantorovich_inequality(Q, num_random_samples=500):
    """
    Verifica numéricamente la Desigualdad de Kantórovich (Luenberger Secc. 8.2):
        [ (x^T x)^2 ] / [ (x^T Q x) (x^T Q^(-1) x) ] >= 4*a*A / (a + A)^2
    
    Comprueba tanto muestras aleatorias como la combinación extremal analítica
    que alcanza la cota inferior con igualdad exacta.
    """
    eigvals, eigvecs = np.linalg.eigh(Q)
    idx_min = np.argmin(eigvals)
    idx_max = np.argmax(eigvals)
    
    a = float(eigvals[idx_min])
    A = float(eigvals[idx_max])
    v_a = eigvecs[:, idx_min]
    v_A = eigvecs[:, idx_max]
    
    # Cota inferior analítica de Kantórovich
    cota_teorica = 4.0 * a * A / ((a + A) ** 2)
    
    Q_inv = np.linalg.inv(Q)
    
    def kantorovich_ratio(v):
        v = np.asarray(v, dtype=float)
        v_norm2 = float(v.T @ v)
        q_val = float(v.T @ Q @ v)
        q_inv_val = float(v.T @ Q_inv @ v)
        return (v_norm2 ** 2) / (q_val * q_inv_val)
    
    # 1. Vector Extremal Canónico: combinación equitativa de los dos autovectores extremos
    v_extremal = (v_a / np.linalg.norm(v_a)) + (v_A / np.linalg.norm(v_A))
    ratio_extremal = kantorovich_ratio(v_extremal)
    
    # 2. Muestreo de vectores aleatorios en R^n
    ratios_aleatorios = []
    np.random.seed(42)
    n = Q.shape[0]
    for _ in range(num_random_samples):
        vec = np.random.randn(n)
        ratios_aleatorios.append(kantorovich_ratio(vec))
        
    min_aleatorio = min(ratios_aleatorios)
    
    return {
        "a": a,
        "A": A,
        "cota_teorica": cota_teorica,
        "ratio_extremal": ratio_extremal,
        "min_aleatorio": min_aleatorio,
        "satisface_desigualdad": (min_aleatorio >= cota_teorica - 1e-12) and np.isclose(ratio_extremal, cota_teorica, atol=1e-10)
    }


def steepest_descent_armijo(f, grad_f, x0, tol=1e-6, max_iter=2000, c1=1e-4, beta=0.5):
    """
    Método de Descenso Más Pronunciado con regla de Armijo y Backtracking
    para funciones no lineales generales en R^n.
    """
    x = np.asarray(x0, dtype=float).copy()
    history = {"x": [x.copy()], "f": [float(f(x))], "grad_norm": []}
    
    for k in range(max_iter):
        g = grad_f(x)
        norm_g = np.linalg.norm(g)
        history["grad_norm"].append(norm_g)
        
        if norm_g < tol:
            break
            
        d = -g  # Dirección de descenso
        alpha = 1.0  # Paso inicial
        f_val = f(x)
        
        # Búsqueda de Armijo (Backtracking)
        while f(x + alpha * d) > f_val + c1 * alpha * float(g.T @ d):
            alpha *= beta
            if alpha < 1e-16:
                break
                
        x = x + alpha * d
        history["x"].append(x.copy())
        history["f"].append(float(f(x)))
        
    return x, history


def newton_multidimensional(f, grad_f, hess_f, x0, tol=1e-8, max_iter=100, damping=True, c1=1e-4, beta=0.5):
    """
    Método de Newton Multidimensional (Luenberger Secc. 8.5):
        Paso de Newton: d_k = - [F(x_k)]^(-1) * grad_f(x_k)
    
    Parámetros:
        damping: Si es True, aplica búsqueda de Armijo a lo largo de d_k.
                 Si es False, aplica el paso puro de Newton con alpha_k = 1.0.
    """
    x = np.asarray(x0, dtype=float).copy()
    history = {
        "x": [x.copy()],
        "f": [float(f(x))],
        "grad_norm": [],
        "step_lengths": [],
        "damping_alphas": []
    }
    
    for k in range(max_iter):
        g = grad_f(x)
        norm_g = np.linalg.norm(g)
        history["grad_norm"].append(norm_g)
        
        if norm_g < tol:
            break
            
        H = hess_f(x)
        
        # Resolver sistema lineal H * d = -g para evitar invertir explícitamente
        try:
            d = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            # Si el hessiano es singular, se aplica corrección mínima de Levenberg
            d = np.linalg.solve(H + 1e-4 * np.eye(len(x)), -g)
            
        alpha = 1.0
        if damping:
            f_val = f(x)
            # Verificar si d es de descenso: g^T d < 0
            if float(g.T @ d) < 0:
                while f(x + alpha * d) > f_val + c1 * alpha * float(g.T @ d):
                    alpha *= beta
                    if alpha < 1e-14:
                        break
            else:
                # Si el hessiano no es def. positivo y d no es de descenso, recurrir al gradiente
                d = -g
                while f(x + alpha * d) > f_val + c1 * alpha * float(g.T @ d):
                    alpha *= beta
                    if alpha < 1e-14:
                        break
                        
        history["damping_alphas"].append(alpha)
        step = alpha * d
        history["step_lengths"].append(np.linalg.norm(step))
        
        x = x + step
        history["x"].append(x.copy())
        history["f"].append(float(f(x)))
        
    return x, history


def main():
    print("=" * 80)
    print("ESFM-IPN :: MÓDULO 14: DESCENSO MÁS PRONUNCIADO, KANTÓROVICH Y NEWTON")
    print("Fiel a Luenberger & Ye (4ta Ed.), Secciones 8.2 y 8.5")
    print("=" * 80)
    
    # -------------------------------------------------------------------------
    # PARTE 1: VERIFICACIÓN MATEMÁTICA DE LA DESIGUALDAD DE KANTÓROVICH
    # -------------------------------------------------------------------------
    print("\n[1] VALIDACIÓN DE LA DESIGUALDAD DE KANTÓROVICH (Luenberger Sec. 8.2)")
    print("-" * 75)
    
    # Matriz mal condicionada con kappa = 100 (a = 2, A = 200)
    Q_test = np.array([
        [101.0, 99.0],
        [99.0, 101.0]
    ])
    # Autovalores: lambda_1 = 2, lambda_2 = 200 => kappa = 100
    res_kant = verify_kantorovich_inequality(Q_test, num_random_samples=1000)
    
    print(f"  • Autovalor mínimo (a):              {res_kant['a']:.4f}")
    print(f"  • Autovalor máximo (A):              {res_kant['A']:.4f}")
    print(f"  • Número de condición kappa:         {res_kant['A']/res_kant['a']:.2f}")
    print(f"  • Cota Inferior Teórica 4aA/(a+A)²:  {res_kant['cota_teorica']:.8f}")
    print(f"  • Cociente de Kantórovich Extremal:  {res_kant['ratio_extremal']:.8f}  (Coincidencia exacta!)")
    print(f"  • Mínimo en 1000 Muestras Random:    {res_kant['min_aleatorio']:.8f}")
    print(f"  • ¿Desigualdad Estrictamente Válida?: {'SÍ (Verificado)' if res_kant['satisface_desigualdad'] else 'FALLÓ'}")
    
    # -------------------------------------------------------------------------
    # PARTE 2: DESCENSO MÁS PRONUNCIADO EN FORMA CUADRÁTICA Y ORTOGONALIDAD
    # -------------------------------------------------------------------------
    print("\n[2] COMPORTAMIENTO DE DESCENSO MÁS PRONUNCIADO Y ORTOGONALIDAD g_{k+1}^T g_k = 0")
    print("-" * 75)
    
    b_test = np.array([0.0, 0.0])
    x0_test = np.array([2.0, 1.0])  # Punto de inicio fuera de los ejes de simetría
    
    x_opt, hist_sd = steepest_descent_quadratic(Q_test, b_test, x0_test, tol=1e-7, max_iter=150)
    
    print(f"  • Punto Inicial x0:             {x0_test}")
    print(f"  • Punto Óptimo Encontrado x*:   [{x_opt[0]:.6e}, {x_opt[1]:.6e}]")
    print(f"  • Iteraciones Requeridas:       {len(hist_sd['grad_norm'])}")
    print(f"  • Tasa Teórica de Kantórovich:  ((kappa-1)/(kappa+1))² = {hist_sd['kantorovich_bound']:.6f}")
    print(f"  • Tasa Máxima Observada:        {max(hist_sd['observed_rates']):.6f}")
    
    print("\n  Primeras 5 iteraciones y comprobación de ortogonalidad exacta:")
    print("  Iter |     Norma Gradiente    |      Error Energía     |  Paso alpha_k  |  Prod g_{k+1}^T g_k")
    print("  " + "-" * 73)
    for i in range(min(5, len(hist_sd["step_sizes"]))):
        print(f"  {i:4d} | {hist_sd['grad_norm'][i]:22.8e} | {hist_sd['energy_error'][i]:22.8e} | {hist_sd['step_sizes'][i]:14.6f} | {hist_sd['orthogonalities'][i]:18.4e}")
    
    # -------------------------------------------------------------------------
    # PARTE 3: IMPACTO DEL CONDICIONAMIENTO KAPPA EN EL NÚMERO DE ITERACIONES
    # -------------------------------------------------------------------------
    print("\n[3] IMPACTO DEL NÚMERO DE CONDICIÓN KAPPA EN LA VELOCIDAD DE CONVERGENCIA")
    print("-" * 75)
    print("  kappa (A/a) | Cota Kantórovich ((k-1)/(k+1))² | Iteraciones SD | Convergencia")
    print("  " + "-" * 73)
    
    kappas = [2, 10, 50, 100, 500]
    for k_val in kappas:
        Q_k = np.array([[float(k_val), 0.0], [0.0, 1.0]])
        _, h_k = steepest_descent_quadratic(Q_k, np.array([0.0, 0.0]), np.array([1.0, 1.0]), tol=1e-6, max_iter=3000)
        cota_k = ((k_val - 1.0) / (k_val + 1.0)) ** 2
        print(f"  {k_val:11d} | {cota_k:31.6f} | {len(h_k['grad_norm']):14d} | {'Muy Rápida' if k_val<=2 else 'Moderada' if k_val<=10 else 'Lenta (Zigzag)'}")
        
    # -------------------------------------------------------------------------
    # PARTE 4: MÉTODO DE NEWTON Y ORDEN CUADRÁTICO LOCAL (p = 2)
    # -------------------------------------------------------------------------
    print("\n[4] MÉTODO DE NEWTON MULTIDIMENSIONAL Y ORDEN CUADRÁTICO LOCAL (Luenberger Sec. 8.5)")
    print("-" * 75)
    
    # Función no cuadrática de prueba: f(x1, x2) = x1⁴ + x1*x2 + (1 + x2)²
    def test_func(x):
        return x[0]**4 + x[0]*x[1] + (1.0 + x[1])**2
        
    def test_grad(x):
        return np.array([
            4.0 * x[0]**3 + x[1],
            x[0] + 2.0 * (1.0 + x[1])
        ])
        
    def test_hess(x):
        return np.array([
            [12.0 * x[0]**2, 1.0],
            [1.0, 2.0]
        ])
        
    x0_newton = np.array([2.0, 2.0])
    x_opt_n, hist_n = newton_multidimensional(test_func, test_grad, test_hess, x0_newton, tol=1e-12, max_iter=30, damping=True)
    
    print(f"  • Punto Inicial x0:             {x0_newton}")
    print(f"  • Mínimo Encontrado x*:         [{x_opt_n[0]:.10f}, {x_opt_n[1]:.10f}]")
    print(f"  • Iteraciones Totales Newton:   {len(hist_n['grad_norm'])}")
    
    print("\n  Evolución del error ||x_k - x*|| y verificación del orden p = 2:")
    print("  Iter |    Norma ||x_k - x*||   |     Norma Gradiente    |  Cociente e_{k+1}/e_k²  |  Régimen")
    print("  " + "-" * 78)
    
    x_exact = x_opt_n
    errors = [np.linalg.norm(xk - x_exact) for xk in hist_n["x"]]
    
    for i in range(len(hist_n["grad_norm"])):
        e_k = errors[i]
        norm_g = hist_n["grad_norm"][i]
        if i > 0 and errors[i-1] > 1e-12:
            ratio_quad = e_k / (errors[i-1] ** 2)
            ratio_str = f"{ratio_quad:22.4f}"
        else:
            ratio_str = "                   N/A"
            
        regimen = "Transitorio" if e_k > 0.1 else ("Cuadrático puro (p=2)" if e_k > 1e-10 else "Tolerancia máquina")
        print(f"  {i:4d} | {e_k:22.8e} | {norm_g:22.8e} | {ratio_str} | {regimen}")

    # -------------------------------------------------------------------------
    # PARTE 5: BENCHMARK DIRECTO: STEEPEST DESCENT VS. NEWTON EN ROSENBROCK
    # -------------------------------------------------------------------------
    print("\n[5] BENCHMARK DIRECTO EN FUNCIÓN DE ROSENBROCK (Valle en Banana)")
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
        
    def rosenbrock_hess(x):
        return np.array([
            [-400.0 * (x[1] - 3.0 * x[0]**2) + 2.0, -400.0 * x[0]],
            [-400.0 * x[0], 200.0]
        ])
        
    x0_rosen = np.array([-1.2, 1.0])
    
    # 1. Steepest Descent con Armijo
    x_sd, hist_rosen_sd = steepest_descent_armijo(rosenbrock, rosenbrock_grad, x0_rosen, tol=1e-5, max_iter=5000)
    
    # 2. Newton Amortiguado
    x_nw, hist_rosen_nw = newton_multidimensional(rosenbrock, rosenbrock_grad, rosenbrock_hess, x0_rosen, tol=1e-5, max_iter=100, damping=True)
    
    print(f"  MÉTODO               | ITERACIONES | ||x_final - x*|| | f(x_final)     | EVALUACIONES CURVATURA")
    print("  " + "-" * 73)
    print(f"  Steepest Descent     | {len(hist_rosen_sd['grad_norm']):11d} | {np.linalg.norm(x_sd - np.array([1.0, 1.0])):16.6e} | {hist_rosen_sd['f'][-1]:14.6e} | 0 (Solo Primer Orden)")
    print(f"  Newton Amortiguado   | {len(hist_rosen_nw['grad_norm']):11d} | {np.linalg.norm(x_nw - np.array([1.0, 1.0])):16.6e} | {hist_rosen_nw['f'][-1]:14.6e} | {len(hist_rosen_nw['grad_norm'])} matrices Hessianas")
    print("=" * 80)
    print("CONCLUSIÓN DOCENTE (Luenberger pp. 235-245):")
    print("  1. Steepest Descent sufre una ralentización extrema por el fenómeno de zig-zag")
    print("     en presencia de anisotropía y curvatura estrecha (tasa acotada por Kantórovich).")
    print("  2. Newton converge en un orden de magnitud dramáticamente menor de iteraciones")
    print("     al incorporar la información de segundo orden, logrando convergencia cuadrática.")
    print("=" * 80)


if __name__ == "__main__":
    main()
