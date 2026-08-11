/* ==============================================================================
   Módulos Interactivos: Visualizadores Matemáticos y Solvers en el Navegador
   Libro Digital Interactivo — Unidad 1: Herramientas para la Optimización
   David G. Luenberger & Yinyu Ye (4ta Edición)
============================================================================== */

window.InteractiveTools = {
  // ============================================================================
  // 1. VISUALIZADOR DE FORMAS CUADRÁTICAS Y SIGNATURA ESPECTRAL (Apéndice A.4)
  // ============================================================================
  initQuadraticFormTool: function() {
    const canvas = document.getElementById('quadCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    
    const sliderA11 = document.getElementById('sliderA11');
    const sliderA12 = document.getElementById('sliderA12');
    const sliderA22 = document.getElementById('sliderA22');
    
    const valA11 = document.getElementById('valA11');
    const valA12 = document.getElementById('valA12');
    const valA22 = document.getElementById('valA22');
    const badgeSig = document.getElementById('signatureBadge');
    const statEig = document.getElementById('eigenvalStats');

    function redraw() {
      const a11 = parseFloat(sliderA11.value);
      const a12 = parseFloat(sliderA12.value);
      const a22 = parseFloat(sliderA22.value);

      valA11.textContent = a11.toFixed(1);
      valA12.textContent = a12.toFixed(1);
      valA22.textContent = a22.toFixed(1);

      // Autovalores de matriz 2x2 simétrica A = [[a11, a12], [a12, a22]]
      const trace = a11 + a22;
      const det = a11 * a22 - a12 * a12;
      const disc = Math.sqrt(Math.max(0, trace * trace - 4 * det));
      const lambda1 = (trace + disc) / 2;
      const lambda2 = (trace - disc) / 2;

      statEig.textContent = `λ₁ = ${lambda1.toFixed(2)}, λ₂ = ${lambda2.toFixed(2)}, det(A) = ${det.toFixed(2)}`;

      // Signatura
      if (lambda1 > 0.05 && lambda2 > 0.05) {
        badgeSig.textContent = 'Definida Positiva (Mínimo Estricto)';
        badgeSig.style.background = '#dbeafe';
        badgeSig.style.color = '#1d4ed8';
      } else if (lambda1 < -0.05 && lambda2 < -0.05) {
        badgeSig.textContent = 'Definida Negativa (Máximo Estricto)';
        badgeSig.style.background = '#fee2e2';
        badgeSig.style.color = '#b91c1c';
      } else if (lambda1 * lambda2 < -0.05) {
        badgeSig.textContent = 'Indefinida (Punto de Ensilladura)';
        badgeSig.style.background = '#fef3c7';
        badgeSig.style.color = '#b45309';
      } else {
        badgeSig.textContent = 'Semidefinida (Autovalor Nulo)';
        badgeSig.style.background = '#f1f5f9';
        badgeSig.style.color = '#475569';
      }

      // Dibujar curvas de nivel en Canvas
      const width = canvas.width;
      const height = canvas.height;
      ctx.clearRect(0, 0, width, height);

      const cx = width / 2;
      const cy = height / 2;
      const scale = 40; // px por unidad

      // Ejes coordenados
      ctx.strokeStyle = 'rgba(150, 150, 150, 0.4)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(0, cy); ctx.lineTo(width, cy);
      ctx.moveTo(cx, 0); ctx.lineTo(cx, height);
      ctx.stroke();

      // Malla de curvas de nivel
      const levels = [0.5, 1.0, 2.0, 4.0, 8.0, -0.5, -1.0, -2.0, -4.0];
      levels.forEach(c => {
        ctx.beginPath();
        ctx.strokeStyle = c > 0 ? 'rgba(37, 99, 235, 0.65)' : 'rgba(239, 68, 68, 0.65)';
        ctx.lineWidth = 1.5;

        let first = true;
        for (let angle = 0; angle <= Math.PI * 2 + 0.1; angle += 0.05) {
          const cos = Math.cos(angle);
          const sin = Math.sin(angle);
          const denom = a11 * cos * cos + 2 * a12 * cos * sin + a22 * sin * sin;

          if (denom * c > 0) {
            const r = Math.sqrt(c / denom);
            if (r > 0 && r < 10) {
              const px = cx + r * cos * scale;
              const py = cy - r * sin * scale;
              if (first) { ctx.moveTo(px, py); first = false; }
              else { ctx.lineTo(px, py); }
            }
          }
        }
        ctx.stroke();
      });
    }

    sliderA11.addEventListener('input', redraw);
    sliderA12.addEventListener('input', redraw);
    sliderA22.addEventListener('input', redraw);
    redraw();
  },

  // ============================================================================
  // 2. SIMULADOR DEL TEOREMA 1 DE SEPARACIÓN DE LUENBERGER (Apéndice B.3)
  // ============================================================================
  initSeparatingHyperplaneTool: function() {
    const canvas = document.getElementById('sepCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const statDelta = document.getElementById('sepDeltaStat');
    const statVerdict = document.getElementById('sepVerdict');

    // Conjunto convexo: Elipse (x1/2)^2 + x2^2 <= 1
    const rx = 2.0;
    const ry = 1.0;
    let y_ext = { x: 2.8, y: 1.8 };
    let isDragging = false;

    const width = canvas.width;
    const height = canvas.height;
    const cx = width / 2 - 30;
    const cy = height / 2;
    const scale = 50; // px por unidad

    function toScreen(p) { return { x: cx + p.x * scale, y: cy - p.y * scale }; }
    function toMath(s) { return { x: (s.x - cx) / scale, y: (cy - s.y) / scale }; }

    // Proyección de mínima distancia sobre la elipse (x0)
    function findProjection(y) {
      let bestX0 = { x: 0, y: 0 };
      let minD = Infinity;
      for (let t = 0; t < Math.PI * 2; t += 0.01) {
        const ex = rx * Math.cos(t);
        const ey = ry * Math.sin(t);
        const d = (ex - y.x) ** 2 + (ey - y.y) ** 2;
        if (d < minD) {
          minD = d;
          bestX0 = { x: ex, y: ey };
        }
      }
      return { x0: bestX0, delta: Math.sqrt(minD) };
    }

    function render() {
      ctx.clearRect(0, 0, width, height);

      // Ejes
      ctx.strokeStyle = 'rgba(150, 150, 150, 0.3)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(0, cy); ctx.lineTo(width, cy);
      ctx.moveTo(cx, 0); ctx.lineTo(cx, height);
      ctx.stroke();

      // Dibujar conjunto convexo C
      ctx.beginPath();
      ctx.ellipse(cx, cy, rx * scale, ry * scale, 0, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(59, 130, 246, 0.18)';
      ctx.fill();
      ctx.strokeStyle = '#2563eb';
      ctx.lineWidth = 2;
      ctx.stroke();

      // Etiqueta conjunto C
      ctx.fillStyle = '#2563eb';
      ctx.font = 'bold 13px Inter, sans-serif';
      ctx.fillText('Conjunto Convexo C', cx - 55, cy + 5);

      const proj = findProjection(y_ext);
      const x0 = proj.x0;
      const delta = proj.delta;

      // Vector a = x0 - y
      const a = { x: x0.x - y_ext.x, y: x0.y - y_ext.y };
      const ay = a.x * y_ext.x + a.y * y_ext.y;
      const ax0 = a.x * x0.x + a.y * x0.y;

      statDelta.textContent = `δ = ||x₀ - y|| = ${delta.toFixed(3)} | a = (${a.x.toFixed(2)}, ${a.y.toFixed(2)})`;
      statVerdict.innerHTML = `<strong>Teorema 1 Verificado:</strong> aᵀy = ${ay.toFixed(2)} < inf aᵀx = aᵀx₀ = ${ax0.toFixed(2)}`;

      const sy = toScreen(y_ext);
      const sx0 = toScreen(x0);

      // Segmento delta
      ctx.strokeStyle = '#ef4444';
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(sy.x, sy.y);
      ctx.lineTo(sx0.x, sx0.y);
      ctx.stroke();
      ctx.setLineDash([]);

      // Hiperplano separador (de apoyo en x0): a^T x = a^T x0
      ctx.strokeStyle = '#10b981';
      ctx.lineWidth = 2;
      ctx.beginPath();
      // Vector direccional de la recta (perpendicular a normal 'a')
      const dir = { x: -a.y, y: a.x };
      const t = 20; // Longitud suficientemente grande
      const p1 = { x: x0.x + t * dir.x, y: x0.y + t * dir.y };
      const p2 = { x: x0.x - t * dir.x, y: x0.y - t * dir.y };
      const sp1 = toScreen(p1);
      const sp2 = toScreen(p2);
      ctx.moveTo(sp1.x, sp1.y);
      ctx.lineTo(sp2.x, sp2.y);
      ctx.stroke();

      // Dibujar punto x0
      ctx.fillStyle = '#10b981';
      ctx.beginPath();
      ctx.arc(sx0.x, sx0.y, 6, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillStyle = '#10b981';
      ctx.fillText('x₀ (Mín. Distancia)', sx0.x + 8, sx0.y - 8);

      // Dibujar punto exterior y (interactivo)
      ctx.fillStyle = '#ef4444';
      ctx.beginPath();
      ctx.arc(sy.x, sy.y, 8, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 2;
      ctx.stroke();
      ctx.fillStyle = '#ef4444';
      ctx.fillText('y (Exterior)', sy.x + 10, sy.y + 4);
    }

    canvas.addEventListener('mousedown', (e) => {
      const rect = canvas.getBoundingClientRect();
      const clickS = { x: e.clientX - rect.left, y: e.clientY - rect.top };
      const sy = toScreen(y_ext);
      if (Math.hypot(clickS.x - sy.x, clickS.y - sy.y) < 15) {
        isDragging = true;
      }
    });

    window.addEventListener('mousemove', (e) => {
      if (!isDragging) return;
      const rect = canvas.getBoundingClientRect();
      const s = { x: e.clientX - rect.left, y: e.clientY - rect.top };
      const m = toMath(s);
      // Validar que permanezca exterior
      if ((m.x / rx) ** 2 + (m.y / ry) ** 2 > 1.05) {
        y_ext = m;
        render();
      }
    });

    window.addEventListener('mouseup', () => { isDragging = false; });
    render();
  },

  // ============================================================================
  // 3. CALCULADORA DE FACTORIZACIÓN LU DE GAUSS (Apéndice C)
  // ============================================================================
  initLUSolverTool: function() {
    const btnSolve = document.getElementById('btnSolveLU');
    if (!btnSolve) return;
    const outputDiv = document.getElementById('luOutput');

    btnSolve.addEventListener('click', () => {
      const a11 = parseFloat(document.getElementById('m_a11').value) || 2;
      const a12 = parseFloat(document.getElementById('m_a12').value) || 1;
      const a13 = parseFloat(document.getElementById('m_a13').value) || 1;
      const a21 = parseFloat(document.getElementById('m_a21').value) || 4;
      const a22 = parseFloat(document.getElementById('m_a22').value) || 3;
      const a23 = parseFloat(document.getElementById('m_a23').value) || 3;
      const a31 = parseFloat(document.getElementById('m_a31').value) || 8;
      const a32 = parseFloat(document.getElementById('m_a32').value) || 7;
      const a33 = parseFloat(document.getElementById('m_a33').value) || 9;

      const b1 = parseFloat(document.getElementById('m_b1').value) || 5;
      const b2 = parseFloat(document.getElementById('m_b2').value) || 13;
      const b3 = parseFloat(document.getElementById('m_b3').value) || 37;

      let A = [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]];
      let L = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      let U = JSON.parse(JSON.stringify(A));

      // Paso 1: Eliminar columna 1
      const m21 = U[1][0] / U[0][0];
      const m31 = U[2][0] / U[0][0];
      L[1][0] = m21;
      L[2][0] = m31;
      for (let j = 0; j < 3; j++) {
        U[1][j] -= m21 * U[0][j];
        U[2][j] -= m31 * U[0][j];
      }

      // Paso 2: Eliminar columna 2
      const m32 = U[2][1] / U[1][1];
      L[2][1] = m32;
      for (let j = 0; j < 3; j++) {
        U[2][j] -= m32 * U[1][j];
      }

      // Solución bifásica: Ly = b
      const y1 = b1;
      const y2 = b2 - L[1][0] * y1;
      const y3 = b3 - L[2][0] * y1 - L[2][1] * y2;

      // Ux = y
      const x3 = y3 / U[2][2];
      const x2 = (y2 - U[1][2] * x3) / U[1][1];
      const x1 = (y1 - U[0][1] * x2 - U[0][2] * x3) / U[0][0];

      outputDiv.innerHTML = `
        <div style="margin-top: 1.25rem; background: var(--bg-tertiary); padding: 1.25rem; border-radius: var(--radius-md); border: 1px solid var(--card-border);">
          <h4 style="color: var(--accent-primary); margin-top:0;">Resultados de la Factorización LU de Luenberger:</h4>
          <p><strong>1. Multiplicadores de Gauss:</strong> $m_{21} = ${m21.toFixed(2)}, \\; m_{31} = ${m31.toFixed(2)}, \\; m_{32} = ${m32.toFixed(2)}$</p>
          <div style="display: flex; gap: 2.5rem; flex-wrap: wrap; margin: 1rem 0;">
            <div>
              <strong>Matriz Triangular Inferior L:</strong><br>
              $\\begin{pmatrix} 1 & 0 & 0 \\\\ ${L[1][0].toFixed(2)} & 1 & 0 \\\\ ${L[2][0].toFixed(2)} & ${L[2][1].toFixed(2)} & 1 \\end{pmatrix}$
            </div>
            <div>
              <strong>Matriz Triangular Superior U:</strong><br>
              $\\begin{pmatrix} ${U[0][0].toFixed(2)} & ${U[0][1].toFixed(2)} & ${U[0][2].toFixed(2)} \\\\ 0 & ${U[1][1].toFixed(2)} & ${U[1][2].toFixed(2)} \\\\ 0 & 0 & ${U[2][2].toFixed(2)} \\end{pmatrix}$
            </div>
          </div>
          <p style="margin-top: 1rem;"><strong>2. Sustitución Progresiva (Hacia Adelante $Ly = b$):</strong> $y = (${y1.toFixed(2)}, \\; ${y2.toFixed(2)}, \\; ${y3.toFixed(2)})^T$</p>
          <p><strong>3. Sustitución Regresiva (Hacia Atrás $Ux = y$):</strong> <span style="font-size:1.1rem; color: #10b981; font-weight:bold;">$x^* = (${x1.toFixed(2)}, \\; ${x2.toFixed(2)}, \\; ${x3.toFixed(2)})^T$</span></p>
        </div>
      `;

      if (window.renderMathInElement) {
        window.renderMathInElement(outputDiv, {
          delimiters: [
            {left: '$$', right: '$$', display: true},
            {left: '$', right: '$', display: false}
          ]
        });
      }
    });
  }
};
