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
  },

  // ============================================================================
  // 4. VISUALIZADOR 2D: TEOREMA FUNDAMENTAL Y PUNTOS EXTREMOS (SBF) (Capítulo 2)
  // ============================================================================
  initLpFundamentalTool: function() {
    const canvas = document.getElementById('lpCanvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    const sliderC1 = document.getElementById('sliderC1');
    const sliderC2 = document.getElementById('sliderC2');
    const sliderB1 = document.getElementById('sliderB1');
    const sliderB2 = document.getElementById('sliderB2');

    const valC1 = document.getElementById('valC1');
    const valC2 = document.getElementById('valC2');
    const valB1 = document.getElementById('valB1');
    const valB2 = document.getElementById('valB2');

    const statVertices = document.getElementById('statLpVertices');
    const statOpt = document.getElementById('statLpOptimal');

    function redraw() {
      const c1 = parseFloat(sliderC1.value);
      const c2 = parseFloat(sliderC2.value);
      const b1 = parseFloat(sliderB1.value);
      const b2 = parseFloat(sliderB2.value);

      valC1.textContent = c1.toFixed(1);
      valC2.textContent = c2.toFixed(1);
      valB1.textContent = b1.toFixed(1);
      valB2.textContent = b2.toFixed(1);

      const width = canvas.width;
      const height = canvas.height;
      ctx.clearRect(0, 0, width, height);

      // Escala: plano cartesiano [0, 8] x [0, 8] con padding
      const xMax = 8.0, yMax = 8.0;
      const padLeft = 45, padBottom = 40, padTop = 20, padRight = 20;
      const plotW = width - padLeft - padRight;
      const plotH = height - padTop - padBottom;

      function toScreen(x, y) {
        return {
          px: padLeft + (x / xMax) * plotW,
          py: height - padBottom - (y / yMax) * plotH
        };
      }

      // 1. Dibujar Grid y Ejes
      ctx.strokeStyle = 'rgba(148, 163, 184, 0.2)';
      ctx.lineWidth = 1;
      for (let x = 0; x <= xMax; x += 1) {
        const p1 = toScreen(x, 0);
        const p2 = toScreen(x, yMax);
        ctx.beginPath();
        ctx.moveTo(p1.px, p1.py);
        ctx.lineTo(p2.px, p2.py);
        ctx.stroke();
      }
      for (let y = 0; y <= yMax; y += 1) {
        const p1 = toScreen(0, y);
        const p2 = toScreen(xMax, y);
        ctx.beginPath();
        ctx.moveTo(p1.px, p1.py);
        ctx.lineTo(p2.px, p2.py);
        ctx.stroke();
      }

      // Ejes coordenados principales
      ctx.strokeStyle = '#64748b';
      ctx.lineWidth = 2;
      const o = toScreen(0, 0);
      const ex = toScreen(xMax, 0);
      const ey = toScreen(0, yMax);

      ctx.beginPath();
      ctx.moveTo(o.px, o.py);
      ctx.lineTo(ex.px, ex.py);
      ctx.stroke();

      ctx.beginPath();
      ctx.moveTo(o.px, o.py);
      ctx.lineTo(ey.px, ey.py);
      ctx.stroke();

      // Etiquetas ejes
      ctx.fillStyle = '#94a3b8';
      ctx.font = '11px JetBrains Mono, monospace';
      for (let x = 2; x <= xMax; x += 2) {
        const p = toScreen(x, 0);
        ctx.fillText(x.toString(), p.px - 4, p.py + 18);
      }
      for (let y = 2; y <= yMax; y += 2) {
        const p = toScreen(0, y);
        ctx.fillText(y.toString(), p.px - 22, p.py + 4);
      }

      // 2. Determinar Vértices Factibles (Puntos Extremos)
      // Restricciones:
      // (R1) x1 + x2 <= b1
      // (R2) x1 + 2*x2 <= b2
      // x1 >= 0, x2 >= 0
      const candidates = [
        { x: 0, y: 0, name: '(0,0)' },
        { x: b1, y: 0, name: `(${b1.toFixed(1)},0)` },
        { x: 0, y: b2 / 2, name: `(0,${(b2/2).toFixed(1)})` }
      ];

      // Intersección de R1 y R2: x1 + x2 = b1, x1 + 2x2 = b2
      // x2 = b2 - b1,  x1 = 2*b1 - b2
      const interX2 = b2 - b1;
      const interX1 = 2 * b1 - b2;
      if (interX1 >= -1e-6 && interX2 >= -1e-6) {
        candidates.push({ x: interX1, y: interX2, name: `(${interX1.toFixed(1)},${interX2.toFixed(1)})` });
      }

      // Filtrar factibles
      const feasibleVertices = [];
      candidates.forEach(pt => {
        const sat1 = pt.x + pt.y <= b1 + 1e-6;
        const sat2 = pt.x + 2 * pt.y <= b2 + 1e-6;
        const satPos = pt.x >= -1e-6 && pt.y >= -1e-6;
        if (sat1 && sat2 && satPos) {
          if (!feasibleVertices.some(v => Math.abs(v.x - pt.x) < 1e-4 && Math.abs(v.y - pt.y) < 1e-4)) {
            feasibleVertices.push(pt);
          }
        }
      });

      // Ordenar vértices angularmente alrededor del centroide para dibujar el polígono
      if (feasibleVertices.length > 2) {
        const cx = feasibleVertices.reduce((sum, v) => sum + v.x, 0) / feasibleVertices.length;
        const cy = feasibleVertices.reduce((sum, v) => sum + v.y, 0) / feasibleVertices.length;
        feasibleVertices.sort((a, b) => Math.atan2(a.y - cy, a.x - cx) - Math.atan2(b.y - cy, b.x - cx));
      }

      // 3. Dibujar Región Factible (Poliedro Convexo)
      if (feasibleVertices.length >= 3) {
        ctx.beginPath();
        const start = toScreen(feasibleVertices[0].x, feasibleVertices[0].y);
        ctx.moveTo(start.px, start.py);
        for (let i = 1; i < feasibleVertices.length; i++) {
          const p = toScreen(feasibleVertices[i].x, feasibleVertices[i].y);
          ctx.lineTo(p.px, p.py);
        }
        ctx.closePath();
        ctx.fillStyle = 'rgba(56, 189, 248, 0.22)';
        ctx.fill();
        ctx.strokeStyle = '#38bdf8';
        ctx.lineWidth = 2;
        ctx.stroke();
      }

      // 4. Dibujar Líneas de Restricción
      // R1: x1 + x2 = b1 => (0, b1) a (b1, 0)
      const r1_p1 = toScreen(0, b1);
      const r1_p2 = toScreen(b1, 0);
      ctx.strokeStyle = '#f59e0b';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(r1_p1.px, r1_p1.py);
      ctx.lineTo(r1_p2.px, r1_p2.py);
      ctx.stroke();

      // R2: x1 + 2x2 = b2 => (0, b2/2) a (b2, 0)
      const r2_p1 = toScreen(0, b2 / 2);
      const r2_p2 = toScreen(b2, 0);
      ctx.strokeStyle = '#ec4899';
      ctx.beginPath();
      ctx.moveTo(r2_p1.px, r2_p1.py);
      ctx.lineTo(r2_p2.px, r2_p2.py);
      ctx.stroke();
      ctx.setLineDash([]);

      // 5. Evaluar Óptimo según z = c1*x1 + c2*x2 (Minimización)
      let optVertex = null;
      let minZ = Infinity;

      feasibleVertices.forEach(v => {
        const z = c1 * v.x + c2 * v.y;
        v.z = z;
        if (z < minZ) {
          minZ = z;
          optVertex = v;
        }
      });

      // 6. Dibujar Vértices y Resaltar Óptimo
      feasibleVertices.forEach(v => {
        const p = toScreen(v.x, v.y);
        const isOpt = (v === optVertex);

        ctx.beginPath();
        ctx.arc(p.px, p.py, isOpt ? 8 : 5, 0, Math.PI * 2);
        ctx.fillStyle = isOpt ? '#10b981' : '#38bdf8';
        ctx.fill();
        ctx.strokeStyle = isOpt ? '#ffffff' : '#0f172a';
        ctx.lineWidth = 2;
        ctx.stroke();

        ctx.fillStyle = isOpt ? '#10b981' : '#e2e8f0';
        ctx.font = isOpt ? 'bold 12px Inter, sans-serif' : '10px JetBrains Mono, monospace';
        ctx.fillText(`${v.name} [z=${v.z.toFixed(1)}]`, p.px + 10, p.py - 6);
      });

      // 7. Curva de Nivel Óptima: c1*x1 + c2*x2 = minZ
      if (optVertex && (Math.abs(c1) > 0.05 || Math.abs(c2) > 0.05)) {
        ctx.strokeStyle = '#10b981';
        ctx.lineWidth = 2;
        ctx.setLineDash([6, 3]);

        // Intersecciones de la recta c1*x + c2*y = minZ con los bordes
        let pLine1, pLine2;
        if (Math.abs(c2) > 0.05) {
          const y0 = (minZ - c1 * 0) / c2;
          const yX = (minZ - c1 * xMax) / c2;
          pLine1 = toScreen(0, y0);
          pLine2 = toScreen(xMax, yX);
        } else {
          const xVal = minZ / c1;
          pLine1 = toScreen(xVal, 0);
          pLine2 = toScreen(xVal, yMax);
        }
        ctx.beginPath();
        ctx.moveTo(pLine1.px, pLine1.py);
        ctx.lineTo(pLine2.px, pLine2.py);
        ctx.stroke();
        ctx.setLineDash([]);
      }

      // 8. Actualizar Estadísticas y Diagnóstico
      if (statVertices) {
        statVertices.textContent = `Vértices Factibles (SBFs): ${feasibleVertices.length} | Combinaciones C(4,2) = 6`;
      }
      if (statOpt && optVertex) {
        statOpt.innerHTML = `Vértice Óptimo (Teorema Fundamental): <strong style="color: #10b981;">x* = (${optVertex.x.toFixed(2)}, ${optVertex.y.toFixed(2)})</strong> con <strong style="color: #38bdf8;">z* = ${minZ.toFixed(2)}</strong>`;
      }
    }

    [sliderC1, sliderC2, sliderB1, sliderB2].forEach(slider => {
      if (slider) slider.addEventListener('input', redraw);
    });

    redraw();
  },

  // ============================================================================
  // 6. VISUALIZADOR INTERACTIVO DE SOPORTE VECTORIAL (SVM) — FIGURA 2.2
  // ============================================================================
  initSvmVisualizerTool: function() {
    const widget = document.getElementById('svmWidget');
    if (!widget) return;

    const btnCanonical = document.getElementById('svmModeCanonical');
    const btnMargins = document.getElementById('svmModeMargins');
    const btnSoft = document.getElementById('svmModeSoft');
    const sliderBeta = document.getElementById('svmBetaSlider');
    const valBeta = document.getElementById('svmBetaVal');
    const btnReset = document.getElementById('svmResetBtn');

    const dynamicGroup = document.getElementById('svmDynamicGroup');
    const dimLayer = document.getElementById('svmDimLayer');
    const softLayer = document.getElementById('svmSoftMarginLayer');
    const svTags = widget.querySelectorAll('.svm-sv-tag');
    const svPulses = widget.querySelectorAll('.svm-sv-pulse');

    const statusText = document.getElementById('svmStatusText');
    const metricText = document.getElementById('svmMetricText');

    // Vector normal unitario perpendicular en coordenadas de pantalla
    const ux = 0.534;
    const uy = -0.845;

    function setMode(mode) {
      [btnCanonical, btnMargins, btnSoft].forEach(btn => {
        if (btn) btn.classList.remove('active');
      });

      if (mode === 'canonical') {
        if (btnCanonical) btnCanonical.classList.add('active');
        if (dimLayer) dimLayer.style.display = 'block';
        if (softLayer) softLayer.style.display = 'none';
        svTags.forEach(el => el.style.display = 'none');
        svPulses.forEach(el => el.style.opacity = '0.25');
        if (statusText) statusText.innerHTML = 'Modo: <strong>Canónica de Luenberger (Figura 2.2)</strong>';
      } else if (mode === 'margins') {
        if (btnMargins) btnMargins.classList.add('active');
        if (dimLayer) dimLayer.style.display = 'block';
        if (softLayer) softLayer.style.display = 'none';
        svTags.forEach(el => el.style.display = 'block');
        svPulses.forEach(el => el.style.opacity = '1');
        if (statusText) statusText.innerHTML = 'Modo: <strong>Vectores de Soporte Activos (λᵢ*, μⱼ* &gt; 0)</strong>';
      } else if (mode === 'soft') {
        if (btnSoft) btnSoft.classList.add('active');
        if (dimLayer) dimLayer.style.display = 'block';
        if (softLayer) softLayer.style.display = 'block';
        svTags.forEach(el => el.style.display = 'block');
        svPulses.forEach(el => el.style.opacity = '0.5');
        if (statusText) statusText.innerHTML = 'Modo: <strong>Margen Suave con Holguras (Soft-Margin LP)</strong>';
      }
      updateDynamicElements();
    }

    function updateDynamicElements() {
      const deltaBeta = parseFloat(sliderBeta ? sliderBeta.value : 0);
      if (valBeta) {
        valBeta.textContent = (deltaBeta >= 0 ? '+' : '') + (deltaBeta * 0.1).toFixed(1);
      }

      // Desplazamiento del hiperplano y márgenes en el plano cartesiano
      const dx = deltaBeta * ux;
      const dy = deltaBeta * uy;
      if (dynamicGroup) {
        dynamicGroup.setAttribute('transform', `translate(${dx.toFixed(2)}, ${dy.toFixed(2)})`);
      }

      // Diagnóstico matemático de separabilidad
      if (metricText) {
        if (Math.abs(deltaBeta) <= 6) {
          metricText.innerHTML = 'Separación: <strong style="color: #10b981;">100% Factible (Margen Máximo)</strong>';
        } else if (deltaBeta > 6) {
          metricText.innerHTML = 'Separación: <strong style="color: #ef4444;">⚠️ Invasión de Margen (+1) por Clase 1</strong>';
        } else {
          metricText.innerHTML = 'Separación: <strong style="color: #f59e0b;">⚠️ Invasión de Margen (-1) por Clase 2</strong>';
        }
      }
    }

    if (btnCanonical) btnCanonical.addEventListener('click', () => setMode('canonical'));
    if (btnMargins) btnMargins.addEventListener('click', () => setMode('margins'));
    if (btnSoft) btnSoft.addEventListener('click', () => setMode('soft'));

    if (sliderBeta) {
      sliderBeta.addEventListener('input', updateDynamicElements);
    }

    if (btnReset) {
      btnReset.addEventListener('click', () => {
        if (sliderBeta) sliderBeta.value = 0;
        setMode('canonical');
      });
    }

    // Inicialización predeterminada
    setMode('canonical');
  },

  // ============================================================================
  // 7. CALCULADORA INTERACTIVA PRIMAL-DUAL Y HOLGURA COMPLEMENTARIA (Capítulo 4)
  // ============================================================================
  initDualityTool: function() {
    const widget = document.getElementById('dualityWidget');
    if (!widget) return;

    const sliderB1 = document.getElementById('sliderDualB1');
    const sliderB2 = document.getElementById('sliderDualB2');
    const sliderC1 = document.getElementById('sliderDualC1');
    const sliderC2 = document.getElementById('sliderDualC2');

    const valB1 = document.getElementById('valDualB1');
    const valB2 = document.getElementById('valDualB2');
    const valC1 = document.getElementById('valDualC1');
    const valC2 = document.getElementById('valDualC2');

    const dispPrimal = document.getElementById('dispPrimalResult');
    const dispDual = document.getElementById('dispDualResult');
    const badgeGap = document.getElementById('badgeDualGap');
    const dispSlackness = document.getElementById('dispSlacknessList');

    function updateDuality() {
      const b1 = parseFloat(sliderB1 ? sliderB1.value : 4);
      const b2 = parseFloat(sliderB2 ? sliderB2.value : 6);
      const c1 = parseFloat(sliderC1 ? sliderC1.value : 6);
      const c2 = parseFloat(sliderC2 ? sliderC2.value : 8);

      if (valB1) valB1.textContent = b1.toFixed(1);
      if (valB2) valB2.textContent = b2.toFixed(1);
      if (valC1) valC1.textContent = c1.toFixed(1);
      if (valC2) valC2.textContent = c2.toFixed(1);

      // Primal Vertices (Minimización)
      // Restricciones: x1 + x2 >= b1, x1 + 2x2 >= b2, x1, x2 >= 0
      const pv = [
        { x1: Math.max(b1, b2), x2: 0, label: 'Eje X₁' },
        { x1: 0, x2: Math.max(b1, b2 / 2), label: 'Eje X₂' }
      ];
      if (b2 / 2 <= b1 && b1 <= b2) {
        pv.push({ x1: 2 * b1 - b2, x2: b2 - b1, label: 'Intersección R1 y R2' });
      }

      let bestZ = Infinity;
      let optX = { x1: 0, x2: 0, label: '' };
      pv.forEach(v => {
        if (v.x1 >= -1e-5 && v.x2 >= -1e-5 &&
            (v.x1 + v.x2) >= b1 - 1e-5 &&
            (v.x1 + 2 * v.x2) >= b2 - 1e-5) {
          const z = c1 * v.x1 + c2 * v.x2;
          if (z < bestZ) {
            bestZ = z;
            optX = v;
          }
        }
      });

      // Dual Vertices (Maximización)
      // Restricciones: y1 + y2 <= c1, y1 + 2y2 <= c2, y1, y2 >= 0
      const dv = [
        { y1: 0, y2: 0, label: 'Origen' },
        { y1: Math.min(c1, c2), y2: 0, label: 'Eje Y₁' },
        { y1: 0, y2: Math.min(c1, c2 / 2), label: 'Eje Y₂' }
      ];
      if (c2 / 2 <= c1 && c1 <= c2) {
        dv.push({ y1: 2 * c1 - c2, y2: c2 - c1, label: 'Intersección D1 y D2' });
      }

      let bestW = -Infinity;
      let optY = { y1: 0, y2: 0, label: '' };
      dv.forEach(v => {
        if (v.y1 >= -1e-5 && v.y2 >= -1e-5 &&
            (v.y1 + v.y2) <= c1 + 1e-5 &&
            (v.y1 + 2 * v.y2) <= c2 + 1e-5) {
          const w = b1 * v.y1 + b2 * v.y2;
          if (w > bestW) {
            bestW = w;
            optY = v;
          }
        }
      });

      // Holguras
      const s1 = optX.x1 + optX.x2 - b1;
      const s2 = optX.x1 + 2 * optX.x2 - b2;
      const e1 = c1 - (optY.y1 + optY.y2);
      const e2 = c2 - (optY.y1 + 2 * optY.y2);

      const gap = Math.abs(bestZ - bestW);

      if (dispPrimal) {
        dispPrimal.innerHTML = `
          <div>• <strong>Vértice Óptimo x*:</strong> (${optX.x1.toFixed(2)}, ${optX.x2.toFixed(2)}) <span style="font-size:0.85rem; color:var(--text-muted);">[${optX.label}]</span></div>
          <div>• <strong>Costo Óptimo z*:</strong> <span style="font-weight:700; color:#2563eb; font-size:1.1rem;">${bestZ.toFixed(2)}</span></div>
          <div>• <strong>Holguras Primales:</strong> s₁ = ${s1.toFixed(2)} ${s1 < 0.01 ? '<span style="color:#10b981;">(Activa)</span>' : '<span style="color:var(--text-muted);">(Inactiva)</span>'}, s₂ = ${s2.toFixed(2)} ${s2 < 0.01 ? '<span style="color:#10b981;">(Activa)</span>' : '<span style="color:var(--text-muted);">(Inactiva)</span>'}</div>
        `;
      }

      if (dispDual) {
        dispDual.innerHTML = `
          <div>• <strong>Vector Dual y*:</strong> (${optY.y1.toFixed(2)}, ${optY.y2.toFixed(2)}) <span style="font-size:0.85rem; color:var(--text-muted);">[${optY.label}]</span></div>
          <div>• <strong>Beneficio Dual w*:</strong> <span style="font-weight:700; color:#10b981; font-size:1.1rem;">${bestW.toFixed(2)}</span></div>
          <div>• <strong>Holguras Duales:</strong> e₁ = ${e1.toFixed(2)}, e₂ = ${e2.toFixed(2)}</div>
        `;
      }

      if (badgeGap) {
        if (gap < 1e-4) {
          badgeGap.style.background = '#10b981';
          badgeGap.textContent = `Δ = 0.000 (Dualidad Fuerte Verificada: z* = w* = ${bestZ.toFixed(2)})`;
        } else {
          badgeGap.style.background = '#ef4444';
          badgeGap.textContent = `Δ = ${gap.toFixed(3)} (Brecha de Dualidad Positiva)`;
        }
      }

      if (dispSlackness) {
        const prodX1 = Math.abs(optX.x1 * e1);
        const prodX2 = Math.abs(optX.x2 * e2);
        const prodY1 = Math.abs(optY.y1 * s1);
        const prodY2 = Math.abs(optY.y2 * s2);

        dispSlackness.innerHTML = `
          <div style="background: var(--bg-card); padding: 0.5rem; border-radius: 4px; border: 1px solid var(--border-color);">
            <strong>x₁ · e₁ = 0:</strong> ${(optX.x1).toFixed(2)} × ${(e1).toFixed(2)} = <strong>${prodX1.toFixed(2)}</strong>
            <span style="color:#10b981; font-weight:bold;"> ✓</span>
          </div>
          <div style="background: var(--bg-card); padding: 0.5rem; border-radius: 4px; border: 1px solid var(--border-color);">
            <strong>x₂ · e₂ = 0:</strong> ${(optX.x2).toFixed(2)} × ${(e2).toFixed(2)} = <strong>${prodX2.toFixed(2)}</strong>
            <span style="color:#10b981; font-weight:bold;"> ✓</span>
          </div>
          <div style="background: var(--bg-card); padding: 0.5rem; border-radius: 4px; border: 1px solid var(--border-color);">
            <strong>y₁ · s₁ = 0:</strong> ${(optY.y1).toFixed(2)} × ${(s1).toFixed(2)} = <strong>${prodY1.toFixed(2)}</strong>
            <span style="color:#10b981; font-weight:bold;"> ✓</span>
          </div>
          <div style="background: var(--bg-card); padding: 0.5rem; border-radius: 4px; border: 1px solid var(--border-color);">
            <strong>y₂ · s₂ = 0:</strong> ${(optY.y2).toFixed(2)} × ${(s2).toFixed(2)} = <strong>${prodY2.toFixed(2)}</strong>
            <span style="color:#10b981; font-weight:bold;"> ✓</span>
          </div>
        `;
      }
    }

    [sliderB1, sliderB2, sliderC1, sliderC2].forEach(sl => {
      if (sl) sl.addEventListener('input', updateDuality);
    });

    updateDuality();
  }
};

