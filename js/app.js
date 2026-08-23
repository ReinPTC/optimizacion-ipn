/* ==============================================================================
   Lógica de Aplicación Principal: Navegación, Búsqueda, Temas y KaTeX
   Libro Digital Interactivo — Unidad 1 (Optimización IPN)
============================================================================== */

function triggerKaTeX() {
  if (typeof renderMathInElement === 'function') {
    try {
      renderMathInElement(document.body, {
        delimiters: [
          { left: '$$', right: '$$', display: true },
          { left: '$', right: '$', display: false },
          { left: '\\(', right: '\\)', display: false },
          { left: '\\[', right: '\\]', display: true }
        ],
        throwOnError: false,
        errorColor: '#ef4444'
      });
    } catch (e) {
      console.warn('[KaTeX] Advertencia de renderizado:', e);
    }
  }
}

// Ejecutar KaTeX tanto en DOMContentLoaded como en window.onload para evitar condiciones de carrera con scripts defer
document.addEventListener('DOMContentLoaded', () => {
  triggerKaTeX();

  // 1. Control de Tema Oscuro / Claro
  const themeToggleBtn = document.getElementById('themeToggleBtn');
  const currentTheme = localStorage.getItem('app-theme') || 'dark';
  document.documentElement.setAttribute('data-theme', currentTheme);
  updateThemeButton(currentTheme);

  if (themeToggleBtn) {
    themeToggleBtn.addEventListener('click', () => {
      const activeTheme = document.documentElement.getAttribute('data-theme');
      const newTheme = activeTheme === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', newTheme);
      localStorage.setItem('app-theme', newTheme);
      updateThemeButton(newTheme);
    });
  }

  function updateThemeButton(theme) {
    if (!themeToggleBtn) return;
    themeToggleBtn.innerHTML = theme === 'dark' 
      ? '☀️ Modo Claro' 
      : '🌙 Modo Oscuro';
  }

  // 2. Menú Móvil & Control de Barra Lateral Responsive
  const mobileMenuBtn = document.getElementById('mobileMenuBtn');
  const sidebarCloseBtn = document.getElementById('sidebarCloseBtn');
  const sidebarOverlay = document.getElementById('sidebarOverlay');
  const sidebar = document.getElementById('sidebar');

  function openSidebar() {
    if (sidebar) sidebar.classList.add('open');
    if (sidebarOverlay) sidebarOverlay.classList.add('active');
    document.body.style.overflow = window.innerWidth <= 1024 ? 'hidden' : '';
  }

  function closeSidebar() {
    if (sidebar) sidebar.classList.remove('open');
    if (sidebarOverlay) sidebarOverlay.classList.remove('active');
    document.body.style.overflow = '';
  }

  if (mobileMenuBtn) {
    mobileMenuBtn.addEventListener('click', () => {
      if (sidebar && sidebar.classList.contains('open')) {
        closeSidebar();
      } else {
        openSidebar();
      }
    });
  }

  if (sidebarCloseBtn) {
    sidebarCloseBtn.addEventListener('click', closeSidebar);
  }

  if (sidebarOverlay) {
    sidebarOverlay.addEventListener('click', closeSidebar);
  }

  // Cerrar la barra lateral en móviles al tocar cualquier enlace de navegación
  document.querySelectorAll('.sidebar-nav a').forEach(link => {
    link.addEventListener('click', () => {
      if (window.innerWidth <= 1024) {
        closeSidebar();
      }
    });
  });

  // Cerrar con tecla Escape
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && sidebar && sidebar.classList.contains('open')) {
      closeSidebar();
    }
  });

  // 3. Búsqueda en Tiempo Real en la Barra Lateral
  const searchInput = document.getElementById('searchInput');
  const navItems = document.querySelectorAll('.nav-item');
  if (searchInput) {
    searchInput.addEventListener('input', (e) => {
      const query = e.target.value.toLowerCase().trim();
      navItems.forEach(item => {
        const text = item.textContent.toLowerCase();
        if (text.includes(query)) {
          item.style.display = 'flex';
        } else {
          item.style.display = 'none';
        }
      });
    });
  }

  // 4. ScrollSpy & Active Link Indicator
  const sections = document.querySelectorAll('section[id], h2[id], h3[id], div[id^="sec-"]');
  const breadcrumbCurrent = document.getElementById('breadcrumbCurrent');

  const observerOptions = {
    root: null,
    rootMargin: '-70px 0px -65% 0px',
    threshold: 0
  };

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        const id = entry.target.getAttribute('id');
        const activeLink = document.querySelector(`.sidebar-nav a[href="#${id}"]`);
        if (activeLink) {
          navItems.forEach(link => link.classList.remove('active'));
          activeLink.classList.add('active');
          if (breadcrumbCurrent) {
            breadcrumbCurrent.textContent = activeLink.textContent.replace(/[🎯📋🐍📚⚡🧮📐]/g, '').trim();
          }
        }
      }
    });
  }, observerOptions);

  sections.forEach(section => observer.observe(section));

  // 5. Botones para Copiar Código
  document.querySelectorAll('.btn-copy').forEach(button => {
    button.addEventListener('click', () => {
      const targetId = button.getAttribute('data-target');
      const codeElement = document.getElementById(targetId);
      if (codeElement) {
        navigator.clipboard.writeText(codeElement.textContent).then(() => {
          const originalText = button.textContent;
          button.textContent = '✓ ¡Copiado!';
          button.style.color = '#34d399';
          setTimeout(() => {
            button.textContent = originalText;
            button.style.color = '';
          }, 2000);
        });
      }
    });
  });

  // 6. Registro Seguro de Service Worker PWA
  if ('serviceWorker' in navigator && window.location.protocol.startsWith('http')) {
    window.addEventListener('load', () => {
      navigator.serviceWorker.register('./sw.js')
        .then((reg) => console.log('[PWA] Service Worker activo:', reg.scope))
        .catch((err) => console.log('[PWA] Info SW:', err.message));
    });
  }

  // 7. Inicialización de Herramientas Interactivas
  if (window.InteractiveTools) {
    window.InteractiveTools.initQuadraticFormTool();
    window.InteractiveTools.initSeparatingHyperplaneTool();
    window.InteractiveTools.initLUSolverTool();
  }
});

// Re-ejecución de KaTeX en evento load
window.addEventListener('load', () => {
  triggerKaTeX();
});
