import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Circle
import os
import shapely.geometry as sg
from shapely.ops import unary_union

# Create assets folder
os.makedirs("assets", exist_ok=True)

# -----------------
# Fig B.1 Convexity
# -----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

# Convex Set
t = np.linspace(0, 2*np.pi, 100)
x_c = 2 * np.cos(t) + 0.5 * np.cos(2*t)
y_c = 2 * np.sin(t) + 0.5 * np.sin(2*t)
ax1.plot(x_c, y_c, 'k-', lw=1.5)

x1, y1 = -1, -0.5
x2, y2 = 1.5, 1.2
ax1.plot([x1, x2], [y1, y2], 'k-', lw=1.5)
ax1.plot(x1, y1, 'ko', markersize=5)
ax1.plot(x2, y2, 'ko', markersize=5)
ax1.text(x1-0.3, y1-0.4, r'$\mathbf{x}_1$', fontsize=14)
ax1.text(x2+0.1, y2+0.2, r'$\mathbf{x}_2$', fontsize=14)
ax1.set_xlim(-3, 3)
ax1.set_ylim(-3, 3)
ax1.axis('off')
ax1.text(0, -3.5, 'convex', fontsize=14, ha='center')

# Nonconvex Set
x_nc = 2 * np.cos(t) - 1.2 * np.cos(2*t)
y_nc = 2 * np.sin(t) - 0.5 * np.sin(2*t)
ax2.plot(x_nc, y_nc, 'k-', lw=1.5)

x1_nc, y1_nc = -1.5, -1
x2_nc, y2_nc = 1, 1.5
ax2.plot([x1_nc, x2_nc], [y1_nc, y2_nc], 'k-', lw=1.5)
ax2.plot(x1_nc, y1_nc, 'ko', markersize=5)
ax2.plot(x2_nc, y2_nc, 'ko', markersize=5)
ax2.text(x1_nc-0.3, y1_nc-0.4, r'$\mathbf{x}_1$', fontsize=14)
ax2.text(x2_nc+0.1, y2_nc+0.2, r'$\mathbf{x}_2$', fontsize=14)
ax2.set_xlim(-4, 4)
ax2.set_ylim(-3, 3)
ax2.axis('off')
ax2.text(0, -3.5, 'nonconvex', fontsize=14, ha='center')

plt.tight_layout()
plt.savefig('assets/fig_b1.png', dpi=300, bbox_inches='tight')
plt.close()

# -----------------
# Fig B.2 Properties
# -----------------
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# i) Scaling (2*C)
ax = axes[0]
circle1 = plt.Circle((1, 1.5), 0.8, fill=False, color='k', lw=1.5)
circle2 = plt.Circle((2, 3), 1.6, fill=False, color='k', lw=1.5)
ax.add_patch(circle1)
ax.add_patch(circle2)
ax.plot(0, 0, 'ko', markersize=5)
ax.text(0, -0.5, '0', fontsize=14, ha='center')
ax.text(1, 1.5, '$C$', fontsize=14, ha='center', va='center')
ax.text(2, 3, '$2 \cdot C$', fontsize=14, ha='center', va='center')
ax.set_xlim(-1, 5)
ax.set_ylim(-1, 6)
ax.axis('off')

# ii) Sum (C + D)
ax = axes[1]
t = np.linspace(0, 2*np.pi, 100)
C_poly = sg.Polygon(np.column_stack((1 + 0.6*np.cos(t), 1.5 + 0.8*np.sin(t))))
D_poly = sg.Polygon(np.column_stack((2 + 0.7*np.cos(t) + 0.2*np.cos(2*t), 0.5 + 0.3*np.sin(t))))

x_c, y_c = C_poly.exterior.xy
x_d, y_d = D_poly.exterior.xy
ax.plot(x_c, y_c, 'k-', lw=1.5)
ax.plot(x_d, y_d, 'k-', lw=1.5)

# Minkowski Sum
sum_points = []
for p1 in np.array(C_poly.exterior.coords):
    for p2 in np.array(D_poly.exterior.coords):
        sum_points.append(p1 + p2)
sum_poly = sg.MultiPoint(sum_points).convex_hull
x_s, y_s = sum_poly.exterior.xy
ax.plot(x_s, y_s, 'k-', lw=1.5)

C_center = np.array([1, 1.5])
D_center = np.array([2, 0.5])
sum_center = C_center + D_center

ax.plot([0, C_center[0]], [0, C_center[1]], 'k-', lw=0.5, alpha=0.5)
ax.plot([0, D_center[0]], [0, D_center[1]], 'k-', lw=0.5, alpha=0.5)
ax.plot([C_center[0], sum_center[0]], [C_center[1], sum_center[1]], 'k-', lw=0.5, alpha=0.5)
ax.plot([D_center[0], sum_center[0]], [D_center[1], sum_center[1]], 'k-', lw=0.5, alpha=0.5)

ax.plot(0, 0, 'ko', markersize=5)
ax.text(0, -0.4, '0', fontsize=14, ha='center')
ax.text(C_center[0], C_center[1], '$C$', fontsize=14, ha='center', va='center')
ax.text(D_center[0], D_center[1], '$D$', fontsize=14, ha='center', va='center')
ax.text(sum_center[0], sum_center[1], '$C + D$', fontsize=14, ha='center', va='center')
ax.set_xlim(-1, 5)
ax.set_ylim(-1, 4)
ax.axis('off')

# iii) Intersection
ax = axes[2]
C_poly = sg.Polygon(np.column_stack((1.5 + 1.2*np.cos(t) + 0.3*np.cos(2*t), 2.0 + 1.5*np.sin(t))))
D_poly = sg.Polygon(np.column_stack((2.5 + 1.5*np.cos(t), 1.5 + 1.2*np.sin(t))))
intersection = C_poly.intersection(D_poly)

x_c, y_c = C_poly.exterior.xy
x_d, y_d = D_poly.exterior.xy
ax.plot(x_c, y_c, 'k-', lw=1.5)
ax.plot(x_d, y_d, 'k-', lw=1.5)

if not intersection.is_empty:
    x_i, y_i = intersection.exterior.xy
    ax.fill(x_i, y_i, color='gray', alpha=0.6)

ax.text(0.8, 2.5, '$C$', fontsize=14)
ax.text(3.5, 2.0, '$D$', fontsize=14)
ax.set_xlim(-0.5, 4.5)
ax.set_ylim(-0.5, 4.5)
ax.axis('off')

plt.tight_layout()
plt.savefig('assets/fig_b2.png', dpi=300, bbox_inches='tight')
plt.close()

# -----------------
# Fig B.3 Cones
# -----------------
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Non-convex cone (intersecting regions)
ax = axes[0]
ax.fill_between([0, 1, 2], [0, 4, 8], [0, 8, 4], color='gray', alpha=0.6, edgecolor='k')
ax.fill_between([0, -1, -2], [0, -4, -8], [0, -8, -4], color='gray', alpha=0.6, edgecolor='k')
ax.plot(0, 0, 'ko', markersize=5)
ax.text(0.5, 0, '0', fontsize=14)
ax.set_xlim(-3, 3)
ax.set_ylim(-4, 4)
ax.set_title('Not convex', y=-0.2)
ax.axis('off')

# Non-convex cone (missing wedge)
ax = axes[1]
th1 = np.linspace(0.1, 1.4, 100)
x1 = np.concatenate(([0], 3*np.cos(th1), [0]))
y1 = np.concatenate(([0], 3*np.sin(th1), [0]))
ax.fill(x1, y1, color='gray', alpha=0.6, edgecolor='k')

th2 = np.linspace(0.6, 0.9, 50)
x2 = np.concatenate(([0], 3*np.cos(th2), [0]))
y2 = np.concatenate(([0], 3*np.sin(th2), [0]))
ax.fill(x2, y2, color='white', edgecolor='k') 

ax.plot(0, 0, 'ko', markersize=5)
ax.text(0.2, -0.2, '0', fontsize=14)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.5, 3.5)
ax.set_title('Not convex', y=-0.2)
ax.axis('off')

# Convex cone
ax = axes[2]
th = np.linspace(0.3, 1.2, 100)
x = np.concatenate(([0], 3*np.cos(th), [0]))
y = np.concatenate(([0], 3*np.sin(th), [0]))
ax.fill(x, y, color='gray', alpha=0.6, edgecolor='k')
ax.plot(0, 0, 'ko', markersize=5)
ax.text(0.2, -0.2, '0', fontsize=14)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.5, 3.5)
ax.set_title('Convex', y=-0.2)
ax.axis('off')

plt.tight_layout()
plt.savefig('assets/fig_b3.png', dpi=300, bbox_inches='tight')
plt.close()
