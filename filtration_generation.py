from itertools import combinations
import numpy as np
from matplotlib.gridspec import  GridSpecFromSubplotSpec
import matplotlib.pyplot as plt

### This py file contains the functions to generate filtration for standard objects.


def generate_sphere_filtration(d, filename):
    """
  Generate a filtration for the d-sphere S^d,
  represented as the boundary of the standard (d+1)-simplex.
  Store it in filename
  """

    with open(filename, 'w') as f:
        n_vertices = d + 2  # d+2 vertices pour S^d
        time = 1.0
        # Vertices (0-simplices)
        for i in range(n_vertices):
            f.write(f"{time} 0 {i}\n")
        time += 1.0

        # Generate all simplices of dimension k
        def generate_simplices(vertices, dim, max_dim):
            nonlocal time
            if dim > max_dim:
                return


            for simplex in combinations(range(n_vertices), dim + 1):
                if dim < d + 1:  # all simplices except the full
                    f.write(f"{time} {dim} " + " ".join(map(str, simplex)) + "\n")
                    time += 1.0

        # Edges (1-simplices)
        for edge in combinations(range(n_vertices), 2):
            f.write(f"{time} 1 " + " ".join(map(str, edge)) + "\n")
            time += 1.0

        # Higher dimensional simplices up to dimension d (boundary of (d+1)-simplex)
        for dim in range(2, d + 1):
            for simplex in combinations(range(n_vertices), dim + 1):
                f.write(f"{time} {dim} " + " ".join(map(str, simplex)) + "\n")
                time += 1.0

def generate_ball_filtration(d, filename):
    """
    Generate filtration d-Ball B^d
    Represented as the boundary of the standard (d+1)-simplex
    """
    with open(filename, 'w') as f:
        n_vertices = d + 1
        time = 1.0
        for i in range(n_vertices):
            f.write(f"{time} 0 {i}\n")
        time += 1.0
        for dim in range(1, d + 1):
            for simplex in combinations(range(n_vertices), dim + 1):
                f.write(f"{time} {dim} " + " ".join(map(str, simplex)) + "\n")
                time += 1.0



def generate_moebius_filtration(filename):
    # layout: bottom 0-1-2, top 3-4-5; twist: (2)~(3)
    with open(filename, 'w') as f:
        t = 1.0
        for v in range(6):                      # 0-simplices
            f.write(f"{t} 0 {v}\n")
        t += 1.0

        edges = {
            (0,1), (1,2), (3,4), (4,5),        # rows
            (0,3), (1,4), (2,5),               # verticals
            (2,3),                              # twist
            (0,4), (1,3), (1,5), (2,4)         # minimal diagonals
        }
        for a,b in sorted(tuple(sorted(e)) for e in edges):
            f.write(f"{t} 1 {a} {b}\n")         # 1-simplices
        t += 1.0

        triangles = [                           # 2-simplices (χ=0, β1=1)
            (0,1,4), (0,4,3),
            (1,2,5), (1,5,4),
            (2,3,4), (2,4,5)
        ]
        for a,b,c in triangles:
            f.write(f"{t} 2 {a} {b} {c}\n")
            t += 1.0


def generate_torus_filtration(filename):
    """
    Generate a filtration for a torus
    with a periodical grid
    """
    with open(filename, 'w') as f:
        time = 1.0


        for v in range(9):
            f.write(f"{time} 0 {v}\n")
        time += 1.0

        h_edges = [
            (0,1), (1,2), (2,0),
            (3,4), (4,5), (5,3),
            (6,7), (7,8), (8,6)
        ]
        for e in h_edges:
            f.write(f"{time} 1 {e[0]} {e[1]}\n")
            time += 1.0

        # Vertical edges
        v_edges = [
            (0,3), (3,6), (6,0),  # colonne 0
            (1,4), (4,7), (7,1),  # colonne 1
            (2,5), (5,8), (8,2)   # colonne 2 = colonne 0
        ]
        for e in v_edges:
            f.write(f"{time} 1 {e[0]} {e[1]}\n")
            time += 1.0

        # Diagonal edges
        d_edges = [
            (0,4), (1,5), (2,3),
            (3,7), (4,8), (5,6),
            (6,1), (7,2), (8,0)
        ]
        for e in d_edges:
            f.write(f"{time} 1 {e[0]} {e[1]}\n")
            time += 1.0


        triangles = [
            (0,1,4), (0,4,3),
            (1,2,5), (1,5,4),
            (2,0,3), (2,3,5),
            (3,4,7), (3,7,6),
            (4,5,8), (4,8,7),
            (5,3,6), (5,6,8),
            (6,7,1), (6,1,0),
            (7,8,2), (7,2,1),
            (8,6,0), (8,0,2)
        ]
        for t in triangles:
            f.write(f"{time} 2 {t[0]} {t[1]} {t[2]}\n")
            time += 1.0



def generate_klein_filtration(filename):
    """
    Generates a filter for the correct Klein bottle.
   Matching:
   (x,0) ~ (x,1)
   (0,y) ~ (1,1-y)
   Uses a 3x3 grid with adapted periodic triangulation
    """
    from itertools import product, combinations

    with open(filename, 'w') as f:
        t = 1.0

        # vertices
        vertices = [(i, j) for i, j in product(range(3), repeat=2)]
        vertex_index = {v: idx for idx, v in enumerate(vertices)}

        # 0-simplices
        for v in vertex_index.values():
            f.write(f"{t} 0 {v}\n")
        t += 1

        # 1-simplices
        def mod3(x): return x % 3
        edges = set()

        for i in range(3):
            for j in range(3):
                v = vertex_index[(i, j)]
                # droite sans twist
                edges.add(tuple(sorted((v, vertex_index[(mod3(i+1), j)]))))
                # haut sans twist
                edges.add(tuple(sorted((v, vertex_index[(i, mod3(j+1))]))))

        for i in range(3):
            for j in range(3):
                v0 = vertex_index[(i, j)]
                v1 = vertex_index[(mod3(i+1), j)]
                v2 = vertex_index[(i, mod3(j+1))]
                edges.add(tuple(sorted((v0, v1))))
                edges.add(tuple(sorted((v0, v2))))
                edges.add(tuple(sorted((v1, v2))))

        for e in edges:
            f.write(f"{t} 1 {e[0]} {e[1]}\n")
        t += 1

        # 2-simplices
        for i in range(3):
            for j in range(3):
                # normal triangle
                v0 = vertex_index[(i, j)]
                v1 = vertex_index[(mod3(i+1), j)]
                v2 = vertex_index[(i, mod3(j+1))]
                v3 = vertex_index[(mod3(i+1), mod3(j+1))]
                f.write(f"{t} 2 {v0} {v1} {v2}\n")
                f.write(f"{t} 2 {v1} {v2} {v3}\n")
        t += 1


def generate_projective_plane_filtration(filename):
    """
    Projective Plane RP^2
    Configuration: Star with center + outer pentagon
    """
    with open(filename, 'w') as f:
        time = 1.0

        # vertices (1,2,3,4,5)
        for v in range(6):
            f.write(f"{time} 0 {v}\n")
        time += 1.0

        # edges
        edges = []
        for i in range(1, 6):
            edges.append((i, i % 5 + 1))  # 1-2, 2-3, 3-4, 4-5, 5-1
        for i in range(1, 6):
            edges.append((0, i))

        edges.append((1, 3))
        edges.append((1, 4))
        edges.append((2, 4))
        edges.append((2, 5))
        edges.append((3, 5))
        for e in edges:
            f.write(f"{time} 1 {e[0]} {e[1]}\n")
            time += 1.0
        triangles = [
            (0, 1, 2),
            (0, 2, 3),
            (0, 3, 4),
            (0, 4, 5),
            (0, 5, 1),
            (1, 2, 4),
            (2, 3, 5),
            (3, 4, 1),
            (4, 5, 2),
            (5, 1, 3)
        ]
        for t in triangles:
            f.write(f"{time} 2 {t[0]} {t[1]} {t[2]}\n")
            time += 1.0
            

### This part of the codes allows to plot the filtrations: 

SHAPES_BARCODE = {
    "Sphere S^2": {
        0: [(1.0, 2.0), (1.0, 3.0), (1.0, 4.0), (1.0, np.inf)],
        1: [(5.0, 8.0), (6.0, 9.0), (7.0, 10.0)],
        2: [(11.0, np.inf)],
    },
    "Torus T^2": {
        0: [(1.0, 2.0), (1.0, 3.0), (1.0, 5.0), (1.0, 6.0),
            (1.0, 8.0), (1.0, 9.0), (1.0, 11.0), (1.0, 12.0), (1.0, np.inf)],
        1: [(4.0, np.inf), (7.0, 34.0), (10.0, 40.0), (13.0, np.inf),
            (14.0, 30.0), (15.0, 36.0), (16.0, 42.0), (17.0, 32.0),
            (18.0, 38.0), (19.0, 44.0), (20.0, 29.0), (21.0, 31.0),
            (22.0, 33.0), (23.0, 35.0), (24.0, 37.0), (25.0, 39.0),
            (26.0, 41.0), (27.0, 43.0), (28.0, 45.0)],
        2: [(46.0, np.inf)],
    },
    "Ball B^2": {
        0: [(1.0, 2.0), (1.0, 3.0), (1.0, np.inf)],
        1: [(4.0, 5.0)],
        2: [],
    },
    "Klein bottle": {
        0: [(1.0, 2.0)]*8 + [(1.0, np.inf)],
        1: [(2.0, 3.0)]*17 + [(2.0, np.inf), (2.0, np.inf)],
        2: [(3.0, np.inf)],
    },
    "Möbius strip": {
        0: [(1.0, 2.0)]*5 + [(1.0, np.inf)],
        1: [(2.0, 3.0), (2.0, 4.0), (2.0, 5.0), (2.0, 6.0),
            (2.0, 7.0), (2.0, 8.0), (2.0, np.inf)],
        2: [],
    },
    "Projective plane $\mathbb{RP}^2$": {
        0: [(1.0, 2.0), (1.0, 3.0), (1.0, 4.0), (1.0, 5.0),
            (1.0, 7.0), (1.0, np.inf)],
        1: [(6.0, 21.0), (8.0, 17.0), (9.0, 18.0), (10.0, 19.0),
            (11.0, 20.0), (12.0, np.inf), (13.0, 24.0), (14.0, 22.0),
            (15.0, 25.0), (16.0, 23.0)],
        2: [(26.0, np.inf)],
    },
}


def geom_sphere_S2():
    ang = np.linspace(0, 2*np.pi, 5)[:-1]
    R = 1.0
    V = np.c_[R*np.cos(ang), R*np.sin(ang)]
    faces = [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]
    edges = [(i,j) for i in range(4) for j in range(i+1,4)]
    return V, edges, faces

def geom_ball_B2():
    V = np.array([[0,0],[1,0],[0.5,0.866]])  # equilateral triangle
    edges = [(0,1),(1,2),(0,2)]
    faces = [(0,1,2)]
    return V, edges, faces

def geom_moebius():
    V = np.array([[0,0],[1,0],[2,0],[0,1],[1,1],[2,1]])
    edges = [(0,1),(1,2),(3,4),(4,5),(0,3),(1,4),(2,5),(2,3),
             (0,4),(1,3),(1,5),(2,4)]
    faces = [(0,1,4),(0,4,3),(1,2,5),(1,5,4),(2,3,4),(2,4,5)]
    return V, edges, faces

def geom_torus():
    V = np.array([[i,j] for j in range(3) for i in range(3)])
    edges = [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),
             (6,7),(7,8),(8,6),(0,3),(3,6),(6,0),
             (1,4),(4,7),(7,1),(2,5),(5,8),(8,2)]
    faces = [(0,1,4),(0,4,3),(1,2,5),(1,5,4),
             (3,4,7),(3,7,6),(4,5,8),(4,8,7)]
    return V, edges, faces

def geom_klein():
    V = np.array([[i,j] for j in range(3) for i in range(3)])
    def mod3(x): return x % 3
    edges = []
    for i in range(3):
        for j in range(3):
            v = 3*j + i
            edges += [(v, 3*j + mod3(i+1)), (v, 3*mod3(j+1) + i)]
    faces = []
    for i in range(3):
        for j in range(3):
            v0 = 3*j + i
            v1 = 3*j + mod3(i+1)
            v2 = 3*mod3(j+1) + i
            v3 = 3*mod3(j+1) + mod3(i+1)
            faces += [(v0,v1,v2),(v1,v2,v3)]
    return V, edges, faces

def geom_projective():
    C = np.array([[0,0]])
    P = np.array([[np.cos(2*np.pi*i/5), np.sin(2*np.pi*i/5)] for i in range(5)])
    V = np.vstack([C,P])
    edges = [(0,i) for i in range(1,6)] + [(i,(i%5)+1) for i in range(1,6)] \
            + [(1,3),(1,4),(2,4),(2,5),(3,5)]
    faces = [(0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,5,1),
             (1,2,4),(2,3,5),(3,4,1),(4,5,2),(5,1,3)]
    return V, edges, faces

geoms = {
    "Sphere S^2": geom_sphere_S2,
    "Ball B^2": geom_ball_B2,
    "Möbius strip": geom_moebius,
    "Torus T^2": geom_torus,
    "Klein bottle": geom_klein,
    "Projective plane $\mathbb{RP}^2$": geom_projective,
}


def betti_numbers(data):
    return {k: sum(np.isinf(e) for _, e in data.get(k, [])) for k in (0,1,2)}

def plot_barcode(ax, data, title):
    color = "navy"
    finite = [e for L in data.values() for (_,e) in L if np.isfinite(e)]
    births = [b for L in data.values() for (b,_) in L]
    xmax = (max(finite) + 0.5) if finite else (max(births)+1.0)
    for d in (0,1,2):
        bars = data.get(d, [])
        for i,(b,e) in enumerate(bars):
            y = d + 0.1 + i*0.07
            ax.plot([b, e if np.isfinite(e) else xmax], [y,y], lw=2, color=color)
            if not np.isfinite(e):
                ax.annotate("", xy=(xmax,y), xytext=(xmax-0.01,y),
                            arrowprops=dict(arrowstyle="->", lw=2, color=color))
    β = betti_numbers(data)
    ax.set_title(f"{title}\n$\\beta_0={β[0]},\\ \\beta_1={β[1]},\\ \\beta_2={β[2]}$", fontsize=9)
    ax.set_yticks([0.5,1.5,2.5]); ax.set_yticklabels([r"$H_0$",r"$H_1$",r"$H_2$"])
    ax.set_xlim(min(births)-0.5, xmax); ax.grid(axis="x", linestyle=":", alpha=0.5)
    ax.set_xlabel("Filtration value")

def plot_persistence_diagram(ax, data, title=""):
    finite = [(b,e) for L in data.values() for (b,e) in L if np.isfinite(e)]
    if finite:
        xs, ys = zip(*finite)
    else:
        xs, ys = [], []
    M = max(xs+ys) if (xs and ys) else 1.0
    ax.scatter(xs, ys, s=20)
    ax.plot([0,M], [0,M], 'k--', lw=1)
    ax.set_aspect('equal'); ax.set_xlim(0,M); ax.set_ylim(0,M)
    ax.set_title("Persistence diagram", fontsize=9)
    ax.set_xlabel("birth"); ax.set_ylabel("death")

def draw_stage(ax, V, edges, faces, stage):
    # stage 1: vertices + small discs
    # stage 2: + edges
    # stage 3: + faces
    x, y = V[:,0], V[:,1]
    ax.set_aspect('equal'); ax.axis('off')

    span = max(x.max()-x.min(), y.max()-y.min())
    r = {1: 0.08*span, 2: 0.18*span, 3: 0.28*span}[stage]
    for (xi,yi) in V:
        circ = plt.Circle((xi,yi), r, fill=False, ls='--', lw=1, alpha=0.6)
        ax.add_patch(circ)
    if stage >= 3:
        for (a,b,c) in faces:
            ax.fill(V[[a,b,c],0], V[[a,b,c],1], alpha=0.25, zorder=0)
    if stage >= 2:
        for (a,b) in edges:
            ax.plot(V[[a,b],0], V[[a,b],1], lw=1.5)
    ax.scatter(x,y, s=30)

def plot_progressive_filtration(ax_parent, V, edges, faces, title):
    ax_parent.set_title("Filtration (stages)", fontsize=9)
    ax_parent.axis('off')
    gs = GridSpecFromSubplotSpec(1, 3, subplot_spec=ax_parent.get_subplotspec(), wspace=0.05)
    for j,stage in enumerate([1,2,3]):
        ax = plt.subplot(gs[0, j])
        draw_stage(ax, V, edges, faces, stage)
        if j>0:

            for spine in ['left']:
                ax.spines[spine].set_visible(True)
                ax.spines[spine].set_linewidth(1)
                ax.spines[spine].set_linestyle('-')
        ax.set_title(f"t = {stage}", fontsize=8)
