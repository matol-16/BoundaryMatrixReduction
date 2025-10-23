import matplotlib.pyplot as plt
from filtration_generation import *
from SparseBoundaryMatrix import SparseBoundaryMatrix
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


# This file allows to test the algorithms on the standard filtrations provided in filtration_generation

EXPECTED_BETTI = {
    'sphere_1d': {0: 1, 1: 1, 2: 0},
    'sphere_2d': {0: 1, 1: 0, 2: 1},
    'sphere_3d': {0: 1, 1: 0, 2: 0, 3: 1},
    'ball_1d': {0: 1, 1: 0},
    'ball_2d': {0: 1, 1: 0, 2: 0},
    'ball_3d': {0: 1, 1: 0, 2: 0, 3: 0},
    'moebius': {0: 1, 1: 1, 2: 0},
    'torus': {0: 1, 1: 2, 2: 1},
    'klein': {0: 1, 1: 2, 2: 1},
    'projective': {0: 1, 1: 1, 2: 1}
}

def count_infinite_intervals(barcode):
    """
    Count infinite intervals per dim
    Return a dictionnary {dimension: count}
    """
    infinite_by_dim = {}
    for (birth, death, dim) in barcode:
        if death == float('inf'):
            infinite_by_dim[dim] = infinite_by_dim.get(dim, 0) + 1
    return infinite_by_dim

def compute_euler_characteristic(betti_numbers):
    """
    Calcul χ = Σ (-1)^p * β_p
    """
    chi = 0
    max_dim = max(betti_numbers.keys()) if betti_numbers else 0
    for p in range(max_dim + 1):
        beta_p = betti_numbers.get(p, 0)
        chi += ((-1) ** p) * beta_p
    return chi

def test_filtration(filename, space_name):
    """
    Test a filtration and compare with the results
    """
    print(f"\n{'='*60}")
    print(f"Test: {space_name}")
    print(f"Fichier: {filename}")
    print(f"{'='*60}")

    # Calcul barcode
    matrix = SparseBoundaryMatrix()
    matrix.from_simplices(filename)
    matrix.gaussian_elimination()
    barcode = matrix.get_barcode()


    # Analysis of results
    print(f"\nBarcode complet ({len(barcode)} intervalles):")
    for (birth, death, dim) in sorted(barcode, key=lambda x: (x[2], x[0])):
        death_str = "inf" if death == float('inf') else f"{death:.1f}"
        print(f"  Dim {dim}: [{birth:.1f}, {death_str})")

    # Count infinites  intervals  (= Betti number)
    betti_computed = count_infinite_intervals(barcode)
    print(f"\nNombres de Betti calculés:")
    max_dim = max(betti_computed.keys()) if betti_computed else 0
    for dim in range(max_dim + 1):
        beta = betti_computed.get(dim, 0)
        print(f"  β_{dim} = {beta}")

    # Compare with theorical values
    if space_name in EXPECTED_BETTI:
        expected = EXPECTED_BETTI[space_name]
        print(f"\nNombres de Betti attendus:")
        for dim in sorted(expected.keys()):
            print(f"  β_{dim} = {expected[dim]}")
            
import numpy as np
import matplotlib.pyplot as plt

shapes = {
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

expected_total = {
    "Sphere S^2": 8,
    "Torus T^2": 29,
    "Ball B^2": 4,
    "Klein bottle": 29,
    "Möbius strip": 13,
    "Projective plane $\mathbb{RP}^2$": 17,
}

def count_intervals(d): return sum(len(v) for v in d.values())
def betti_numbers(d): return {k: sum(np.isinf(e) for _, e in d.get(k, [])) for k in (0, 1, 2)}

# Plot persistent barcodes
def plot_barcode(ax, intervals, title, xlim=None):
    dims = [0, 1, 2]
    births = [b for d in dims for (b, _) in intervals.get(d, [])]
    finite_deaths = [e for d in dims for (_, e) in intervals.get(d, []) if np.isfinite(e)]
    xmin = min(births) - 0.5 if births else 0
    xmax_auto = max(finite_deaths) + 0.5 if finite_deaths else max(births) + 1
    xmin, xmax = xlim if xlim else (xmin, xmax_auto)
    color = "navy"

    for d in dims:
        bars = intervals.get(d, [])
        n = max(1, len(bars))
        step = 0.9 / n
        for i, (b, e) in enumerate(bars):
            y = d + (i + 0.5) * step
            ax.plot([b, e if np.isfinite(e) else xmax], [y, y], lw=2, color=color)
            if not np.isfinite(e):
                ax.annotate("", xy=(xmax, y), xytext=(xmax - 0.01, y),
                            arrowprops=dict(arrowstyle="->", lw=2, color=color))
    β = betti_numbers(intervals)
    ax.set_title(f"{title}\n$\\beta_0={β[0]},\\ \\beta_1={β[1]},\\ \\beta_2={β[2]}$", fontsize=11)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.1, 3.1)
    ax.set_yticks([0.45, 1.45, 2.45])
    ax.set_yticklabels([r"$H_0$", r"$H_1$", r"$H_2$"])
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    ax.set_xlabel("Filtration value")

            
            
if __name__ == "__main__":
    
    generate_filtrations = False # Set True if not already computed
    
    plot_filtrations= True #Set true to plot the filtrations
    
    plot_barcodes = True # True if you intend to plot the barcodes
    
    if generate_filtrations:
        # ---------------------------
        #  Generate filtrations
        # ---------------------------
        print("=== Generate filtrations ===")
        generate_sphere_filtration(1, "output/sphere_1d.txt")
        generate_sphere_filtration(2, "output/sphere_2d.txt")
        generate_ball_filtration(2, "output/ball_2d.txt")
        generate_torus_filtration("output/torus.txt")
        generate_klein_filtration("output/klein.txt")
        generate_moebius_filtration("output/moebius.txt")
        generate_projective_plane_filtration("output/projective.txt")

    # ---------------------------
    #  tests
    # ---------------------------
    test_filtration("output/sphere_1d.txt", "sphere_1d")
    test_filtration("output/sphere_2d.txt", "sphere_2d")
    test_filtration("output/ball_2d.txt", "ball_2d")
    test_filtration("output/torus.txt", "torus")
    test_filtration("output/klein.txt", "klein")
    test_filtration("output/moebius.txt", "moebius")
    test_filtration("output/projective.txt", "projective")
    
    if plot_filtrations:
        order = [
            "Sphere S^2",
            "Torus T^2",
            "Ball B^2",
            "Klein bottle",
            "Möbius strip",
            "Projective plane $\mathbb{RP}^2$",
        ]

        fig = plt.figure(figsize=(14, 18))
        outer = GridSpec(len(order), 3, figure=fig, wspace=0.25, hspace=0.4)

        for i, name in enumerate(order):

            V, E, F = geoms[name]()

            ax_filtr = fig.add_subplot(outer[i, 0])
            plot_progressive_filtration(ax_filtr, V, E, F, title=name)

            ax_bar = fig.add_subplot(outer[i, 1])
            plot_barcode(ax_bar, SHAPES_BARCODE[name], title=name)

            ax_diag = fig.add_subplot(outer[i, 2])
            plot_persistence_diagram(ax_diag, SHAPES_BARCODE[name])


            fig.text(0.01, ax_bar.get_position().y1 + 0.01, name, fontsize=11, weight='bold')
        plt.savefig("plots/filtration_plots.png")
        plt.show()
    
    if plot_barcodes:
        for name, data in shapes.items():
            assert count_intervals(data) == expected_total[name]

        fig, axes = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
        axes = axes.ravel()
        order = [
            "Sphere S^2",
            "Torus T^2",
            "Ball B^2",
            "Klein bottle",
            "Möbius strip",
            "Projective plane $\mathbb{RP}^2$",
        ]

        for ax, name in zip(axes, order):
            plot_barcode(ax, shapes[name], name)
        
        plt.savefig("plots/barcodes.png")
        plt.show()
        