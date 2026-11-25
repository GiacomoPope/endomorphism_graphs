"""
The type I graph (g = 4: type I Weil structures and the middle steps between them) next to the supersingular
2-isogeny graph of elliptic curves over F_{p^2}, prime by prime, from the saved graphs in graph_data/
(sage -python storage.py p --type-one).

For each p and each graph we print:

    - n      : the number of vertices
    - mass   : sum 1 / |Aut| (for E/F_{p^2} this is (p - 1) / 24)
    - d      : the degree (3 and 30)
    - lambda : the largest non-trivial |eigenvalue| of the adjacency matrix, also divided by the Ramanujan bound
               2 sqrt(d - 1) (at most 1 for a Ramanujan graph)
    - diam   : the diameter
    - F_p    : the number of vertices defined over F_p
    - loops  : the fraction of kernels from a vertex back to itself

The adjacency matrices are not symmetric (vertices with extra automorphisms have fewer kernels back), but they are
reversible, m(v, w) |Aut w| = m(w, v) |Aut v|, so the eigenvalues are real.

Usage: sage -python compare_type_one.py
"""

import glob
import json
import os
import re

import numpy
from sage.all import GF, EllipticCurve, Graph

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "graph_data")


def statistics(n, edges, mass, fp):
    """
    edges: (v, w, multiplicity) of a regular directed multigraph on range(n)
    """
    A = numpy.zeros((n, n))
    for v, w, m in edges:
        A[v, w] += m
    d = int(round(A.sum(axis=1)[0]))
    eigenvalues = sorted(numpy.linalg.eigvals(A).real, key=abs, reverse=True)
    second = abs(eigenvalues[1]) if n > 1 else 0.0
    U = Graph([list(range(n)), [(v, w) for v, w, m in edges if v != w]], format="vertices_and_edges")
    return {
        "n": n,
        "mass": mass,
        "d": d,
        "lambda": second,
        "ratio": second / (2 * (d - 1) ** 0.5),
        "diam": U.diameter() if n > 1 else 0,
        "F_p": fp,
        "loops": sum(m for v, w, m in edges if v == w) / (n * d),
    }


def elliptic(p):
    """
    The supersingular 2-isogeny graph of E/F_{p^2}, from Sage
    """
    F = GF(p**2, "z")
    G = EllipticCurve(F, [1, 0]).isogeny_ell_graph(2, directed=True, label_by_j=True)
    vertices = G.vertices(sort=True)
    index = {j: k for k, j in enumerate(vertices)}
    edges = [(index[a], index[b], 1) for a, b, _ in G.edges(sort=False)]
    aut = lambda j: 6 if F(j) == 0 else (4 if F(j) == 1728 else 2)
    from sage.all import QQ

    mass = sum(QQ(1) / aut(j) for j in vertices)
    fp = sum(1 for j in vertices if F(j) ** p == F(j))
    return statistics(len(vertices), edges, mass, fp)


def type_one(path):
    data = json.load(open(path))
    from sage.all import QQ

    mass = QQ(data["mass"])
    fp = sum(1 for v in data["vertices"] if v["over_Fp"])
    return statistics(len(data["vertices"]), data["edges"], mass, fp)


def row(name, s):
    return (f"  {name:10s} n = {s['n']:5d}  mass = {str(s['mass']):>16s} ({float(s['mass']):9.4f})  d = {s['d']:3d}  "
            f"lambda = {s['lambda']:7.2f} ({s['ratio']:.2f} x Ramanujan)  diam = {s['diam']}  "
            f"F_p = {s['F_p']:4d} ({s['F_p'] / s['n']:.0%})  loops = {s['loops']:.1%}")


if __name__ == "__main__":
    files = glob.glob(os.path.join(DATA, "lattice_g4_p*_ell2_type1.json"))
    for path in sorted(files, key=lambda f: int(re.search(r"_p(\d+)_", f).group(1))):
        p = int(re.search(r"_p(\d+)_", path).group(1))
        E, M = elliptic(p), type_one(path)
        print(f"p = {p}")
        print(row("E/F_p^2", E))
        print(row("type I", M))
        print(f"  type I / E: vertices {M['n'] / E['n']:.2f}, mass {float(M['mass'] / E['mass']):.4f}, "
              f"mass / p^3 = {float(M['mass']) / p**3:.3e}")
