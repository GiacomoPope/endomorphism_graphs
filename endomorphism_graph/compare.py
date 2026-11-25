"""
The superspecial isogeny graph from the endomorphism ring of E0, compared with independent computations:

    - g = 1, any ell: Sage's E0.isogeny_ell_graph(ell, directed=True),
    - g = 2, ell = 2: the Richelot / Velu / gluing graph of dimension_two/richelot.py,
    - g = 2, 3, ell odd: no independent graph, but the vertices must be the same isomorphism classes as for ell = 2.

Checks, for each p:
    - the directed multigraphs are isomorphic, by an isomorphism respecting |Aut| (and, for g = 2, sending products
      to products),
    - sum 1 / |Aut| matches the mass formula,
    - for g = 2, the number of products is n(n + 1) / 2 for n supersingular j-invariants,
    - for odd ell: the same vertices as the 2-isogeny graph, prod (ell^i + 1) kernels out of every vertex, and every
      edge balanced by its dual: m(v, w) |Aut w| = m(w, v) |Aut v|.

Usage: sage -python compare.py p [p ...] [--ell 3]
"""

import os
import sys
from collections import Counter

from sage.all import GF, EllipticCurve

from graph import isogeny_graph, mass, mass_formula
from quaternions import EndomorphismRing
from torsion import lagrangian_count

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dimension_two"))
import richelot  # noqa: E402


def isomorphic_with_colours(G, colour_G, H, colour_H):
    """
    Whether there is an isomorphism G -> H sending each vertex to one of the same colour: compare canonical forms
    computed with the vertices partitioned by colour (in the same order on both sides)
    """
    def canonical(X, colour):
        classes = {}
        for v in X.vertices(sort=False):
            classes.setdefault(colour(v), []).append(v)
        keys = sorted(classes)
        return keys, [len(classes[k]) for k in keys], X.canonical_label(partition=[classes[k] for k in keys]).edges(sort=True)

    return canonical(G, colour_G) == canonical(H, colour_H)


def compare_dimension_one(p, ell=2):
    vertices, G = isogeny_graph(EndomorphismRing(p), 1, ell)
    E0 = EllipticCurve(GF(p**2, "z"), [1, 0])
    H = E0.isogeny_ell_graph(ell, directed=True, label_by_j=True)
    F = GF(p**2, "z")
    sage_aut = lambda j: richelot.elliptic_automorphism_count(F(j))

    print(f"p = {p}, g = 1, ell = {ell}: {G.num_verts()} vertices")
    print(f"  isomorphic to Sage (with |Aut|)  : {isomorphic_with_colours(G, vertices.automorphism_count, H, sage_aut)}")
    print(f"  mass = (p - 1) / 24             : {mass(vertices) == mass_formula(p, 1)}")


def compare_dimension_two(p):
    vertices, G = isogeny_graph(EndomorphismRing(p), 2)
    r_vertices, H = richelot.isogeny_graph(p)
    products = [v for v in range(len(vertices.lattices)) if vertices.is_decomposable(v)]
    n = richelot.supersingular_count(p)

    ours = lambda v: (vertices.is_decomposable(v), vertices.automorphism_count(v))
    theirs = lambda v: (isinstance(r_vertices[v], richelot.Product), richelot.automorphism_count(r_vertices[v]))

    print(f"p = {p}, g = 2: {G.num_verts()} vertices, {len(products)} products")
    print(f"  isomorphic to Richelot (with products, |Aut|) : {isomorphic_with_colours(G, ours, H, theirs)}")
    print(f"  mass = (p - 1)(p^2 + 1) / 5760                : {mass(vertices) == mass_formula(p, 2)}")
    print(f"  products = n(n + 1) / 2 = {n * (n + 1) // 2:<20}: {len(products) == n * (n + 1) // 2}")


def compare_odd_ell(p, g, ell):
    ring = EndomorphismRing(p)
    vertices, G = isogeny_graph(ring, g, ell)
    vertices_2, _ = isogeny_graph(ring, g, 2)
    images = [vertices_2.index(L) for L in vertices.lattices]
    same = not any(new for _, new, _ in images) and len({w for w, _, _ in images}) == len(vertices_2.lattices)
    same = same and len(vertices.lattices) == len(vertices_2.lattices)
    m = Counter((v, w) for v, w, _ in G.edges(sort=False))
    aut = vertices.automorphism_count
    regular = all(G.out_degree(v) == lagrangian_count(g, ell) for v in G.vertices(sort=False))
    balanced = all(m[(v, w)] * aut(w) == m[(w, v)] * aut(v) for v, w in m)

    print(f"p = {p}, g = {g}, ell = {ell}: {G.num_verts()} vertices")
    print(f"  same vertices as ell = 2                 : {same}")
    print(f"  {f'{lagrangian_count(g, ell)} kernels out of every vertex':<41}: {regular}")
    print(f"  m(v, w) |Aut w| = m(w, v) |Aut v|        : {balanced}")
    print(f"  mass formula                             : {mass(vertices) == mass_formula(p, g)}")


if __name__ == "__main__":
    args = sys.argv[1:]
    ell = 2
    if "--ell" in args:
        k = args.index("--ell")
        ell = int(args[k + 1])
        args = args[:k] + args[k + 2:]
    for p in map(int, args):
        compare_dimension_one(p, ell)
        if ell == 2:
            compare_dimension_two(p)
        else:
            compare_odd_ell(p, 2, ell)
