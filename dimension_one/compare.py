"""
The supersingular 2-isogeny graph over F_{p^2}, computed three ways and compared:

    - Sage: E0.isogeny_ell_graph(2, directed=True), vertices labelled by j-invariant,
    - modular polynomials: a BFS on j-invariants using the roots of Phi2 (modular_polynomial_graph.py),
    - endomorphism ring: the g = 1 case of the generic lattice search in ../endomorphism_graph/, a BFS on right
      ideals of O = End(E0) (lattices of rank 4) through the two torsion.

Checks:
    - the modular polynomial graph has exactly Sage's labelled edges (both are labelled by j-invariants),
    - the endomorphism ring graph is isomorphic to Sage's as a directed multigraph, by an isomorphism which also
      respects |Aut| (so j = 0 and j = 1728 are sent to ideals with the right unit groups),
    - both satisfy Eichler's mass formula  sum 1 / |Aut| = (p - 1) / 24.

Usage: sage -python compare.py p [p ...]
"""

import os
import sys

from sage.all import GF, QQ, EllipticCurve

import modular_polynomial_graph

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "endomorphism_graph"))
from graph import isogeny_graph, mass  # noqa: E402  (the generic search, here with g = 1)
from quaternions import EndomorphismRing  # noqa: E402


def sage_graph(p):
    """
    Sage's directed 2-isogeny multigraph of the supersingular curves over F_{p^2}, vertices labelled by j
    """
    E0 = EllipticCurve(GF(p**2, "z"), [1, 0])
    return E0.isogeny_ell_graph(2, directed=True, label_by_j=True)


def labelled_edges(G):
    """
    The sorted list of edges (with multiplicity) as pairs of strings
    """
    return sorted((str(a), str(b)) for a, b, _ in G.edges(sort=False))


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


def compare(p):
    eichler = QQ(p - 1) / 24
    H = sage_graph(p)

    # Modular polynomials: same labelled graph as Sage
    js, M = modular_polynomial_graph.isogeny_graph(p)
    m_mass = modular_polynomial_graph.mass(js)

    # Endomorphism ring (the generic search with g = 1): isomorphic to Sage, respecting |Aut|
    vertices, Q = isogeny_graph(EndomorphismRing(p), 1)
    q_mass = mass(vertices)
    sage_aut = lambda j: modular_polynomial_graph.automorphism_count(GF(p**2, "z")(j))
    same_graph = isomorphic_with_colours(Q, vertices.automorphism_count, H, sage_aut)

    print(f"p = {p}: {H.num_verts()} vertices, {H.num_edges()} edges")
    print(f"  modular polynomial = Sage (labelled edges)        : {labelled_edges(M) == labelled_edges(H)}")
    print(f"  endomorphism ring  = Sage (isomorphic, with |Aut|) : {same_graph}")
    print(f"  {f'mass = (p - 1)/24 = {eichler}':<50}: modular polynomial {m_mass == eichler}, "
          f"endomorphism ring {q_mass == eichler}")


if __name__ == "__main__":
    for p in map(int, sys.argv[1:]):
        compare(p)
