"""
The supersingular 2-isogeny graph over F_{p^2} from the classical modular polynomials.
The main reason to do this is if you squint hard enough what we do with Richelot in
dim two will look similar.

For the graph we build, each vertex is a j-invarant and the edge is some 2-isogeny
(which we don't know explicitly). We find vertices instead by computing roots of the
ell = 2 modular polynomial which satisfies Phi_ell(j1, j2) = 0 iff curves E1 and E2
are ell-isogenous.

For the two isogeny graph, we use

    Phi2(X, Y) = X^3 + Y^3 - X^2 Y^2 + 1488 (X^2 Y + X Y^2) - 162000 (X^2 + Y^2) + 40773375 X Y
                 + 8748000000 (X + Y) - 157464000000000.

and so substituting in X = j1, the three remaining roots are the three 2-isogenous curves.
Additionally, we can optimise this by substituting in j1, j2 for ell-isogenous curves and
the other two roots are the non-backtracking verticies which makes things even faster.

So, the pattern is:

    1. Pick a starting curve and compute the three outgoing edges with `cubic_neighbours()`
    2. From each new vertex remember the previous vertex and find the remaining two with `quadratic_neighbours()`
    3. Perform a BFS on this graph storing the vertex (j-invarant) and the edge which is the tuple of two verticies
    4. From this BFS data we make the SageMath graph with DiGraph
"""

import time

from sage.all import GF, QQ, DiGraph, PolynomialRing


def phi2(X, Y):
    """
    The classical modular polynomial of level two
    """
    return (
        X**3
        + Y**3
        - X**2 * Y**2
        + 1488 * (X**2 * Y + X * Y**2)
        - 162000 * (X**2 + Y**2)
        + 40773375 * X * Y
        + 8748000000 * (X + Y)
        - 157464000000000
    )


def cubic_neighbours(j):
    """
    The three roots of Phi2(j, Y), with multiplicity
    """
    R = PolynomialRing(j.parent(), "Y")
    f = phi2(j, R.gen())
    roots = [r for r, m in f.roots() for _ in range(m)]
    assert len(roots) == 3, "supersingular j-invariants have all their neighbours in F_{p^2}"
    return roots


def quadratic_neighbours(j, j_prev):
    """
    The two roots of Phi2(j, Y) / (Y - j_prev) = Y^2 + b Y + c, with multiplicity
    """
    b = -(j**2) + 1488 * j + j_prev - 162000
    c = (
        j_prev**2
        - j**2 * j_prev
        + 1488 * (j**2 + j * j_prev)
        + 40773375 * j
        - 162000 * j_prev
        + 8748000000
    )
    d = (b**2 - 4 * c).sqrt()
    return [(-b + d) / 2, (-b - d) / 2]


def neighbours(j, j_prev=None):
    """
    The codomains of the three 2-isogenies out of a curve with j-invariant j.
    If we know a neighbour j_prev already, we only need to solve a quadratic.
    """
    if j_prev is None:
        return cubic_neighbours(j)
    return [j_prev] + quadratic_neighbours(j, j_prev)


# ============================================================
# Graph search
# ============================================================


def isogeny_graph(p, verbose=False):
    """
    Breadth first search from j = 1728 (E0 : y^2 = x^3 + x) over F_{p^2}.

    Returns the list of j-invariants and the directed multigraph with one
    edge for every subgroup of order two.
    """
    F = GF(p**2, "z")
    start = F(1728)

    parent = {start: None}  # the vertex each j-invariant was discovered from
    edges = []
    queue, done = [start], set()
    t0 = time.time()
    while queue:
        j = queue.pop(0)
        if j in done:
            continue
        done.add(j)
        for j_new in neighbours(j, parent[j]):
            assert phi2(j, j_new) == 0
            edges.append((j, j_new))
            if j_new not in parent:
                parent[j_new] = j
                queue.append(j_new)
        if verbose:
            print(f"  processed {len(done)} / found {len(parent)}  ({time.time() - t0:.1f}s)")

    return list(parent), DiGraph(edges, multiedges=True, loops=True)


def automorphism_count(j):
    """
    |Aut(E)| for a supersingular curve with j-invariant j (p > 3): 6 for j = 0, 4 for j = 1728 and 2 (just -1)
    otherwise. All of these automorphisms are defined over F_{p^2}.
    """
    if j == 0:
        return 6
    if j == 1728:
        return 4
    return 2


def mass(vertices):
    """
    sum over the j-invariants of 1 / |Aut|, to compare with Eichler's mass formula (p - 1) / 24
    """
    return sum(QQ(1) / automorphism_count(j) for j in vertices)
