"""
The superspecial (2,2)-isogeny graph over F_{p^2}, computed with genus two curves and products of elliptic curves.

There are two kinds of vertex, each with 15 outgoing (2,2)-isogenies:

    - a Jacobian J(C), C : y^2 = c prod(x - r_k), stored by its six roots r_k, all in F_{p^2}.
      The 15 isogenies correspond to the 15 ways of splitting the six roots into three pairs,
      f = G1 G2 G3 with quadratics G_k (Richelot).

      Write G_k = g_k0 + g_k1 x + g_k2 x^2 and let delta = det(g_kl).
        - If delta != 0 the codomain is the Jacobian of y^2 = delta H1 H2 H3, where H1 = G2' G3 - G2 G3'.
        - If delta = 0 the codomain is a product of elliptic curves E1 x E2 (split).

    - a product E1 x E2, each E : y^2 = (x - r1)(x - r2)(x - r3) given by its three roots. The 15 isogenies are
        - 9 kernels <P> x <Q>: a 2-isogeny on each factor (velu), giving another product,
        - 6 kernels {(P, psi(P))} for the bijections psi : E1[2] -> E2[2] (glue). If psi comes from an isomorphism
          E1 -> E2 the quotient is E1 x E2 again, otherwise it is a Jacobian.

Gluing and splitting use closed formulae from my unpublished note.

To label vertices we do the following:
    - For a Product, the label is the j-invariant of each of the factors
    - For a Jacobian, we use the Igusa-Clebsch invariants (I2 : I4 : I6 : I10), but implementing this
      I discovered there's edge cases for when I2 or I4 vanishes. I followed arXiv:2111.07270 to deal
      with this, but I think this should probably be added to SageMath...

The graph is then made by:

    1. Start from E0 x E0, E0 : y^2 = x^3 + x, as in dimension one
    2. For each vertex compute its 15 neighbours with `neighbours()`
    3. Perform a BFS on this graph storing each vertex by its key and the edges between keys
    4. From this BFS data we make the SageMath graph with DiGraph

Usage: sage -python richelot.py p [p ...]
"""

import sys
import time
from collections import Counter
from itertools import combinations, permutations

from sage.all import GF, QQ, DiGraph, EllipticCurve

# ============================================================
# Vertices and their keys
# ============================================================

class Jacobian:
    """
    J(C) for C : y^2 = c prod(x - r_k), given by its six Weierstrass points r_k (all in the base field). The leading
    coefficient c only changes the twist, so we never need it.

    A Weierstrass point may arrive at infinity (None), e.g. from a Richelot codomain of degree five; we then move all
    six by x -> 1 / (x - t), which does not change the curve, so that the stored roots are always finite.
    """

    def __init__(self, roots):
        if any(r is None for r in roots):
            F = next(r for r in roots if r is not None).parent()
            t = next(t for t in F if t not in roots)
            roots = [F(0) if r is None else 1 / (r - t) for r in roots]
        self.roots = tuple(roots)
        self.key = jacobian_key(self.roots)


class Product:
    """
    E1 x E2, each E : y^2 = (x - r1)(x - r2)(x - r3) given by its three roots
    """
    def __init__(self, roots1, roots2):
        self.roots = (tuple(roots1), tuple(roots2))
        self.key = product_key(*self.roots)



def pairings(roots):
    """
    The 15 ways to split six roots into three unordered pairs
    """
    if not roots:
        return [[]]
    first, rest = roots[0], roots[1:]
    out = []
    for k, partner in enumerate(rest):
        for others in pairings(rest[:k] + rest[k + 1 :]):
            out.append([(first, partner)] + others)
    return out


# ============================================================
# Fast dim two invariants for p > 5
# ============================================================

# index form of the 15 pairings, used by jacobian_key
ROOT_PAIRINGS = pairings(list(range(6)))

# The classical root formulae for the Igusa-Clebsch invariants sum over: the 15 pairings, and the 10 ways to split
# the six roots into two triples (with the 6 matchings between the triples for I6)
ROOT_TRIPLES = [(t, tuple(k for k in range(6) if k not in t)) for t in combinations(range(6), 3) if 0 in t]


def igusa_clebsch_invariants(roots):
    """
    The Igusa-Clebsch invariants (I2, I4, I6, I10) of y^2 = c prod(x - r_k), from the six roots. Writing
    (ij) = (r_i - r_j)^2:

        I2  = sum over the 15 pairings             (12)(34)(56),
        I4  = sum over the 10 splits into triples  (12)(23)(31) (45)(56)(64),
        I6  = sum over the 60 splits and matchings (12)(23)(31) (45)(56)(64) (14)(25)(36),
        I10 = prod over all pairs (ij).

    These are Sage's HyperellipticCurve(f).igusa_clebsch_invariants() rescaled by t = c^2 / 16 in each weight
    (I2, I4, I6, I10) -> (t I2, t^2 I4, t^3 I6, t^5 I10), so they define the same point of weighted projective space.
    We compute them ourselves because we always know the roots which means things are faster.
    """
    r = roots
    d = {}
    for i, j in combinations(range(6), 2):
        d[(i, j)] = d[(j, i)] = (r[i] - r[j]) ** 2
    I2 = sum(d[P[0]] * d[P[1]] * d[P[2]] for P in ROOT_PAIRINGS)
    I4 = I6 = 0
    for (a, b, e), (u, v, w) in ROOT_TRIPLES:
        triangles = d[(a, b)] * d[(b, e)] * d[(e, a)] * d[(u, v)] * d[(v, w)] * d[(w, u)]
        I4 += triangles
        for x, y, z in permutations((u, v, w)):
            I6 += triangles * d[(a, x)] * d[(b, y)] * d[(e, z)]
    I10 = 1
    for i, j in combinations(range(6), 2):
        I10 *= d[(i, j)]
    return I2, I4, I6, I10


def weighted_point_key(I2, I4, I6, I10):
    """
    We need absolute invariants for the key when we construct the graph, but apparently this isn't handled
    by absolute invariants in Sage because of edge cases when I2 or I4 is zero.

    I found the reference, "Genus 2 Curves in Small Characteristic" (arXiv:2111.07270, Section 1.2), which uses a
    result of Cardona, Quer, Nart and Pujolas for p > 5 such that we get absolute invairants from the following

        I2 != 0          : (I2^5 / I10, I2^3 I4 / I10, I2^2 I6 / I10),
        I2 = 0, I4 != 0  : (0, I4^5 / I10^2, I4 I6 / I10),
        I2 = I4 = 0      : (0, 0, I6^5 / I10^3).
    """
    assert I10 != 0, "repeated root"
    if I2 != 0:
        return (I2**5 / I10, I2**3 * I4 / I10, I2**2 * I6 / I10)
    if I4 != 0:
        return (0, I4**5 / I10**2, I4 * I6 / I10)
    return (0, 0, I6**5 / I10**3)


def jacobian_key(roots):
    """
    A hashable key for the curve y^2 = c prod(x - r_k): equal keys iff the curves are isomorphic over the algebraic
    closure
    """
    return ("jacobian",) + weighted_point_key(*igusa_clebsch_invariants(roots))


def product_key(roots1, roots2):
    """
    E1 x E2 is determined (as a principally polarised surface) by the unordered pair of j-invariants
    """
    return ("product",) + tuple(sorted([j_invariant(roots1), j_invariant(roots2)], key=str))


def j_invariant(roots):
    """
    The j-invariant of y^2 = (x - r1)(x - r2)(x - r3):

        j = 256 (r1^2 + r2^2 + r3^2 - r1 r2 - r2 r3 - r3 r1)^3 / delta^2,   delta = (r1 - r2)(r2 - r3)(r3 - r1).

    (Moving r1 to 0 gives the Legendre form with lambda = (r3 - r1) / (r2 - r1) and
    j = 256 (lambda^2 - lambda + 1)^3 / (lambda^2 (lambda - 1)^2); clearing denominators gives the formula above.)
    """
    r1, r2, r3 = roots
    delta = (r1 - r2) * (r2 - r3) * (r3 - r1)
    return 256 * (r1**2 + r2**2 + r3**2 - r1 * r2 - r2 * r3 - r3 * r1) ** 3 / delta**2


# ============================================================
# Neighbours of a Jacobian: Richelot isogenies
# ============================================================


def quadratic(r, s):
    """
    Coefficients (A, B, C) of A x^2 + B x + C = (x - r)(x - s), or of x - r when s is at infinity (None)
    """
    if s is None:
        return (0, 1, -r)
    return (1, -(r + s), r * s)


def quadratic_roots(A, B, C):
    """
    The two roots of A x^2 + B x + C (None for a root at infinity when A = 0)
    """
    if A == 0:
        return (-C / B, None)
    D = (B**2 - 4 * A * C).sqrt()
    return ((-B + D) / (2 * A), (-B - D) / (2 * A))


def richelot_codomain(pairs):
    """
    The Richelot isogeny for the splitting f = G1 G2 G3 given by three pairs of roots. Returns the three pairs of
    roots of H1, H2, H3 (the codomain is y^2 = delta H1 H2 H3, and these pairs give the dual isogeny), or None when
    delta = det(coefficients of the G_k) = 0.

    For G_j = A_j x^2 + B_j x + C_j, H1 = G2' G3 - G2 G3' has coefficients given by 2 x 2 minors:

        H1 = (A2 B3 - A3 B2) x^2 + 2 (A2 C3 - A3 C2) x + (B2 C3 - B3 C2)   (and cyclically).
    """
    G = [quadratic(r, s) for r, s in pairs]
    (A1, B1, C1), (A2, B2, C2), (A3, B3, C3) = G
    delta = A1 * (B2 * C3 - B3 * C2) - B1 * (A2 * C3 - A3 * C2) + C1 * (A2 * B3 - A3 * B2)
    if delta == 0:
        return None
    out = []
    for j, k in [(1, 2), (2, 0), (0, 1)]:
        (Aj, Bj, Cj), (Ak, Bk, Ck) = G[j], G[k]
        out.append(quadratic_roots(Aj * Bk - Ak * Bj, 2 * (Aj * Ck - Ak * Cj), Bj * Ck - Bk * Cj))
    return out


def split(pairs):
    """
    The codomain E1 x E2 of a Richelot isogeny with delta = 0, given the three pairs of roots.
    The closed formulae are from my unpublished note on gluing and splitting.
    """
    quadratics = [(-(r + s), r * s) for r, s in pairs]
    (a1, b1), (a2, b2), (a3, _) = quadratics
    if a1 != a2:
        D = ((a1 - a2) * (a1 * b2 - b1 * a2) + (b1 - b2) ** 2).sqrt()
        e = [D * (2 * (b1 - b2) - a * (a1 - a2)) for a, _ in quadratics]
        return even_split([2 * D**2 - t for t in e], [2 * D**2 + t for t in e])
    assert a3 == a1, "a_1 = a_2 forces a_3 = a_1"
    one = a1.parent()(1)
    return even_split([one, one, one], [b - a1**2 / 4 for _, b in quadratics])


def even_split(beta, gamma):
    """
    E1 : Y^2 = prod(beta_i X + gamma_i) and E2 : V^2 = prod(beta_i + gamma_i U), as the product of the curves with
    roots -gamma_i / beta_i and -beta_i / gamma_i (scaling y^2 = c g(x) to y^2 = g(x) only changes the twist)
    """
    assert all(t != 0 for t in beta + gamma), "a root at infinity"
    return Product([-g / b for b, g in zip(beta, gamma)], [-b / g for b, g in zip(beta, gamma)])


def richelot_neighbours(J):
    """
    The 15 codomains of Richelot isogenies out of J(C)
    """
    out = []
    for pairs in pairings(list(J.roots)):
        H = richelot_codomain(pairs)
        out.append(Jacobian([r for pair in H for r in pair]) if H is not None else split(pairs))
    return out


# ============================================================
# Neighbours of a product: Diagonal and gluing
# ============================================================


def velu(roots, k):
    """
    The roots of E / <(r_k, 0)> for E : y^2 = (x - r1)(x - r2)(x - r3).

    Moving r_k to 0, E : y^2 = x (x - u)(x - v) and the quotient by (0, 0) is y^2 = x (x^2 + 2 (u + v) x + (u - v)^2),
    whose roots are 0 and -(u + v) -+ 2 sqrt(u v) = -(sqrt(u) -+ sqrt(v))^2.
    """
    r = roots[k]
    u, v = [s - r for n, s in enumerate(roots) if n != k]
    w = (u * v).sqrt()
    return (r.parent()(0), -(u + v) + 2 * w, -(u + v) - 2 * w)


def glue(a, b):
    """
    The gluing (E1 x E2) / {(P, psi(P))} for E1 : y^2 = prod(x - a_i), E2 : y^2 = prod(x - b_i)
    The closed formulae are from my unpublished note on gluing and splitting.
    """
    a1, a2, a3 = a
    b1, b2, b3 = b
    RD = b1 * (a2 - a3) + b2 * (a3 - a1) + b3 * (a1 - a2)
    if RD == 0:
        return None
    # s1 = -delta_a / RD only scales H, so for the vertex we only need the six roots +- sqrt(alpha_i)
    roots = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        t = ((b[j] - b[k]) / (a[j] - a[k])).sqrt()
        roots += [t, -t]
    return roots


def product_neighbours(P):
    """
    The 15 codomains out of E1 x E2: 9 products of 2-isogenies, then the 6 gluings
    """
    roots1, roots2 = P.roots
    out = []
    for k in range(3):
        for l in range(3):
            out.append(Product(velu(roots1, k), velu(roots2, l)))
    for b in permutations(roots2):
        H = glue(roots1, b)
        # (E1 x E2) / {(P, phi(P))} for an isomorphism phi is E1 x E2 again
        out.append(Jacobian(H) if H is not None else P)
    return out


# ============================================================
# Graph search
# ============================================================


def neighbours(V):
    return product_neighbours(V) if isinstance(V, Product) else richelot_neighbours(V)


def isogeny_graph(p, verbose=False):
    """
    Breadth first search from E0 x E0. Returns a dictionary from keys to vertices and the directed multigraph with
    one edge for every (2,2)-isogeny (15 out of every vertex).
    """
    F = GF(p**2, "z")
    i = F(-1).sqrt()
    E0_roots = (F(0), i, -i)
    start = Product(E0_roots, E0_roots)

    vertices = {start.key: start}
    edges = []
    queue, done = [start.key], set()
    t0 = time.time()
    while queue:
        v = queue.pop(0)
        if v in done:
            continue
        done.add(v)
        for W in neighbours(vertices[v]):
            if W.key not in vertices:
                vertices[W.key] = W
                queue.append(W.key)
            edges.append((v, W.key))
        if verbose:
            print(f"  processed {len(done)} / found {len(vertices)}  ({time.time() - t0:.1f}s)")

    return vertices, DiGraph(edges, multiedges=True, loops=True)


# ============================================================
# Checks
# ============================================================


def supersingular_count(p):
    """
    Number of supersingular j-invariants over F_{p^2}
    """
    E0 = EllipticCurve(GF(p**2, "z"), [1, 0])
    return E0.isogeny_ell_graph(2, directed=True).num_verts()


def check_duals(vertices):
    """
    Richelot isogenies come in dual pairs: from the codomain y^2 = delta H1 H2 H3, the splitting into H1, H2, H3 leads
    back to the original curve. Returns whether this holds for every Jacobian to Jacobian edge.
    """
    for V in vertices.values():
        if isinstance(V, Product):
            continue
        for pairs in pairings(list(V.roots)):
            H = richelot_codomain(pairs)
            if H is None:
                continue
            back = richelot_codomain(H)
            if back is None or Jacobian([r for pair in back for r in pair]).key != V.key:
                return False
    return True


# ============================================================
# Automorphism groups and the mass formula
# ============================================================


def elliptic_automorphism_count(j):
    """
    |Aut(E)| for a supersingular elliptic curve (p > 3): 6 for j = 0, 4 for j = 1728, 2 otherwise
    """
    return 6 if j == 0 else 4 if j == 1728 else 2


def jacobian_automorphism_count(roots):
    """
    WARNING: this was LLM written, I don't understand it yet

    |Aut(C)| = |Aut(J(C), principal polarisation)| (Torelli) for C with (finite) Weierstrass points `roots`.

    Every automorphism of C induces a Mobius transformation of the x-line permuting the six Weierstrass points, and
    each such transformation lifts to exactly two automorphisms (differing by the hyperelliptic involution). A Mobius
    transformation is determined by the images of three points, so we try all 6 * 5 * 4 = 120 choices.

    A Mobius map x -> (a x + b) / (c x + d) is stored as the tuple (a, b, c, d) (up to scaling), and we compose
    S_images^-1 S_points, where S_t sends t1, t2, t3 to 0, 1, infinity.
    """
    points = list(roots)
    point_set = set(points)

    def to_standard(t1, t2, t3):
        return (t2 - t3, -t1 * (t2 - t3), t2 - t1, -t3 * (t2 - t1))

    def inverse_times(M, N):
        # M^-1 N, using the adjugate of M (scaling does not matter)
        a, b, c, d = M
        e, f, g, h = N
        return (d * e - b * g, d * f - b * h, -c * e + a * g, -c * f + a * h)

    start = to_standard(*points[:3])
    count = 0
    for images in permutations(points, 3):
        a, b, c, d = inverse_times(to_standard(*images), start)
        moved = set()
        for t in points[3:]:  # the first three points go to `images` by construction
            denominator = c * t + d
            if denominator == 0:
                break
            moved.add((a * t + b) / denominator)
        else:
            count += moved | set(images) == point_set
    return 2 * count


def product_automorphism_count(j1, j2):
    """
    |Aut(E1 x E2)| with its product polarisation: |Aut(E1)| |Aut(E2)|, twice that for E x E (which also has the swap)
    """
    count = elliptic_automorphism_count(j1) * elliptic_automorphism_count(j2)
    return 2 * count if j1 == j2 else count


def automorphism_count(V):
    """
    |Aut| of a vertex
    """
    if isinstance(V, Jacobian):
        return jacobian_automorphism_count(V.roots)
    return product_automorphism_count(*V.key[1:])


def mass(aut):
    """
    sum over the vertices of 1 / |Aut|, to compare with the Hashimoto-Ibukiyama mass formula (p - 1)(p^2 + 1) / 5760.
    `aut` maps each vertex key to |Aut|.
    """
    return sum(QQ(1) / a for a in aut.values())


def dual_edges_balance(aut, G):
    """
    Every isogeny A -> B has a dual B -> A, and the kernels A -> B are the isogenies up to Aut(B), so the numbers of
    edges satisfy  m(A, B) |Aut(B)| = m(B, A) |Aut(A)|  for every pair of vertices.
    """
    m = Counter(G.edges(labels=False, sort=False))
    return all(m[(a, b)] * aut[b] == m[(b, a)] * aut[a] for a, b in m)


def report(p):
    if p <= 5:
        print("p must be greater than 5, limitation based on the invariants paper I found")
        return

    vertices, G = isogeny_graph(p)
    products = [v for v, V in vertices.items() if isinstance(V, Product)]
    squares = [v for v in products if v[1] == v[2]]
    n = supersingular_count(p)
    target = QQ((p - 1) * (p**2 + 1)) / 5760
    aut = {v: automorphism_count(V) for v, V in vertices.items()}
    m = mass(aut)
    print(f"p = {p}")
    print(f"  {len(vertices)} vertices: {len(vertices) - len(products)} Jacobians, {len(products)} products, n = {n}")
    checks = [
        ("every vertex has 15 outgoing edges", set(G.out_degree()) == {15}),
        (f"products = n(n + 1) / 2 = {n * (n + 1) // 2}", len(products) == n * (n + 1) // 2),
        (f"squares E x E = n = {n}", len(squares) == n),
        (f"mass sum 1/|Aut| = {m} = (p - 1)(p^2 + 1) / 5760", m == target),
        ("m(A, B) |Aut B| = m(B, A) |Aut A| for all edges", dual_edges_balance(aut, G)),
        ("Richelot duals lead back", check_duals(vertices)),
    ]
    for label, ok in checks:
        print(f"  {label:<52}: {ok}")


if __name__ == "__main__":
    for p in map(int, sys.argv[1:]):
        report(p)
