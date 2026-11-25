"""
A canonical form for Hermitian O-lattices: a complete isomorphism invariant that can be used as a dictionary key or,
hashed, as the name of a vertex.

Two vertices L1, L2 are isomorphic iff there is a right O-linear isometry L1 -> L2. We describe every lattice in a way
that only depends on (L, h, right O-action), so that isomorphic lattices get the same description. This follows the
method of Dutour Sikirić, Haensch, Voight and van Woerden [DHVW20], adapted to the Hermitian setting, with one
change to keep the graphs small. The canonical labelling is done by bliss [JK07], [JK11] through SageMath. (Keys refer
to REFERENCES.md.)

    1. Take a characteristic set of short vectors V: go through the shells {x in L : h(x, x) = n} for n = 1, 2, ...
       and keep a shell when it enlarges the space spanned (over Q, as a right O-module) by the shells kept so far,
       until V spans L (x) Q. Any O-linear isometry maps V onto V.
    2. Make the arc-coloured digraph on V with colour h(v, w) in O (dropping arcs with h(v, w) = 0) and vertices
       coloured by h(v, v). Its automorphisms are the O-linear isometries of L (x) Q that permute V.
    3. Canonically label this digraph with bliss, which orders V canonically up to an automorphism of the digraph. The
       vectors v b_t (v in V in canonical order, b_t the fixed Z-basis of O) give a basis B of L (x) Q: the first 4g
       independent ones. The Gram matrices of B are then canonical.
    4. L itself is a lattice in L (x) Q, written in B-coordinates. It is only determined up to the automorphisms of the
       digraph, so we apply all of them and keep the smallest Hermite normal form.

The key is (the forms in the basis B, that smallest Hermite normal form). Equal keys give an O-linear isometry.

With an extra structure J (an automorphism, e.g. a Z[i]-action) the arcs are also coloured by h(J v, w), so the
automorphisms of the digraph are the isometries commuting with J, and the key describes the pair (L, J). (V is
J-stable since J is an isometry.) This needs no qfauto, which at large p overflows PARI's stack.

Why not ask V to span L itself, as in the paper? Then V needs every shell up to the norm where L is generated, and a
few lattices (those close to E0^g) have thousands of such vectors. The digraph has |V|^2 arcs, and naming one fourfold
at p = 23 took 8 GB that way. Spanning L (x) Q needs far fewer vectors (for that fourfold 32 instead of 3304); the
price is step 4, a loop over a group that is at most Aut of the sublattice spanned by V.
"""

from sage.all import QQ, ZZ, DiGraph, block_matrix, matrix, pari



# ============================================================
# The characteristic vector set
# ============================================================


def right_action_matrices(L):
    """
    The integer matrices of x -> x b_t on the coordinates (rows) of L, for the Z-basis b_t of O
    """
    ambient = L.ambient
    return [(L.M * ambient.right(b) * L.M_inverse).change_ring(ZZ) for b in ambient.ring.O.basis()]


def characteristic_vectors(L):
    """
    The shells {x in L : h(x, x) = n} (both signs, rows in coordinates of L) that enlarge the O-span of the shells
    before them, until they span L (x) Q
    """
    G = pari(L.form())  # b(x, x) = Trd h(x, x) = 2 h(x, x)
    rank = 4 * L.ambient.g
    actions = right_action_matrices(L)
    kept, span_rank, n = [], 0, 0
    while span_rank < rank:
        n += 1
        count, _, vectors = G.qfminim(2 * n)  # vectors of norm <= 2n (up to sign); keep those of norm exactly 2n
        if not int(count):
            continue
        shell = [v for v in matrix(ZZ, vectors.mattranspose()).rows() if v * L.form() * v == 2 * n]
        if not shell:
            continue
        candidate = kept + shell + [-v for v in shell]
        V = matrix(ZZ, candidate)
        new_rank = block_matrix([[V * A] for A in actions]).rank()
        if new_rank > span_rank:
            kept, span_rank = candidate, new_rank
    return matrix(ZZ, kept)


# ============================================================
# Canonical form
# ============================================================


def vector_digraph(V, H):
    """
    The digraph on range(len(V)) with an arc v -> w coloured by the coordinates of h(v, w) whenever it is non-zero,
    and the partition of the vertices by h(v, v), in increasing order
    """
    grams = [V * Ht * V.transpose() for Ht in H]
    N = V.nrows()
    arcs = []
    for a in range(N):
        rows = [Gm.row(a) for Gm in grams]
        for b in range(N):
            if a != b:
                label = tuple(int(r[b]) for r in rows)
                if any(label):
                    arcs.append((a, b, label))
    norms = {}
    for a in range(N):
        norms.setdefault(tuple(int(Gm[a, a]) for Gm in grams), []).append(a)
    partition = [norms[k] for k in sorted(norms)]
    return DiGraph([range(N), arcs], format="vertices_and_edges", loops=False, multiedges=False), partition


def spanning_rows(V, order, actions):
    """
    The first 4g linearly independent vectors among v b_t, for v in V in the given order: a basis of L (x) Q, as rows
    in coordinates of L, together with the (index into V, t) each comes from
    """
    rows, labels = [], []
    for a in order:
        for t, A in enumerate(actions):
            candidate = rows + [V.row(a) * A]
            if matrix(ZZ, candidate).rank() == len(candidate):
                rows, labels = candidate, labels + [(a, t)]
    return matrix(ZZ, rows), labels


def smallest_hermite_form(X):
    """
    The Hermite normal form of the rational lattice spanned by the rows of X, as a hashable tuple
    """
    d = X.denominator()
    H = (d * X).change_ring(ZZ).echelon_form(include_zero_rows=False)
    return (int(d),) + tuple(int(x) for x in H.list())


def canonical_key(L, structure=None):
    """
    A hashable key which is equal for two lattices iff they are isomorphic (O-linear isometric); with structure = R
    (a matrix on flat coordinates of an automorphism J), for the pairs (L, J) up to isometries commuting with J
    """
    V = characteristic_vectors(L)
    H = L.hermitian_forms()
    if structure is not None:  # h(x J, y) as well
        A = L.M * structure * L.M_inverse
        H = H + [A * Ht for Ht in H]
    D, partition = vector_digraph(V, H)
    _, certificate = D.canonical_label(partition=partition, edge_labels=True, certificate=True, algorithm="bliss")
    order = sorted(range(V.nrows()), key=lambda a: certificate[a])
    actions = right_action_matrices(L)
    B, labels = spanning_rows(V, order, actions)
    B_inverse = B.change_ring(QQ).inverse()

    # The Gram data of B only depends on the canonical order, so it is canonical
    forms = [B * F * B.transpose() for F in H]

    # L in B-coordinates is the row span of B^-1. An automorphism s of the digraph moves B to s(B) (the vectors
    # v_{s(a)} b_t), and moves L to the lattice with coordinates B^-1 s(B) B^-1. Keep the smallest Hermite form.
    group = D.automorphism_group(partition=partition, edge_labels=True, algorithm="bliss")
    best = None
    for s in group:
        moved = matrix(ZZ, [V.row(s(a)) * actions[t] for a, t in labels])
        candidate = smallest_hermite_form(B_inverse * moved * B_inverse)
        if best is None or candidate < best:
            best = candidate
    return tuple(tuple(F.list()) for F in forms) + (best,)
