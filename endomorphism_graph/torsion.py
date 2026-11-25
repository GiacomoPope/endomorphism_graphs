"""
The ell-torsion A[ell] and the (ell, ..., ell)-isogenies out of a vertex, for a Hermitian lattice L = Hom(E0, A) of
rank g and a prime ell != p.

Everything is linear algebra on L / ell L:

    - Evaluating homomorphisms on E0[ell] gives L / ell L -> Hom(E0[ell], A[ell]), and A[ell] = (L / ell L) e, where e
      is a rank one idempotent of O / ell O = M_2(F_ell) ("evaluation at one point of E0[ell]"). It is a
      2g-dimensional subspace of F_ell^4g.
    - For a subgroup K of A[ell], the quotient A' = A / K has

          Hom(E0, A') = L_K = {f in L : f(E0[ell]) in K} = {f in L : f a e in K mod ell L for all a in O},

      because f factors through the dual isogeny A' -> A exactly when it sends E0[ell] into K.
    - The isogeny A -> A / K respects principal polarisations (it is an (ell, ..., ell)-isogeny) exactly when K is a
      Lagrangian: a maximal (dimension g) subspace on which the Weil pairing vanishes. L_K then carries the form
      divided by ell, and is again unimodular. In dimension g there are prod_{i=1}^g (ell^i + 1) Lagrangians.

The Weil pairing. For u = x e and v = y e in A[ell] the value h(x e, y e) = conj(e) h(x, y) e = (1 - e) h(x, y) e lies
in (1 - e) M_2(F_ell) e modulo ell, which is one dimensional (EndomorphismRing.weil_functional reads off its
coordinate). We can take x = u, y = v (lifted to L), since u e = u. It is alternating because h(u, u) is an integer
and (1 - e) e = 0. (For ell = 2: omega(u, v) = 1 iff h(u, v) is not in 2O.)

Sources. This is Kneser's neighbour method [Kne57] for quaternion Hermitian lattices, as developed by Bachoc [Bac95]
and, in exactly this symplectic form, by Chisholm [Chi14]. Over F_ell the quaternion algebra splits, O / ell O =
M_2(F_ell), and the Morita equivalence V -> e_11 V, with the symplectic form h_12 read off from one entry of h
[Chi14, Section 2.4, Prop. 36, and Section 3.6, Prop. 88], turns L / ell L into a symplectic F_ell-space of dimension
2g. Our A[ell] = (L / ell L) e and the coordinate of (1 - e) h e are this construction for right modules. The
ell-neighbours of L are then in bijection with the Lagrangians of that space, the neighbour being the preimage of the
Lagrangian plus ell L [Chi14, Lemma 83 and (3.5.1)]. That these neighbours are the lattices Hom(E0, A / K) of the
(ell, ..., ell)-isogenies, so that the neighbour graph is the isogeny graph, is the setting of Jordan and Zaytman
[JZ26]. Keys like [Chi14] refer to REFERENCES.md at the top of the repository.

The frame. A symplectic basis of A[ell] identifies it with the standard F_ell^2g, whose Lagrangians we list once per
(g, ell). Everything per vertex (kernels, Aut-orbits, dual kernels) is then done on small tuples of integers, and a
Lagrangian is an index into the standard list.
"""

from itertools import combinations, product

from sage.all import GF, ZZ, block_matrix, matrix, vector

from lattice import HermitianLattice


# ============================================================
# Torsion and the Weil pairing
# ============================================================


def torsion(L, ell):
    """
    Return A[ell] = (L / ell L) e as a subspace of F_ell^4g, together with the four maps x -> x a e (a in the basis of
    O) which record how an element of L acts on E0[ell]
    """
    ambient = L.ambient
    e = ambient.ring.idempotent(ell)
    A = L.mod_matrix(ambient.right(e), ell).column_space()
    evaluation_maps = [L.mod_matrix(ambient.right(a * e), ell) for a in ambient.ring.O.basis()]
    return A, evaluation_maps


def weil_pairing_matrix(L, A, ell):
    """
    The Weil pairing on A[ell] in the basis of A: the coordinate of (1 - e) h(u, v) e, which is h(u, v) itself mod ell.
    This is the symplectic form h_12 of the Morita equivalence [Chi14, Lemma 34, Prop. 88].
    """
    F = GF(ell)
    t, c = L.ambient.ring.weil_functional(ell)
    U = A.basis_matrix().change_ring(ZZ)
    values = (U * L.hermitian_forms()[t] * U.transpose()).change_ring(F)
    W = c * values
    assert W.transpose() == -W and W.det() != 0, "not a symplectic form"
    return W


# ============================================================
# The lattice of an isogeny with kernel K, and its dual
# ============================================================


def kernel_lattice(L, K, evaluation_maps, ell):
    """
    L_K = {f in L : f a e in K mod ell L for all a in O} = Hom(E0, A / K), with the form divided by ell. K is given
    by a matrix whose rows span it.

    Up to the scaling of the form, this is the ell-neighbour of L attached to the Lagrangian K: the preimage of K plus
    ell L [Chi14, Lemma 83, (3.5.1)].
    """
    F = GF(ell)
    annihilator = K.right_kernel_matrix()  # linear forms vanishing on K
    conditions = block_matrix(F, [[annihilator * M] for M in evaluation_maps])
    kept = conditions.right_kernel().basis_matrix().change_ring(ZZ)  # coordinates (mod ell) of the f we keep
    return HermitianLattice(L.ambient, (kept * L.M).stack(ell * L.M), L.scale * ell)


def dual_kernel(L, N, ell):
    """
    For N = L_K, the kernel K' in N[ell] of the dual isogeny A / K -> A, as a matrix whose rows span it (coordinates
    of N). The composition is multiplication by ell, so N_{K'} = ell L (with the form divided by ell^2, isometric to
    L), and K' = (ell L / ell N) e: the image of ell L in N[ell].
    """
    ell_L = (ell * L.M * N.M_inverse).change_ring(ZZ).change_ring(GF(ell))  # ell L in coordinates of N, mod ell
    e_map = N.mod_matrix(N.ambient.right(N.ambient.ring.idempotent(ell)), ell)  # acts on columns
    return (e_map * ell_L.transpose()).transpose().row_space().basis_matrix()


# ============================================================
# Standard Lagrangians
# ============================================================
#
# In the standard symplectic space F_ell^2g, with basis e_1, ..., e_g, f_1, ..., f_g and omega(e_i, f_j) = delta_ij,
# every Lagrangian is transversal to one of the 2^g coordinate Lagrangians W_J = <f_i (i not in J), -e_j (j in J)>,
# and the Lagrangians transversal to W_J are exactly the graphs {u + S u} of symmetric matrices S from the
# complementary U_J = <e_i (i not in J), f_j (j in J)> to W_J. So we run through all J and symmetric S, and keep one
# copy of each (by reduced row echelon form). Vectors are tuples of integers mod ell.
#
# Chisholm lists the same Lagrangians for g = 2 by their row echelon forms and counts (p^2 + 1)(p + 1) of them
# [Chi14, Example 25, Example 85]; the chart by symmetric matrices is the standard way to do this for every g.


def rref(rows, ell):
    """
    The reduced row echelon form (non-zero rows, as a tuple of tuples) of a list of vectors mod ell
    """
    rows = [list(r) for r in rows]
    n = len(rows[0]) if rows else 0
    out, r = [], 0
    for c in range(n):
        pivot = next((k for k in range(r, len(rows)) if rows[k][c] % ell), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inverse = pow(rows[r][c], -1, ell)
        rows[r] = [(x * inverse) % ell for x in rows[r]]
        for k in range(len(rows)):
            if k != r and rows[k][c] % ell:
                factor = rows[k][c]
                rows[k] = [(x - factor * y) % ell for x, y in zip(rows[k], rows[r])]
        r += 1
    return tuple(tuple(row) for row in rows[:r])


def lagrangian_count(g, ell):
    """
    The number of Lagrangians of a symplectic F_ell^2g: prod_{i=1}^g (ell^i + 1)
    """
    out = 1
    for i in range(1, g + 1):
        out *= ell**i + 1
    return out


def apply_map(images, v, ell):
    """
    The image of the vector v under the linear map sending basis vector r to images[r]
    """
    out = [0] * len(images[0])
    for x, image in zip(v, images):
        if x:
            out = [(a + x * b) % ell for a, b in zip(out, image)]
    return tuple(out)


class StandardLagrangians:
    """
    The Lagrangians of the standard symplectic F_ell^2g, numbered 0, 1, ..., each stored by its reduced row echelon
    basis
    """

    _cache = {}

    @classmethod
    def of(cls, g, ell):
        if (g, ell) not in cls._cache:
            cls._cache[(g, ell)] = cls(g, ell)
        return cls._cache[(g, ell)]

    def __init__(self, g, ell):
        self.g, self.ell = g, ell
        unit = lambda r: tuple(int(k == r) for k in range(2 * g))
        e = [unit(i) for i in range(g)]
        f = [unit(g + i) for i in range(g)]
        neg = lambda v: tuple((-x) % ell for x in v)
        found = set()
        positions = [(i, j) for i in range(g) for j in range(i, g)]  # the entries of a symmetric g x g matrix
        for size in range(g + 1):
            for J in combinations(range(g), size):
                # the chart: U_J and W_J, with omega(U_J[i], W_J[k]) = delta_ik
                U = [f[i] if i in J else e[i] for i in range(g)]
                W = [neg(e[i]) if i in J else f[i] for i in range(g)]
                # the graph of every symmetric S : U_J -> W_J, rows U[i] + sum_k S[i][k] W[k]
                for entries in product(range(ell), repeat=len(positions)):
                    S = [[0] * g for _ in range(g)]
                    for (i, j), x in zip(positions, entries):
                        S[i][j] = S[j][i] = x
                    rows = [tuple((U[i][c] + sum(S[i][k] * W[k][c] for k in range(g))) % ell for c in range(2 * g))
                            for i in range(g)]
                    found.add(rref(rows, ell))
        self.bases = sorted(found)
        assert len(self.bases) == lagrangian_count(g, ell), "wrong number of Lagrangians"
        self.index = {B: k for k, B in enumerate(self.bases)}

    def __len__(self):
        return len(self.bases)

    def image(self, k, images):
        """
        The index of the image of Lagrangian k under the linear map sending basis vector r to images[r]
        """
        return self.span_index([apply_map(images, v, self.ell) for v in self.bases[k]])

    def span_index(self, vectors):
        """
        The index of the Lagrangian spanned by the given vectors
        """
        return self.index[rref(vectors, self.ell)]


# ============================================================
# The torsion of a vertex in a symplectic frame
# ============================================================


def symplectic_basis(W):
    """
    A change of basis T (rows) with T W T^t = the standard form (vector r paired with r + g, omega = 1), by
    symplectic Gram-Schmidt (as in [Chi14, Lemma 27])
    """
    n = W.nrows()
    pair = lambda u, v: u * W * v
    rest = list(matrix.identity(W.base_ring(), n).rows())
    es, fs = [], []
    while rest:
        a = rest.pop(0)
        if a == 0:
            continue
        k = next(k for k, b in enumerate(rest) if pair(a, b) != 0)
        b = rest.pop(k)
        b = b / pair(a, b)  # now omega(a, b) = 1, omega(b, a) = -1
        es.append(a)
        fs.append(b)
        rest = [x - pair(x, b) * a + pair(x, a) * b for x in rest]
    return matrix(W.base_ring(), es + fs)


class Frame:
    """
    A[ell] of a vertex L with a symplectic basis (the rows of self.basis, vectors of F_ell^4g), so that Lagrangians of
    A[ell] are indices into StandardLagrangians.of(g, ell)
    """

    def __init__(self, L, ell):
        self.L, self.ell = L, ell
        self.g = L.ambient.g
        self.lagrangians = StandardLagrangians.of(self.g, ell)
        self.A, self.evaluation_maps = torsion(L, ell)
        assert self.A.dimension() == 2 * self.g
        T = symplectic_basis(weil_pairing_matrix(L, self.A, ell))
        self.basis = T * self.A.basis_matrix()
        self._pivots = self.basis.pivots()
        self._solve = self.basis.matrix_from_columns(self._pivots).inverse()  # pivot entries -> frame coordinates

    def coordinates(self, y):
        """
        The frame coordinates of a vector y of A[ell], as a tuple of integers
        """
        c = vector(self.basis.base_ring(), [y[j] for j in self._pivots]) * self._solve
        assert c * self.basis == y, "vector not in A[ell]"
        return tuple(int(x) for x in c)

    def vector(self, c):
        """
        The vector of F_ell^4g with frame coordinates c
        """
        return vector(self.basis.base_ring(), c) * self.basis

    def kernel_lattice(self, k):
        """
        L_K for the Lagrangian with index k: the neighbour along that edge
        """
        K = matrix(self.basis.base_ring(), [self.vector(c) for c in self.lagrangians.bases[k]])
        return kernel_lattice(self.L, K, self.evaluation_maps, self.ell)

    def index_of(self, K):
        """
        The index of the Lagrangian spanned by the rows of K
        """
        return self.lagrangians.span_index([self.coordinates(y) for y in K.rows()])

    def action(self, s):
        """
        An automorphism s of L (an integer matrix on coordinate column vectors) on A[ell], as images of the frame
        basis. It preserves A[ell] and the Weil pairing because it is O-linear and an isometry.
        """
        S = s.change_ring(self.basis.base_ring())
        return [self.coordinates(S * self.basis.row(r)) for r in range(2 * self.g)]

    def is_scalar(self, images):
        """
        Whether the map with these images of the frame basis is a scalar (it then fixes every subspace)
        """
        return all(images[r] == tuple(images[0][0] * int(k == r) for k in range(2 * self.g)) for r in range(2 * self.g))

    def stable_lagrangians(self, J):
        """
        The indices of the Lagrangians mapped to themselves by J (an automorphism, on coordinate column vectors)
        """
        images = self.action(J)
        return [k for k in range(len(self.lagrangians)) if self.lagrangians.image(k, images) == k]

    def orbits(self, generators, among=None):
        """
        The orbits of the group generated by the automorphisms on the Lagrangians (or on the subset among, which they
        must preserve), as lists of indices. L_K and L_{s(K)} are isomorphic (s maps one onto the other), so one
        kernel per orbit is enough. The scalars (in particular -1) fix every subspace, so they are dropped.
        """
        n = len(self.lagrangians)
        indices = range(n) if among is None else among
        actions = [a for a in (self.action(s) for s in generators) if not self.is_scalar(a)]
        if not actions:
            return [[k] for k in indices]
        orbit_of, orbits = [None] * n, []
        for k in indices:
            if orbit_of[k] is not None:
                continue
            orbit, frontier = [k], [k]
            orbit_of[k] = len(orbits)
            while frontier:
                new = []
                for m in frontier:
                    for a in actions:
                        t = self.lagrangians.image(m, a)
                        if orbit_of[t] is None:
                            orbit_of[t] = len(orbits)
                            orbit.append(t)
                            new.append(t)
                frontier = new
            orbits.append(orbit)
        return orbits


# ============================================================
# Neighbours
# ============================================================


def neighbours(L, ell=2):
    """
    The lattices L_K for the Lagrangians K of A[ell]: one for each (ell, ..., ell)-isogeny out of the vertex
    """
    frame = Frame(L, ell)
    return [frame.kernel_lattice(k) for k in range(len(frame.lagrangians))]
