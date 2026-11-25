"""
Hermitian O-lattices of rank g: the vertices of the superspecial isogeny graph in dimension g.

A principally polarised superspecial abelian variety A of dimension g is recorded by the lattice

    L = Hom(E0, A),   a right O-submodule of B^g (rank 4g over Z),

together with the Hermitian form coming from the principal polarisation of A. For E0^g with the product polarisation
this is L = O^g with h(x, y) = sum conj(x_k) y_k. Every vertex we reach from E0^g by (ell, ..., ell)-isogenies is a
sublattice of O^g with the form scaled down by a power of ell (the "scale" n = ell^s),

    h_n(x, y) = n^-1 sum_k conj(x_k) y_k,

and it is principally polarised exactly when h_n takes values in O on L (unimodular). For g = 1 this is the right
ideal I = Hom(E0, E) with h(x, x) = Nrd(x) / Nrd(I).

Sources (keys refer to REFERENCES.md). For g >= 2 every superspecial A is isomorphic to E0^g as an unpolarised
variety, and principal polarisations on it correspond to the classes of unimodular quaternion Hermitian forms
[IKO86, Thm 2.10], so these lattices see every vertex. The equivalence A -> Hom(E, A) between abelian varieties
isogenous to a power of E and modules over End(E) is [JKPRST18].

Coordinates. We write an element x = (x_1, ..., x_g) of B^g in the Z-basis (b_1, b_2, b_3, b_4) of O in each slot,
which gives a row vector of 4g rationals ("flat coordinates"). Then everything is linear algebra over Z:

    - a lattice is a 4g x 4g rational matrix M whose rows are a Z-basis,
    - right multiplication x -> x a is a 4g x 4g matrix R_a acting on rows: flat(x a) = flat(x) R_a,
    - the trace form b(x, y) = Trd h_n(x, y) has Gram matrix M G0 M^T / n.

Coordinates of L itself (in its basis M) are used for everything that only depends on the lattice: Gram matrices, the
maps mod ell, and the automorphisms PARI returns, which act on coordinate column vectors.
"""

from sage.all import GF, QQ, ZZ, MatrixSpace, block_diagonal_matrix, identity_matrix, lcm, matrix, vector


# ============================================================
# The ambient space B^g
# ============================================================


class Ambient:
    """
    B^g in flat coordinates, with the matrices we need precomputed
    """

    def __init__(self, ring, g):
        self.ring, self.g = ring, g
        basis = ring.O.basis()
        to_O = matrix(QQ, [list(b) for b in basis]).inverse()  # B-coordinates -> O-coordinates
        self.O_coordinates = lambda q: vector(QQ, list(q)) * to_O

        # one slot: Trd(conj(b_r) b_c) and the O-coordinates of conj(b_r) b_c, as 4 x 4 matrices
        trace = matrix(ZZ, 4, 4, lambda r, c: (basis[r].conjugate() * basis[c]).reduced_trace())
        values = [matrix(ZZ, 4, 4, lambda r, c, t=t: self.O_coordinates(basis[r].conjugate() * basis[c])[t])
                  for t in range(4)]
        self.G0 = block_diagonal_matrix([trace] * g)              # Gram matrix of Trd h on O^g
        self.C0 = [block_diagonal_matrix([C] * g) for C in values]  # h itself, one matrix per O-coordinate
        self._right = {}

    def right(self, a):
        """
        The 4g x 4g matrix of x -> x a on flat coordinates (the same 4 x 4 block in every slot). We only ever use a
        handful of elements a (i, j, e and a e), so the matrices are cached.
        """
        key = tuple(a)
        if key not in self._right:
            block = matrix(QQ, [self.O_coordinates(b * a) for b in self.ring.O.basis()])
            self._right[key] = block_diagonal_matrix([block] * self.g)
        return self._right[key]

    def start(self):
        """
        O^g with the product polarisation: the lattice of E0^g
        """
        return HermitianLattice(self, identity_matrix(QQ, 4 * self.g), 1)


# ============================================================
# Lattices
# ============================================================


class HermitianLattice:
    """
    The Z-span of the rows of a rational matrix (flat coordinates), with the Hermitian form h_n = h / n (n = scale)
    """

    def __init__(self, ambient, rows, scale):
        self.ambient, self.scale = ambient, scale
        n = 4 * ambient.g
        d = lcm([x.denominator() for x in rows.list()])
        H = (d * rows).change_ring(ZZ).hermite_form(include_zero_rows=False)
        assert H.nrows() == n, "not a full rank lattice"
        self.M = H / d
        self.M_inverse = self.M.inverse()

    # ---------- forms ----------

    def form(self, R=None):
        """
        The integer matrix of b(x R, y) = Trd h_n(x R, y) on the basis (R = None: the trace form itself)
        """
        A = self.M if R is None else self.M * R
        return (A * self.ambient.G0 * self.M.transpose() / self.scale).change_ring(ZZ)

    def right_action_form(self, a):
        """
        b(x a, y): together with the trace form these force isometries to be O-linear
        """
        return self.form(self.ambient.right(a))

    def hermitian_forms(self):
        """
        Four integer matrices H_t with h_n(x, y) = sum_t (x H_t y^T) b_t, for x, y in coordinates of L (the four
        O-coordinates of the Hermitian form)
        """
        return [(self.M * C * self.M.transpose() / self.scale).change_ring(ZZ) for C in self.ambient.C0]

    def is_principally_polarised(self):
        """
        Whether h_n is a unimodular Hermitian form on L: it takes values in O, and L has the volume of O^g (the same
        determinant of the trace form). Integrality alone is not enough: a sublattice of a principally polarised
        lattice is still integral.
        """
        n = 4 * self.ambient.g
        integral = all((self.M * C * self.M.transpose() / self.scale) in MatrixSpace(ZZ, n) for C in self.ambient.C0)
        return integral and self.form().det() == self.ambient.G0.det()

    # ---------- maps on coordinates ----------

    def column_matrix(self, R):
        """
        The integer matrix of x -> x R (R on flat coordinates, preserving L) on coordinate column vectors of L, the
        convention of qfauto's automorphisms
        """
        A = self.M * R * self.M_inverse
        assert A.denominator() == 1, "R does not preserve L"
        return A.change_ring(ZZ).transpose()

    def mod_matrix(self, R, ell):
        """
        The F_ell-matrix of x -> x R on L / ell L in the basis of L (columns are images of basis vectors)
        """
        return self.column_matrix(R).change_ring(GF(ell))
