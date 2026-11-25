"""
The endomorphism ring of the starting curve

    E0 : y^2 = x^3 + x over F_{p^2}, for p = 3 mod 4.

E0 is supersingular and End(E0) is a maximal order O in the quaternion algebra

    B = (-1, -p)_Q  with  i^2 = -1, j^2 = -p, k = ij,

where i is the automorphism (x, y) -> (-x, sqrt(-1) y) and j is the Frobenius endomorphism. Everything in the
graph search is done inside B, so the curves themselves never appear.
"""

from sage.all import GF, QQ, QuaternionAlgebra, cartesian_product, matrix, vector


class EndomorphismRing:
    """
    The maximal order O = End(E0) together with a rank one idempotent e of O / 2O = M_2(F_2) (and, on request, of
    O / ell O = M_2(F_ell) for any prime ell != p).

    The idempotent is how we see E0[ell] inside O: the left O-module (O / ell O) e is two dimensional over F_ell and
    O acts on it through End(E0[ell]), so it is a model of E0[ell]. More generally, for a right ideal I = Hom(E0, E)
    we get E[ell] = (I / ell I) e: an isogeny f is recorded by where it sends one point of E0[ell].
    """

    def __init__(self, p):
        assert p % 4 == 3, "E0 : y^2 = x^3 + x is supersingular for p = 3 mod 4"
        self.p = p
        self.B = QuaternionAlgebra(QQ, -1, -p)
        self.i, self.j, self.k = self.B.gens()
        self.O = self.B.maximal_order()
        self._idempotents = {}
        self._weil_functionals = {}
        self.e = self.idempotent(2)

    def unit_ideal(self):
        """
        The right ideal O = Hom(E0, E0): the starting vertex
        """
        return self.O.right_ideal(self.O.basis())

    def idempotent(self, ell):
        """
        Find e in O whose image in O / ell O = M_2(F_ell) is idempotent of rank one:
        e^2 = e mod ell O, with neither e nor 1 - e in ell O
        """
        if ell in self._idempotents:
            return self._idempotents[ell]
        assert ell != self.p, "O / pO is not a matrix algebra"
        O = self.O
        for c in cartesian_product([range(ell)] * 4):
            e = sum(ci * b for ci, b in zip(c, O.basis()))
            if e / ell in O or (1 - e) / ell in O:
                continue
            if (e * e - e) / ell in O:
                self._idempotents[ell] = e
                return e
        raise ValueError("no rank one idempotent found")

    def coordinates(self, q):
        """
        The coordinates of q in the Z-basis of O
        """
        basis = matrix(QQ, [list(b) for b in self.O.basis()])
        return vector(QQ, list(q)) * basis.inverse()

    def weil_functional(self, ell):
        """
        The Weil pairing on E[ell] = (I / ell I) e takes values in (1 - e)(O / ell O) e, which is one dimensional over
        F_ell (the (2, 1) entry of a 2 x 2 matrix, when e = diag(1, 0)). Return (t, c): for z in this line, given by
        its coordinates in the basis of O mod ell, its value in F_ell is c * z[t].
        """
        if ell in self._weil_functionals:
            return self._weil_functionals[ell]
        e, F = self.idempotent(ell), GF(ell)
        line = [vector(F, self.coordinates((1 - e) * b * e)) for b in self.O.basis()]
        w = next(z for z in line if z != 0)
        assert all(matrix(F, [w, z]).rank() == 1 for z in line if z != 0), "not a line"
        t = next(t for t in range(4) if w[t] != 0)
        self._weil_functionals[ell] = (t, 1 / w[t])
        return self._weil_functionals[ell]
