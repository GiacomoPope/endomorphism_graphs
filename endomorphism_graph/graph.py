"""
Isomorphism classes and the breadth first search of the superspecial (ell, ..., ell)-isogeny graph in dimension g.

Two vertices L1, L2 are isomorphic principally polarised abelian varieties iff there is a right O-linear isometry
L1 -> L2. In dimension one this is Sage's is_right_equivalent; in general there is no Sage built in, so we use PARI's
qfisom on the trace form b(x, y) = Trd h(x, y) together with the forms b(x i, y) and b(x j, y), which force the
isometry to commute with the right O-action. The automorphism group of a vertex is the matching qfauto group.

Sources. Encoding the quaternionic structure by extra Z-bilinear forms, so that a quaternion Hermitian isometry test
becomes one for Z-lattices, is due to Greenberg and Voight [GV14, Lemma 6.2] and is used the same way by Chisholm
[Chi14, Section 3.8.1]. PARI's qfisom and qfauto implement the algorithm of Plesken and Souvignier [PS97]. The search
as a whole (neighbours, isometry testing, and the mass formula as a check that every class was found) is Chisholm's
algorithm [Chi14, Section 3.8]. The mass is

    sum over the vertices of 1 / |Aut|  =  prod_{k=1}^g |zeta(1 - 2k)| / 2 * (p^k + (-1)^k)

([Eke87], [HI80]; keys refer to REFERENCES.md): (p - 1) / 24 for g = 1 (Eichler) and (p - 1)(p^2 + 1) / 5760 for
g = 2.

So, the pattern is:

    1. Start from O^g, the lattice of E0^g with the product polarisation
    2. For each vertex compute its neighbours L_K, one for every Lagrangian K of A[ell] (torsion.py)
    3. Perform a BFS, keeping one lattice per isomorphism class (VertexSet) and recording an edge for every K. Two
       things cut down the lattices we build: only one K per orbit of the automorphism group of the vertex, and no K
       whose edge we already know as the dual of an edge found earlier
    4. From this BFS data we make the SageMath graph with DiGraph
"""

import time
from collections import deque

from sage.all import GF, QQ, ZZ, DiGraph, bernoulli, identity_matrix, matrix, pari, prod

from lattice import Ambient
from structures import conjugation_orbit, kernel_classifier, weil_class_type
from torsion import Frame, dual_kernel

# ============================================================
# Isomorphism classes
# ============================================================


def isometry_forms(L):
    """
    The trace form together with the forms encoding the right O-action [GV14, Lemma 6.2], [Chi14, Lemma 97]
    """
    ring = L.ambient.ring
    return [L.form(), L.right_action_form(ring.i), L.right_action_form(ring.j)]


def cheap_invariant(forms):
    """
    Number of vectors of each small norm of the trace form (its theta series): equal for isomorphic vertices
    """
    return tuple(int(c) for c in pari(forms[0]).qfrep(8, 1))


class VertexSet:
    """
    Representatives of the isomorphism classes found so far, bucketed by a cheap invariant (the theta series of the
    trace form). Only lattices in the same bucket are compared with qfisom, and for every representative we keep
    PARI's qfisominit precomputation, which makes each comparison 10 - 20 times faster (about 1 ms in g = 2, 3). The
    automorphism group (qfauto) is only computed for new classes.
    """

    def __init__(self):
        self.lattices = []
        self.forms = []
        self.orders = []      # |Aut| of each class
        self.generators = []  # generators of Aut, acting on coordinate column vectors
        self.buckets = {}
        self._isom_data = []  # qfisominit of each representative

    def index(self, L):
        """
        Return (index of the class of L, whether it is new, U) where U is an isometry from L to the stored
        representative, as an integer matrix on coordinate column vectors (the identity if L is new)
        """
        forms = isometry_forms(L)
        key = cheap_invariant(forms)
        for idx in self.buckets.get(key, []):
            if self._isom_data[idx] is None:
                self._isom_data[idx] = pari(self.forms[idx]).qfisominit()
            S = self._isom_data[idx].qfisom(pari(forms))  # forms[idx] = S^T forms S: S maps the representative to L
            if S != 0:
                return idx, False, matrix(ZZ, S).inverse().change_ring(ZZ)
        return self.add(L, forms, key), True, identity_matrix(ZZ, forms[0].nrows())

    def add(self, L, forms, key):
        """
        Store L as the representative of a new class, with its automorphism group
        """
        idx = len(self.lattices)
        order, generators = pari(forms).qfauto()
        G = forms[0]
        self.lattices.append(L)
        self.forms.append(forms)
        self.orders.append(int(order))
        self.generators.append([matrix(ZZ, s) for s in generators])
        assert all(s.transpose() * G * s == G for s in self.generators[-1])
        self._isom_data.append(None)  # qfisominit, computed the first time another lattice lands in this bucket
        self.buckets.setdefault(key, []).append(idx)
        return idx

    def automorphism_count(self, idx):
        """
        |Aut| of the vertex: the number of O-linear isometries of its lattice
        """
        return self.orders[idx]

    def automorphism_group(self, idx):
        """
        All O-linear isometries, as integer matrices acting on coordinate column vectors (closure of qfauto's
        generators)
        """
        gens = self.generators[idx]
        n = self.forms[idx][0].nrows()
        one = identity_matrix(ZZ, n)
        group = {tuple(one.list()): one}
        frontier = [one]
        while frontier:
            new = []
            for a in frontier:
                for g in gens:
                    b = a * g
                    if tuple(b.list()) not in group:
                        group[tuple(b.list())] = b
                        new.append(b)
            frontier = new
        assert len(group) == self.orders[idx]
        return list(group.values())

    def decomposition(self, idx):
        """
        The dimensions of the indecomposable factors of the vertex, as a sorted tuple: (g,) for an indecomposable
        variety (in g = 2: a Jacobian), (1, 1) for E1 x E2, (1, 1, 2) for E1 x E2 x (a surface), ...

        An automorphism s with s^2 = 1 whose +1 and -1 eigenlattices together span L gives an orthogonal splitting
        (h(x, y) = h(s x, s y) = -h(x, y) between them). The splitting of a definite lattice into indecomposables is
        unique (Eichler), so every such splitting is a union of the indecomposable factors, and the factor on which
        s = -1 and the rest = +1 is one of them. Refining by all such s therefore gives exactly the factors.
        """
        n = self.forms[idx][0].nrows()
        one = identity_matrix(ZZ, n)
        pieces = [one.row_space()]
        for s in self.automorphism_group(idx):
            if s == one or s == -one or s * s != one:
                continue
            plus = (s - one).right_kernel()
            minus = (s + one).right_kernel()
            if plus.rank() + minus.rank() != n or abs(plus.basis_matrix().stack(minus.basis_matrix()).det()) != 1:
                continue
            refined = []
            for P in pieces:
                refined += [Q for Q in (P.intersection(plus), P.intersection(minus)) if Q.rank()]
            pieces = refined
        return tuple(sorted(P.rank() // 4 for P in pieces))

    def is_decomposable(self, idx):
        """
        Whether the vertex is a product of lower dimensional principally polarised varieties
        """
        return len(self.decomposition(idx)) > 1


# ============================================================
# Isomorphism classes of vertices with extra structure
# ============================================================


class StructuredVertexSet:
    """
    Pairs (A, J) of a vertex and an automorphism J (given by a fixed matrix R on flat coordinates, which acts on every
    lattice we reach), up to isomorphisms commuting with J. Same interface as VertexSet.

    Giving PARI the extra form b(x J, y) makes qfauto and qfisom much slower (and at larger p overflow its stack), so
    we use the underlying varieties instead: (A, J) and (A', J') are isomorphic iff there is an isomorphism
    phi : A' -> A (found by VertexSet) with phi J' phi^-1 conjugate to J in Aut(A). For every variety we store the
    conjugacy classes of the J met so far, with a conjugator for each element (structures.conjugation_orbit), and the
    centraliser of J (the automorphisms of (A, J)) comes from the Schreier generators of that orbit.
    """

    def __init__(self, structure):
        self.structure = structure
        self.varieties = VertexSet()
        self.lattices = []    # the first lattice found for each pair, with its own J = structure
        self.orders = []      # |Aut(A, J)|
        self.generators = []  # generators of Aut(A, J) on coordinate column vectors of self.lattices[v]
        self.variety = []     # the index of A in self.varieties
        self._to_variety = [] # U0: an isometry from self.lattices[v] to the representative of A
        self._classes = {}    # (variety, J as a tuple) -> (pair, t) with J = t J0 t^-1, J0 the pair's J

    def index(self, L):
        """
        Return (index of the class of (L, J), whether it is new, U) with U an isometry from L to the stored lattice
        that commutes with J
        """
        f, _, U = self.varieties.index(L)
        J = L.column_matrix(self.structure)
        K = U * J * U.inverse()  # J moved to the representative of A
        found = self._classes.get((f, tuple(K.list())))
        if found is not None:
            v, t = found
            return v, False, self._to_variety[v].inverse() * t.inverse() * U
        return self.add(L, f, U, K), True, identity_matrix(ZZ, J.nrows())

    def add(self, L, f, U, K):
        """
        A new pair: record the conjugacy class of K under Aut(A), and the centraliser of K moved back to L
        """
        v = len(self.lattices)
        conjugator, schreier = conjugation_orbit(K, self.varieties.generators[f])
        for key, t in conjugator.items():
            self._classes[(f, key)] = (v, t)
        n, U_inverse = K.nrows(), U.inverse().change_ring(ZZ)
        self.lattices.append(L)
        self.variety.append(f)
        self._to_variety.append(U)
        self.orders.append(self.varieties.orders[f] // len(conjugator))
        self.generators.append([U_inverse * matrix(ZZ, n, n, list(z)) * U for z in schreier])
        J = L.column_matrix(self.structure)
        assert all(z * J == J * z for z in self.generators[-1])
        return v

    def automorphism_count(self, idx):
        return self.orders[idx]


# ============================================================
# Graph search
# ============================================================


def kernels_to_follow(F, structure, type_one):
    """
    Which Lagrangians of the frame F to follow, and a label for each: all of them, unlabelled, without structure;
    with a structure J only the J-stable ones, labelled by |K cap A[1 + i]| (constant on Aut(A, J)-orbits); and for
    the type I graph only the middle steps (label 4)
    """
    if structure is None:
        return None, lambda k: None
    J = F.L.column_matrix(structure)
    among, classify = F.stable_lagrangians(J), kernel_classifier(F, J)
    if type_one:
        among = [k for k in among if classify(k) == 4]
    return among, classify


def isogeny_graph(ring, g, ell=2, structure=None, type_one=False, verbose=False):
    """
    Breadth first search from E0^g through (ell, ..., ell)-isogenies. Returns the vertex set and the directed
    multigraph with one edge for every Lagrangian kernel of A[ell], labelled as in kernels_to_follow.

    With structure = R (the matrix on flat coordinates of an automorphism J of E0^g, e.g. a Z[i]-action from
    structures.py), the vertices are pairs (A, J) up to isomorphisms commuting with J, and the edges are only the
    J-stable Lagrangians: the connected component of (E0^g, J) in the graph of varieties with this structure. J is a
    fixed map on B^g, so it acts on every vertex we reach.

    With type_one=True as well (g = 4, ell = 2, J of Weil type), only the type I graph: type I vertices, and the
    middle steps between them (kernels K with |K cap A[1 + i]| = 4, which factor through a middle variety) whose
    target is again of type I. The type of a target is read off from the Weil pairing on its A[2] before it is looked
    up, so the steps leaving type I cost no isomorphism test.

    Two things cut down the lattices we build:

        - Orbits. Only one kernel per Aut-orbit is computed, and its edge is repeated once for every kernel in the
          orbit.
        - Duals. Every edge v -> w along K gives the dual edge w -> v along K' = dual_kernel. We move K' to the
          representative of w with the isometry found by VertexSet.index and remember it in known[w], so when w is
          processed the orbit of K' needs no lattice computation at all. This roughly halves the work.
    """
    vertices = VertexSet() if structure is None else StructuredVertexSet(structure)
    start, _, _ = vertices.index(Ambient(ring, g).start())

    frames = {}  # the frames of the vertices found but not yet processed (they are dropped once processed)

    def frame(v):
        if v not in frames:
            frames[v] = Frame(vertices.lattices[v], ell)
        return frames[v]

    known = {}  # known[w][k] = v: the Lagrangian with index k at w is the kernel of an edge w -> v
    edges = []
    queue, done = deque([start]), set()
    computed = skipped = 0
    t0 = time.time()
    while queue:
        v = queue.popleft()
        done.add(v)
        F = frame(v)
        known_here = known.pop(v, {})
        among, classify = kernels_to_follow(F, structure, type_one)
        orbits = F.orbits(vertices.generators[v], among)
        for orbit in orbits:
            label = classify(orbit[0])

            # the dual of an edge found earlier: nothing to compute
            targets = {known_here[k] for k in orbit if k in known_here}
            if targets:
                assert len(targets) == 1, "isomorphic kernels lead to different vertices"
                edges += [(v, targets.pop(), label)] * len(orbit)
                skipped += 1
                continue

            N = F.kernel_lattice(orbit[0])
            if type_one and weil_class_type(N, N.column_matrix(structure)) != "I":
                continue  # this step leaves the type I vertices
            w, new, U = vertices.index(N)
            computed += 1
            edges += [(v, w, label)] * len(orbit)
            if new:
                queue.append(w)

            # remember the dual edge w -> v, in the frame of w
            if w not in done or w == v:
                dual = (U.change_ring(GF(ell)) * dual_kernel(F.L, N, ell).transpose()).transpose()
                k = frame(w).index_of(dual)
                (known_here if w == v else known.setdefault(w, {}))[k] = v
        frames.pop(v)
        if verbose:
            print(f"  processed {len(done)} / found {len(vertices.lattices)}  ({len(orbits)} orbits, "
                  f"{computed} lattices computed, {skipped} known duals, {time.time() - t0:.1f}s)")

    return vertices, DiGraph(edges, multiedges=True, loops=True)


# ============================================================
# Mass formula
# ============================================================


def mass(vertices):
    """
    sum over the vertices of 1 / |Aut|
    """
    return sum(QQ(1) / vertices.automorphism_count(v) for v in range(len(vertices.lattices)))


def mass_formula(p, g):
    """
    prod_{k=1}^g |zeta(1 - 2k)| / 2 * (p^k + (-1)^k), using zeta(1 - 2k) = -B_2k / 2k
    """
    return prod(abs(bernoulli(2 * k)) / (4 * k) * (p**k + (-1) ** k) for k in range(1, g + 1))
