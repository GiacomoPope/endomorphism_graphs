# Superspecial isogeny graphs from quaternions

A SageMath implementation of the $(\ell, \ldots, \ell)$-isogeny graphs of principally polarised superspecial abelian varieties of any dimension $g$, computed with quaternion Hermitian lattices. It works for any $g$, $\ell$ and $p$, provided everything is small (very small).

The method is essentially that of Sarah Chisholm's PhD thesis, *Algorithmic enumeration of quaternionic lattices* (University of Calgary, 2014, [doi:10.11575/PRISM/24766](https://doi.org/10.11575/PRISM/24766)), which came with Magma code. If you use this code, or ideas from it, please cite **Sarah Chisholm** and the other sources in [REFERENCES.md](REFERENCES.md). They are the foundation for everything here. The connection to the thesis was found by Eda Kirimli.

The intention is a fun toy for experimenting with, and visualising, these graphs. For $g = 1$ and $2$, [Benjamin Wesolowski's Gallery](https://bweso.com/gallery.html) is a more interactive and detailed exploration. This code is slower and clunkier, but it also works for $g = 3$ and $4$.

**How this code was written.** Most of the code, and the first drafts of the docstrings, were generated with an LLM (Claude), starting from my own code for dimensions one and two. The sources were found and checked afterwards (see REFERENCES.md). Treat everything here as unverified until the code and the theory have been checked carefully.

Contributions are welcome:

1. proof reading the code, to make sure the docstrings are theoretically correct;
2. making sure the state of the art is implemented and properly cited;
3. opening issues for features that would be useful to you.

Happy graphing!

## How it works

Keys like [Chi14] refer to [REFERENCES.md](REFERENCES.md).

- **Vertices are lattices.** Let $E_0 : y^2 = x^3 + x$ with $\mathcal{O} = \mathrm{End}(E_0)$, a maximal order in the quaternion algebra $B$ ramified at $p$ and $\infty$. A principally polarised superspecial abelian variety $A$ of dimension $g$ is recorded by the lattice $L = \mathrm{Hom}(E_0, A) \subset B^g$, with the Hermitian form coming from its polarisation. The form is unimodular exactly when the polarisation is principal [IKO86, Thm 2.10], [JKPRST18].
- **Edges are neighbours.** Modulo $\ell$, $\mathcal{O}/\ell\mathcal{O} \cong M_2(\mathbb{F}_\ell)$. The Morita equivalence turns $L/\ell L$ into a symplectic space of dimension $2g$: our $A[\ell] = (L/\ell L)e$ for a rank one idempotent $e$, with the Weil pairing read off from $h$. The $\ell$-neighbours of $L$ correspond to the Lagrangians of that space. This is Kneser's neighbour method [Kne57], for quaternionic lattices [Bac95], in the symplectic form of [Chi14, Sections 2.4 and 3.5-3.7]. The neighbour attached to a Lagrangian $K$ is $\mathrm{Hom}(E_0, A/K)$, so the neighbour graph is the isogeny graph [JZ26].
- **Isomorphism testing.** Two vertices are isomorphic when there is an $\mathcal{O}$-linear isometry between their lattices. Following [GV14, Lemma 6.2] and [Chi14, Section 3.8.1], this is tested on $\mathbb{Z}$-lattices with extra forms that force $\mathcal{O}$-linearity, using PARI's `qfisom` and `qfauto` [PS97]. Vertex names come from a canonical form [DHVW20] computed with bliss [JK07].
- **Checks.** The sum of $1/|\mathrm{Aut}|$ over the vertices must equal the mass formula [Eke87], [HI80], as in [Chi14, Section 3.8.2]. For $g = 1$ the graph is compared with Sage's, and for $g = 2$ with the Richelot graph.

What we add is mostly engineering:
- one kernel per orbit of the automorphism group of a vertex;
- never recomputing the dual of a known edge;
- the restriction to Weil-type fourfolds;
- the painter.

### Vertices without edges: stopping at the mass

If only the vertices are wanted, there is a much cheaper search [Chi14, Sections 3.7-3.8]. Explore neighbours, at one or several primes $\ell$, until the sum of $1/|\mathrm{Aut}|$ over the classes found equals the mass formula, then stop. No edges need to be recorded, and most kernels never have to be built.

This is not implemented here, and it gives something up. In the full search, reaching the mass is a *check*: the graph closes up, every edge has its dual, the degrees and the comparison with Richelot come out right, and the mass agrees as well. Stopping at the mass makes the formula the *proof* that the list is complete. So it assumes that the construction is right (neighbours stay in the genus, the isometry test and the automorphism groups are correct) and that the formula counts exactly the set being enumerated. Where the formula is only conjectural, for example for the type I graph below, a mass-stopped search cannot be used as evidence for it.

## Running it

Everything is plain Python importing from SageMath (tested with Sage 10.9). The canonical vertex names use the `bliss` package. If your Sage does not have it, pass `--no-names` when saving graphs. The painter only needs Python 3 and a browser. All commands are run from the top of the repository.

### Dimension one (`dimension_one/`)

Dimension one is the only place where we can compare against Sage directly. There are three ways to compute the graph:

1. Sage, using the `E.ell_isogeny_graph` method.
2. The $\ell = 2$ modular polynomial for a draft of the BFS.
3. The $g = 1$ case of the generic code in `endomorphism_graph/`.

```
sage -python dimension_one/compare.py 7 11 43 83     # the three graphs agree, and the mass formula holds
```

### Dimension two (`dimension_two/`)

Dimension two is not covered by Sage, but the graph can be computed in two independent ways: with Richelot isogenies [CDS20], [FS22] (edge cases following [Zob21]), and as the $g = 2$ case of `endomorphism_graph/`. The two graphs agree. The Richelot code also includes my gluing and splitting formulas for elliptic curves which are currently unpublished.

```
sage -python dimension_two/richelot.py 43 83     # mass formula, dual edges, products and squares checks
```

### Any dimension (`endomorphism_graph/`)

| file | what it does |
|---|---|
| `quaternions.py` | $\mathcal{O} = \mathrm{End}(E_0)$, its idempotents mod $\ell$ and the Weil pairing functional |
| `lattice.py` | $B^g$ in coordinates and the Hermitian lattices |
| `torsion.py` | $A[\ell]$, the Weil pairing, the Lagrangians, the neighbours $L_K$ and the dual kernels |
| `graph.py` | isomorphism classes, the breadth first search, the decomposition into factors and the mass formula |
| `structures.py` | Frobenius (which vertices are defined over $\mathbb{F}_p$) and the $\mathbb{Z}[i]$-actions of Weil type (types I and II) |
| `canonical.py` | a canonical form for the lattices, which names every vertex independently of the search |
| `storage.py` | saving and loading graphs as JSON in `graph_data/` |
| `compare.py` | checks against Sage ($g = 1$), Richelot ($g = 2$) and the mass formula |
| `compare_type_one.py` | the type I graph next to the supersingular elliptic curve graph |

```
sage -python endomorphism_graph/compare.py 7 11 43 83            # ell = 2: Sage, Richelot, mass formula, products
sage -python endomorphism_graph/compare.py 11 43 83 --ell 3      # ell = 3

sage -python endomorphism_graph/storage.py 11 3                  # search g = 3, p = 11, ell = 2 and save it
sage -python endomorphism_graph/storage.py 43 2 3                # g = 2, p = 43, ell = 3
sage -python endomorphism_graph/storage.py 19 --weil             # g = 4, p = 19: only the Weil type fourfolds
sage -python endomorphism_graph/storage.py 19 --type-one         # g = 4, p = 19: only the type I graph
sage -python endomorphism_graph/storage.py --update graph_data/lattice_g4_p7_ell2.json   # recompute descriptions
```

Rough timings on a laptop, $\ell = 2$:

| graph | vertices | time |
|---|---|---|
| $g = 2$, $p = 131$ | 1008 | 40 s |
| $g = 3$, $p = 19$ | 109 | 30 s |
| $g = 4$, $p = 7$ | 21 | 1 min |
| $g = 4$, $p = 11$ | 215 | 20 min |
| $g = 4$, $p = 11$, Weil type only | 45 | 30 s |
| $g = 4$, $p = 19$, Weil type only | 182 | 2 min |

#### Weil type fourfolds

A Weil structure on a fourfold is an automorphism $J$ with $J^2 = -1$ that acts on the tangent space with eigenvalues $i, i, -i, -i$. Each one has a type, I or II, according to whether $x \mapsto \omega(x, Jx)$ vanishes on $A[2]$. In the graph of pairs $(A, J)$, a type I vertex has 151 $J$-stable Lagrangians and a type II vertex has 39. The *type I graph* consists of the type I vertices and the "middle steps" between them, the $J$-stable Lagrangians $K$ with $|K \cap A[1+i]| = 4$.

`storage.py p --weil` searches only the Weil-type graph, starting from $E_0^4$ with $J$ the rotation $(x_1, x_2, x_3, x_4) \mapsto (-x_2, x_1, -x_4, x_3)$. It is much cheaper than the full graph. For $p = 7$ and $11$ the full graph has also been computed, and the Weil-only search reaches every Weil structure it contains.

The type I/II split, the type I graph and anything observed about them are our own computations, not results from the literature.

### Saved graphs (`graph_data/`)

There is one JSON file per graph, `lattice_g{g}_p{p}_ell{ell}.json`, with the suffix `_weil` or `_type1` for the Weil-type graphs. Each file contains the edges $(v, w, \text{multiplicity})$ and, for each vertex:

- its canonical name;
- $|\mathrm{Aut}|$;
- its decomposition into indecomposable factors;
- its Frobenius twist, and whether it is defined over $\mathbb{F}_p$;
- its Weil structures;
- the lattice itself.

`storage.load_graph` reads a file back without any search.

### Pictures (`graph_painter/`)

```
python3 graph_painter/paint.py         # writes graph_painter/graphs.html from everything in graph_data/
open graph_painter/graphs.html
```

The page draws every saved graph:

- the colour of a vertex is its decomposition (indecomposable, $E \times$ threefold, surface $\times$ surface, ...);
- a ring marks a Weil structure, solid for type I and dashed for type II only;
- a dot marks a vertex defined over $\mathbb{F}_p$;
- hovering over an edge shows its multiplicities in both directions, and clicking a vertex labels all of its edges;
- for the Weil-type graphs, the "type I graph" checkbox dims everything of type II and keeps only the middle steps.
