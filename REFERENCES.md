# References

Every entry below was checked against an online record (arXiv, a DOI via Crossref, zbMATH, Numdam, or the project's own website) in November 2026. Entries marked **check** have a detail that could not be confirmed and should be verified by hand before being cited elsewhere.

The keys (for example `[Chi14]`) are the ones used in the docstrings and the READMEs.

## Where the method comes from

The graph search in `endomorphism_graph/` is not a new algorithm. Its pieces come from the following sources.

| What the code does | Where it comes from |
|---|---|
| A superspecial principally polarised abelian variety A is recorded by the quaternion Hermitian lattice Hom(E0, A), and principal polarisations correspond to unimodular forms | [IKO86, Thm 2.10]; the equivalence of categories A ↦ Hom(E, A) is [JKPRST18] |
| (ℓ, ..., ℓ)-isogenies correspond to ℓ-neighbours of these lattices, and the isogeny graph to a Brandt-matrix (Hecke) graph | [JZ26]; see also [Jor21]. **check:** we have not yet compared their setup with ours line by line |
| ℓ-neighbours of a quaternion Hermitian lattice correspond to Lagrangians of a symplectic space over F_ℓ, via the Morita equivalence O/ℓO = M_2(F_ℓ) (our A[ℓ] = (L/ℓL)e and the Weil pairing read off from h) | Kneser's neighbour method [Kne57], for quaternionic lattices [Bac95], in exactly this symplectic form in [Chi14, §2.4, §3.5–3.7, Lemma 83] |
| Testing isometry of quaternion Hermitian lattices with auxiliary Z-bilinear forms that force O-linearity | [GV14, Lemma 6.2], [Chi14, §3.8.1] |
| Isometry testing and automorphism groups of Z-lattices (PARI's `qfisom`, `qfauto`) | [PS97], implemented in [PARI] |
| Using the mass formula to check that every class has been found | [Chi14, §3.8.2], [GV14] |
| The mass Σ 1/\|Aut\| of superspecial principally polarised abelian varieties | [Eke87] for all g, [HI80] for g = 2; see also [Yu08], [Yu10]. **check:** theorem numbers in [Eke87] |
| The canonical form of a lattice (characteristic vectors, then a canonical labelling of a coloured graph) | [DHVW20]; graph canonical labelling with bliss [JK07], [JK11] |
| The comparison graph in dimension two (Richelot isogenies) | [CDS20], [FS22]; edge cases in [Zob21] |

The following are, as far as we know, not in the literature: the type I / type II split of Weil structures on superspecial fourfolds, and the type I graph. They are computational observations, checked only for the small primes in `graph_data/`. If you know of prior work on them, please tell us.

### A note on a published table

[Chi14, §3.8.3] tabulates the class numbers of the principal genus of O × O for every prime discriminant below 100. For p = 43 and p = 83 it lists 52 and 271 classes, with masses 1295/96 and 28249/288. Our lattice search and our independent Richelot computation both find 55 and 296 classes, with the same masses. The class number formula of [IKO86, Thm 3.3, Cor. 2.12], in the form restated in [KT20, §2] and [Kud21, Thm 2.5], also gives 55 and 296, with 10 and 36 of them products of elliptic curves. We have not investigated the cause of the difference.

## Lattices and the neighbour method

- **[Chi14]** S. Chisholm, *Algorithmic enumeration of quaternionic lattices*, PhD thesis, University of Calgary, 2014. Supervised by M. Bauer and M. Greenberg. DOI [10.11575/PRISM/24766](https://doi.org/10.11575/PRISM/24766), handle [11023/1920](https://prism.ucalgary.ca/handle/11023/1920).
- **[Kne57]** M. Kneser, Klassenzahlen definiter quadratischer Formen, *Arch. Math.* 8 (1957), no. 4, 241–250. DOI [10.1007/BF01898782](https://doi.org/10.1007/BF01898782).
- **[Bac95]** C. Bachoc, Voisinage au sens de Kneser pour les réseaux quaternioniens, *Comment. Math. Helv.* 70 (1995), no. 1, 350–374. DOI [10.1007/BF02566012](https://doi.org/10.1007/BF02566012).
- **[Sch98]** A. Schiemann, Classification of Hermitian forms with the neighbour method, *J. Symbolic Comput.* 26 (1998), no. 4, 487–508. DOI [10.1006/jsco.1998.0225](https://doi.org/10.1006/jsco.1998.0225).
- **[GV14]** M. Greenberg, J. Voight, Lattice methods for algebraic modular forms on classical groups, in *Computations with Modular Forms*, Contrib. Math. Comput. Sci., Springer, 2014, 147–179. DOI [10.1007/978-3-319-03847-6_6](https://doi.org/10.1007/978-3-319-03847-6_6), arXiv:[1209.2460](https://arxiv.org/abs/1209.2460). **check:** series volume number.
- **[PS97]** W. Plesken, B. Souvignier, Computing isometries of lattices, *J. Symbolic Comput.* 24 (1997), no. 3–4, 327–334. DOI [10.1006/jsco.1996.0130](https://doi.org/10.1006/jsco.1996.0130).
- **[Shi63]** G. Shimura, Arithmetic of alternating forms and quaternion hermitian forms, *J. Math. Soc. Japan* 15 (1963), no. 1, 33–65. DOI [10.2969/jmsj/01510033](https://doi.org/10.2969/jmsj/01510033).
- **[Iya69]** K. Iyanaga, Class numbers of definite Hermitian forms, *J. Math. Soc. Japan* 21 (1969), no. 3, 359–374. DOI [10.2969/jmsj/02130359](https://doi.org/10.2969/jmsj/02130359). Unimodular Hermitian forms over Z[i].

## Canonical forms

- **[DHVW20]** M. Dutour Sikirić, A. Haensch, J. Voight, W. P. J. van Woerden, A canonical form for positive definite matrices, *ANTS XIV*, Open Book Series 4 (2020), 179–195. DOI [10.2140/obs.2020.4.179](https://doi.org/10.2140/obs.2020.4.179), arXiv:[2004.14022](https://arxiv.org/abs/2004.14022).
- **[JK07]** T. Junttila, P. Kaski, Engineering an efficient canonical labeling tool for large and sparse graphs, *Proc. ALENEX 2007*, SIAM, 135–149. DOI [10.1137/1.9781611972870.13](https://doi.org/10.1137/1.9781611972870.13).
- **[JK11]** T. Junttila, P. Kaski, Conflict propagation and component recursion for canonical labeling, *TAPAS 2011*, LNCS, Springer, 151–162. DOI [10.1007/978-3-642-19754-3_16](https://doi.org/10.1007/978-3-642-19754-3_16). **check:** LNCS volume number.

## Superspecial abelian varieties and mass formulas

- **[IKO86]** T. Ibukiyama, T. Katsura, F. Oort, Supersingular curves of genus two and class numbers, *Compositio Math.* 57 (1986), no. 2, 127–152. [Numdam](https://numdam.org/item/CM_1986__57_2_127_0/). Thm 2.10 (p. 144), Cor. 2.12 (p. 145), Thm 3.3 (pp. 150–151).
- **[JKPRST18]** B. W. Jordan, A. G. Keeton, B. Poonen, E. M. Rains, N. Shepherd-Barron, J. T. Tate, Abelian varieties isogenous to a power of an elliptic curve, *Compositio Math.* 154 (2018), 934–959. DOI [10.1112/S0010437X17007990](https://doi.org/10.1112/S0010437X17007990), arXiv:[1602.06237](https://arxiv.org/abs/1602.06237).
- **[Eke87]** T. Ekedahl, On supersingular curves and abelian varieties, *Math. Scand.* 60 (1987), 151–178. DOI [10.7146/math.scand.a-12178](https://doi.org/10.7146/math.scand.a-12178).
- **[HI80]** K. Hashimoto, T. Ibukiyama, On class numbers of positive definite binary quaternion hermitian forms, *J. Fac. Sci. Univ. Tokyo Sect. IA Math.* 27 (1980), 549–601; II, 28 (1981), no. 3, 695–699; III, 30 (1983), 393–401. **check:** part II is dated 1982 in some bibliographies.
- **[Yu08]** C.-F. Yu, An exact geometric mass formula, *Int. Math. Res. Not. IMRN* 2008, Art. ID rnn113. DOI [10.1093/imrn/rnn113](https://doi.org/10.1093/imrn/rnn113), arXiv:[0706.3474](https://arxiv.org/abs/0706.3474).
- **[Yu10]** C.-F. Yu, Simple mass formulas on Shimura varieties of PEL-type, *Forum Math.* 22 (2010), no. 3, 565–582. DOI [10.1515/FORUM.2010.030](https://doi.org/10.1515/FORUM.2010.030), arXiv:[math/0603451](https://arxiv.org/abs/math/0603451).
- **[GHY01]** W. T. Gan, J. P. Hanke, J.-K. Yu, On an exact mass formula of Shimura, *Duke Math. J.* 107 (2001), no. 1, 103–133. DOI [10.1215/S0012-7094-01-10716-3](https://doi.org/10.1215/S0012-7094-01-10716-3).
- **[GY00]** W. T. Gan, J.-K. Yu, Group schemes and local densities, *Duke Math. J.* 105 (2000), no. 3, 497–524. DOI [10.1215/S0012-7094-00-10535-2](https://doi.org/10.1215/S0012-7094-00-10535-2).
- **[IKY25]** T. Ibukiyama, V. Karemaker, C.-F. Yu, When is a polarised abelian variety determined by its p-divisible group?, *Trans. Amer. Math. Soc. Ser. B* 12 (2025), 65–111. DOI [10.1090/btran/222](https://doi.org/10.1090/btran/222), arXiv:[2205.13180](https://arxiv.org/abs/2205.13180).
- **[KYY21]** V. Karemaker, F. Yobuko, C.-F. Yu, Mass formula and Oort's conjecture for supersingular abelian threefolds, *Adv. Math.* 386 (2021), 107812. DOI [10.1016/j.aim.2021.107812](https://doi.org/10.1016/j.aim.2021.107812), arXiv:[2002.10960](https://arxiv.org/abs/2002.10960).
- **[KT20]** T. Katsura, K. Takashima, Counting Richelot isogenies between superspecial abelian surfaces, *ANTS XIV* (2020), arXiv:[2003.00633](https://arxiv.org/abs/2003.00633). Used here only for its restatement of [IKO86, Thm 3.3].
- **[Kud21]** M. Kudo, Counting isomorphism classes of superspecial curves, arXiv:[2106.12409](https://arxiv.org/abs/2106.12409). Used here only for Thm 2.5.

## Isogeny graphs

- **[JZ26]** B. W. Jordan, Y. Zaytman, Isogeny graphs of superspecial abelian varieties and Brandt matrices, *J. Math. Soc. Japan* 78 (2026), no. 4, 1035–1074. DOI [10.2969/jmsj/91629162](https://doi.org/10.2969/jmsj/91629162), arXiv:[2005.09031](https://arxiv.org/abs/2005.09031).
- **[Jor21]** B. W. Jordan, Isogeny graphs of superspecial abelian varieties, arXiv:[2101.01579](https://arxiv.org/abs/2101.01579) (lecture notes).
- **[JZ22]** B. W. Jordan, Y. Zaytman, Isogeny complexes of superspecial abelian varieties, arXiv:[2205.07383](https://arxiv.org/abs/2205.07383).
- **[ATY24]** Y. Aikawa, R. Tanaka, T. Yamauchi, Isogeny graphs on superspecial abelian varieties: eigenvalues and connection to Bruhat–Tits buildings, *Canad. J. Math.* 76 (2024), no. 6, 1891–1916. DOI [10.4153/S0008414X23000676](https://doi.org/10.4153/S0008414X23000676), arXiv:[2201.04293](https://arxiv.org/abs/2201.04293). Relevant to the spectra of our graphs.
- **[Piz90]** A. K. Pizer, Ramanujan graphs and Hecke operators, *Bull. Amer. Math. Soc. (N.S.)* 23 (1990), no. 1, 127–137. DOI [10.1090/S0273-0979-1990-15918-X](https://doi.org/10.1090/S0273-0979-1990-15918-X).
- **[Mes86]** J.-F. Mestre, La méthode des graphes. Exemples et applications, *Proc. Internat. Conf. on Class Numbers and Fundamental Units of Algebraic Number Fields* (Katata, 1986), 217–242. zbMATH 0621.14021. **check:** French title, not shown in the zbMATH record.
- **[CDS20]** W. Castryck, T. Decru, B. Smith, Hash functions from superspecial genus-2 curves using Richelot isogenies, *J. Math. Cryptol.* 14 (2020), no. 1, 268–292. DOI [10.1515/jmc-2019-0021](https://doi.org/10.1515/jmc-2019-0021), arXiv:[1903.06451](https://arxiv.org/abs/1903.06451).
- **[FS22]** E. Florit, B. Smith, Automorphisms and isogeny graphs of abelian varieties, with applications to the superspecial Richelot isogeny graph, in *Arithmetic, Geometry, Cryptography, and Coding Theory 2021*, Contemp. Math. 779, AMS, 2022, 103–132. DOI [10.1090/conm/779/15672](https://doi.org/10.1090/conm/779/15672), arXiv:[2101.00919](https://arxiv.org/abs/2101.00919).
- **[Zob21]** L. Zobernig, Genus 2 curves in small characteristic, arXiv:[2111.07270](https://arxiv.org/abs/2111.07270).
- **[OOKYN25]** R. Ohashi, H. Onuki, M. Kudo, R. Yoshizumi, K. Nuida, Listing superspecial curves of genus three using Richelot isogeny graphs, *Res. Number Theory* 11 (2025), Art. 76. DOI [10.1007/s40993-025-00657-3](https://doi.org/10.1007/s40993-025-00657-3), arXiv:[2401.10500](https://arxiv.org/abs/2401.10500). The dimension three graph, computed with curves.

## Software

- **[Sage]** The Sage Developers, *SageMath, the Sage Mathematics Software System* (Version 10.9), https://www.sagemath.org.
- **[PARI]** The PARI Group, *PARI/GP*, Univ. Bordeaux, https://pari.math.u-bordeaux.fr/. Used through SageMath; cite the version you run.
- bliss [JK07], [JK11], used through SageMath's `canonical_label` and `automorphism_group`.
