"""
Saving and loading isogeny graphs, so they only have to be searched once.

The search is expensive (every neighbour is a lattice to build and identify) but its result is small: for each vertex
a 4g x 4g lattice basis and a few numbers, and the edges with multiplicities. We store this as plain JSON in
graph_data/ at the top of the repository, which graph_painter/ reads to draw the graphs.

Every vertex gets a name which only depends on its isomorphism class: a short hash of its canonical form (see
canonical.py), so the same vertex has the same name in every run, independent of the order of the search. In the
Weil type graphs a vertex is a pair (A, J) and its name is that of the fourfold A followed by /0, /1, ... for the
different Weil structures on it.

So, the pattern is:

    1. Search the graph with graph.isogeny_graph
    2. For each vertex record |Aut|, its decomposition type, its canonical name, its Frobenius twist (so whether it is
       defined over F_p), its Weil structures (see structures.py) and its lattice
    3. Write the vertices and the edges (as v, w, multiplicity) to graph_data/lattice_g{g}_p{p}_ell{ell}.json, with a
       suffix _weil or _type1 for the Weil type graphs
    4. load_graph reads the file back into a SageMath DiGraph (and the lattices, if wanted) without any search

Usage:

    sage -python storage.py p g [ell]                  the full graph (ell = 2 by default)
    sage -python storage.py p --weil                   only the Weil type fourfolds around (E0^4, J)
    sage -python storage.py p --type-one               only the type I graph: type I vertices and the middle steps
    sage -python storage.py --rename FILE [FILE ...]   recompute the names in saved graphs
    sage -python storage.py --update FILE [FILE ...]   recompute everything except the edges in saved (full) graphs

Add --no-names to skip the names, which are the slowest part of the description.
"""

import argparse
import hashlib
import json
import os
import time
from collections import Counter

import numpy
from sage.all import QQ, ZZ, DiGraph, matrix

from canonical import canonical_key
from graph import VertexSet, isogeny_graph, mass, mass_formula
from lattice import Ambient, HermitianLattice
from quaternions import EndomorphismRing
from structures import frobenius_map, rotation_structure, weil_class_type, weil_structures

DATA_DIRECTORY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "graph_data")


def data_path(p, g, ell=2, structure=None):
    suffix = f"_{structure}" if structure else ""
    return os.path.join(DATA_DIRECTORY, f"lattice_g{g}_p{p}_ell{ell}{suffix}.json")


# ============================================================
# Pieces of a record
# ============================================================


def canonical_name(L):
    """
    A 12 character hexadecimal name which only depends on the isomorphism class of L (a hash of its canonical form)
    """
    return hashlib.sha256(repr(canonical_key(L)).encode()).hexdigest()[:12]


def lattice_record(L):
    """
    The lattice as integers: rows of d * M, the denominator d and the scale n of the form h / n
    """
    d = L.M.denominator()
    return {"scale": int(L.scale), "denominator": int(d), "basis": [[int(x) for x in row] for row in (d * L.M).rows()]}


def edge_list(G):
    """
    The edges of the multigraph as [v, w, multiplicity]
    """
    multiplicities = Counter((v, w) for v, w, _ in G.edges(sort=False))
    return [[int(v), int(w), int(m)] for (v, w), m in sorted(multiplicities.items())]


def spectrum(G, n):
    """
    The eigenvalues of the adjacency matrix, largest first (real, since the matrix is reversible:
    m(v, w) |Aut w| = m(w, v) |Aut v|)
    """
    A = numpy.zeros((n, n))
    for v, w, _ in G.edges(sort=False):
        A[v, w] += 1
    return [round(float(x), 4) for x in sorted(numpy.linalg.eigvals(A).real, reverse=True)]


def write_json(data, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        json.dump(data, f, separators=(",", ":"))
    os.replace(path + ".tmp", path)  # never leave a half written file behind


def progress(verbose, v, total, t0):
    if verbose and (v + 1) % 50 == 0:
        print(f"  described {v + 1} / {total} vertices ({time.time() - t0:.1f}s)")


# ============================================================
# Saving
# ============================================================


def describe_vertices(vertices, names=True, verbose=False):
    """
    The JSON record of every vertex of a full graph
    """
    frobenius = frobenius_map(vertices)
    records, t0 = [], time.time()
    for v, L in enumerate(vertices.lattices):
        records.append({
            "id": v,
            "name": canonical_name(L) if names else None,
            "aut": int(vertices.automorphism_count(v)),
            "decomposition": [int(d) for d in vertices.decomposition(v)],
            "frobenius": int(frobenius[v]),
            "over_Fp": frobenius[v] == v,
            "weil": weil_structures(vertices, v),
            "lattice": lattice_record(L),
        })
        progress(verbose, v, len(vertices.lattices), t0)
    return records


def describe_weil_vertices(vertices, J, names=True, verbose=False):
    """
    The JSON record of every vertex (A, J) of a Weil type graph: its automorphisms are those commuting with J, "weil"
    is the type (I or II) of its J, and "fourfold" says which vertices have the same underlying fourfold A
    """
    fourfolds = vertices.varieties
    frobenius = frobenius_map(vertices)
    fourfold_names, count = {}, Counter()
    records, t0 = [], time.time()
    for v, L in enumerate(vertices.lattices):
        f = vertices.variety[v]
        if f not in fourfold_names:
            fourfold_names[f] = canonical_name(fourfolds.lattices[f]) if names else f"A{f}"
        records.append({
            "id": v,
            "name": f"{fourfold_names[f]}/{count[f]}",  # the fourfold, and which of its Weil structures
            "aut": int(vertices.automorphism_count(v)),
            "decomposition": [int(d) for d in fourfolds.decomposition(f)],
            "frobenius": int(frobenius[v]),
            "over_Fp": frobenius[v] == v,
            "weil": [weil_class_type(L, L.column_matrix(J))],
            "fourfold": int(f),
            "lattice": lattice_record(L),
        })
        count[f] += 1
        progress(verbose, v, len(vertices.lattices), t0)
    return records


def save_graph(p, g, ell, vertices, G, names=True, verbose=False):
    """
    Write a full graph from isogeny_graph to JSON
    """
    total = mass(vertices)
    data = {
        "method": "lattice",
        "p": int(p),
        "g": int(g),
        "ell": int(ell),
        "start": 0,
        "mass": str(total),
        "mass_formula_ok": bool(total == mass_formula(p, g)),
        "vertices": describe_vertices(vertices, names, verbose),
        "edges": edge_list(G),
    }
    path = data_path(p, g, ell)
    write_json(data, path)
    return path


def save_weil_graph(p, vertices, G, J, structure, names=True, verbose=False):
    """
    Write a Weil type graph (structure "weil": all J-stable kernels; "type1": the type I graph) to JSON. The edges
    also carry their kernel classes |K cap A[1 + i]| (4 for the middle steps).
    """
    classes = Counter((v, w, size) for v, w, size in G.edges(sort=False))
    data = {
        "method": "lattice",
        "structure": structure,
        "p": int(p),
        "g": 4,
        "ell": 2,
        "start": 0,
        "mass": str(mass(vertices)),
        "fourfolds": len(vertices.varieties.lattices),
        "vertices": describe_weil_vertices(vertices, J, names, verbose),
        "edges": edge_list(G),
        "kernel_classes": [[int(v), int(w), int(s), int(m)] for (v, w, s), m in sorted(classes.items())],
        "eigenvalues": spectrum(G, len(vertices.lattices)),
    }
    path = data_path(p, 4, 2, structure)
    write_json(data, path)
    return path


# ============================================================
# Loading and updating
# ============================================================


def load_graph(p, g, ell=2, path=None, lattices=False):
    """
    Read a saved graph: returns (data, G) with data the JSON dictionary and G the directed multigraph. With
    lattices=True, also rebuild the HermitianLattice of every vertex as data["lattices"].
    """
    with open(path or data_path(p, g, ell)) as f:
        data = json.load(f)
    edges = [(v, w) for v, w, m in data["edges"] for _ in range(m)]
    G = DiGraph([list(range(len(data["vertices"]))), edges], format="vertices_and_edges", multiedges=True, loops=True)
    if lattices:
        ambient = Ambient(EndomorphismRing(data["p"]), data["g"])
        data["lattices"] = [
            HermitianLattice(ambient, matrix(ZZ, r["lattice"]["basis"]).change_ring(QQ) / r["lattice"]["denominator"],
                             r["lattice"]["scale"])
            for r in data["vertices"]
        ]
    return data, G


def rename_file(path):
    """
    Recompute the canonical names in a saved graph (for example after a change to canonical.py)
    """
    data, _ = load_graph(None, None, path=path, lattices=True)
    lattices = data.pop("lattices")
    if data.get("structure"):  # vertices (A, J): the name of the fourfold A, then /k
        names = {}
        for r, L in zip(data["vertices"], lattices):
            if r["fourfold"] not in names:
                names[r["fourfold"]] = canonical_name(L)
            r["name"] = names[r["fourfold"]] + "/" + r["name"].split("/")[1]
    else:
        for r, L in zip(data["vertices"], lattices):
            r["name"] = canonical_name(L)
    write_json(data, path)


def update_file(path, names=True, verbose=False):
    """
    Recompute the vertex descriptions of a saved full graph, without searching it again
    """
    data, _ = load_graph(None, None, path=path, lattices=True)
    vertices = VertexSet()
    for L in data.pop("lattices"):
        _, new, _ = vertices.index(L)
        assert new, "two saved vertices are isomorphic"
    data["vertices"] = describe_vertices(vertices, names, verbose)
    write_json(data, path)


def cached_isogeny_graph(p, g, ell=2, names=True, verbose=False):
    """
    Load the graph from graph_data/ if it is there, otherwise search for it and save it
    """
    if not os.path.exists(data_path(p, g, ell)):
        vertices, G = isogeny_graph(EndomorphismRing(p), g, ell, verbose=verbose)
        save_graph(p, g, ell, vertices, G, names=names, verbose=verbose)
    return load_graph(p, g, ell)


# ============================================================
# Command line
# ============================================================


def main():
    parser = argparse.ArgumentParser(description="Search, save and update superspecial isogeny graphs")
    parser.add_argument("arguments", nargs="*", help="p g [ell]; or p with --weil / --type-one; or files")
    parser.add_argument("--weil", action="store_true", help="only the Weil type fourfolds around (E0^4, J)")
    parser.add_argument("--type-one", action="store_true", help="only the type I graph")
    parser.add_argument("--rename", action="store_true", help="recompute the names in the given files")
    parser.add_argument("--update", action="store_true", help="recompute the descriptions in the given files")
    parser.add_argument("--no-names", action="store_true", help="do not compute canonical names")
    args = parser.parse_args()
    names = not args.no_names
    t0 = time.time()

    if args.rename or args.update:
        for path in args.arguments:
            if args.rename:
                rename_file(path)
            else:
                update_file(path, names, verbose=True)
            print(f"{'renamed' if args.rename else 'updated'} {path} ({time.time() - t0:.1f}s)")
        return

    if args.weil or args.type_one:
        p = int(args.arguments[0])
        ring = EndomorphismRing(p)
        J = rotation_structure(Ambient(ring, 4))
        vertices, G = isogeny_graph(ring, 4, 2, structure=J, type_one=args.type_one, verbose=True)
        t1 = time.time()
        path = save_weil_graph(p, vertices, G, J, "type1" if args.type_one else "weil", names, verbose=True)
    else:
        p, g, ell = (list(map(int, args.arguments)) + [2])[:3]
        vertices, G = isogeny_graph(EndomorphismRing(p), g, ell, verbose=True)
        t1 = time.time()
        path = save_graph(p, g, ell, vertices, G, names, verbose=True)
    print(f"p = {p}: {G.num_verts()} vertices, search {t1 - t0:.1f}s, description {time.time() - t1:.1f}s")
    print(f"saved to {os.path.normpath(path)}")


if __name__ == "__main__":
    main()
