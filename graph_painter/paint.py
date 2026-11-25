"""
Paint the isogeny graphs saved in graph_data/ as one interactive web page.

Reads every graph_data/*.json written by endomorphism_graph/storage.py, works out what the picture needs (an
undirected simple graph with both edge multiplicities, distances from the start vertex, summary statistics) and embeds
it into template.html. The output is a single file, graph_painter/graphs.html, which opens in any browser (it only
loads d3 from a CDN).

So, the pattern is:

    1. Load each saved graph (no SageMath needed: this is plain Python)
    2. Merge the directed edges v -> w and w -> v into one link carrying both multiplicities (and, for the Weil type
       graphs, both numbers of middle steps)
    3. Compute distances from E0^g by breadth first search, and the statistics shown in the side panel
    4. Write the data into the template

Usage: python3 paint.py [graph_data/lattice_g2_p83_ell2.json ...]    (default: every file in graph_data/)
"""

import glob
import json
import os
import sys
from collections import Counter, deque

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIRECTORY = os.path.join(HERE, "..", "graph_data")
TEMPLATE = os.path.join(HERE, "template.html")
OUTPUT = os.path.join(HERE, "graphs.html")

TITLE_SUFFIX = {None: "", "weil": ", Weil type only", "type1": ", type I graph"}
ORDER = {None: 0, "weil": 1, "type1": 2}


# ============================================================
# Small graph helpers
# ============================================================


def undirected(n, pairs):
    """
    The adjacency sets of the undirected simple graph on range(n) with an edge for every (v, w), v != w
    """
    adjacency = [set() for _ in range(n)]
    for v, w in pairs:
        if v != w:
            adjacency[v].add(w)
            adjacency[w].add(v)
    return adjacency


def distances(adjacency, start):
    """
    Breadth first search distances from the start vertex (None if unreachable)
    """
    dist = [None] * len(adjacency)
    dist[start] = 0
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for w in adjacency[v]:
            if dist[w] is None:
                dist[w] = dist[v] + 1
                queue.append(w)
    return dist


def component_count(adjacency, among):
    """
    The number of connected components of the graph induced on the vertices among
    """
    seen, count = set(), 0
    for s in among:
        if s in seen:
            continue
        count += 1
        stack = [s]
        seen.add(s)
        while stack:
            v = stack.pop()
            for w in adjacency[v]:
                if w in among and w not in seen:
                    seen.add(w)
                    stack.append(w)
    return count


def decomposition_label(decomposition, g):
    """
    A readable name for a decomposition type, e.g. (1, 1, 2) -> "E x E x surface"
    """
    if len(decomposition) == 1:
        return {1: "elliptic curve", 2: "Jacobian (indecomposable)"}.get(g, f"indecomposable {g}-fold")
    return " x ".join({1: "E", 2: "surface", 3: "threefold"}.get(d, f"{d}-fold") for d in decomposition)


# ============================================================
# Graph data for the picture
# ============================================================


def middle_steps(data):
    """
    {(v, w): number of middle steps (|K cap A[1 + i]| = 4) from v to w}, between type I vertices (Weil type graphs)
    """
    weil_type = {r["id"]: (r.get("weil") or [None])[0] for r in data["vertices"]}
    return {(v, w): m for v, w, size, m in data.get("kernel_classes", [])
            if size == 4 and weil_type[v] == "I" and weil_type[w] == "I"}


def links(multiplicity, middle):
    """
    One link per unordered pair {a, b} (a < b) with an edge either way: [a, b, m(a -> b), m(b -> a), then the same
    for the middle steps]
    """
    pairs = sorted({(min(v, w), max(v, w)) for v, w in multiplicity if v != w})
    both = lambda counts, a, b: [counts.get((a, b), 0), counts.get((b, a), 0)]
    return [[a, b] + both(multiplicity, a, b) + both(middle, a, b) for a, b in pairs]


def vertex_records(data, multiplicity, middle, dist):
    """
    What the page needs to know about each vertex
    """
    return [{
        "id": r["id"],
        "name": r["name"] or f"v{r['id']}",
        "aut": r["aut"],
        "type": decomposition_label(r["decomposition"], data["g"]),
        "distance": dist[r["id"]],
        "loops": multiplicity.get((r["id"], r["id"]), 0),
        "middleLoops": middle.get((r["id"], r["id"]), 0),
        "frobenius": r.get("frobenius"),
        "overFp": bool(r.get("over_Fp")),
        "weil": ", ".join(r.get("weil") or []) or "none",
    } for r in data["vertices"]]


def statistics(data, vertices, multiplicity, middle, dist):
    """
    The table in the side panel
    """
    structure, g, start = data.get("structure"), data["g"], data["start"]
    out_degree = Counter()
    for (v, w), m in multiplicity.items():
        out_degree[v] += m
    stats = {
        "vertices": len(vertices),
        "kernels at E0^g" if structure else "kernels per vertex": out_degree[start],
        "neighbour pairs": len(links(multiplicity, {})),
        "vertices with loops": sum(1 for x in vertices if x["loops"]),
        "defined over F_p": sum(1 for x in vertices if x["overFp"]),
        "distance from E0^g (max)": max(d for d in dist if d is not None),
        "mass (sum 1/|Aut|)": data["mass"],
    }
    if "mass_formula_ok" in data:  # there is no mass formula for the Weil type graphs
        stats["equals the mass formula"] = "yes" if data["mass_formula_ok"] else "no"

    if structure == "weil":
        type_one = {x["id"] for x in vertices if x["weil"] == "I"}
        stats["type I vertices"] = len(type_one)
        if middle:
            stats["middle steps per type I vertex"] = sum(m for (v, w), m in middle.items() if v == start)
            stats["components of the type I graph"] = component_count(undirected(len(vertices), middle), type_one)
    if structure == "type1":
        stats.pop("kernels at E0^g")
        stats["middle steps per vertex"] = out_degree[start]
        stats["components"] = component_count(undirected(len(vertices), multiplicity), set(range(len(vertices))))
    if structure:
        stats["underlying fourfolds"] = data["fourfolds"]
    elif g % 2 == 0:
        stats["with a Weil structure"] = sum(1 for x in vertices if x["weil"] != "none")

    if "eigenvalues" in data:
        eigenvalues = data["eigenvalues"]
        stats["largest eigenvalue"] = eigenvalues[0]
        stats["second |eigenvalue|"] = max((abs(x) for x in eigenvalues[1:]), default=0)
        degrees = set(out_degree[v] for v in range(len(vertices)))
        if len(degrees) == 1:  # regular: compare with the Ramanujan bound
            stats["Ramanujan bound 2 sqrt(d - 1)"] = round(2 * (degrees.pop() - 1) ** 0.5, 2)
    return stats


def prepare(data):
    """
    The dictionary used by template.html for one saved graph
    """
    multiplicity = {(v, w): m for v, w, m in data["edges"]}
    middle = middle_steps(data)
    dist = distances(undirected(len(data["vertices"]), multiplicity), data["start"])
    vertices = vertex_records(data, multiplicity, middle, dist)
    structure = data.get("structure")
    return {
        "title": f"g = {data['g']}, p = {data['p']}, ℓ = {data['ell']}" + TITLE_SUFFIX[structure],
        "sort": (ORDER[structure], data["ell"], data["g"], data["p"]),
        "structure": structure,
        "hasClasses": "kernel_classes" in data,
        "g": data["g"],
        "start": data["start"],
        "vertices": vertices,
        "links": links(multiplicity, middle),
        "stats": statistics(data, vertices, multiplicity, middle, dist),
        "types": sorted(Counter(x["type"] for x in vertices).items(), key=lambda t: (-len(t[0]), t[0])),
    }


# ============================================================
# Main
# ============================================================


def main(paths):
    graphs = []
    for path in paths:
        with open(path) as f:
            graphs.append(prepare(json.load(f)))
    graphs.sort(key=lambda G: G["sort"])
    for G in graphs:
        del G["sort"]
    with open(TEMPLATE) as f:
        template = f.read()
    page = template.replace("/*GRAPH_DATA*/[]", json.dumps(graphs, separators=(",", ":")))
    with open(OUTPUT, "w") as f:
        f.write(page)
    print(f"painted {len(graphs)} graphs into {os.path.normpath(OUTPUT)}")


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(glob.glob(os.path.join(DATA_DIRECTORY, "*.json"))))
