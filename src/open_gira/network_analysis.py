"""
Network + spatial metrics for spatial energy transmission/distribution grids.
"""

import logging
import time
from functools import wraps

import numpy as np
import pandas as pd
import networkx as nx
import networkit as nk


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:
            elapsed = time.perf_counter() - start
            logging.info(f"{func.__name__} took {elapsed:.4f} seconds")

    return wrapper


@timer
def _nx_graph(
    edges, nodes, from_col="from_id", to_col="to_id",
    geometry_col="geometry", id_col="id"
):
    """Create a networkx graph from DataFrames."""
    G = nx.Graph()

    for _, row in nodes.iterrows():
        geom = row[geometry_col]
        G.add_node(row[id_col], pos=(geom.x, geom.y))

    lengths = edges.to_crs("EPSG:3857").geometry.length.values

    for (_, row), length_m in zip(edges.iterrows(), lengths):
        u, v = row[from_col], row[to_col]
        if u == v:
            continue
        if G.has_edge(u, v):
            if length_m < G[u][v].get("length_m", np.inf):
                G[u][v]["length_m"] = length_m
            continue
        G.add_edge(u, v, length_m=length_m, edge_id=row.get(id_col, None))

    return G


@timer
def _nk_graph(G, weight=None):
    """Convert a networkx graph to a NetworKit graph."""
    node_list = list(G.nodes())
    id_to_idx = {nid: i for i, nid in enumerate(node_list)}
    nkG = nk.Graph(len(node_list), weighted=(weight is not None), directed=False)
    for u, v, data in G.edges(data=True):
        w = data.get(weight, 1.0) if weight else 1.0
        nkG.addEdge(id_to_idx[u], id_to_idx[v], w)
    return nkG, node_list


@timer
def basic_topology(G):
    N, E = G.number_of_nodes(), G.number_of_edges()
    degrees = np.array([d for _, d in G.degree()])
    return {
        "n_nodes": N,
        "n_edges": E,
        "avg_degree": degrees.mean() if N else np.nan,
        "degree_std": degrees.std() if N else np.nan,
        "max_degree": degrees.max() if N else np.nan,
        "degree_assortativity": nx.degree_assortativity_coefficient(G) if E else np.nan,
    }


@timer
def _nk_topo_stats(
        Gc, n_threads, rel_tol=0.01, batch_size=50,
        min_samples=100, max_samples=2000, seed=42
    ):
    nk.setNumberOfThreads(n_threads)
    nkG, _ = _nk_graph(Gc, weight="length_m")
    n = nkG.numberOfNodes()

    logging.info("Calculate graph 'diameter' (longest shortest path)")
    diam = nk.distance.Diameter(nkG, algo=nk.distance.DiameterAlgo.EXACT)
    diam.run()
    diameter = diam.getDiameter()[0]

    rng = np.random.default_rng(seed)
    order = rng.permutation(n)

    per_source_means = []
    inv_means = []
    n_used = 0
    idx = 0
    logging.info(f"Sample avg. shortest path until error < {rel_tol * 100}%")
    while idx < min(max_samples, n):
        batch_end = min(idx + batch_size, max_samples, n)
        for s in order[idx:batch_end]:
            bfs = nk.distance.BFS(nkG, int(s), storePaths=False)
            bfs.run()
            d = np.array(bfs.getDistances())
            d = np.delete(d, s)  # drop self-distance
            per_source_means.append(d.mean())
            inv_means.append(np.mean(1.0 / d))
            n_used += 1
        idx = batch_end
        logging.info(batch_end)

        if n_used >= min_samples:
            arr = np.array(per_source_means)
            se = arr.std(ddof=1) / np.sqrt(n_used)
            if se / arr.mean() < rel_tol:
                break
    logging.info("Converged")

    avg_dist = float(np.mean(per_source_means))
    avg_inv = float(np.mean(inv_means))

    return {
        "avg_shortest_path_len": avg_dist,
        "diameter": diameter,
        "global_efficiency": avg_inv,
        "topo_stats_n_samples": n_used,
    }


@timer
def _nk_circuity(Gc, pos, pairs, n_threads):
    nkG, node_list = _nk_graph(Gc, weight="length_m")
    id_to_idx = {nid: i for i, nid in enumerate(node_list)}

    sources = sorted({id_to_idx[u] for u, v in pairs})
    spsp = nk.distance.SPSP(nkG, sources)
    spsp.run()

    circuities = []
    for u, v in pairs:
        graph_dist = spsp.getDistance(id_to_idx[u], id_to_idx[v])
        if not np.isfinite(graph_dist):
            continue
        lon1, lat1 = pos[u]
        lon2, lat2 = pos[v]
        euclid_dist_m = haversine_m(lat1, lon1, lat2, lon2)
        if euclid_dist_m > 0:
            circuities.append(graph_dist / euclid_dist_m)

    return float(np.mean(circuities)) if circuities else np.nan


@timer
def connectivity_metrics(G, n_workers=1):
    N = G.number_of_nodes()
    logging.info("Find connected components")
    components = list(nx.connected_components(G))
    largest_cc_nodes = max(components, key=len)
    logging.info("Copy largest component")
    Gc = G.subgraph(largest_cc_nodes).copy()  # work on giant component for path-based stats

    out = {
        "n_components": len(components),
        "largest_cc_fraction": len(largest_cc_nodes) / N if N else np.nan,
    }

    if Gc.number_of_nodes() > 1:
        logging.info("Large network: use networkit")
        out.update(_nk_topo_stats(Gc, n_workers))
        # Algebraic connectivity, on giant component only
        logging.info("Algebraic connectivity")
        # N.B. tracemin_lu 100x faster than default (tracemin_pcg) for gridfinder type networks
        out["algebraic_connectivity"] = nx.algebraic_connectivity(Gc, method="tracemin_lu", tol=1e-6)
        out["algebraic_connectivity_norm"] = out["algebraic_connectivity"] / Gc.number_of_nodes()
    else:
        out.update({k: np.nan for k in
                    ["avg_shortest_path_len", "diameter", "global_efficiency",
                     "local_efficiency", "algebraic_connectivity",
                     "algebraic_connectivity_norm"]})

    return out


@timer
def meshedness_metrics(G):
    N, E = G.number_of_nodes(), G.number_of_edges()
    n_components = nx.number_connected_components(G)

    mu = E - N + n_components  # Cyclomatic number, generalised for multiple components
    alpha = mu / (2 * N - 5 * n_components) if N > 2 else np.nan
    beta = E / N if N else np.nan
    gamma = E / (3 * (N - 2 * n_components)) if N > 2 * n_components else np.nan

    return {
        "cyclomatic_number": mu,
        "alpha_index_meshedness": alpha,
        "beta_index": beta,
        "gamma_index": gamma,
    }


@timer
def centrality_metrics(G, k_sample=500, n_workers=1):
    N = G.number_of_nodes()
    k = min(k_sample, N) if N > k_sample else None

    logging.info("Large network, use networkit")
    nk.setNumberOfThreads(n_workers)
    logging.info(f"{n_workers} threads")
    nkG, _ = _nk_graph(G, weight="length_m")
    logging.info(f"Sampling betweenness k={k_sample} times")
    eb = nk.centrality.EstimateBetweenness(nkG, k, normalized=True, parallel=True)
    eb.run()
    logging.info("Sampling done")
    bc_vals = np.array(eb.scores())

    return {
        "betweenness_mean": bc_vals.mean() if N else np.nan,
        "betweenness_max": bc_vals.max() if N else np.nan,
        "betweenness_centralization": (bc_vals.max() - bc_vals.mean()) if N else np.nan,
        "betweenness_gini": gini(bc_vals) if N else np.nan,
    }


@timer
def clustering_metrics(G):
    return {
        "avg_clustering": nx.average_clustering(G),
        "transitivity": nx.transitivity(G),
    }


@timer
def community_metrics(G):
    if G.number_of_edges() == 0:
        return {"modularity": np.nan, "n_communities": np.nan}
    communities = nx.community.louvain_communities(G, seed=42)
    Q = nx.community.modularity(G, communities)
    return {"modularity": Q, "n_communities": len(communities)}


@timer
def spatial_metrics(G, sample_pairs=2000, seed=42, n_workers=1):
    """Requires 'pos' (lon, lat) node attrs and 'length_m' edge attrs."""
    pos = nx.get_node_attributes(G, "pos")
    total_len_m = sum(d.get("length_m", 0) or 0 for _, _, d in G.edges(data=True))

    # Circuity: graph distance / straight-line distance, sampled over pairs
    # within the giant component
    largest_cc = max(nx.connected_components(G), key=len)
    Gc = G.subgraph(largest_cc).copy()
    cc_nodes = list(Gc.nodes())
    rng = np.random.default_rng(seed)

    if len(cc_nodes) > 1:
        n_pairs = min(sample_pairs, len(cc_nodes) * (len(cc_nodes) - 1) // 2)
        pairs = set()
        tries = 0
        while len(pairs) < n_pairs and tries < n_pairs * 5:
            tries += 1
            u, v = rng.choice(cc_nodes, size=2, replace=False)
            pairs.add((u, v))
        pairs = list(pairs)

        avg_circuity = _nk_circuity(Gc, pos, pairs, n_workers)
    else:
        avg_circuity = np.nan

    return {
        "total_length_km": total_len_m / 1000,
        "avg_circuity": avg_circuity,
    }


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return np.nan
    cum = np.cumsum(x)
    return (n + 1 - 2 * (cum.sum() / cum[-1])) / n


def haversine_m(lat1, lon1, lat2, lon2):
    R = 6_371_000
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlambda / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


@timer
def compute_all_metrics(
        edges, nodes, from_col="from_id", to_col="to_id", geometry_col="geometry",
        id_col="id", n_workers=1,
    ) -> pd.Series:

    G = _nx_graph(edges, nodes, from_col, to_col, geometry_col, id_col)

    metrics = {"n_isolated_nodes": nx.number_of_isolates(G)}
    metrics.update(basic_topology(G))
    metrics.update(connectivity_metrics(G, n_workers=n_workers))
    metrics.update(meshedness_metrics(G))
    metrics.update(centrality_metrics(G, n_workers=n_workers))
    metrics.update(clustering_metrics(G))
    metrics.update(community_metrics(G))
    metrics.update(spatial_metrics(G, n_workers=n_workers))

    return pd.Series(metrics)
