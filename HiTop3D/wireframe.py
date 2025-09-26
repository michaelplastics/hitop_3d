#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wireframe.py — PyGEL3D MSLS skeletonizer (path-safe, stable kNN)

Usage (any of these work):
  • MATLAB wrapper sets module globals before main():
        m.PATH_INPUT_OBJ, m.PATH_SKELETON_OBJ, m.PATH_CYLINDER_OBJ
  • env vars:
        PATH_INPUT_OBJ=... PATH_SKELETON_OBJ=... PATH_CYLINDER_OBJ=... python wireframe.py
  • CLI:
        python wireframe.py <input_mesh> [out_skeleton_obj] [out_cylinder_obj]
"""

import os, sys, math
from pathlib import Path
import numpy as np
from scipy.spatial import KDTree
from pygel3d import hmesh, graph

# ────────────────────────── Path resolution ──────────────────────────
def _resolve(name, default=None):
    if name in globals() and globals()[name] is not None:
        return globals()[name]
    env = os.environ.get(name)
    if env:
        return env
    if name == "PATH_INPUT_OBJ" and len(sys.argv) > 1:      return sys.argv[1]
    if name == "PATH_SKELETON_OBJ" and len(sys.argv) > 2:   return sys.argv[2]
    if name == "PATH_CYLINDER_OBJ" and len(sys.argv) > 3:   return sys.argv[3]
    return default

DEFAULT_IN  = None
DEFAULT_SK  = "skeleton.obj"
DEFAULT_CYL = "cylinders.obj"

PATH_INPUT_OBJ    = _resolve("PATH_INPUT_OBJ",    DEFAULT_IN)
PATH_SKELETON_OBJ = _resolve("PATH_SKELETON_OBJ", DEFAULT_SK)
PATH_CYLINDER_OBJ = _resolve("PATH_CYLINDER_OBJ", DEFAULT_CYL)

BANNER = "─" * 59
def mark(msg): print(f"\n{BANNER}\n{msg}\n{BANNER}", flush=True)

# ─────────────────── Mesh loading & validation ───────────────────────
def load_manifold(path: str) -> hmesh.Manifold:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(p)
    ext = p.suffix.lower()
    if   ext == ".obj":  m = hmesh.obj_load(str(p))
    elif ext == ".off":  m = hmesh.off_load(str(p))
    elif ext == ".ply":  m = hmesh.ply_load(str(p))
    else:                m = hmesh.load(str(p))
    if m is None:
        raise IOError(f"✗  Could not parse '{p.name}' — try triangulating it first.")
    if not hmesh.valid(m):
        raise ValueError(f"✗  '{p.name}' is not a manifold mesh.")
    hmesh.triangulate(m)   # idempotent; ensures triangles
    return m

# ───────────────────── Point sampling (stable) ───────────────────────
def tri_area(a, b, c): return 0.5 * np.linalg.norm(np.cross(b - a, c - a))

def surf_samples(m, n, rng=None):
    rng = rng or np.random.default_rng(42)
    V   = m.positions()
    F   = [list(m.circulate_face(fid, 'v')) for fid in m.faces()]
    A   = np.array([tri_area(V[a], V[b], V[c]) for a, b, c in F], dtype=np.float64)
    acc = np.cumsum(A / A.sum())
    out = np.empty((n, 3), dtype=np.float64)
    for i in range(n):
        a, b, c = F[np.searchsorted(acc, rng.random())]
        u, v = rng.random(), rng.random()
        if u + v > 1: u, v = 1 - u, 1 - v
        out[i] = V[a] + u * (V[b] - V[a]) + v * (V[c] - V[a])
    return out

def inner_jitter(m, n, eps=0.003, rng=None):
    rng = rng or np.random.default_rng(123)
    V = np.asarray(m.positions())
    idx = rng.integers(0, len(V), size=n)
    jit = np.empty((n,3), dtype=np.float64)
    for i, ii in enumerate(idx):
        jit[i] = V[ii] - eps * m.vertex_normal(int(ii))
    return jit

# ─────────────────────── k-NN graph (robust) ─────────────────────────
def _sanitize_points(P):
    P = np.ascontiguousarray(P, dtype=np.float64)
    P = P[np.isfinite(P).all(axis=1)]
    if len(P) == 0:
        return P
    bbox = P.max(0) - P.min(0)
    tol  = (np.linalg.norm(bbox) / 1e5) or 1e-9
    cell = tol
    keys = np.floor((P - P.min(0)) / cell + 0.5).astype(np.int64)
    view = keys.view([('', keys.dtype)]*3)
    _, uniq_idx = np.unique(view, return_index=True)
    return P[np.sort(uniq_idx)]

def build_knn_graph(points, k=10):
    P = _sanitize_points(points)
    if len(P) < 4:
        raise ValueError("Too few points after sanitization.")
    k   = max(3, min(int(k), len(P) - 1))
    g   = graph.Graph()
    ids = [g.add_node(p) for p in P]
    kdt = KDTree(P)
    for i, p in enumerate(P):
        dists, nbrs = kdt.query(p, k + 1)   # includes self
        for j, d in zip(nbrs[1:], dists[1:]):
            if i == j or d <= 1e-12: 
                continue
            g.connect_nodes(ids[i], ids[j])
    return g

def graph_edges(g):
    for u in g.nodes():
        for v in g.neighbors(u, 'n'):
            if u < v:
                yield u, v

# ───────────────── Cylinders (skip tiny segments) ────────────────────
def make_cyl(p0, p1, r, sides=14):
    axis = np.asarray(p1, dtype=float) - np.asarray(p0, dtype=float)
    h = float(np.linalg.norm(axis))
    if h == 0:  return [], []
    axis /= h
    tmp  = np.array([1,0,0]) if abs(axis[0]) < .9 else np.array([0,1,0])
    u    = np.cross(axis, tmp); u /= np.linalg.norm(u)
    v    = np.cross(axis, u)
    verts = [ list(p0 + r*(math.cos(2*math.pi*i/sides)*u + math.sin(2*math.pi*i/sides)*v)) for i in range(sides) ]
    verts += [ list(p1 + r*(math.cos(2*math.pi*i/sides)*u + math.sin(2*math.pi*i/sides)*v)) for i in range(sides) ]
    faces = [ [i, (i+1)%sides, sides+(i+1)%sides, sides+i] for i in range(sides) ]
    return verts, faces

def save_cylinders(path, skel, dist, fallback_r=0.03, sides=14):
    verts, faces = [], []
    P = skel.positions()
    for u, v in graph_edges(skel):
        seg = np.asarray(P[v]) - np.asarray(P[u])
        if np.dot(seg, seg) < 1e-12:
            continue
        r = abs(dist.signed_distance([P[u]])[0]) or fallback_r
        V, F = make_cyl(P[u], P[v], r, sides)
        off  = len(verts)
        verts.extend(V)
        faces.extend([[off+idx for idx in f] for f in F])
    with open(path, "w") as fp:
        for p in verts:  fp.write(f"v {p[0]} {p[1]} {p[2]}\n")
        for f in faces:  fp.write("f " + " ".join(str(i+1) for i in f) + "\n")
    print(f"[✓] cylinders  →  {path}")

# ───────────────────────────── Pipeline ───────────────────────────────
def main():
    print(f"\nINPUT : {PATH_INPUT_OBJ}")
    print(f"SKELE : {PATH_SKELETON_OBJ}")
    print(f"CYL   : {PATH_CYLINDER_OBJ}")

    if not PATH_INPUT_OBJ:
        raise SystemExit("PATH_INPUT_OBJ not set (wrapper/env/cli).")
    p_in = Path(PATH_INPUT_OBJ).expanduser().resolve()
    if not p_in.exists():
        raise FileNotFoundError(f"Input mesh not found: {p_in}")

    mark("Loading mesh")
    mesh = load_manifold(str(p_in))

    # compute bbox/extent from vertices (Manifold has no bbox_min/max)
    V_all = np.asarray(mesh.positions())
    bb_min = V_all.min(axis=0)
    bb_max = V_all.max(axis=0)
    extent = float(np.linalg.norm(bb_max - bb_min))

    # adaptive, reproducible sampling
    rng   = np.random.default_rng(42)
    Nsurf = int(np.clip(extent * 150.0, 5000, 40000))
    Njit  = max(1000, Nsurf // 3)

    mark("Sampling points")
    pts = np.vstack([surf_samples(mesh, Nsurf, rng=rng),
                     inner_jitter(mesh, Njit, rng=rng)])

    mark("Building k-NN graph")
    G = build_knn_graph(pts, k=10)

    mark("Skeletonizing (MSLS)")
    skel = graph.MSLS_skeleton(G)

    # save line skeleton
    P = skel.positions()
    with open(PATH_SKELETON_OBJ, "w") as fp:
        for p in P: fp.write(f"v {p[0]} {p[1]} {p[2]}\n")
        for u, v in graph_edges(skel):
            fp.write(f"l {u+1} {v+1}\n")
    print(f"[✓] skeleton   →  {PATH_SKELETON_OBJ}")

    # radii + cylinders
    mark("Computing radii + writing cylinders")
    dist = hmesh.MeshDistance(mesh)
    save_cylinders(PATH_CYLINDER_OBJ, skel, dist)

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("\n‼️  Python raised:", exc)
        sys.exit(1)
