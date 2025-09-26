
INPUT_OBJ            = "minW3dmgcg_SMreMF_capped_skel_ls.obj" 
ANGLE_THRESHOLD_DEG  = 120.0              

import json, math, pathlib

##############################
def load_obj(path):
    verts, edges = [], []
    with open(path) as fp:
        for ln in fp:
            if ln.startswith("v "):
                _, x, y, z = ln.split()
                verts.append((float(x), float(y), float(z)))
            elif ln.startswith("l "):
                ids = list(map(int, ln.split()[1:]))
                for a, b in zip(ids, ids[1:]):         # support multi-point line
                    edges.append((a-1, b-1))          # OBJ → 0-based
    return verts, edges

##############################
def vec(a, b):  
    return (b[0]-a[0], b[1]-a[1], b[2]-a[2])

def is_sharp(prev, mid, nxt, cos_thresh):
    u, v = vec(mid, prev), vec(mid, nxt)
    nu = math.sqrt(sum(t*t for t in u)); nv = math.sqrt(sum(t*t for t in v))
    if nu*nv == 0:       
        return False
    cos = (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]) / (nu*nv)
    return cos > cos_thresh

def build_adj(nv, edges):
    adj = [set() for _ in range(nv)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    return adj

##############################
def find_nodes(V, adj, cos_thresh):
    nodes = set()
    for v, nbrs in enumerate(adj):
        deg = len(nbrs)
        if deg != 2: #ends/junctions
            nodes.add(v)
        else:                      
            a, b = tuple(nbrs)
            if is_sharp(V[a], V[v], V[b], cos_thresh):
                nodes.add(v)
    return nodes

##############################
def extract_branches(V, E, angle_thresh_deg):
    cos_thresh = math.cos(math.radians(angle_thresh_deg))
    adj   = build_adj(len(V), E)
    nodes = find_nodes(V, adj, cos_thresh)

    visited = set()                 # undirected edges
    branches = []

    def mark(a, b):
        visited.add(tuple(sorted((a, b))))

    # walk paths from nodes
    for start in nodes:
        for nb in adj[start]:
            if tuple(sorted((start, nb))) in visited:
                continue
            branch = [start]
            prev, cur = start, nb
            mark(prev, cur)
            while True:
                branch.append(cur)
                if cur in nodes and cur != start:
                    break
                nxt = [k for k in adj[cur] if k != prev and
                       tuple(sorted((cur, k))) not in visited]
                if not nxt:                     # reached tip
                    break
                nxt = nxt[0]
                mark(cur, nxt)
                prev, cur = cur, nxt
            branches.append(branch)

    # handle loops with no nodes
    for a, b in E:
        key = tuple(sorted((a, b)))
        if key in visited:
            continue
        loop = [a]
        prev, cur = a, b
        mark(prev, cur)
        while True:
            nxt = [k for k in adj[cur] if k != prev and
                   tuple(sorted((cur, k))) not in visited]
            if not nxt:
                break
            nxt = nxt[0]
            mark(cur, nxt)
            prev, cur = cur, nxt
            loop.append(prev)
            if cur == loop[0]:
                break
        branches.append(loop)

    return branches

##############################
def save_branch(folder, idx, vlist, verts):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"seg{idx}.obj"

    remap, uniq = {}, []
    for vid in vlist:
        if vid not in remap:
            remap[vid] = len(uniq) + 1
            uniq.append(vid)

    with open(path, "w") as fp:
        for gid in uniq:
            x, y, z = verts[gid]
            fp.write(f"v {x} {y} {z}\n")
        for a, b in zip(vlist, vlist[1:]):
            fp.write(f"l {remap[a]} {remap[b]}\n")

##############################
def main():
    obj_path = pathlib.Path(INPUT_OBJ).expanduser()
    if not obj_path.exists():
        raise FileNotFoundError(obj_path)

    V, E = load_obj(obj_path)
    branches = extract_branches(V, E, ANGLE_THRESHOLD_DEG)

    json_path = obj_path.with_suffix(".branches.json")
    with open(json_path, "w") as fp:
        json.dump({"angle_threshold": ANGLE_THRESHOLD_DEG,
                   "branches": [{"verts": br} for br in branches]}, fp, indent=2)

    seg_dir = obj_path.with_suffix("").with_name(f"{obj_path.stem}_segments")
    for i, br in enumerate(branches):
        save_branch(seg_dir, i, br, V)

    print(f"✓ {len(branches)} branches → folder '{seg_dir}/' and {json_path}")

if __name__ == "__main__":
    main()