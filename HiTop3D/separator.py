import pathlib
from pygel3d import hmesh, graph
import trimesh          

PATH_MESH = "minW3dmgcg_SMreMF_capped.stl"
BASENAME  = pathlib.Path(PATH_MESH).with_suffix('').name + "_skel"

##############################

def load_mesh_any(path: str) -> hmesh.Manifold:
    ext = pathlib.Path(path).suffix.lower()
    if ext in (".obj", ".off", ".ply", ".x3d"):
        m = (hmesh.obj_load(path) or hmesh.off_load(path) or
             hmesh.ply_load(path) or hmesh.x3d_load(path))
    elif ext == ".stl":
        # read with trimesh, send triangles to PyGEL
        tm = trimesh.load(path, force='mesh')
        if not tm.is_watertight:
            raise RuntimeError("mesh is not watertight")
        m = hmesh.Manifold.from_triangles(tm.vertices, tm.faces)
    else:
        raise ValueError("unsupported extension.")
    if m is None or not hmesh.valid(m) or not hmesh.closed(m):
        raise RuntimeError("mesh is faulty")
    hmesh.triangulate(m)
    return m

# nodes from mesh vertices
def mesh_to_graph(m: hmesh.Manifold) -> graph.Graph:
    g = graph.Graph()
    for p in m.positions():
        g.add_node(p)
    for fid in m.faces():
        vs = list(m.circulate_face(fid, 'v'))
        g.connect_nodes(vs[0], vs[1])
        g.connect_nodes(vs[1], vs[2])
        g.connect_nodes(vs[2], vs[0])
    return g


def write_obj(fn: str, G: graph.Graph):
    P = G.positions()
    with open(fn, "w") as fp:
        for x, y, z in P:
            fp.write(f"v {x} {y} {z}\n")
        for n in G.nodes():
            for nb in G.neighbors(n, 'n'):
                if n < nb:
                    fp.write(f"l {n+1} {nb+1}\n")

##############################

def main():
    mesh = load_mesh_any(PATH_MESH)
    g_src = mesh_to_graph(mesh)
    print(f"mesh_to_graph nodes: {len(g_src.nodes())}, edges ≈ {sum(len(g_src.neighbors(n,'n')) for n in g_src.nodes())//2}")

    # local separator
    sk_ls = graph.LS_skeleton(g_src)
    write_obj(BASENAME + "_ls.obj", sk_ls)

    # multi scale local separator
    # grow_thresh = 64  # tweak?
    # sk_msls = graph.MSLS_skeleton(g_src, grow_thresh=grow_thresh)
    # write_obj(BASENAME + "_msls.obj", sk_msls)

    # front separator (bad on arm)
    # ----------------------------------------------------------------
    # try:
    #     P = [g_src.positions()[i] for i in g_src.nodes()]  # python list of xyz
    #     colors = np.ascontiguousarray(P, dtype=np.float64)  # (N,3) float64
    #     sk_fr = graph.front_skeleton(g_src, colors, intervals=120)
    #     write_obj(BASENAME + "_front.obj", sk_fr)
    # except Exception as e:
    #     print(f"front separator failed {e.__class__.__name__}: {e}")

    print("created", BASENAME + "_ls.obj")
    #print("created", BASENAME + "_msls.obj")

if __name__ == "__main__":
    main()