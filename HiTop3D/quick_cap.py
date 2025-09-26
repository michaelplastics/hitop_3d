#!/usr/bin/env python3
"""
quick_cap.py  ––  cap open borders (“box cuts”) in an STL using PyMeshFix

usage
    python quick_cap.py in_mesh.stl  [out_mesh.stl]

requirements
    pip install pymeshfix meshio
"""
import sys, pathlib, meshio, pymeshfix as mf

# ---------------------------------------------------------------- CLI ----
if len(sys.argv) < 2:
    sys.exit("usage: python quick_cap.py in_mesh.stl [out_mesh.stl]")

infile  = pathlib.Path(sys.argv[1])
outfile = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 \
         else infile.with_name(infile.stem + "_capped.stl")

# ---------------------------------------------------------------- load ---
mesh  = meshio.read(infile)
verts = mesh.points
faces = mesh.cells_dict["triangle"]

# ---------------------------------------------------------------- fix ----
fixer = mf.MeshFix(verts, faces)

if hasattr(fixer, "fill_small_boundaries"):          # newer PyMeshFix ≥0.16
    fixer.fill_small_boundaries()                    # cap holes only
else:                                                # older builds
    # try full repair with widest keyword set, then fall back as needed
    repaired = False
    for kwargs in (
        dict(verbose=False, joincomp=True,
             remove_smallest_components=False, remove_degenerate=False),
        dict(verbose=False, joincomp=True,
             remove_smallest_components=False),
        dict(verbose=False, joincomp=True),           # minimal
        dict()                                       # as last resort
    ):
        try:
            fixer.repair(**kwargs)
            repaired = True
            break
        except TypeError:
            continue
    if not repaired:
        sys.exit("PyMeshFix.repair() failed – version mismatch?")

v_fix, f_fix = fixer.v, fixer.f

# ---------------------------------------------------------------- save ---
meshio.write(outfile, meshio.Mesh(v_fix, [("triangle", f_fix)]))
print(f"Capped STL written to  {outfile}")
