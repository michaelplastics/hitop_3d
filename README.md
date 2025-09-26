# HiTop3D Pipeline

A pipeline to integrate user-annotation of 3D voxel models into 3D topology optimization. Utilizes skeletonization, curvature/connectivity analysis, and distance mapping to create distinct and intuitively selectable segments. The workflow starts in MATLAB with `runTOpipeline(...)` and organizes outputs into consistent subfolders for later analysis.

## Installation

Clone the repo:

```
git clone https://github.com/yourusername/hitop-3d.git
cd hitop-3d
```

## Usage

To run from MATLAB:

```matlab
minW3dmgcg_SMreMF(nelx,nely,nelz,penal,rmin,ft,nl,cgtol,cgmax,compconst,filename)
```

To run Demo:

```matlab
runTOpipeline(48,24,24,3,sqrt(3),1,4,1e0,100,3330,'minW3dmgcg_SMreMF')
```

Arguments correspond to:


* ```nelx, nely, nelz``` → number of elements in x, y, z.
* ```volfrac``` → volume fraction (constraint).
* ```penal``` → penalization power.
* ```rmin``` → filter radius.
* ```maxiter``` → maximum iterations.
* ```e0``` → Young’s modulus.
* ```nu``` → Poisson’s ratio.
* ```filename``` → base name for outputs.


Optional flags:

* ```'SkipExisting'``` → default true: skip steps if artifacts exist.
* ```'Force'``` → default false: if true, overwrite/redo.
* ```'Visualize'``` → default true: open viewer + write JSON.
* ```'Python'``` → pick interpreter.
* ```'AutoPip'``` → default false: auto-install missing py deps.
* ```'PipArgs'``` → extra args for pip when AutoPip=true.


## Output Structure

After completion, results are placed in two subfolders (inside the run directory):

```
<filename>/
├─ annotated_voxel_data/
│  ├─ <filename>_voxBin.mat
│  ├─ *.mat                   
│  └─ <filename>_selected_clusters.json
└─ skeletonization_data/
   ├─ <filename>.stl
   ├─ <filename>_capped.stl
   ├─ <filename>_capped_skel_ls.obj
   ├─ <filename>_capped_skel_ls_segments/ 
   └─ *.json (except *_selected_clusters.json)

```

To access user-annotated voxel data, find ```annotated_voxel_data/<base>_selected_clusters.json```.

## Script Overview

MATLAB:

* ```runTOpipeline.m``` → Orchestrates the entire pipeline; handles options, Python preflight, step-skips, and output organization.
* ```minW3dmgcg_SMreMF.m``` → 3D minimum-volume TO solver (Andreassen/Amir line), writes ```<base>_voxBin.mat```.
* ```voxBin2stl.m``` → Converts binarized voxel field (```*_voxBin.mat```) to triangulated surface STL.
* ```viewSkeletonVsVoxel.m``` → Side-by-side viewer; overlays segments vs voxels, writes ```<base>_selected_clusters.json```.

Python:

* ```quick_cap.py``` → Mesh IO + repair/capping (uses meshio, pymeshfix).
* ```separator.py``` → Level-set skeletonization from capped surface; outputs ```*_skel_ls.obj```.
* ```segmentor.py``` → Graph/branch extraction & per-segment exports (creates *_segments/).

## Requirements

MATLAB R202x or newer

Python 3.10+

numpy, scipy, matplotlib, pygel3d
