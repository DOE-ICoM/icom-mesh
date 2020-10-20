## `Delaware-local-3d`

A "boundary-constrained" Delaware-centric configuration, in a local quasi-2d stereographic plane, with the geometry "lifted" onto the 3d DEM. Various watershed and river polyline constraints are imposed.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/delaware_local_3d"
    conda deactivate

will build the mesh in the `../delaware_local_3d` directory. Output is written to `../delaware_local_3d/out`. 
