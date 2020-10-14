## `Delaware-local-2d`

A "boundary-constrained" Delaware-centric configuration, in a local 2d stereographic plance. Various watershed and river polyline constraints are imposed.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/delaware_local_2d"
    conda deactivate

will build the mesh in the `../delaware_local_2d` directory. Output is written to `../delaware_local_2d/out`. 
