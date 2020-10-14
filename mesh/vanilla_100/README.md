## `Vanilla-100`

A "vanilla" quasi-uniform (100km) global configuration. Nothing fancy, no variable resolution, streams, etc, just a partitioned land/ocean mesh.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/vanilla_100"
    conda deactivate

will build the mesh in the `../vanilla_100` directory. Output is written to `../vanilla_100/out`. No data dependencies.
