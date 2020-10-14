## `Vanilla-50`

A "vanilla" quasi-uniform (50km) global configuration. Nothing fancy, no variable resolution, streams, etc, just a partitioned land/ocean mesh.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/vanilla_50"
    conda deactivate

will build the mesh in the `../vanilla_50` directory. Output is written to `../vanilla_50/out`. No data dependencies.
