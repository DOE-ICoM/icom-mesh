## `Vanilla-50`

A "vanilla" quasi-uniform (50km) global configuration. Nothing fancy, no variable resolution, streams, etc, just a partitioned land/ocean mesh.

    conda activate icom_mesh_plus
    python meshify.py --mesh_path="mesh/vanilla_50"
    conda deactivate

will build the mesh in the `../vanilla_50` directory. Output in `../vanilla_50/out`. No data dependencies.