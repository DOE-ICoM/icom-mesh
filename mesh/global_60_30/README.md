## `Global-60-30`

The "eddy-closure" global configuration: smooth 60-to-30km variation with latitude.

    conda activate icom_mesh_plus
    python meshify.py --mesh_path="mesh/global_60_30"
    conda deactivate

will build the mesh in the `../global_60_30` directory. Output in `../global_60_30/out`. No data dependencies.