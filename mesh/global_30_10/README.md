## `Global-30-10`

The "Rossby-radius" global configuration: smooth 30-to-10km variation with latitude.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/global_30_10"
    conda deactivate

will build the mesh in the `../global_30_10` directory. Output is written to `../global_30_10/out`. No data dependencies.
