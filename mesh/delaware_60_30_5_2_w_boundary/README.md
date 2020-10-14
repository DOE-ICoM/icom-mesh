## `Delaware-60-30-5-2-w-boundary`

A "boundary-constrained" Delaware-centric configuration, adding various watershed and river polyline constraints to the `Delaware-60-30-5-2` configuration.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/delaware_60_30_5_2_w_boundary"
    conda deactivate

will build the mesh in the `../delaware_60_30_5_2_w_boundary` directory. Output is written to `../delaware_60_30_5_2_w_boundary/out`. Requires the `../delaware_60_30_5_2` configuration to first be run, to be used as initial conditions. 
