## `Delaware-60-30-5-2`

A simple variable-resolution Delaware-centric configuration, combining the `60-to-30km` global ocean spacing, `45km` global land spacing, `30km` North Atlantic ocean spacing, `5km` spacing in Delaware-adjacent watersheds, and `2km` spacing at the coastal fringe.

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/delaware_60_30_5_2"
    conda deactivate

will build the mesh in the `../delaware_60_30_5_2` directory. Output in `../delaware_60_30_5_2/out`. Requires various elevation, coastline, watershed and ocean-halo datasets.
