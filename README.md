## `ICoM-MESH`

The "unified" multi-model meshing workflow.

### `Quickstart`

`ICoM-MESH` requires a number of dependencies be installed, including: <a href="https://github.com/dengwirda/jigsaw-python">`JIGSAW`</a> (the underlying mesh generator), <a href="https://github.com/MPAS-Dev/MPAS-Tools">`MPAS-Tools`</a> (to manipulate `MPAS`-format mesh data-structures), as well as various `Python` packages and utilities and data files. Once installed, meshes can be built via:

    clone/download + unpack this repository
    python extract.py
    python meshify.py --mesh-path="PATH-TO-USER-MESH-DIRECTORY"
    
Dependencies can either be built directly from src, or installed using the <a href="https://anaconda.org/conda-forge/">`conda`</a> package manager (highly recommended!):

    conda create -n icom_mesh_plus python=3.7 netCDF4 scipy geojson githubrelease mpas_tools jigsaw jigsawpy

Data files can be downloaded from the <a href="https://github.com/DOE-ICoM/icom-mesh-data">`ICoM-MESH-DATA`</a> repository using the `extract.py` utility. Assets are unpacked into the local `../data` directory.

Each time you want to use `ICoM-MESH` first activate the environment using: `conda activate icom_mesh_plus`

Once activated, the various command-line utilities and `Python` packages will be available to the `ICoM-MESH` workflows contained in this repository. For example:

    conda activate icom_mesh_plus
    python meshify.py --mesh-path="mesh/vanilla_100"
    conda deactivate

These commands will build the "vanilla-100" mesh configuration defined in the `../mesh/vanilla_100` directory.

### `Mesh Build`

`ICoM-MESH` defines various mesh configurations, with each `../mesh/USER-MESH-DIRECTORY` containing the user-defined script `compose.py` used to design a particular grid. Running:

    python meshify.py --mesh-path="PATH-TO-USER-MESH-DIRECTORY"
    
initiates the mesh build process, roughly consisting of the following steps:

    The local ../USER-MESH-DIRECTORY/compose.py script is called to define: the geometry, 
    initial conditions, mesh spacing pattern, etc.
    JIGSAW is called to build the triangular mesh.
    MPAS mesh-tools is used to assemble MPAS-format data / file structures.
    Output is exported for the various ESM, FVCOM and ATS model components.

Output is written to `../mesh/USER-MESH-DIRECTORY/out/` with `../mesh/USER-MESH-DIRECTORY/tmp/` used for scratch storage. Various `*.vtk` output is also generated, allowing meshes to be visualised with, for example, <a href=https://www.paraview.org/>Paraview</a>.

### `Mesh Setup`

Mesh configurations are defined using a local `compose.py` script that should contain the following functions:

    compose.setgeom: defines the geometry for the mesh: ellipsoid shape, boundaries, 
                     coastlines, stream networks, etc.
    compose.setspac: defines the spacing pattern for the mesh: the variation in 
                     target edge length throughout the domain.
    compose.setinit: defines the initial conditions: (optional) an initialising mesh 
                     object on which subsequent refinement / optimisation operations 
                     are based.
    compose.setopts: defines any runtime meshing options.

This information is used by `JIGSAW` / `MPAS mesh-tools` to build a specific mesh configuration.

To define a new configuration, create a new `../mesh/USER-MESH-DIRECTORY` containing a `compose.py` script and `../tmp` and `../out` subdirectories. Customise the routines in `compose.py` to define your mesh! 



