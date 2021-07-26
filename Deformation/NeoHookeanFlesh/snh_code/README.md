This code implements the paper [Stable Neo-Hookean Flesh Simulation](https://dl.acm.org/citation.cfm?id=3180491)
on hexahedral and tetrahedral meshes. This implementation includes a core library
named ```cubesim``` that provides the necessary 3D mesh data structures to
implement a FEM-based Lagrangian simulation, code to implement the Stable
Neo-Hookean material model, and Newton solvers to minimize the strain energy
over these meshes.

Authors: [Breannan Smith](https://breannansmith.com), [Fernando de Goes](http://fernandodegoes.org), [Theodore Kim](http://www.tkim.graphics)


Dependencies
------------

* [Eigen](http://eigen.tuxfamily.org/): Required. A linear algebra library.
* [OpenMP](https://www.openmp.org): Optional. Multi-threading support.


Build
-----

This code requires the [CMake](https://cmake.org) build system. From a clean
download, to build the code:

1. Download and extract Eigen with the provided script:

        ./download_eigen.sh

2. Create and enter a build directory:

        mkdir build
        cd build

3. Configure the build:

        cmake ..
        # With OpenMP
        cmake -DUSE_OPENMP=ON ..

4. Run the build:

        # Add -j#threads if desired
        make


Run
---

After building the code, executables will be available in the ```hexcli``` and
```tetcli``` subdirectories.

To run a simulation and save obj meshes to a directory:

        # After following the build instructions
        cd tetcli
        mkdir output
        # Parameters: resolution, model, mu, lambda, directory
        ./tetcli 10 stable_neo_hookean 1.0 10.0 output
