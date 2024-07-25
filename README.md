# FESA

FESA: Field-Driven Enriched Segmentation Algorithm

FESA is an algorithm for partitioning a planar or volumetric domain, discretized through a mesh of any kind and with a scalar field defined on it, into a given number of regions in which the field is homogeneous. The resulting segmentation is called “enriched” because for each segment we are able to produce information and statistics that can be useful in applicative scenarios.

# How to get it 

Please, use --recursive when cloning this repository:

```
git clone --recursive https://github.com/TommasoSorgente/FESA
```

In the following, please consider `${REPO_ROOT}` variable as the folder where this `README.md` lies. 
FESA can be built by running CMake, and the executable will be available in the `${REPO_ROOT}/build` folder.

```
cd ${REPO_ROOT}
mkdir build
cd build
cmake ..
make
```

# Content of the repository

 - `data` folder containing the data used in the paper.
 - `external` folder containing links to external libraries.
 - `merge_regions` post-processing routine for merging regions into larger domains (used for the _Liguria_ example).
 - `src` folder containing the source code.
 - `CMakeLists.txt` CMake file to compile the code.
 - `main.cpp` sample program that reads input parameters, partitions the mesh and optionally displays it in a graphical interface.
 - `parameters.run` file containing all the parameters involved in the algorithm.
