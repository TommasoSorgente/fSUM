# FESA

FESA: _Field-driven Enriched Segmentation Algorithm_ is an algorithm for partitioning a planar or volumetric domain, discretized through a mesh of any kind and with a scalar field defined on it, into a given number of regions in which the field is homogeneous. The resulting segmentation is called “enriched” because for each segment we are able to produce information and statistics that can be useful in applicative scenarios.

## How to get it 

Please use --recursive when cloning this repository:

```
git clone --recursive git@github.com:TommasoSorgente/FESA.git
```

FESA can be built through CMake by running the following instructions.
The executable will be available in `${REPO_ROOT}/build`, being `${REPO_ROOT}` the folder where this `README.md` lies. 

```
cd ${REPO_ROOT}
mkdir build
cd build
cmake ..
make
```

## Content of the repository

 - `data` folder containing the data used in the paper.
 - `external` folder containing links to external libraries.
 - `merge_regions` post-processing routine for merging regions into larger domains (used for the _Liguria_ example).
 - `src` folder containing the source code.
 - `CMakeLists.txt` CMake file to compile the code.
 - `main.cpp` sample program that reads input parameters, partitions the mesh and optionally displays it in a graphical interface.
 - `parameters.run` file containing all the parameters involved in the algorithm.
