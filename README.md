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

 - `data` folder containing the data used in the paper;
 - `external` folder containing links to external libraries;
 - `src` folder containing the source code;
 - `segmentation` main routine for segmenting a scalar field associated to a mesh;
 - `merge_subdomains` post-processing routine for merging multiple segmentations into a single one;
 - `misclassification` post-processing routine for computing the misclassification of a segmentation;
 - `optimization` post-processing script for computing multiple segmentations with different _epsilon_ parameters, and choosing the optimal one.

The three routines are independent, and each of them needs to be compiled through CMake by running:
```
cd ${ROUTINE_NAME}
mkdir build
cd build
cmake ..
make
```
After that, launch the `launcher.sh` Bash script to automatically reproduce the paper results.
