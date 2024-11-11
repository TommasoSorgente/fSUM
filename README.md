# FESA

FESA: _Field-driven Enriched Segmentation Algorithm_ is an algorithm for partitioning a planar or volumetric domain, discretized through a mesh of any kind and with a scalar field defined on it, into a given number of regions in which the field is homogeneous.

## How to get it 

Please use --recursive when cloning this repository:

```
git clone --recursive git@github.com:TommasoSorgente/FESA.git
```

## Content of the repository

 - `data` meshes and scalar fields used in the paper;
 - `external` links to external libraries;
 - `src` the source code;
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
After that, launch the `launcher.sh` Bash script in each folder to automatically reproduce the paper results.
