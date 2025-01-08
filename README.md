# fSUM

fSUM: _field Segmentation of Unstructured Meshes_ is an algorithm for partitioning a planar or volumetric domain, discretized through a mesh of any kind and with a scalar field defined on it, into a given number of regions in which the field is homogeneous.

## How to get it 

Please use --recursive when cloning this repository:

```
git clone --recursive git@github.com:TommasoSorgente/fSUM.git
```

## Content of the repository

 - `data` meshes and scalar fields used in the paper;
 - `external` links to external libraries;
 - `src` the source code;
 - `segmentation` main routine for segmenting a scalar field associated to a mesh;
 - `merge_subdomains` post-processing routine for merging multiple segmentations into a single one;
 - `misclassification` post-processing routine for computing the misclassification of a segmentation;
 - `optimization` post-processing script for computing multiple segmentations with different _epsilon_ parameters, and choosing the optimal one.

## Build the source code

To streamline and automate the build process, two scripts are provided:

 - `build.sh`: Designed for Unix-like systems, including Linux and macOS.
 - `build.ps1`: Tailored for Windows systems.

Select the appropriate script for your operating system and execute it from its directory using the command line.

## Reproduce paper results

After that, launch the `launcher.sh` Bash script in each folder to automatically reproduce the paper results.
The following scripts are available:

 - `segmentation/launcher_sassello.sh` main example from Section 2, 3;
 - `merge_subdomains/launcher_arenzano.sh` example of merging two subdomains, from Section 3.4;
 - `misclassification/launcher.sh` example of computing the misclassification from Section 5.1;
 - `optimization/launcher.sh` example of optimizing the epsilon value, from Section 5.2;
 - `merge_subdomains/launcher_liguria.sh` environmental geochemistry application from Section 6.1;
 - `segmentation/launcher_harbor.sh` marine science application from Section 6.2;
 - `segmentation/launcher_groundwater.sh` groundwater modeling application from Section 6.3;


