# FESA
FESA: Field-Driven Enriched Segmentation Algorithm

FESA is an algorithm for partitioning a planar or volumetric domain, discretized through a mesh of any kind and with a scalar field defined on it, into a given number of regions in which the field is homogeneous. The resulting segmentation is called “enriched” because for each segment we are able to produce information and statistics that can be useful in applicative scenarios.

**How to get it:** 
 - clone with '--recursive' to get the cinolib library as well
 - install the library 'shapelib' and specify its path in 'CMakeLists.txt'
 - compile with cmake

**Content of the repository:**
 - _CMakeLists.txt_: CMake file to compile the code.
 - _parameters.run_: file containing all the parameters involved in the algorithm.
 - _main.cpp_: sample program that reads input parameters, partitions the mesh and optionally displays it in a graphical interface.
 - _src_: folder containing the source code.
 - _data_: folder containing the data used in the paper.

