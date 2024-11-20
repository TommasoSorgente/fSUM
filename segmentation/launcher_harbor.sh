#!/bin/bash

segmentation="build/segmentation"
#segmentation="../segmentation/build/Qt_6_8_0_for_macOS-Debug/segmentation"
#segmentation="../segmentation/build/Desktop_Qt_6_8_0-Debug/segmentation"

dim=3
mesh_path="../data/Genova_Harbor/mesh.vtk"
field_path="../data/Genova_Harbor/field.csv"
field_type=1
FGLOBAL=0
fieldG_path="."

n_regions=5
isoval_type=2
isoval_vals=(0 25 50 75 95 100)
DENOISE=1
ISOCONTOURS=0

ANALYZE=1
CLEAN=1
SMOOTH=1
clean_thresh=0.1
n_iter=50

out_path="out/"
out_level=1
GUI=1
VERBOSE=0

## LEGEND
##
## segmentation: path to the executable file
##
## INPUT
## dim:          dimension of the problem (2 | 3)
## mesh_path:    path of the input mesh
## field_path:   path of the input field
## field_type:   input field defined on cells or vertices (1 | 2)
## FGLOBAL:      use of the global field (0 | 1)
## fieldG_path:  path of the input global field
##
## SEGMENTATION
## n_regions:    number of regions to be computed in the domain (in [1,inf])
## isoval_type:  (1) equispaced, (2) percentiles, (3) assigned
## isoval_vals:  isovalues (percentiles or explicit)
## DENOISE:      use the denoised field for isoregions (0 | 1)
## ISOCONTOURS:  compute isocontours/isosurfaces (0 | 1)
##
## POSTPROCESSING
## ANALYZE:       analyze regions (0 | 1)
## CLEAN:         cleaning of small regions (0 | 1)
## SMOOTH:        smoothing of the boundaries (0 | 1)
## clean_thresh:  percentual min size of the regions (in [0,100])
## n_iter:        max number of iterations (in [0,inf])
##
## OUTPUT
## out_path:     path do a directory where to save all outputs
## out_level:	 output level: (0) domain, (1) regions, (2) subregions
## GUI:          launch graphical interface (0 | 1)
## VERBOSE:      print debug information (0 | 1)

##########################################################################################

values_string="-m "$mesh_path" -f "$field_path" -g "$fieldG_path" -d "$dim" -t "$field_type" -r "$n_regions" -i "$isoval_type" -n "$n_iter" -e "$clean_thresh" -o "$out_path" -l "$out_level

isoval_string=""
for v in ${isoval_vals[@]}; do
  isoval_string=$isoval_string" -v "$v""
done

switch_string=""
if [ $FGLOBAL = 1 ]; then
    switch_string=$switch_string"-G "
fi
if [ $DENOISE = 1 ]; then
    switch_string=$switch_string"-D "
fi
if [ $ISOCONTOURS = 1 ]; then
    switch_string=$switch_string"-I "
fi
if [ $ANALYZE = 1 ]; then
    switch_string=$switch_string"-A "
fi
if [ $CLEAN = 1 ]; then
    switch_string=$switch_string"-C "
fi
if [ $SMOOTH = 1 ]; then
    switch_string=$switch_string"-S "
fi
if [ $GUI = 1 ]; then
    switch_string=$switch_string"-U "
fi
if [ $VERBOSE = 1 ]; then
    switch_string=$switch_string"-V "
fi

echo "launcher segmentation:" $segmentation $values_string $isoval_string $switch_string
echo ""

./$segmentation $values_string $isoval_string $switch_string
