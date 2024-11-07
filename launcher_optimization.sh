#!/bin/bash

FESA="build/Qt_6_8_0_for_macOS-Release/FESA"
#FESA="build/Desktop_Qt_6_8_0-Release/FESA"

dim=2
mesh_path="data/sec_rot/mesh.obj"
field_path="data/sec_rot/mean.csv"
field_type=1
FGLOBAL=0
fieldG_path="data/Liguria_tri/cr_mean_global.csv"
out_path="out/"

n_regions=8
isoval_type=1
isoval_vals=(0 25 50 75 95 100)
DENOISE=1
ISOCONTOURS=0

ANALYZE=0
CLEAN=0
SMOOTH=0
#clean_thresh=0.1
n_iter=50
SIGMA=0

gui=1
verbose=0

threshold=(1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10)
src="./out/sec_rot/global_stats.txt"
dst="./out/optimization/sec_rot_global_stats.txt"

##########################################################################################

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
if [ $SIGMA = 1 ]; then
    switch_string=$switch_string"-M "
fi
if [ $gui = 1 ]; then
    switch_string=$switch_string"-U "
fi
if [ $verbose = 1 ]; then
    switch_string=$switch_string"-V "
fi

# Remove the dst file if it exists
if [ -f $dst ]; then
    rm $dst
fi
# Add the header from the src file to the dst file
head -n 1 $src > $dst

for t in ${threshold[@]}; do
    clean_thresh=$t
    values_string="-m "$mesh_path" -f "$field_path" -g "$fieldG_path" -o "$out_path" -d "$dim" -t "$field_type" -r "$n_regions" -i "$isoval_type" -n "$n_iter" -e "$clean_thresh
    #
    echo "launcher:" $values_string $isoval_string $switch_string
    echo ""
    #./${FESA} $values_string $isoval_string $switch_string
    #
    # Skip the first line and append the rest to the output file
    tail -n +2 $src >> $dst
    echo ""
done

