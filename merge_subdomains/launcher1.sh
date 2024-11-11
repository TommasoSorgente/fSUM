#!/bin/bash

merge="build/merge_subdomains"
#merge="build/Qt_6_8_0_for_macOS-Debug/merge_subdomains"
#merge="build/Desktop_Qt_6_8_0-Debug/merge_subdomains"

data_folder="../data/Liguria_tri/"
domains=("10_ARENZANO 11_SASSELLO")
fieldG_path="../data/Liguria_tri/field_global.csv"
isovalues_path="out/segmentation/"$(echo $domains | awk '{print $1}')"/isovalues.txt"
output_path="out/10_ARENZANO+11_SASSELLO"

##########################################################################################

segmentation="../segmentation/build/segmentation"
#segmentation="../segmentation/build/Qt_6_8_0_for_macOS-Debug/segmentation"
#segmentation="../segmentation/build/Desktop_Qt_6_8_0-Debug/segmentation"

dim=2
mesh_path=""
field_path=""
field_type=1
FGLOBAL=1
fieldG_path="../data/Liguria_tri/field_global.csv"

n_regions=5
isoval_type=2
isoval_vals=(0 25 50 75 95 100)
DENOISE=1
ISOCONTOURS=0

ANALYZE=1
CLEAN=1
SMOOTH=1
clean_thresh=2
n_iter=50

out_path="../merge_subdomains/out/segmentation/"
out_level=3
gui=0
verbose=0

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
if [ $gui = 1 ]; then
    switch_string=$switch_string"-U "
fi
if [ $verbose = 1 ]; then
    switch_string=$switch_string"-V "
fi

for r in ${domains[@]}; do
    mesh_path=${data_folder}${r}"/mesh.obj"
    field_path=${data_folder}${r}"/field.csv"
    values_string="-m "$mesh_path" -f "$field_path" -g "$fieldG_path" -d "$dim" -t "$field_type" -r "$n_regions" -i "$isoval_type" -n "$n_iter" -e "$clean_thresh" -o "$out_path" -l "$out_level
    #
    echo ""
    echo "launcher segmentation:" $segmentation $values_string $isoval_string $switch_string
    echo ""
    ./$segmentation $values_string $isoval_string $switch_string
done

##########################################################################################

echo ""
echo "launcher merge:" $merge $domains $fieldG_path $isovalues_path $output_path
echo ""
./$merge "$domains" $fieldG_path $isovalues_path $output_path

