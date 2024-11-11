#!/bin/bash

misclassification="build/misclassification"
#misclassification="build/Qt_6_8_0_for_macOS-Debug/misclassification"
#misclassification="build/Desktop_Qt_6_8_0-Debug/misclassification"

data1_path="../data/Liguria_tri/11_SASSELLO/"
data2_path="out/segmentation/11_SASSELLO/"

mesh_path=$data1_path"mesh.obj"
field_m_sigma=$data1_path"field_m_sigma.csv"
field_p_sigma=$data1_path"field_p_sigma.csv"
cells_data_path=$data2_path"domain/domain_cells_data.csv"
isovalues_path=$data2_path"isovalues.txt"
output_path="out/11_SASSELLO"
misclass_gui=1

##########################################################################################

segmentation="../segmentation/build/segmentation"
#segmentation="../segmentation/build/Qt_6_8_0_for_macOS-Debug/segmentation"
#segmentation="../segmentation/build/Desktop_Qt_6_8_0-Debug/segmentation"

dim=2
mesh_path="../data/Liguria_tri/11_SASSELLO/mesh.obj"
field_path="../data/Liguria_tri/11_SASSELLO/field.csv"
field_type=1
FGLOBAL=0
fieldG_path="../data/Liguria_tri/field_global.csv"

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

out_path="../misclassification/out/segmentation/"
out_level=1
GUI=0
VERBOSE=0

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

echo ""
echo "launcher segmentation:" $segmentation $values_string $isoval_string $switch_string
echo ""

./$segmentation $values_string $isoval_string $switch_string

##########################################################################################

echo ""
echo "launcher misclassification:" $misclassification $mesh_path $cells_data_path $isovalues_path $output_path $misclass_gui $field_m_sigma $field_p_sigma
echo ""
./$misclassification $mesh_path $cells_data_path $isovalues_path $output_path $misclass_gui $field_m_sigma $field_p_sigma

