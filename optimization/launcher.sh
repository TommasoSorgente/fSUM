#!/bin/bash

threshold=(0.001 0.5)
#threshold=(0.001 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5)
src="out/misclassification/global_misclassification.txt"
dst="out/optimization.txt"

##########################################################################################

misclassification="../misclassification/build/misclassification"
#misclassification="../misclassification/build/Qt_6_8_0_for_macOS-Debug/misclassification"
#misclassification="../misclassification/build/Desktop_Qt_6_8_0-Debug/misclassification"

data1_path="../data/Liguria_tri/11_SASSELLO/"
data2_path="../optimization/out/segmentation/11_SASSELLO/"

mesh_path=$data1_path"mesh.obj"
field_m_sigma=$data1_path"field_m_sigma.csv"
field_p_sigma=$data1_path"field_p_sigma.csv"
cells_data_path=$data2_path"domain/domain_cells_data.csv"
isovalues_path=$data2_path"isovalues.txt"
output_path="../optimization/out/misclassification"
misclass_gui=0

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
clean_thresh=""
n_iter=50

out_path="../optimization/out/segmentation/"
out_level=1
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

# Remove the dst file if it exists
if [ -f $dst ]; then
    rm $dst
fi
# Add the header from the src file to the dst file
echo -n "# clean_thresh," >> $dst
head -n 1 $src | tr '#' ' ' >> $dst

for t in ${threshold[@]}; do
    clean_thresh=$t
    values_string="-m "$mesh_path" -f "$field_path" -g "$fieldG_path" -d "$dim" -t "$field_type" -r "$n_regions" -i "$isoval_type" -n "$n_iter" -e "$clean_thresh" -o "$out_path" -l "$out_level
    #
    echo ""
    echo "launcher segmentation:" $segmentation $values_string $isoval_string $switch_string
    echo ""
    ./$segmentation $values_string $isoval_string $switch_string
    #
    echo ""
    echo "launcher misclassification:" $misclassification $mesh_path $cells_data_path $isovalues_path $output_path $misclass_gui $field_m_sigma $field_p_sigma
    echo ""
    ./$misclassification $mesh_path $cells_data_path $isovalues_path $output_path $misclass_gui $field_m_sigma $field_p_sigma

    # Skip the first line and append the rest to the output file
    echo -n "$t, " >> $dst
    tail -n +2 $src | tr '\n' ' ' >> $dst
    echo >> $dst
    echo ""
done

