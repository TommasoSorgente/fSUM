#!/bin/bash

FESA="build/Qt_6_8_0_for_macOS-Release/FESA"
#FESA="build/Desktop_Qt_6_8_0-Release/FESA"

dim=2
#mesh_path="data/sec_rot/mesh.obj"
#field_path="data/sec_rot/mean.csv"
field_type=1
FGLOBAL=1
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
clean_thresh=0.1
n_iter=50
SIGMA=0

gui=1
verbose=0

data_folder="data/Liguria_tri/"
regions=("1_F5TERRE" 
	 "2_FMAGRA" 
	 "3_FPETRONIO" 
	 "4_FPORTOFINO" 
	 "5_FENTELLA"
	 "6_FAVETO" 
	 "7_FPADANO" 
	 "8_FBISAGNO" 
	 "9_FPOLCEVERA" 
	 "10_FARENZANO" 
	 "11_FSASSELLO" 
	 "12_FBORMIDE" 
	 "13_FSAVONESE"
	 "14_FIMPERIESE"
	 )

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

for r in ${regions[@]}; do
    mesh_path=${data_folder}${r}"/mesh_tri.obj"
    field_path=${data_folder}${r}"/cr_mean.csv"
    values_string="-m "$mesh_path" -f "$field_path" -g "$fieldG_path" -o "$out_path" -d "$dim" -t "$field_type" -r "$n_regions" -i "$isoval_type" -n "$n_iter" -e "$clean_thresh
    #
    echo "launcher:" $values_string $isoval_string $switch_string
    echo ""
    ./${FESA} $values_string $isoval_string $switch_string
done


