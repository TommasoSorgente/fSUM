#!/bin/bash

# Automatic routine for launching FESA on all Ligurian regions.
# Please specify the path to the regions folder and to the FESA executable. 

data_folder="./"
#FESA="../../build/FESA"
#FESA="../../build/Qt_6_7_2_for_macOS-Release/FESA"
FESA="../../build/Desktop_Qt_6_7_2-Debug/FESA"

regions=("cr_F5TERRE" 
	 "cr_FARENZANO" 
	 "cr_FAVETO" 
	 "cr_FBISAGNO" 
	 "cr_FBORMIDE" 
	 "cr_FENTELLA" 
	 "cr_FIMPERIESE" 
	 "cr_FMAGRA" 
	 "cr_FPADANO" 
	 "cr_FPETRONIO" 
	 "cr_FPOLCEVERA" 
	 "cr_FPORTOFINO" 
	 "cr_FSASSELLO" 
	 "cr_FSAVONESE"
	 )

mesh_name="mesh_grid.obj"
field_name="cr_mean"

for r in ${regions[@]}; do
    mesh=${data_folder}${r}"/"${mesh_name}
    field_local=${data_folder}${r}"/"${field_name}".csv"
    field_global=${data_folder}${field_name}"_global.csv"
    echo ${FESA} ${mesh} ${field_local} ${field_global}
    #
    ./${FESA} ${mesh} ${field_local} ${field_global}
    echo ""
done
