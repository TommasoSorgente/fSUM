#!/bin/bash

data_folder="data/Liguria_tri/"
FESA="build/Qt_6_8_0_for_macOS-Release/FESA"
#FESA="build/Desktop_Qt_6_8_0-Release/FESA"

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

mesh_name="mesh_tri.obj"
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
