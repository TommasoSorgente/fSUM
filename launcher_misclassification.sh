#!/bin/bash

#exe="misclassification/build/Qt_6_8_0_for_macOS-Release/misclassification"
exe="misclassification/build/Desktop_Qt_6_8_0-Debug/misclassification"

mesh_path="../data/Liguria_tri/11_FSASSELLO/mesh_tri.obj"
cells_data_path="../out/11_FSASSELLO/domain/domain_cells_data.csv"
isovalues_path="../out/11_FSASSELLO/isovalues.txt"
field_m_sigma="../data/Liguria_tri/11_FSASSELLO/cr_mean_m_sigma.csv"
field_p_sigma="../data/Liguria_tri/11_FSASSELLO/cr_mean_p_sigma.csv"

##########################################################################################

echo "launcher:" $exe $mesh_path $cells_data_path $isovalues_path $field_m_sigma $field_p_sigma
echo ""
./$exe $mesh_path $cells_data_path $isovalues_path $field_m_sigma $field_p_sigma

