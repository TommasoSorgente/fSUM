#!/bin/bash

data_folder="data/Liguria_tri/"
FESA="build/Qt_6_8_0_for_macOS-Release/FESA"
#FESA="build/Desktop_Qt_6_8_0-Release/FESA"

region=("11_FSASSELLO")
threshold=($(seq 0.0 0.01 1.0))

src="./out/"$region"/global_stats.txt"
dst="./out/optimization/"$region"_global_stats.txt"
echo $src $dst

# Remove the dst file if it exists
if [ -f $dst ]; then
    rm $dst
fi
# Add the header from the src file to the dst file
head -n 1 $src > $dst

for t in ${threshold[@]}; do
    mesh=$data_folder$region"/mesh_tri.obj"
    echo $FESA $mesh $t
    #
    ./$FESA $mesh $t
    #
    # Skip the first line and append the rest to the output file
    tail -n +2 $src >> $dst
    echo ""
done
