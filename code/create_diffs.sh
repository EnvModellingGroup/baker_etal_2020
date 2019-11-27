#!/bin/bash

dir=$1

# Variable names to pull out of the PVTUs
declare -a species=("bib" 
                      "brown_shrimp" 
                      "catworm"
                      "dover_sole"
                      "lugworm"
                      "mud_shrimp"
                      "mud_snail"
                      "periwinkle"
                      "pink_clam"
                      "poor_cod"
                      "sand_digger_shrimp"
                      "sand_goby"
                      "velvet_swimming_crab"
                      "whelk"
                     )

for sp in "${species[@]}"
do
    echo "Dealing with $sp"
    gdal_calc.py -A "${dir}/${sp}_pd.tif" -B "${dir}/${sp}_barrage.tif" --calc="A-B" --outfile="${dir}/${sp}_diff.tif" --NoDataValue=0.0 --overwrite
done

