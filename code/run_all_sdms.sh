#!/bin/bash

for i in {1..10}
do
    #create output dir
    mkdir -p "severn_sdm/run_$i"
    mkdir -p "severn_sdm_web/run_$i"
    # run basic SDM
    python3 run_sdms.py -n 1 ../data/species_layers.csv severn_barrage.csv "severn_sdm/run_$i"
    # run food web SDM
    python3 run_sdms_web.py ../data/species_layers_web.csv severn_barrage_linked.csv "severn_sdm_web/run_$i" "severn_sdm/run_$i"
    # clear up
done
