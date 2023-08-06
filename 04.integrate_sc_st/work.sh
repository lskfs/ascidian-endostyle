#!/bin/bash

cwd=$PWD

index=1
for rds_file in ../04.single_cell_clustering/*/*.rds; do
    ln -s $rds_file sc${index}.rds
    index=$(($index+1))
    echo $index
done

index=1
for rds_file in ../02.stereo-seq_clustering/*/*.rds; do
    ln -s $rds_file st${index}.rds
    index=$(($index+1))
done

Rscript integrate.sc-st.R -i ./ -o ./ -r 0.8 -s st_sc

