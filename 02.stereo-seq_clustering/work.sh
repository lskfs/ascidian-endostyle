#!/bin/bash

cwd=$PWD

### seurat clustering for each slice
mkdir $cwd/01.each_slice
for gem_file in *.cellbin.gem; do
    prefix=$(basename $gem_file .cellbin.gem)
    Rscript Stereo_seurat.cellbin.R -i $gem_file -o $cwd/01.each_slice/$prefix -s $prefix --topn 30 --resolution 0.5 --minCount 25 --maxCount 7500
done

### integrate multiple sections stereo-seq data
mkdir $cwd/02.integrate
for rds_file in $cwd/01.each_slice/*/*.rds; do 
    prefix=$(basename $rds_file .rds)
    ln -s $rds_file $cwd/02.integrate/${prefix}.rds
done
Rscript seurat_integrate.R -i $cwd/02.integrate -o $cwd/02.integrate --resolution 0.8

