
### cell segmentation using cellprofiler
for tif_file in *.ssDNA.tif; do
    prefix=$(basename $tif_file .ssDNA.tif)
    python3 cellsegment.py --image $tif_file --cpi styela_clava.cppipe --prefix $prefix
done

### manully detect the tissue region as well as dense region using Fiji
# '${prefix}.en.tif' dense region masked tif
# '${prefix}.tissue.tif' tissue region masked tif

### also you need affine matrix from registration step
### you should have finished this step using TrakEM2
### PLS following https://www.ini.uzh.ch/~acardona/trakem2_manual.html#registration
# '${prefix}.registration.affineR.txt' affine matrix

### mix segmente with cell and bin
for tif_file in *.ssDNA.tif; do
    prefix=$(basename $tif_file .ssDNA.tif)
    python3 mix_segmentation.py -c ./${prefix}_mask.txt -b ./${prefix}.en.tif -t ./${prefix}.tissue.tif -i ./$tif_file -g ./${prefix}.gem -a ./${prefix}.registration.affineR.txt -p ${prefix}_mix
    python3 centroid.py ${prefix}_mix.mask.txt ${prefix}_mix
    python3 gem_mask.py -g ./${prefix}.gem -m ${prefix}_mix.mask.txt -p ${prefix}_mix
done

### after finishing these steps, you should obtain cell segmentated gem file


