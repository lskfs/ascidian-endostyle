> [!IMPORTANT]
> Manual processes are involved so no script exists for some steps. But the final results were included in our deposited data at CNGBdb.

# Main Steps
- Registration: alignment of ssDNA staining image and spatial color-coding images.
- Segementation: detection of cell/bin boundary.
- Aggregation: group DNB expression data into cell/bin level.

# How to run
> [!NOTE]
> `${prefix}` defined the file name prefix of each sample section

### Registration

Manully detect the tissue region as well as dense region using Fiji, from which two mask image could be obtained:

- `${prefix}.en.tif` dense region masked tif 
- `${prefix}.tissue.tif` tissue region masked tif

Also, you need affine matrix from registration step which
you should have finished using TrakEM2. Please following the [document](https://www.ini.uzh.ch/~acardona/trakem2_manual.html#registration). After this step you should obtain one matrix file:

- `${prefix}.registration.affineR.txt` affine matrix

### Segementation
We first detected cell objects using [CellProfiler](https://cellprofiler.org/) based on the ssDNA nucleotide images.

The input file `styela_clava.cppipe` came from a pre-tested configure using local CellProfiler. 

```shell
$ python3 cellsegment.py --image ${prefix}.ssDNA.tif --cpi styela_clava.cppipe --prefix $prefix
```
The main file generated in this step is the follow one: 
- `${prefix}_mask.txt` cell mask matrix file

### Aggregation

In this step, the cell mask matrix is first correctted by combine the mask of dense region, after which the dense region is segmented using square bins while other tissue region with cells.
```shell
$ python3 mix_segmentation.py \
    -c ./${prefix}_mask.txt \
    -b ./${prefix}.en.tif \
    -t ./${prefix}.tissue.tif \
    -i ./${prefix}.ssDNA.tif \
    -g ./${prefix}.gem \
    -a ./${prefix}.registration.affineR.txt \
    -p ${prefix}_mix
```

Here the centroid coordinates of each cell/bin were calculated.
```shell
$ python3 centroid.py ${prefix}_mix.mask.txt ${prefix}_mix
```

finally using the correctted mask to extract valid DNB and further aggregated the DNBs into cells.
```shell
$ python3 gem_mask.py -g ./${prefix}.gem -m ${prefix}_mix.mask.txt -p ${prefix}_mix
```

After finishing these steps, you should obtain cell segmentated gem file and could convert these GEM into anndata using package [stereopy](https://github.com/STOmics/Stereopy).


