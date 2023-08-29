
# How to run

define the work directory and slice prefix which to calculate

```shell
cwd=$PWD
prefix='C3-2'
```

**step 1:** extract alignment records from original bam file based on the cell mask, from which you will obtain resultant bam file named `${prefix}.bam` including tissue covered region.
```shell
python3 pre-velocyto.py \
    --out ./${prefix}.bam \
    --mask ../01.stereo-seq_cellseg/${prefix}_mix.mask.txt \
    --bam ./DP8400015160TR_C3.bam \ ### original aligned bam file from SAW pipeline
    --offsetx 13 \ ### offsetx of original GEM file
    --offsety 4209 \ ### offsety of original GEM file
    --shiftx 2072 \ ### shiftx calculated in 01.stereo-seq_cellseg step
    --shifty 4732 ### shifty calculated in 01.stereo-seq_cellseg step
```

**step 2:** run velocyto based on the bam file and gene structure annotation to generate splice/un-splice matrix in loom  format
```shell
python3 velocyto run \
    ./${prefix}.bam \
    ./updated-20211028.MT.renamed.gtf \
    --outputfolder ./ \
    --samtools-threads 1 \
    --samtools-memory 3000 \
    -t int
```

**step 3:** run dynamo using the loom file
```shell
python3 dynamo_pipe.single.py 
```
