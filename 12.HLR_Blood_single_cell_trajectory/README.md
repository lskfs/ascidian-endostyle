
# How to run

## *Pre-process of sequencing reads

The sequencing reads from HLR samples were downsampled to eliminate influence of different sequencing deepth 
```shell
$ sh 01.downsampling.sh
```

Then we used the SeekSoulTools to generate gene expression matrix from raw reads.
```shell
$ sh 01.seeksoultools.sh
```

## *Seurat clustering of each sample

We first clustering the blood and HLR samples seperately using Seurat.
```shell
$ sh 02.work-Seurat.sh
```

## *Integration of blood and HLR data

Then the processed seurat objects from blood and HLR samples were integrated using Seurat.
```shell
$ Rscript 03.Integration.R
```

## *Trajectory inference using monocle2

Clusters C2, C9 and C11 that we inferred as HSC-related cells were extracted for further analysis
```shell
$ Rscript 04.SeuratSubset.R
```

Potential differentiation trajectory was inferred using monocle2
```shell
$ Rscript 05.Monocle2.R
```

