
# How to run

## *Pre-process and integrate zebrafish datasets

> The zebrafish datasets were collectted from [Farnsworth et al.](https://www.sciencedirect.com/science/article/pii/S0012160619304919), [Qian et al.](https://link.springer.com/article/10.1007/s00018-022-04410-2) and [Gillotay et al.](https://www.embopress.org/doi/full/10.15252/embr.202050612).

we first integrated zebrafish datasets using seurat
```shell
$ Rscript zebrafish.integrate.R
```
![integrated umap](./img/umap.tissue.png)

the result shows that *Central Nervous System* cells occupied large part of the data, so we further filtered the developmental dataset to keep cells with clear sampling time recorded.

```shell
$ Rscript zebrafish.subset.R
``` 

this will generate the final zebrafish dataset which shall be used in our comparison analysis, including files `zebrafish.final.rds` and `zebrafish.final.loom`


## *SAMap comparison for thyroid like cell and hair cell-like cell

**step 1:** *pseudo-metacell* aggregation and data subsetting

```shell
$ Rscript cal_metacell.R
```
This will generate pseudo-meta cell expression matrix data for endostyle TLC and HCLC, as well as zebrafish pharyngeal endoderm related groups and hair cell types.

**step 2:** SAMap perform

- for endostyle TLC and zebrafish thyroid
```shell
$ python3 performSAMap.py \
    -species1 sc \
    -species2 Da \
    -fn1 sc.TLC.h5ad \
    -fn2 zf.thyroid.h5ad \
    -celltype1 inte_anno \
    -celltype2 Tissue \
    -NUMITERS 3 \
    -cpu 1 \
    -mappingtable TLC.mappingtable.txt \
    -enrichedGene TLC.table.txt

$ awk -F '\t' 'BEGIN{print "source\ttarget\tValue"}{if(NR==1){source=$2}else{print source"\t"$1"\t"$2}}' TLC.mappingtable.txt | awk -F '\t' '$1!=$2' > TLC.for_plot.txt
```

- for endostyle HCLC and zebrafish hair cells
```shell
$ python3 performSAMap.py \
    -species1 sc \
    -species2 Da \
    -fn1 sc.HCLC.h5ad \
    -fn2 zf.cochlea.h5ad \
    -celltype1 inte_anno \
    -celltype2 Tissue \
    -NUMITERS 3 \
    -cpu 1 \
    -mappingtable HCLC.mappingtable.txt \
    -enrichedGene HCLC.table.txt

$ awk -F '\t' 'BEGIN{print "source\ttarget\tValue"}{if(NR==1){source=$2}else{print source"\t"$1"\t"$2}}' HCLC.mappingtable.txt | awk -F '\t' '$1!=$2' > HCLC.for_plot.txt
```

**step 3:** sankey plot

- for TLC cells
```shell
$ Rscript sankeyPlot.R -i TLC.for_plot.txt
```
![TLC sankey](./img/TLC.sankey.png)

- for HCLC cells
```shell
$ Rscript sankeyPlot.R -i HCLC.for_plot.txt
```
![HCLC sankey](./img/HCLC.sankey.png)


## *Additional intepretion

As the result of cross-species comparison largely depend on quality of the dataset we can access, it is necessary to perform with reasonalble filteration process.

In our case, both the zebrafish thyroid data and cochlea data somewhat co-embedded with nervous and immune cells from the zebrafish developmental dataset, as shown in below figure, which might caused by residual sampling.
![splitted umap for thyroid and cochlea cell types](./img/umap.split.png)

So we use diverse filtering criteria for different purpose to avoid this problem. And the final used cell groups in each task was labeled in their code or description seperately.
