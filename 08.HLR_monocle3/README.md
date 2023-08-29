
# How to run

**step 1:** extract sinus region cells

```shell
$ python3 get_sinus_target_cells.py
```

**step 2:** get potential stem cell index from the single cell data, which recorded in *`stemcell.index.txt`*

**step 3:** run monocle3 to infer trajectory

```shell
$ Rscript monocle3.R
```

