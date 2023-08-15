
cwd=$PWD

# single cell seurat pipeline for each sample

Rscript single-cell-seurat-MT-211101.R \
-i $cwd/Sc-En-1/outs/filtered_feature_bc_matrix \
-R $cwd/Sc-En-1/outs/raw_feature_bc_matrix \
-o $cwd/Sc_En_1_output \
-P Sc-En-1  \
-M $cwd/mito-gene-list.txt \
-r 1.0 -c 0 -C 5000 -f 300 -F 1200 -d 15 -t 10 -m 30 -D 0.075

Rscript single-cell-seurat-MT-211101.R \
-i $cwd/Sc-En-2/outs/filtered_feature_bc_matrix \
-R $cwd/Sc-En-2/outs/raw_feature_bc_matrix \
-o $cwd/Sc_En_2_output \
-P Sc-En-2  \
-M $cwd/mito-gene-list.txt \
-r 1.0 -c 0 -C 5000 -f 300 -F 1200 -d 15 -t 10 -m 45 -D 0.075

Rscript single-cell-seurat-MT-211101.R \
-i $cwd/Sc-En-3/outs/filtered_feature_bc_matrix \
-R $cwd/Sc-En-3/outs/raw_feature_bc_matrix \
-o $cwd/Sc_En_3_output \
-P Sc-En-3  \
-M $cwd/mito-gene-list.txt \
-r 1.0 -c 0 -C 5000 -f 300 -F 1200 -d 15 -t 10 -m 40 -D 0.075

# integrate three filited single cell datasets
Rscript integration.singlecell.SCT.R \
	-d path \
	-b Sc-mito30.45.40-inte \
	-o $cwd/220118-Sc-mito30.45.40/ \
	-r 0.3 -z 2 -D 10

