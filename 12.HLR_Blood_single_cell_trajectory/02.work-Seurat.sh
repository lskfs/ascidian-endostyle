Rscript single-cell-seurat.R \
	-i stcr-HLRonly-downsampling/step3/filtered_feature_bc_matrix \
	-R stcr-HLRonly-downsampling/step3/raw_feature_bc_matrix \
	-o HLRdownsp_output \
	-P HLRdownsp \
	-M mito-gene-list.txt \
	-r 0.6 -c 100 -C 3000 -f 50 -F 1000 -d 15 -t 10 -m 40 


Rscript single-cell-seurat.R \
	-i stcr-blood-2/step3/filtered_feature_bc_matrix \
	-R stcr-blood-2/step3/raw_feature_bc_matrix \
	-o stcr-blood-2_output \
	-P stcr-blood-2  \
	-M mito-gene-list.txt \
	-r 0.6 -c 0 -C 1000 -f 10 -F 500 -d 15 -t 10 -m 50 
	