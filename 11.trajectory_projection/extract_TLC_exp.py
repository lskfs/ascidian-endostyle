
import pandas as pd
import anndata
import scanpy as sc

trajectory_genes = pd.read_csv('styela_gene_trajectory_class.txt', sep='\t', header=0)
trajectory_genes = trajectory_genes[['class', 'styela']].drop_duplicates()
genes = trajectory_genes['styela'].unique()

adata = anndata.read('../data/Styela_clava.anno.h5ad')
adata.X = adata.layers['data']
adata = adata[adata.obs['short_inte_anno'] == 'TLC']

sc.pp.scale(adata)

genes = [x for x in genes if x in adata.var_names]
exp = adata[:, list(genes)].to_df()
exp = exp.T.reset_index().rename(columns={'index': 'styela'})
exp = exp.merge(trajectory_genes, how='left', on=['styela'])
exp.to_csv('exp.txt', sep='\t', index=False)

