
import pandas as pd
import anndata

### extract sinus region cells
adata = anndata.read('../data/Styela_clava.anno.h5ad')
target_celltypes = ['BC', 'CPC', 'CSC', 'LLF', 'M', 'MLLC', 'NC', 'NCSC', 'PLLC', 'SCA', 'SCB', 'SLLC']
adata = adata[adata.obs['short_inte_anno'].isin(target_celltypes)]

sc_adata = adata[adata.obs['Batch'].isin(['sc1', 'sc2', 'sc3'])]
sc_obs = sc_adata.obs.copy()

st_adata = adata[~adata.obs['Batch'].isin(['sc1', 'sc2', 'sc3'])]
st_adata.obs['label'] = st_adata.obs['label'].astype(str).astype(int)
st_adata = st_adata[st_adata.obs['label'] > 0]
st_obs = st_adata.obs.copy()

obs = pd.concat([sc_obs, st_obs])
obs.to_csv('selectted_cells.obs.txt', sep='\t')
