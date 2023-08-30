
import pandas as pd

df = pd.read_csv('./genSmoothCurves_mat.scaled.txt', sep='\t', header=0, index_col=0)
df = df.T
maxPT = df.idxmax().to_frame(name='class').reset_index().rename(columns={'index': 'gene_name'})

mapping = pd.read_csv('./data/styela_to_zebraname.orth.txt', sep='\t', header=0)

maxPT = maxPT.merge(mapping, how='left', on=['gene_name'])
maxPT = maxPT.dropna().drop_duplicates()
maxPT.to_csv('styela_gene_trajectory_class.txt', sep='\t', index=False)



