
import pandas as pd

df = pd.read_csv('../monocle2/genSmoothCurves_mat.scaled.txt', sep='\t', header=0, index_col=0)
df = df.T
maxPT = df.idxmax().to_frame(name='class').reset_index().rename(columns={'index': 'gene_name'})

mapping = pd.read_csv('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/13.sonic/styela_to_zebraname.txt', 
        sep='\t', header=0)
"""
header = ['styela', 'gene_name', 'ident', 'len', 'mismatch', 'gap', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
mapping = pd.read_csv('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/00.data/styela_clava2zebrafish.blast/sc2da.renamed.filtered.out', 
        sep='\t', header=None, names=header)
mapping = mapping[mapping['ident'] >= 60]
mapping = mapping.sort_values('ident', ascending=False).groupby(['styela', 'gene_name']).head(1)
mapping = mapping[['styela', 'gene_name']].drop_duplicates()
"""
maxPT = maxPT.merge(mapping, how='left', on=['gene_name'])
maxPT = maxPT.dropna().drop_duplicates()
maxPT.to_csv('styela_gene_trajectory_class.txt', sep='\t', index=False)



