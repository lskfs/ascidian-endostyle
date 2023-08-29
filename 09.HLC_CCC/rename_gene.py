
import sys
from collections import defaultdict
import pandas as pd
import anndata

data = defaultdict(list)
r = open('styela_human.ortholog.tsv')
next(r)
for line in r:
    lst = line.split()
    styela = lst[4].split(',')
    human = lst[5].split(',')
    if len(styela) > 1:
        continue
    for g in human:
        g = g.split('.')[0]
        data['styela'].append(styela[0])
        data['human'].append(g)
r.close()
name_mapping = pd.DataFrame.from_dict(data, orient='columns')

gene2name = pd.read_csv('./human.gene_name.txt', sep='\t', header=None, names=['human', 'gene_name'])
name_mapping = name_mapping.merge(gene2name, how='left', on=['human'])
name_mapping = name_mapping[~name_mapping['gene_name'].isna()]
name_mapping = name_mapping[['styela', 'gene_name']].drop_duplicates(subset='styela')
name_mapping = name_mapping.to_dict('tight')['data']
name_mapping = dict(name_mapping)
#name_mapping = name_mapping[['styela', 'gene_name']].drop_duplicates(subset='styela').set_index('styela')

adata = anndata.read('../data/Styela_clava.anno.h5ad')
adata = adata[:, list(name_mapping.keys())]
new_gene_names = [name_mapping[x] if name_mapping.get(x) else x for x in adata.var_names]
adata.var_names = new_gene_names

adata.var_names_make_unique()
adata.write('Styela_clava.renamed_by_human.h5ad')


