
import argparse as arg
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 6})

import anndata
import scanpy

parser = arg.ArgumentParser( description = 'spearman correlation calculation')
parser.add_argument('-input', action='store', help='TF matches list')
parser.add_argument('-figname', action='store', help='output pdf figure filename')
parser.add_argument('-type', choices=['single', 'multi'], action='store', default='single', help='which type of the input TF list, "single" for one-to-one and "multi" for one-to-many')
parser.add_argument('-strategy', choices=['mean', 'copy'], action='store', default='mean', help='"mean" for averaging one-to-many, "copy" for making copies for one-to-many')
args = parser.parse_args()

### load homologous TF list
suffix = args.suffix
tf = pd.read_csv(args.input, sep='\t', header=0)
tf = tf[['Symbol', 'Styela_clava.pep', 'is_single_copy']].drop_duplicates()
tf = tf.rename(columns={'Symbol':'zf', 'Styela_clava.pep': 'sc'})
tf['sc'] = tf['sc'].str.replace('_', '-', regex=False)

### load h5ad data
zf = anndata.read_loom('zebrafish.loom')
zf.layers['data'] = zf.X
zf.raw = zf
zf.X = zf.layers['counts']
zf_gene_subset, number_per_gene = scanpy.pp.filter_genes(zf, min_cells=3, inplace=False)
zf.X = zf.layers['data']
zf_genes = zf[:, zf_gene_subset].var_names

sc = anndata.read_loom('styela_clava.loom')
sc.layers['data'] = sc.X
sc.raw = sc
sc.X = sc.layers['counts']
sc_gene_subset, number_per_gene = scanpy.pp.filter_genes(sc, min_cells=3, inplace=False)
sc.X = sc.layers['data']
sc_genes = sc[:, sc_gene_subset].var_names

### filtering tf by expression
tf = tf[tf['zf'].isin(zf_genes) & tf['sc'].isin(sc_genes)]

### subset h5ad data
zf_genes = tf['zf'].unique().tolist()
zf = zf[:, zf_genes]

sc_genes = tf['sc'].unique().tolist()
sc = sc[:, sc_genes]

### drop undefined cell types
excluded = ['Undefined cell A', 'Undefined cell B', 'Undefined cell C']
sc_celltypes = [x for x in sc.obs['inte_anno.202302'].unique() if x not in excluded]
sc = sc[sc.obs['inte_anno.202302'].isin(sc_celltypes)]

### function defination for average expression
def grouped_obs_mean(adata, group_key, layer=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
            np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
            columns=list(grouped.groups.keys()),
            index=adata.var_names
            )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out

### function defination for process multiple copy homologous genes
def multicopy_process(df1, df2, mapping):
    """
    for one-to-many homologous relationship, make a copy for the 'one'
    """
    
    copy_idx = lambda x: [x[:(i+1)].count(v) for i, v in enumerate(x)]

    mapping['zf.idx'] = copy_idx(mapping['zf'].to_list())
    mapping['sc.idx'] = copy_idx(mapping['sc'].to_list())

    df1 = df1.loc[mapping['zf'].values]
    df1.index = df1.index + '.' + mapping['zf.idx'].astype(str)

    df2 = df2.loc[mapping['sc'].values]
    df2.index = df2.index + '.' + mapping['sc.idx'].astype(str)

    mapping['zf'] = mapping['zf'] + '.' + mapping['zf.idx'].astype(str)
    mapping['sc'] = mapping['sc'] + '.' + mapping['sc.idx'].astype(str)

    df2_df1 = dict(mapping[['sc', 'zf']].to_dict('tight')['data'])
    df2 = df2.rename(index=df2_df1)
    
    data = df2.join(df1, how='outer')
    return data

def multimean_process(df1, df2, mapping):
    """
    for one-to-many homologous relationship, calculate mean expression for 'many'
    """
    
    name_mapping = mapping[['zf', 'sc']].drop_duplicates().set_index('zf')

    name_mapping = name_mapping[name_mapping['sc'].isin(df2.index.values)]
    
    df1 = df1.merge(name_mapping, how='left', left_index=True, right_index=True)
    df1 = df1.groupby('sc').mean()

    ordered_index = name_mapping['sc'].unique()
    df1 = df1.loc[ordered_index]
    df2 = df2.loc[ordered_index]

    data = df2.join(df1, how='outer')
    return data

zf_mean = grouped_obs_mean(zf, 'Tissue')
sc_mean = grouped_obs_mean(sc, 'inte_anno.202302')

if args.strategy == 'copy':
    data = multicopy_process(zf_mean, sc_mean, tf)
elif args.strategy == 'mean':
    data = multimean_process(zf_mean, sc_mean, tf)

### function defination of plotting
def corr_barplot(df, xlabel='zebrafish', outfile=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(dpi=600, figsize=(2, 5))
    else:
        fig = plt.gcf()
        ax = ax
    df.plot.barh(color='#e85588', ax=ax, 
            xlabel=xlabel, 
            ylabel='Correlation coefficient',
            xlim=(0, 0.25),
            )
    ax.xaxis.set_ticks_position('top')
    if outfile is not None:
        plt.tight_layout()
        fig.savefig(outfile)
    return ax

### figure plotting
fig = plt.figure(layout='tight', figsize=(8, 4), dpi=600)
gs = fig.add_gridspec(nrows=2, ncols=4, hspace=0.1, wspace=0.1)

### correlation cofficient between earlyTissue in zebrafish and sc endostyle
##### for embryo lineages of zebrafish
early_lineage = [
                 'Hypoblast', 'Hypochord', 'Lateral plate mesoderm', 'Lens placode', 
                 'Mesenchyme-related_arch or fin bud', 'Neural crest', 'Notochord',
                 'Paraxial mesoderm', 'Pectoral fin bud', 'Pharyngeal arch', 'Pharyngeal endoderm',
                 'Placode', 'Tailbud', 'Visceral mesoderm'
                ]
sc_to_early = data[early_lineage + sc_celltypes]
early_corr = sc_to_early.corr('spearman')
early_corr = early_corr.loc[sc_celltypes, early_lineage]
early_corr = early_corr.mean(axis=0).to_frame(name='corr').sort_values(by='corr')
ax0 = fig.add_subplot(gs[:, 0])
ax0 = corr_barplot(early_corr, xlabel='Embryo lineage of zebrafish', ax=ax0)

##### for tissues of zebrafish
later_tissue = [
                'Blood', 'Blood vessel', 'Cochlea', 'Connective tissue', 'Integument',
                'Intestine', 'Kidney', 'Liver', 'Muscle', 'Myeloid lineage', 'Nervous system',
                'Pancreas', 'Pineal gland', 'Retina', 'Spleen', 'Thymus', 'Thyroid'
               ]
sc_to_later = data[later_tissue + sc_celltypes]
later_corr = sc_to_later.corr('spearman')
later_corr = later_corr.loc[sc_celltypes, later_tissue]
later_corr = later_corr.mean(axis=0).to_frame(name='corr').sort_values(by='corr')
ax1 = fig.add_subplot(gs[:, 1])
ax1 = corr_barplot(later_corr, xlabel='Tissue of zebrafish', ax=ax1)

import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

early_top = early_corr.index.values.tolist()[::-1]
later_top = later_corr.index.values.tolist()[::-1]
data = data[early_top + later_top + sc_celltypes]
cmap = LinearSegmentedColormap.from_list('', ['#84c9b8', '#ffffff', '#e85588'])

##### for endostyle celltype to embryo lineages
sc_to_early_top = data[early_top + sc_celltypes]
et_corr = sc_to_early_top.corr('spearman')
et_corr = et_corr.loc[sc_celltypes, early_top]
ax2 = fig.add_subplot(gs[0, 2:])
sns.heatmap(et_corr.T, linecolor='#f5f5f5', linewidths=0.8, cmap=cmap, annot=False, ax=ax2, yticklabels=1)

##### for endostyle celltype to tissues
sc_to_later_top = data[later_top + sc_celltypes]
lt_corr = sc_to_later_top.corr('spearman')
lt_corr = lt_corr.loc[sc_celltypes, later_top]
ax3 = fig.add_subplot(gs[1, 2:], sharex=ax2)
sns.heatmap(lt_corr.T, linecolor='#f5f5f5', linewidths=0.8, cmap=cmap, annot=False, ax=ax3, yticklabels=1)

### save figure
fig.savefig(f'{args.figname}.pdf')







