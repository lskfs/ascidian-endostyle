import argparse as arg
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy
import pandas as pd
import numpy as np
import scipy

# 
parser = arg.ArgumentParser(formatter_class=arg.RawDescriptionHelpFormatter, description = '''This programe is used to calculate Spearman's correlation coefficient between zebrafish tissues of early and later development stage and sc endostyle different cell types.''')
parser.add_argument('-genes', action='store', help='Input gene list.')
parser.add_argument('-early', action='store', help='Output early heatmap.')
parser.add_argument('-later', action='store', help='Output later heatmap.')
parser.add_argument('-zfDot', action='store', help='Output early dotplot.')
parser.add_argument('-scDot', action='store', help='Output later dotplot.')
args = parser.parse_args()

# 
zf=scanpy.read('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/09.SAMap/04.correlationsHeatmap/DRTFGeneFilter/00.zf.h5ad')
zf.var_names_make_unique()
sc=scanpy.read('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/09.SAMap/04.correlationsHeatmap/DRTFGeneFilter/Styela_clava.anno.h5ad')
sc.var_names_make_unique()
tf=pd.read_csv(args.genes, sep='\t', header=None)
tf=tf.rename({0:'zf', 1:'sc'}, axis=1)

zf.var.index=zf.var['features']
zfbool=[i in zf.var.index for i in tf['zf']]
scbool=[i in sc.var.index for i in tf['sc']]
allbool=[ (i+j)==2 for (i, j) in zip(zfbool, scbool)]
tf=tf[allbool]
print(tf.shape)


zf=zf[: , tf['zf'].to_list()]
sc=sc[: , tf['sc'].to_list()]

earlyTissue=["Mesenchyme-related_arch or fin bud", "Notochord", "Pectoral fin bud", "Pharyngeal endoderm", "Placode", "Tailbud", "Visceral mesoderm"]

laterTissue=["Myeloid lineage", "Thyroid", "Integument", "Cochlea", "Intestine", "Thymus", "Retina", "Connective tissue"]

zf=zf[zf.obs[[i in (earlyTissue + laterTissue) for i in zf.obs.Tissue]].index, :]

# tf dotplot
plt.close()
plt.figure(figsize=(15,15))
scanpy.pl.dotplot(zf, tf['zf'], groupby='Tissue', swap_axes=True, dot_max=1, dot_min=0.2, standard_scale='var')
plt.savefig(args.zfDot, bbox_inches='tight', dpi=300)

plt.close()
plt.figure(figsize=(15,15))
scanpy.pl.dotplot(sc, tf['sc'], groupby='inte_anno', swap_axes=True, dot_max=0.9, dot_min=0.1, standard_scale='var')
plt.savefig(args.scDot, bbox_inches='tight', dpi=300)

zfTissue=pd.DataFrame()
for i in zf.obs['Tissue'].unique():
      tmp=zf.to_df().loc[zf[zf.obs['Tissue']==i].obs.index.to_list(), :]
      zfTissue[i]=np.mean(tmp, axis=0)

scTissue=pd.DataFrame()
for i in sc.obs['inte_anno'].unique():
      tmp=sc.to_df().loc[sc[sc.obs['inte_anno']==i].obs.index.to_list(), :]
      scTissue[i]=np.mean(tmp, axis=0)

scTissue.index=tf['zf']
scTissue.index.name='feature'

df=zfTissue.join(scTissue, how='outer')

# correlation with early tissue
early=df.drop(laterTissue, axis=1)
corr=early.corr('spearman')
corr.drop(earlyTissue, axis=0, inplace=True)
corr.drop(scTissue.columns.to_list(), axis=1, inplace=True)
corr.drop(['Styelin cell A', 'Styelin cell B', 'Undefined cell A', 'Undefined cell B', 'Undefined cell C'], axis=0, inplace=True)

my_colormap = LinearSegmentedColormap.from_list("", ["#84c9b8", "#ffffff", "#e85588"])
plt.close()
plt.figure(figsize=(15,15))
ax1=plt.gca()
sns.heatmap(corr.T, linecolor='#F5F5F5', linewidths=1, cmap=my_colormap, annot=False, ax=ax1)
# add mark of signigicance
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#
#columnD=list(corr.columns)
#indexD=list(corr.index)
#for m in columnD:
#    for n in indexD:
#        c=scipy.stats.spearmanr(early[m], early[n])[0]
#        p=scipy.stats.spearmanr(early[m], early[n])[1]
#        if p<0.01:
#            ax1.text(columnD.index(m)+0.5, indexD.index(n)+0.5, '**', ha='center', color='black')
#        elif p<0.05:
#            ax1.text(columnD.index(m)+0.5, indexD.index(n)+0.5, '*', ha='center', color='black')
#
plt.savefig(args.early, bbox_inches='tight', dpi=600)

# correlation with later tissue
later=df.drop(earlyTissue, axis=1)

#later=later.apply(zscore, axis=0)
corr=later.corr('spearman')
corr.drop(laterTissue, axis=0, inplace=True)
corr.drop(scTissue.columns.to_list(), axis=1, inplace=True)
corr.drop(['Styelin cell A', 'Styelin cell B', 'Undefined cell A', 'Undefined cell B', 'Undefined cell C'], axis=0, inplace=True)

my_colormap = LinearSegmentedColormap.from_list("", ["#84c9b8", "#ffffff", "#e85588"])
plt.close()
plt.figure(figsize=(15,15))
ax1=plt.gca()
sns.heatmap(corr.T, linecolor='#F5F5F5', linewidths=1, cmap=my_colormap, annot=False, ax=ax1)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

## add mark of signigicance
#columnD=list(corr.columns)
#indexD=list(corr.index)
#
#for m in columnD:
#    for n in indexD:
#        c=scipy.stats.spearmanr(later[m], later[n])[0]
#        p=scipy.stats.spearmanr(later[m], later[n])[1]
#        if p<0.01:
#            ax1.text(columnD.index(m)+0.5, indexD.index(n)+0.5, '**', ha='center', color='black')
#        elif p<0.05:
#            ax1.text(columnD.index(m)+0.5, indexD.index(n)+0.5, '*', ha='center', color='black')
plt.savefig(args.later, bbox_inches='tight', dpi=600)

print('------'+ args.genes+ '------')
