
import sys
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd

import anndata
import scanpy as sc
import dynamo as dyn

batch = sys.argv[1]

adata = anndata.read('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/16.sinus_trajectory/dynamo/Styela_clava.anno.velocyto.h5ad')
color_dict = dict(adata.obs[['inte_anno', 'Color']].to_dict('tight')['data'])

centroid_meta = []
centroid_path = '/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/15.registration_segmentation/02.mix_segmentation.v2'
for stat_file in Path(centroid_path).glob('*.centroid.csv'):
    xy = pd.read_csv(stat_file, sep='\t', header=0)
    xy['Batch'] = stat_file.stem.replace('.centroid', '')
    xy['CellID'] = xy['Batch'] + ':' + xy['label'].astype(int).astype(str)
    centroid_meta.append(xy)
centroid_meta = pd.concat(centroid_meta)
centroid_meta = centroid_meta.set_index('CellID')
centroid_meta = centroid_meta[['x', 'y']]
adata.obs = adata.obs.merge(centroid_meta, how='left', left_index=True, right_index=True)
adata.obsm['spatial'] = adata.obs[['x', 'y']].values
adata.obsm['X_spatial'] = adata.obs[['x', 'y']].values
adata = adata[adata.obs['Batch'] == batch]
print(adata.obs)

selectted_cells = pd.read_csv('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/16.sinus_trajectory/selectted_cells.obs.txt', 
        sep='\t', header=0, index_col=0)
selectted_cells = selectted_cells[selectted_cells['Batch'] == batch]
selectted_cells['CellID'] = selectted_cells['Batch'] + ':' + selectted_cells['label'].astype(int).astype(str)
selectted_cells = [x for x in selectted_cells['CellID'].values if x in adata.obs.index.values]
adata = adata[selectted_cells]

dyn.pl.basic_stats(adata, 
        save_show_or_return='save',
        save_kwargs={'prefix': 'basic', 'ext': 'png'},
        )
dyn.pl.show_fraction(adata, 
        save_show_or_return='save',
        save_kwargs={'prefix': 'fraction', 'ext': 'png'},
        )

dyn.pp.recipe_monocle(adata, num_dim=30, normalized=False, 
        #genes_to_append=markers,
        exprs_frac_for_gene_exclusion=0.005, 
        fc_kwargs={'min_expr_genes_s':0, 'min_expr_genes_u':0},
        fg_kwargs={'min_cell_s':5, 'min_cell_u':5, 'shared_count': 5},
        )

dyn.tl.dynamics(adata, cores=1, 
        filter_gene_mode='final', 
        est_method='negbin', 
        )

sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata)
min_dist = 0.05
sc.tl.umap(adata, min_dist=min_dist)
dyn.pl.umap(adata, color=['ntr', 'inte_anno'], 
        ncols=2, 
        pointsize=0.3,
        save_show_or_return='save',
        save_kwargs={'prefix': 'umap_scanpy', 'ext': 'png'},
        )

adata.obsm['X_umap_corrected'] = adata.obsm['X_umap'].copy()
#dyn.tl.reduceDimension(adata, basis='umap', enforce=True)

dyn.pl.umap(adata, color=['ntr', 'inte_anno'], 
        ncols=2, 
        pointsize=0.3,
        save_show_or_return='save',
        save_kwargs={'prefix': 'umap', 'ext': 'png'},
        )

fig, axes = plt.subplots(1, 2)
adata.var.gamma.hist(ax=axes[0])
adata.var.gamma_r2.hist(ax=axes[1])
fig.savefig('hist.png')

# get -5 from 'hist.png'
dyn.tl.cell_velocities(adata, method='pearson', 
        other_kernels_dict={'transform': 'sqrt'}, 
        min_r2=-5
        )
dyn.tl.cell_wise_confidence(adata)

tgn = adata.var.use_for_transition.sum()
print(f'detect {tgn} genes used for transition', flush=True)

#dyn.tl.cell_wise_confidence(adata)
dyn.pl.cell_wise_vectors(adata, color=['inte_anno'], basis='umap', show_legend='on data', 
        #quiver_length=6, quiver_size=6, pointsize=0.1, 
        save_show_or_return='save',
        save_kwargs={'prefix': f'cellwise_umap', 'ext': 'png'},
        )

dyn.pl.streamline_plot(adata, basis='umap', show_legend='on data', color=['inte_anno', 'ntr'],
        s_kwargs_dict={'alpha':1, 'rasterized': False},
        inset_dict={'aspect':'equal'},
        quiver_length=6, quiver_size=6, pointsize=0.1, 
        frontier=True, density=0.5,
        show_arrowed_spines=True, save_show_or_return='save',
        save_kwargs={'prefix': f'streamline_umap', 'ext': 'svg'},
        )

print('hankai: go start spatial dim')
dyn.tl.cell_velocities(adata, basis='spatial')
print('hankai: go start spatial dim plot')
dyn.pl.streamline_plot(adata, basis='spatial', show_legend='on data', color=['inte_anno'],
        color_key=color_dict,
        inset_dict={'aspect':'equal'},
        s_kwargs_dict={'alpha':1, 'rasterized': False},
        show_arrowed_spines=False, pointsize=0.1, 
        frontier=True, density=0.5,
        save_show_or_return='save',
        save_kwargs={'prefix': f'streamline_spatial', 'ext': 'svg'},
        )
print('hankai: finish spatial dim')

dyn.vf.VectorField(adata, basis='umap', M=1000)
dyn.pl.plot_energy(adata, basis='umap', save_show_or_return='save', 
        save_kwargs={'prefix': f'energy', 'ext': 'png'})
dyn.pl.topography(adata, basis='umap', background='white', 
        color=['ntr', 'inte_anno'], 
        inset_dict={'aspect':'equal'},
        s_kwargs_dict={'alpha':1, 'rasterized': False},
        streamline_color='black', show_legend='on data', frontier=True, 
        save_show_or_return='save', save_kwargs={'prefix': f'topography', 'ext': 'svg'},
        )

dyn.ext.ddhodge(adata, basis='umap')
dyn.pl.streamline_plot(adata, basis='umap', show_legend='on data', 
        color=['inte_anno'], color_key=color_dict, 
        alpha=1,
        inset_dict={'aspect':'equal'},
        s_kwargs_dict={'alpha':1, 'rasterized': False},
        frontier=True, density=0.5,
        quiver_length=6, quiver_size=6, pointsize=0.1, 
        show_arrowed_spines=True, save_show_or_return='save',
        save_kwargs={'prefix': f'stream_on_clusters', 'ext': 'svg'},
        )
dyn.pl.streamline_plot(adata, basis='umap', show_legend='on data', 
        color=['umap_ddhodge_potential'], 
        inset_dict={'aspect':'equal'},
        s_kwargs_dict={'alpha':1, 'rasterized': False},
        quiver_length=6, quiver_size=6, pointsize=0.1, 
        frontier=True, density=0.5,
        show_arrowed_spines=True, save_show_or_return='save',
        save_kwargs={'prefix': f'ddhodge', 'ext': 'svg'},
        )

"""
# Differential geometry analysis
dyn.tl.cell_velocities(adata, basis='pca')
dyn.vf.VectorField(adata, basis='pca', M=100)

dyn.vf.rank_velocity_genes(adata, groups='inte_anno', vkey='velocity_S')

dyn.vf.acceleration(adata, basis='pca')
dyn.vf.rank_acceleration_genes(adata, groups='inte_anno', akey='acceleration', prefix_store='rank')

dyn.vf.curvature(adata, basis='pca')
dyn.vf.rank_curvature_genes(adata, groups='inte_anno')

dyn.vf.divergence(adata, basis='pca')
#dyn.vf.rank_divergence_genes(adata, groups='inte_anno', akey='divergence', prefix_store='rank')

dyn.vf.speed(adata, basis='pca')
dyn.vf.curl(adata, basis='umap')

dyn.export_rank_xlsx(adata, path='rank_info.xlsx')

dyn.vf.cluster_field(adata, basis='pca',
        features=['speed', 'potential', 'divergence', 'acceleration', 'curvature']
        )

dyn.pl.umap(adata, color=['field_leiden', 'acceleration_pca', 'curvature_pca', 'divergence_pca', 'speed_pca', 'curl_umap'], 
        save_show_or_return='save',
        save_kwargs={'prefix': f'geometry', 'ext': 'png'},
        )
"""

adata.obs.to_csv(f'dynamo_obs.{batch}.txt', sep='\t')
adata.var.to_csv(f'dynamo_var.{batch}.txt', sep='\t')

#subset = adata[adata.obs['inte_anno'].isin(['c1', 'c31', 'c20'])]
transition_genes = adata.var[adata.var['use_for_transition']].index
df = dyn.pl.kinetic_heatmap(adata,
        genes=transition_genes,
        tkey='umap_ddhodge_potential',
        gene_order_method='maximum',
        mode='pseudotime',
        cell_group='inte_anno',
        standard_scale=1,
        log=False,
        show_colorbar=True,
        yticklabels=True,
        figsize=(10, 15),
        color_map='RdYlBu_r',
        save_show_or_return='save',
        save_kwargs={'prefix': 'transition_genes_kinetic_heatmap', 'ext': 'pdf'}
        )

"""
marker_list = pd.read_csv('../sct.integrated_AllMarkers.xls', sep='\t', header=0)
marker_list = marker_list[marker_list['cluster'].isin([1, 31, 20])]
marker_list = marker_list[marker_list['avg_log2FC'] >= 1]
marker_list = marker_list['gene'].unique()
genes = set(transition_genes).intersection(set(marker_list))
dyn.pl.kinetic_heatmap(subset,
        genes=genes,
        tkey='umap_ddhodge_potential',
        gene_order_method='maximum',
        mode='pseudotime',
        cell_group='inte_anno',
        standard_scale=1,
        log=False,
        show_colorbar=True,
        yticklabels=True,
        figsize=(10, 15),
        color_map='RdYlBu_r',
        save_show_or_return='save',
        save_kwargs={'prefix': 'tgmg_kinetic_heatmap', 'ext': 'pdf'}
        )
"""

dyn.cleanup(adata)
adata.write(f'dynamo_final.{batch}.h5ad')


