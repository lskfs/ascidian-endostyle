import scanpy as sc
import hotspot

import numpy as np
import mplscience
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import pickle

sys.path.append('/dellfsqd2/ST_OCEAN/USER/hankai/Project/10.smed/code/site-packages')
from WBRtools import get_cmap

def performHotSpot(adata, sliceId):
    adata=adata[adata.obs['Batch']==sliceId]
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs["total_counts"] = np.asarray(adata.X.sum(1)).ravel()
    adata.layers["csc_counts"] = adata.X.tocsc()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000) 
    # Create the Hotspot object and the neighborhood graph 
    hs = hotspot.Hotspot(
        adata, 
        layer_key="csc_counts",
        model='normal', 
        latent_obsm_key="spatial", 
        umi_counts_obs_key="total_counts",
    )
    
    hs.create_knn_graph(
        weighted_graph=False, n_neighbors=300,
    
    )
    '''
    weighted_graph: bool
        Whether or not to create a weighted graph
    n_neighbors: int
        Neighborhood size
    neighborhood_factor: float
        Used when creating a weighted graph. Sets how quickly weights decay
        relative to the distances within the neighborhood. The weight for
        a cell with a distance d will decay as exp(-d/D) where D is the distance
        to the `n_neighbors`/`neighborhood_factor`-th neighbor.
    approx_neighbors: bool
        Use approximate nearest neighbors or exact scikit-learn neighbors. Only
        when hotspot initialized with `latent`.
    '''

    hs_results = hs.compute_autocorrelations(jobs=4)
    
    # Select the genes with significant spatial autocorrelation
    hs_genes = hs_results.index[hs_results.FDR < 0.05]
    
    # Compute pair-wise local correlations between these genes
    lcz = hs.compute_local_correlations(hs_genes, jobs=4)
    
    modules = hs.create_modules(
        min_gene_threshold=20, core_only=False, fdr_threshold=0.05
    )
    '''
    min_gene_threshold: int
        Controls how small modules can be. Increase if there are too many
        modules being formed. Decrease if substructre is not being captured
    core_only: bool
        Whether or not to assign ambiguous genes to a module or leave unassigned
    fdr_threshold: float
        Correlation theshold at which to stop assigning genes to modules
    ''' 
    
    with open ("hs.txt", 'wb') as f:
        pickle.dump(hs, f)

    zcmap = get_cmap(['#84c9b8', '#ffffff', '#e85588'])

    plt.figure()
    hs.plot_local_correlations(vmin=-5, vmax=5, z_cmap=zcmap)
    plt.savefig('{}_hotspot_correlations.pdf'.format(sliceId))
    plt.figure()
    hs.plot_local_correlations(vmin=-5, vmax=5)
    plt.savefig('{}_hotspot_correlations.tiff'.format(sliceId))
    
    results = hs.results.join(hs.modules).sort_values('Z', ascending=False).to_csv('{}_result.txt'.format(sliceId))
    
    module_scores=hs.calculate_module_scores()
    
    module_cols=[]
    for c in module_scores.columns:
        key=f"Module {c}"
        adata.obs[key]=module_scores[c]
        module_cols.append(key)
   
    adata.write_h5ad('{}.h5ad'.format(sliceId), compression='gzip') 
    
    plt.figure()
    with mplscience.style_context():
        sc.pl.spatial(adata, color=module_cols, frameon=False, vmin="p0", vmax="p99", spot_size=30)
    plt.savefig('{}_hotspot_moudle.pdf'.format(sliceId))
    
    plt.figure()
    with mplscience.style_context():
        sc.pl.spatial(adata, color=module_cols, frameon=False, vmin="p0", vmax="p99", spot_size=30)
    plt.savefig('{}_hotspot_moudle.tiff'.format(sliceId))

adata=sc.read("/dellfsqd2/ST_OCEAN/USER/wangrui21/04.project/01.styela_clava/03.hotspot/01.all_slide/Styela_clava.h5ad")

sliceList=['C3-2']

for i in sliceList:
    performHotSpot(adata, i)  
