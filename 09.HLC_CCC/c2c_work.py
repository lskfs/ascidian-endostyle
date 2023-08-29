import pandas as pd
import anndata
import cell2cell as c2c

adata = anndata.read('./Styela_clava.renamed_by_human.h5ad')
colors = dict(adata.obs[['short_inte_anno.202302', 'Color']].drop_duplicates().to_dict('tight')['data'])

lr_pairs = pd.read_csv('./Human-2020-Jin-LR-pairs.csv')
lr_pairs = lr_pairs.astype(str)
meta = adata.obs.copy()

interactions = c2c.analysis.SingleCellInteractions(
        rnaseq_data=adata.to_df().T,
        ppi_data=lr_pairs,
        metadata=meta,
        interaction_columns=('ligand_symbol', 'receptor_symbol'),
        communication_score='expression_thresholding',
        expression_threshold=0.0001, # values after aggregation
        cci_score='bray_curtis',
        cci_type='undirected',
        aggregation_method='nn_cell_fraction',
        barcode_col='index',
        celltype_col='short_inte_anno.202302',
        complex_sep='&',
        verbose=False
        )

interactions.compute_pairwise_communication_scores()
interactions.compute_pairwise_cci_scores()
cci_pvals = interactions.permute_cell_labels(
        evaluation='interactions', 
        permutations=10, 
        fdr_correction=False,
        verbose=True
        )

ccc_pvals = interactions.permute_cell_labels(
        evaluation='communication',
        permutations=10,
        fdr_correction=False,
        verbose=True
        )

group_meta = pd.DataFrame(columns=['Celltype', 'Group'])
group_meta['Celltype'] = meta['short_inte_anno.202302'].unique().tolist()
group_meta['Group'] = group_meta['Celltype']

interaction_clustermap = c2c.plotting.clustermap_ccc(
        interactions,
        metric='jaccard',
        method='complete',
        metadata=group_meta,
        sample_col='Celltype',
        group_col='Group',
        colors=colors,
        row_fontsize=14,
        title='Active ligand-receptor pairs for interacting cells',
        filename=None,
        cell_labels=('SENDER-CELLS', 'RECEIVER-CELLS'),
        **{'figsize' : (10,9),}
        )

# Add a legend to know the groups of the sender and receiver cells:
l1 = c2c.plotting.generate_legend(
        color_dict=colors,
        loc='center left',
        bbox_to_anchor=(20, -2), # Indicated where to include it
        ncol=1, fancybox=True,
        shadow=True,
        title='Groups',
        fontsize=14,
        )
interaction_clustermap.savefig('haircell/clustermap.png')

sender_cells = ['HC', 'IRFC', 'P+SC', 'P-SC']
receiver_cells = ['HC', 'IRFC', 'P+SC', 'P-SC']

ligands = lr_pairs['ligand'].unique()
receptors = lr_pairs['receptor'].unique()

c2c.plotting.circos_plot(
        interaction_space=interactions,
        sender_cells=sender_cells,
        receiver_cells=receiver_cells,
        ligands=ligands,
        receptors=receptors,
        excluded_score=0,
        metadata=group_meta,
        sample_col='Celltype',
        group_col='Group',
        colors=colors,
        fontsize=15,
        ligand_label_color='orange',
        receptor_label_color='brown',
        filename='circos_plot.pdf',
        legend=True,
        )


from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('', ['#0a40a3', '#ffffff', '#8d033c'])

fig = c2c.plotting.dot_plot(
        interactions,
        evaluation='communication',
        significance = 0.11,
        figsize=(16, 9),
        cmap=cmap,
        senders=sender_cells,
        receivers=receiver_cells
        )
fig.savefig('./dot_plot.pdf')

cm = c2c.plotting.clustermap_cci(
        interactions,
        method='complete',
        metadata=group_meta,
        sample_col="Celltype",
        group_col="Group",
        colors=colors,
        title='CCI scores for cell-types',
        cmap='Blues'
        )
# Add a legend to know the groups of the sender and receiver cells:
l1 = c2c.plotting.generate_legend(
        color_dict=colors,
        loc='center left',
        bbox_to_anchor=(20, -2), # Indicated where to include it
        ncol=1, fancybox=True,
        shadow=True,
        title='Groups',
        fontsize=14,
        )
cm.savefig('./clustermap_cci.png')

fig = c2c.plotting.dot_plot(
        interactions,
        evaluation='interactions',
        significance = 0.11,
        figsize=(16, 9),
        cmap='Blues',
        )
fig.savefig('./significance_dot_plot.png')


















