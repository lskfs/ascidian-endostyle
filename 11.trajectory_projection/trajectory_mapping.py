
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

##### create pseudotime colormap
norm = mpl.colors.Normalize(vmin=1, vmax=100)
cell_cmap = plt.get_cmap('plasma')
cell_color_arr = np.array([mpl.colors.to_hex(cell_cmap(norm(x))) for x in range(1, 101)])
gene_cmap = sns.color_palette("blend:#7AB,#EDA", as_cmap=True)
gene_color_arr = np.array([mpl.colors.to_hex(gene_cmap(norm(x))) for x in range(1, 101)])
class_arr = [x for x in range(1, 101)]
color_pt = pd.DataFrame({'cell.ptclass.color': cell_color_arr, 'gene.ptclass.color': gene_color_arr, 'ptclass': class_arr})

##### read in expression matrix of trajectory genes
exp = pd.read_csv('exp.txt', sep='\t', header=0)
exp['class'] = exp['class'].str.replace('X', '').astype(int)

##### order cells by trajectory genes
cells = [x for x in exp.columns if x not in ['styela', 'class']]
maxPT_index = exp[cells].idxmax(axis='index').values
maxPT = exp.loc[maxPT_index, ['class']]['class'].values
maxPT_genes = exp.loc[maxPT_index, ['styela']]['styela'].values
cell2pt = pd.DataFrame(data={'cells': cells, 'cell.ptclass': maxPT, 'cell.ptgene': maxPT_genes})

##### reform the expression matrix to dataframe for plot
exp = exp.melt(id_vars=['styela', 'class'], var_name='cells', value_name='exp')
exp = exp.rename(columns={'class': 'gene.ptclass'})
exp = exp.merge(cell2pt, how='left', on=['cells'])
exp.loc[exp['cells'].str.contains('.sc'), 'cell.batch'] = 'SC'
exp.loc[exp['cells'].str.contains('.C3-'), 'cell.batch'] = 'ST'
exp['cell.ptclass.cellnumber'] = exp.groupby('cell.ptclass')['cells'].transform('nunique')

##### add color attribution for ptclass
for attr in ['cell', 'gene']:
    colors = color_pt.copy().rename(columns={
        'ptclass': f'{attr}.ptclass' 
        })
    colors = colors[[f'{attr}.ptclass', f'{attr}.ptclass.color']]
    exp = exp.merge(colors, how='left', on=[f'{attr}.ptclass'])
    del colors
exp.to_csv('exp.trajectory_ordered.txt', sep='\t', header=True, index=False)

import sys
sys.exit()
def plot_cells(exp, color=None, outfile=None, cmap='plasma'):
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(tight_layout=True, figsize=(12, 5), dpi=600)
    gs = GridSpec(4, 1, figure=fig)
    ax1 = fig.add_subplot(gs[:-1, :])
    exp = exp[exp['exp'] >= 1] ##### filter out dot with low expression (to reduce figure size)
    exp = exp.sort_values(by=['cell.ptclass'])
    sns.scatterplot(data=exp, x='cells', y='exp', size='exp', sizes=(0, 50), 
            hue=color, palette=cmap, ax=ax1, linewidth=0, edgecolors='none')
    ax1.get_xaxis().set_visible(False)

    ax2 = fig.add_subplot(gs[-1, :])
    exp['tmpy'] = 1
    sns.scatterplot(data=exp, x='cell.ptclass', y='tmpy', size='cell.ptclass.cellnumber', #sizes=(0, 150), 
            hue='cell.ptclass', palette='plasma', ax=ax2, linewidth=0, edgecolors='none')
    ax2.get_yaxis().set_visible(False)

    fig.savefig(outfile)
    return

plot_cells(exp, color='gene.ptclass', outfile='pt.str.pdf', cmap='blend:#7AB,#EDA') 
plot_cells(exp, color='cell.batch', outfile='batch.pdf', cmap='tab10')




