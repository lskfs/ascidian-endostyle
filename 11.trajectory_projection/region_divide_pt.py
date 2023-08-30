
import numpy as np
import pandas as pd
import anndata
import cv2

mask = cv2.imread('./data/C3-2.TLC_region_mask.tif', 0)
mask = np.isin(mask, [0])
adata = anndata.read('../data/Styela_clava.anno.h5ad') 
adata = adata[~adata.obs['Batch'].isin(['sc1', 'sc2', 'sc3'])]
adata.obs[['x', 'y']] = adata.obsm['spatial']
adata = adata[adata.obs['Batch'] == 'C3-2']
cellseg = adata.uns['seg_spatial']['C3-2']['seg']
cellseg = cellseg[mask]
unique, counts = np.unique(cellseg, return_counts=True)
half_tlc = [str(u) for u, c in zip(unique, counts) if c >= 12]

obs = adata.obs.copy()
obs.loc[(obs['short_inte_anno.202302'] == 'TLC') & obs['label'].isin(half_tlc), 'part'] = 'red'
obs.loc[(obs['short_inte_anno.202302'] == 'TLC') & (~obs['label'].isin(half_tlc)), 'part'] = 'blue'
obs.loc[(obs['part'] == 'blue') & (obs['y'] >= 480), 'part'] = 'red'
obs.loc[(obs['part'] == 'blue') & (obs['y'] <= 380), 'part'] = None
obs = obs[['x', 'y', 'part', 'Batch']]
obs = obs.fillna(value={'part': '#eeeeee'})

def get_perp_vec(u1, u2, direction=1):
    """Return the unit vector perpendicular to the vector u2-u1."""
    x1, y1 = u1
    x2, y2 = u2
    vx, vy = x2-x1, y2-y1
    v = np.linalg.norm((vx, vy))
    wx, wy = -vy/v * direction, vx/v * direction
    return wx, wy

def get_av_vec(u1, u2):
    """Return the average unit vector between u1 and u2."""
    u1x, u1y = u1
    u2x, u2y = u2
    dx, dy = u1x + u2x, u1y + u2y
    dlen = np.linalg.norm((dx,dy))
    return dx/dlen, dy/dlen

def add_ticks(x, y, ax, interval=5, label=None, direction=1, 
        tick_length=1, tick_width=0.1, tick_color='k'):
    
    left_size = len(x) % interval
    step_size = len(x) // interval
    if left_size > 0:
        step_size += 1
    else:
        pass
    idx = range(0, len(x), step_size)

    for j,i in enumerate(idx):
        if i == 0:
            # The first tick is perpendicular to the line between the
            # first two points
            tx, ty = get_perp_vec((x[0], y[0]), (x[1], y[1]),
                                  direction)
        elif i == len(x)-1:
            # The last tick is perpendicular to the line between the
            # last two points
            tx, ty = get_perp_vec((x[-2], y[-2]), (x[-1], y[-1]),
                                  direction)
        else:
            # General tick marks bisect the incoming and outgoing line
            # segments
            u1 = get_perp_vec((x[i-1], y[i-1]), (x[i], y[i]),
                              direction)
            u2 = get_perp_vec((x[i], y[i]), (x[i+1], y[i+1]),
                              direction)
            tx, ty = get_av_vec(u1, u2)
        tx, ty = tick_length * tx, tick_length * ty
        this_tick, = ax.plot((x[i],x[i]+tx), (y[i],y[i]+ty), 
                color=tick_color, lw=tick_width, zorder=0)

        if label:
            this_ticklabel = ax.text(x[i]+tx*2, y[i]+ty*2, label[j],
                    ha='center', va='center', clip_on=True, size=1)
    return

from sklearn.neighbors import NearestNeighbors
def neighbor_project(ps, line, interval=10):
    
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(line)
    distances, indices = nbrs.kneighbors(ps, n_neighbors=1, return_distance=True)
    
    ##### indices is the nearest projection of points on line
    ##### so indices represent the distance of TLC cell to region 6
    ##### scale distance to range (0-10, as the same as the ticks)
    left_size = len(line) % interval
    step_size = len(line) // interval
    if left_size > 0:
        step_size += 1
    else:
        pass
    unit = 1 / step_size
    indices = np.concatenate(indices) * unit
    return indices

batches = obs['Batch'].unique()
colnum = len(batches)
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, colnum, figsize=(4, 4), dpi=600)
for index, batch in enumerate(batches):
    data = obs[obs['Batch'] == batch]
    xs = data['x'].values
    ys = data['y'].values
    cs = data['part'].values
    
    if len(batches) > 1:
        ax = axes[index]
    else:
        ax = axes
    sc = ax.scatter(xs, ys, c=cs, s=6, linewidth=0, edgecolors='none', zorder=-10)
    
    label = [str(x) for x in range(10)]
    top = data[data['part'] == 'red'].copy()
    x_top = top['x'].values[::-1]
    y_top = top['y'].values[::-1]
    p = np.poly1d(np.polyfit(x_top, y_top, deg=7))
    xp = np.linspace(x_top.max(), x_top.min(), 100)
    yp = p(xp)
    ax.plot(xp, p(xp), '-', lw=0.5, zorder=10)
    add_ticks(xp, yp, ax, interval=10, tick_length=2, direction=-1, 
            tick_width=0.2, label=label)
    ps = top[['x', 'y']].values
    line = np.stack((xp, yp), axis=1)
    ps_dis = neighbor_project(ps, line)
    top['dis'] = ps_dis
    
    bottom = data[data['part'] == 'blue'].copy()
    part1 = bottom[bottom['y'] >= 435].copy()
    x_part1 = part1['x'].values
    y_part1 = part1['y'].values
    p1 = np.poly1d(np.polyfit(x_part1, y_part1, deg=5))
    xp1 = np.linspace(x_part1.max(), x_part1.min(), 100)
    yp1 = p1(xp1)
    #ax.plot(x_part1, y_part1, '.', xp1, p1(xp1), '-')
    part2 = bottom[bottom['y'] < 435].copy()
    x_part2 = part2['x'].values
    y_part2 = part2['y'].values
    p2 = np.poly1d(np.polyfit(x_part2, y_part2, deg=5))
    xp2 = np.linspace(x_part2.min(), x_part2.max(), 100)
    yp2 = p2(xp2)
    #ax.plot(x_part2, y_part2, '.', xp2, p2(xp2), '-')
    yp_middle = np.linspace(yp1[-1], yp2[0], 10)[1:-1]
    xp_middle = np.linspace(xp1[-1], xp2[0], 10)[1:-1]
    
    xp = np.concatenate([xp1, xp_middle, xp2])
    yp  = np.concatenate([yp1, yp_middle, yp2])
    ax.plot(xp, yp, '-', lw=0.5, zorder=10)
    add_ticks(xp, yp, ax, interval=10, tick_length=2, direction=-1, 
            tick_width=0.2, label=label)
    ps = bottom[['x', 'y']].values
    line = np.stack((xp, yp), axis=1)
    ps_dis = neighbor_project(ps, line)
    bottom['dis'] = ps_dis

    ax.set_aspect('equal')

plt.tight_layout()
fig.savefig('thyroid_region_division.pdf')

results = pd.concat([top, bottom])
pt_meta = pd.read_csv('exp.trajectory_ordered.txt', sep='\t', header=0)
pt_meta = pt_meta[['cells', 'cell.ptclass']].drop_duplicates()
pt_meta = pt_meta.rename(columns={'cells': 'index'}).set_index('index')
results = results.merge(pt_meta, how='left', left_index=True, right_index=True)

#hist_bins = [i * 2 for i in range(6)]
hist_bins = [i * 3 for i in range(4)]
results['bins'] = pd.cut(results['dis'], bins=hist_bins)
results = results[['part', 'bins', 'cell.ptclass']]

import seaborn as sns
fig, ax = plt.subplots(figsize=(4, 4), dpi=600)
sns.boxplot(data=results, x='bins', y='cell.ptclass', hue='part', ax=ax)
fig.savefig('thyroid_region_boxplot.pdf')



