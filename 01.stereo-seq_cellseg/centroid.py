
import os
import sys

from pathlib import Path
import numpy as np
import pandas as pd
from skimage.measure import regionprops, regionprops_table

def calculate_centroid(label_image):

    label_image = np.loadtxt(label_image, dtype=int)
    zero_mask = np.isin(label_image, [0])
    label_image[zero_mask] = label_image.max() + 1

    unique, unique_indices, unique_inverse = np.unique(label_image, return_index=True, return_inverse=True)

    relabed_image = unique_inverse.reshape(label_image.shape)
    min_mask = np.isin(relabed_image, [0])
    relabed_image[min_mask] = relabed_image.max() + 1
    relabed_image[zero_mask] = 0

    unique[-1] = 0
    label_map = pd.DataFrame(dict(orig_label=unique, label=unique_inverse[unique_indices]))

    props = regionprops_table(
            relabed_image,
            properties=('label', 'centroid')
            )

    df = pd.DataFrame(props)
    df = df.merge(label_map, how='left', on='label')
    df = df[['orig_label', 'centroid-1', 'centroid-0']]
    df = df.rename(
            columns={
                'orig_label':'label', 
                'centroid-0':'y',
                'centroid-1':'x'
                }
            )
    df['label'] = df['label'].astype(int)
    return df

if __name__ == '__main__':
    
    label_image = sys.argv[1]
    prefix = sys.argv[2]

    data = calculate_centroid(label_image)
    data.to_csv(f'{prefix}.centroid.csv', sep='\t', index=False, header=True)


