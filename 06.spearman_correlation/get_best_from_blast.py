
import argparse as arg
import numpy as np
import pandas as pd
from kneed import KneeLocator

parser = arg.ArgumentParser( description = 'Extract homologous hit from blast results.')
parser.add_argument('-best', action='store_true', help='use only best hit.')
parser.add_argument('-out', action='store', help='output file name')
args = parser.parse_args()

class PairWise:
    def __init__(self, row):
        self.query = row.query
        self.subject = row.subject

        if self.query.startswith(('evm.', 'MSTRG.')):
            zf = self.subject
            sc = self.query
        else:
            zf = self.query
            sc = self.subject
        self.pair = (zf, sc)

    def __hash__(self):
        return hash(self.pair)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __str__(self):
        zf, sc = self.pair
        return f'{zf}\t{sc}'

    def __repr__(self):
        return f'PairWise({self.query}, {self.subject})'

    def is_reverse(self, other):
        return (self.query == other.subject) and (self.subject == other.query)

def filtering(m8, iden_thres=None, top=True):
    if iden_thres is None:
        iden_thres = 0
    
    m8 = m8[m8['identity'] >= iden_thres]
    m8 = m8.sort_values('identity', ascending=False).groupby(['query', 'subject']).head(1)
    
    if top:
        m8 = m8.groupby('query').head(1)
        m8 = m8.groupby('subject').head(1)
    
    m8['pair'] = m8.apply(PairWise, axis=1)
    m8['pair.hash'] = m8['pair'].apply(lambda x: hash(x))
    m8 = m8.sort_values('pair.hash')
    m8 = m8[m8.groupby('pair.hash')['pair.hash'].transform('count') == 2]

    m8['zf'] = m8.apply(lambda x: x['pair'].pair[0], axis=1)
    m8['sc'] = m8.apply(lambda x: x['pair'].pair[1], axis=1)
    m8 = m8.sort_values('identity', ascending=False)
    m8 = m8.drop_duplicates(subset=['sc', 'zf'], keep='first')

    #del m8['pair']
    #del m8['pair.hash']
    return m8

def locating_knee(group, on=None):
    if group.shape[0] <= 4:
        return group
    
    group = group.sort_values('identity', ascending=False).head(5)

    if on is None:
        on = ['identity', 'e-value', 'bitscore']
    
    knee_index = []
    for colname in on:
        if colname in ['e-value']:
            ascending = True
            curve = 'convex'
            direction = 'increasing'
            sensitivity = 1
        else:
            ascending = False
            curve = 'concave'
            direction = 'decreasing'
            sensitivity = 100
        data = group.sort_values(colname, ascending=ascending)
        x = list(range(data.shape[0]))
        y = data[colname].values
        kneedle = KneeLocator(x, y, S=sensitivity, 
                curve=curve, direction=direction, 
                )
        if kneedle.knee is None:
            continue
        idx = data.iloc[[kneedle.knee]].index[0]
        knee_index.append(idx)
    knee_iden = group.loc[knee_index]['identity'].max()
    group = group[group['identity'] >= knee_iden]
    return group

def rep_infering(m8, evo_copy_num=4, iden_thres=50):
    if iden_thres:
        m8 = m8[m8['identity'] >= iden_thres]
    m8 = m8.sort_values('identity', ascending=False).groupby('zf').head(1)
    m8 = m8.groupby('sc', as_index=False).apply(locating_knee)
    m8['is_single_copy'] = True
    m8.loc[m8.groupby('sc')['zf'].transform('count') != 1, 'is_single_copy'] = False
    return m8

### read in reciprocal blast results
colnames = ['query', 'subject', 'identity', 'alnlen', 'mismatch', 'gap', 
            'q.start', 'q.end', 's.start', 's.end', 'e-value', 'bitscore']
sc2zf = pd.read_csv('./reciprocal_blast/sc2zf.txt', sep='\t', header=None, names=colnames)
zf2sc = pd.read_csv('./reciprocal_blast/zf2sc.txt', sep='\t', header=None, names=colnames)
ortholog = pd.concat([sc2zf, zf2sc])

### filtering blast result by reciprocal best match
if args.best:
    ortholog = filtering(ortholog, top=True)
    ortholog['is_single_copy'] = True
    suffix = 'single'
else:
    ortholog = filtering(ortholog, top=False)
    ortholog = rep_infering(ortholog, iden_thres=None, evo_copy_num=4)
    suffix = 'multi'

### convert unique match to pd.DataFrame
ortholog = ortholog[['zf', 'sc', 'is_single_copy']].rename(columns={'zf': 'Danio_rerio.pep', 'sc': 'Styela_clava.pep'})

ortholog.to_csv(args.out, sep='\t', index=False)


