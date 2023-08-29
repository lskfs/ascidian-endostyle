
import argparse as arg
import pandas as pd

parser = arg.ArgumentParser( description = 'Extract TF matches from homologous.')
parser.add_argument('-input', action='store', help='TF list file')
parser.add_argument('-out', action='store', help='output file name')
args = parser.parse_args()

### read in reciprocal matches
ortholog = pd.read_csv(args.input, sep='\t', header=0)

### read in TF list
tf = pd.read_csv('zebrafish.TF.txt', sep='\t', header=0)

### extract zebrafish TF record
tf_orth = ortholog.merge(tf, how='left', left_on='Danio_rerio.pep', right_on='Symbol')
tf_orth = tf_orth[tf_orth['Symbol'].notna()]

### write out
tf_orth.to_csv(args.out, sep='\t', index=False)


