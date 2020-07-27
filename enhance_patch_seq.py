from __future__ import print_function, division
import pandas as pd
import numpy as np
import enhance
import os
fpath_expanded = os.path.expanduser('/Users/tangxin/Dropbox (Harvard University)/Harvard_research&course/Sig_seq/sig_seq_translation/Enhance_imputation/enhance/data/pbmc-4k_expression.tsv.gz')

matrix = pd.read_csv(fpath_expanded, index_col=0, sep='\t')
gene_csv='/Users/tangxin/Dropbox (Harvard University)/Harvard_research&course/Sig_seq/sig_seq_translation/NN_based_translation/data/scala2019/patch-seq-counts.csv'
data = pd.read_csv(gene_csv, index_col=0, sep='\t')
counts = data.values[:, 1:].transpose().astype(float)
genes = data.values[:, 0]
cells = np.array(data.columns[1:])
print('Number of cells: {}\nNumber of genes: {}'.format(counts.shape[0], counts.shape[1]))
