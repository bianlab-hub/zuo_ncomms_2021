from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-white')

import multiprocess as mp
import numpy as np
import pandas as pd
import bioframe
import cooltools
import cooler
from cooltools.eigdecomp import cooler_cis_eig

mm10 = bioframe.fetch_chromsizes('mm10')
chromsizes = bioframe.fetch_chromsizes('mm10')
chromosomes = list(chromsizes.index)

binsize=10000
bins = cooler.binnify(mm10, binsize)
fasta_records = bioframe.load_fasta('/data05/genomes/mm10_20chr.fa')
bins['GC'] = bioframe.tools.frac_gc(bins, fasta_records)
bins.head()

import fnmatch
import os

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*_10kb.cool'):
        clr=cooler.Cooler(file)
        cond=file.split('.')[0]
        lam, eigs = cooler_cis_eig(clr, bins,n_eigs=3, phasing_track_col='GC', sort_metric='var_explained')
        # Save text files
        lam.to_csv(f'./{cond}.eigs.cis.lam.txt', sep='\t')
        eigs.to_csv(f'./{cond}.eigs.cis.vecs.txt', sep='\t', index=False)
