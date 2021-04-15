import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocess as mp
import bioframe
import cooler
import itertools
import click
import cooltools
import cooltools.eigdecomp
import cooltools.expected
import cooltools.saddle
from dask.distributed import Client, LocalCluster
from scipy.linalg import toeplitz

stages=["sertoli", "spg", "PreL", "L", "Z", "P", "D","MII"]

for stage in stages:
  coolfile= ''.join([stage, '_500kb.cool'])
  c = cooler.Cooler(coolfile)
  expected=''.join([stage, '_500kb_cis_expected.tsv'])
  expected1=pd.read_csv(expected, sep='\t')
  expected2={k: x.values for k, x in expected1.groupby('chrom')['balanced.avg']}
  for chr in c.chromnames:
    obs_mat1 = c.matrix().fetch(chr)
    exp_mat1 = toeplitz(expected2[chr][:obs_mat1.shape[0]])
    output1 =obs_mat1/exp_mat1
    filename=''.join([chr, '_cis_obs_over_whole_genome_exp.matrix'])
    np.savetxt(filename, output1, delimiter="\t")

for stage in stages:
    name =''.join(['chr1_',stage,'_cis_obs_over_whole_genome_exp.matrix'])
    mat=pd.read_csv(name, sep='\t', header=None)
    mat1=mat.iloc[1789:1949,1789:1949]
    row_chrom= 'chr1'
    col_chrom= 'chr1'
    scale='log2'
    out=''.join(['chr1_',stage,'_log2_obs_over_whole_genome_exp_50kb_region2.pdf'])
    dpi= 300
    colormap='bwr'
    zmin=-1
    zmax=1
    plt.figure(figsize=(10,10))
    plt.gcf().canvas.set_window_title("Contact matrix".format())
    plt.title("")
    plt.imshow(np.log2(mat1), interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
    plt.ylabel("{} coordinate".format(row_chrom))
    plt.xlabel("{} coordinate".format(col_chrom))
    cb = plt.colorbar()
    cb.set_label({"linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
     }[scale])
    plt.savefig(out, dpi=dpi, format='pdf')

