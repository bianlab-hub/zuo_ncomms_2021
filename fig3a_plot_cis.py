# all by all heatmap

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

stages=['L','Z','P','D','PreL','Sertoli','spg','MII']

# cis 
for stage in stages:
  coolfile= ''.join([stage, '_50kb.cool'])
  c = cooler.Cooler(coolfile)
  for chr in chroms:
    mat = c.matrix().fetch(chr)
    row_chrom= chr
    col_chrom= chr
    scale='log2'
    out=''.join([chr,'_', stage,'_log2_obs_50kb.pdf'])
    dpi= 300
    colormap='YlOrRd'
    zmin= -5
    zmax= -1
    plt.figure(figsize=(10,10))
    plt.gcf().canvas.set_window_title("Contact matrix".format())
    plt.title("")
    plt.imshow(np.log10(mat), interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
    plt.ylabel("{} coordinate".format(row_chrom))
    plt.xlabel("{} coordinate".format(col_chrom))
    cb = plt.colorbar()
    cb.set_label({"linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
     }[scale])
    plt.savefig(out, dpi=dpi, format='pdf') 