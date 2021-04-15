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
import scipy.ndimage

stages=["sertoli", "spg", "PreL", "L", "Z", "P", "D", "Sun1"]

pileupshape=500
plotmin=-1
plotmax=1

for stage in stages:
    coolfile=''.join([stage, '_combined_500kb.cool'])
    transexpfile=''.join([stage, '_combined_500kb_trans_expected.tsv'])
    c = cooler.Cooler(coolfile)
    chroms=c.chromnames
    chroms.remove('chrM')
    chroms.remove('chrX')
    chroms.remove('chrY')
    chrompairs=list(itertools.combinations(chroms,2))
    trans_expected=pd.read_csv(transexpfile, sep='\t')
    trans_expected_2 = {k: x.values for k, x in trans_expected.groupby(['chrom1', 'chrom2'])['balanced.avg']}
    pileup=np.zeros((pileupshape, pileupshape))
    n=0
    for pair in chrompairs:
        trans_obs_mat = c.matrix().fetch(pair[0], pair[1])
        t1=trans_obs_mat/ trans_expected_2 [pair]
        t1[np.isnan(t1)]=0
        newt1 = scipy.ndimage.interpolation.zoom(input=t1, zoom=(pileupshape/(t1.shape[0]), pileupshape/(t1.shape[1])), order = 2)    
        pileup+=newt1
        n+=1  
        
    average=pileup/n
    oepileup=''.join([stage, '_average_trans_oe_500x500.matrix'])
    np.savetxt(oepileup, average, delimiter="\t")

#plot heatmap
    mat=np.log2(average)
    row_chrom='chrA'
    col_chrom='chrB'
    scale='log2'
    out=''.join([stage, '_average_trans_obs_exp_500x500_',str(plotmin),'_',str(plotmax),'.pdf'])
    dpi= 300
    colormap='bwr'
    zmin=plotmin
    zmax=plotmax

    plt.figure(figsize=(10,10))
    plt.gcf().canvas.set_window_title("Contact matrix".format())
    plt.title(stage)
    plt.imshow(mat, interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
    plt.ylabel("{} coordinate".format(row_chrom))
    plt.xlabel("{} coordinate".format(col_chrom))
    cb = plt.colorbar()
    cb.set_label({"linear": "relative contact frequency",
        "log2": "log 2 ( relative contact frequency )",
        "log10": "log 10 ( relative contact frequency )",
        }[scale])
    plt.savefig(out, dpi=dpi, format='pdf')   

