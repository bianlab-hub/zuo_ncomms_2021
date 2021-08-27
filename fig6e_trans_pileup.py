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

coolfile='Sertoli_50kb.cool'
transexpfile='Sertoli_50kb_trans_expected.tsv'
stage='Sertoli_50kb'

coolfile='spg_50kb.cool'
transexpfile='spg_50kb_trans_expected.tsv'
stage='spg_50kb'

coolfile='L_50kb.cool'
transexpfile='L_50kb_trans_expected.tsv'
stage='L_50kb'

coolfile='Z_50kb.cool'
transexpfile='Z_50kb_trans_expected.tsv'
stage='Z_50kb'

coolfile='P_50kb.cool'
transexpfile='P_50kb_trans_expected.tsv'
stage='P_50kb'

coolfile='D_50kb.cool'
transexpfile='D_50kb_trans_expected.tsv'
stage='D_50kb'

coolfile='MII_50kb.cool'
transexpfile='MII_50kb_trans_expected.tsv'
stage='MII_50kb'

pileupshape=500
c = cooler.Cooler(coolfile)

#generate chrom pairs for all autosomes
chroms=c.chromnames
chroms.remove('chrX')
chroms.remove('chrM')
chroms.remove('chrY')

chrompairs=list(itertools.combinations(chroms,2))

##############################################################################################
#make trans obs pileup
pileup=np.zeros((pileupshape, pileupshape))
n=0
for pair in chrompairs:
    t1 = c.matrix().fetch(pair[0], pair[1])
    t1[np.isnan(t1)]=0
    newt1 = scipy.ndimage.interpolation.zoom(input=t1, zoom=(pileupshape/(t1.shape[0]), pileupshape/(t1.shape[1])), order = 2)    
    pileup+=newt1
    n+=1     

# seperate
obsaverage = pileup/n
oepileup=''.join([stage, '_average_trans_obs_500x500.matrix'])
np.savetxt(oepileup, obsaverage, delimiter="\t")

#plot heatmap
mat=np.log10(obsaverage)
row_chrom='chrA'
col_chrom='chrB'
scale='log10'
out=''.join([stage, '_average_trans_obs_500x500.pdf'])
dpi= 300
colormap='YlOrRd'
zmin=np.nanquantile(mat,0.05)
zmax=np.nanquantile(mat,0.95)

plt.figure(figsize=(10,10))
plt.gcf().canvas.set_window_title("Contact matrix".format())
plt.title("")
plt.imshow(mat, interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
plt.ylabel("{} coordinate".format(row_chrom))
plt.xlabel("{} coordinate".format(col_chrom))
cb = plt.colorbar()
cb.set_label({"linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
    }[scale])
plt.savefig(out, dpi=dpi, format='pdf')   

###############################################################################################
#load trans expected values
trans_expected=pd.read_csv(transexpfile, sep='\t')
trans_expected_2 = {k: x.values for k, x in trans_expected.groupby(['chrom1', 'chrom2'])['balanced.avg']}

#calculate average trans obs/exp
pileup=np.zeros((pileupshape, pileupshape))
n=0
for pair in chrompairs:
    trans_obs_mat = c.matrix().fetch(pair[0], pair[1])
    t1=trans_obs_mat/ trans_expected_2 [pair]
    t1[np.isnan(t1)]=0
    newt1 = scipy.ndimage.interpolation.zoom(input=t1, zoom=(pileupshape/(t1.shape[0]), pileupshape/(t1.shape[1])), order = 2)    
    pileup+=newt1
    n+=1    


# seperate
average=pileup/n
oepileup=''.join([stage, '_average_trans_obs_exp_500x500.matrix'])
np.savetxt(oepileup, average, delimiter="\t")

#plot heatmap
mat=np.log2(average)
row_chrom='chrA'
col_chrom='chrB'
scale='log2'
out=''.join([stage, '_average_trans_obs_exp_500x500.pdf'])
dpi= 300
colormap='bwr'
zmin=-1
zmax=1

plt.figure(figsize=(10,10))
plt.gcf().canvas.set_window_title("Contact matrix".format())
plt.title("")
plt.imshow(mat, interpolation="none",vmin=zmin,vmax=zmax, cmap=colormap)
plt.ylabel("{} coordinate".format(row_chrom))
plt.xlabel("{} coordinate".format(col_chrom))
cb = plt.colorbar()
cb.set_label({"linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
    }[scale])

plt.savefig(out, dpi=dpi, format='pdf')   
