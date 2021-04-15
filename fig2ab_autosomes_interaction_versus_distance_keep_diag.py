import numpy as np
import pandas as pd

import pybedtools

import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt

import cooler
import cooltools
import cooltools.expected
import cooltools.snipping
import multiprocess

stages=['L','Z','P','D','PreL','Sertoli','spg','MII']
#stages=['L1','Z1','P1','D1','PreL1','Sertoli1','spg1','MII1','L2','Z2','P2','D2','PreL2','Sertoli2','spg2','MII2']

for stage in stages:
    coolfile='/data01/meiosis_Y/hiclib/10kb_cool/' + stage + '_10kb.cool'
    c = cooler.Cooler(coolfile)
    autosomes=c.chromnames
    autosomes.remove('chrX')
    autosomes.remove('chrY')
    
    cis_count_sum={chr: (np.nansum(c.matrix(balance=False).fetch(chr)) + np.nansum(c.matrix(balance=False).fetch(chr).diagonal()))/2 for chr in autosomes }
    total_auto_cis = sum(cis_count_sum.values())
    
    total_x_cis = (np.nansum(c.matrix(balance=False).fetch('chrX')) + np.nansum(c.matrix(balance=False).fetch('chrX').diagonal()))/2
    cisexpfile='/data01/meiosis_Y/hiclib/10kb_cool/expected/' + stage + '_10kb_cis_expected.tsv'
    expected_df=pd.read_csv(cisexpfile, sep='\t')
    
    x_df=expected_df[expected_df['chrom']=='chrX']
    x_df = x_df.groupby(('diag')).aggregate({'n_valid':'sum','count.sum':'sum','balanced.sum':'sum'}).reset_index()
    x_df['distances']=x_df['diag']*10000
    x_df['frequency']=x_df['count.sum']/total_x_cis
    x_outfile=stage+'_x_interaction_versus_distance_keep_diag.tsv'
    x_df.to_csv(x_outfile,sep='\t', index=False, na_rep='nan')
    
    
    autosome_df=expected_df[expected_df['chrom'].isin(autosomes)]
    autosome_df = autosome_df.groupby(('diag')).aggregate({'n_valid':'sum','count.sum':'sum','balanced.sum':'sum'}).reset_index()    
    autosome_df['distances']=autosome_df['diag']*10000
    autosome_df['frequency']=autosome_df['count.sum']/total_auto_cis        
    autosome_outfile=stage+'_autosomes_interaction_versus_distance_keep_diag.tsv'
    autosome_df.to_csv(autosome_outfile,sep='\t', index=False, na_rep='nan')  