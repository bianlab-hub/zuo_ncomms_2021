# fetch the sites and interaction frequency
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

def cis_diag_regions_stats (c, regions, outfile):
    c1, s1, n, cmed, bmed, csum, bsum, baverage= ([] for i in range(10))
    for chr in chrs:
        exec ("bsum_%s = ([] for i in range(10)) "%chr)
        exec ("bmed_%s = ([] for i in range(10)) "%chr)
    for region in regions:
        c1.append(region.split(':')[0])
        s1.append(region)      
        mat_bal=c.matrix(balance=True).fetch(region)
        mat_raw=c.matrix(balance=False).fetch(region)
        n_valid= np.count_nonzero(~np.isnan(mat_bal))
        balance_sum = np.nansum(mat_bal)
        n.append(n_valid)
        cmed.append(np.nanmedian(mat_raw))
        bmed.append(np.nanmedian(mat_bal))
        csum.append(np.nansum(mat_raw))
        bsum.append(balance_sum)
        baverage.append(balance_sum/n_valid)
        chrs=c.chromnames
        chrs.remove('chrM')
        chrs.remove('chrY')
        for chr in chrs:
            mat_chr_bal=c.matrix(balance=True).fetch(region,chr)
            chr_balance_sum = np.nansum(mat_chr_bal)
            globals()[''.join("bsum_" + chr)].append(chr_balance_sum)
            globals()[''.join("bmed_", + chr)].append(np.nanmedian(mat_chr_bal))
    df = pd.DataFrame()
    df['chrom'] = c1
    df['region'] = s1
    df['n_valid'] = n
    df['count.sum'] = csum
    df['count.median'] = cmed
    df['balanced.sum'] = bsum
    df['balanced.median'] = bmed
    df['balanced.average'] = baverage
    for chr in chrs:
        df[''.join("balanced.sum_" + chr)] = globals()[''.join("bsum_" + chr)]
        df[''.join("balanced.median_" + chr)] = globals()[''.join("bmed_" + chr)]
    df.to_csv(outfile,sep='\t', index=False, na_rep='nan')


co = open('CO_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_co_1Mb_chr.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_co_1Mb_chr.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_co_1Mb_chr.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_co_1Mb_chr.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_co_1Mb_chr.tsv')

co = open('NCO_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_nco_1Mb_chr.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_nco_1Mb_chr2.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_nco_1Mb_chr2.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_nco_1Mb_chr.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_nco_1Mb_chr.tsv')

co = open('shuffled_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_shuffled_1Mb_chr.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_shuffled_1Mb_chr.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_shuffled_1Mb_chr.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_shuffled_1Mb_chr.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_shuffled_1Mb_chr.tsv')

########
def cis_diag_regions_stats (c, regions, outfile):
    c1, bmed, bsum = ([] for i in range(3))
    for region in regions:
        c1 = region.split(':')[0]
        obs_mat = c.matrix(balance=True).fetch(region)
        mat = obs_mat[51:151,101]
        balance_sum = np.nansum(mat)
        bmed.append(np.nanmedian(mat))
        bsum.append(balance_sum)
    df = pd.DataFrame()
    df['balanced.sum'] = bsum
    df['balanced.median'] = bmed
    df.to_csv(outfile,sep='\t', index=False, na_rep='nan')


co = open('CO_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_co_500kb_cen2.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_co_500kb_cen2.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_co_500kb_cen2.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_co_500kb_cen2.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_co_500kb_cen2.tsv')

co = open('NCO_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_nco_500kb_cen2.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_nco_500kb_cen2.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_nco_500kb_cen2.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_nco_500kb_cen2.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_nco_500kb_cen2.tsv')

co = open('shuffled_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_shuffled_500kb_cen2.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_shuffled_500kb_cen2.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_shuffled_500kb_cen2.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_shuffled_500kb_cen2.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_shuffled_500kb_cen2.tsv')

co = open('dsb_co_dmc1_300_10kb_bin_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_co_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_co_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_co_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_co_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_co_dmc1_500kb_cen2.tsv')

co = open('dsb_nco_dmc1_300_10kb_bin_plus_minus_1Mb_merge.bed','r')
co = co.readlines()
c = cooler.Cooler('L_10kb.cool') 
cis_diag_regions_stats (c, co, 'L_nco_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('Z_10kb.cool') 
cis_diag_regions_stats (c, co, 'Z_nco_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('PreL_10kb.cool') 
cis_diag_regions_stats (c, co, 'PreL_nco_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('Sertoli_10kb.cool') 
cis_diag_regions_stats (c, co, 'Sertoli_nco_dmc1_500kb_cen2.tsv')
c = cooler.Cooler('P_10kb.cool') 
cis_diag_regions_stats (c, co, 'P_nco_dmc1_500kb_cen2.tsv')

