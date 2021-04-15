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

#Fetch the interactions from trans heatmap
stages=["sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm"]
dist_to_end=5000000
     

def trans_region_stats (c, regions, outfile):
    c1, s1, e1, c2, s2, e2, n, cmed, bmed, csum, bsum, baverage = ([] for i in range(12))
    for transregion in regions:
        c1.append(transregion[0][0])
        s1.append(transregion[0][1])
        e1.append(transregion[0][2])
        c2.append(transregion[1][0])
        s2.append(transregion[1][1])
        e2.append(transregion[1][2])        
       
        mat_bal=c.matrix(balance=True).fetch(transregion[0], transregion[1])
        mat_raw=c.matrix(balance=False).fetch(transregion[0], transregion[1])

        n_valid= np.count_nonzero(~np.isnan(mat_bal))
        balance_sum = np.nansum(mat_bal)
        n.append(n_valid)
        cmed.append(np.nanmedian(mat_raw))
        bmed.append(np.nanmedian(mat_bal))
        csum.append(np.nansum(mat_raw))
        bsum.append(balance_sum)
        baverage.append(balance_sum/n_valid)
        
    df = pd.DataFrame()
    df['chrom1'] = c1
    df['start1'] = s1
    df['end1'] = e1
    df['chrom2'] = c2
    df['start2'] = s2
    df['end2'] = e2
    df['n_valid'] = n
    df['count.sum'] = csum
    df['count.median'] = cmed
    df['balanced.sum'] = bsum
    df['balanced.median'] = bmed
    df['balanced.average'] = baverage
    
    df.to_csv(outfile,sep='\t', index=False, na_rep='nan')

    
#    
def cis_offdiag_regions_stats (c, regions, outfile):
    c1, s1, e1, c2, s2, e2, n, cmed, bmed, csum, bsum, baverage = ([] for i in range(12))
    for region in regions:
        c1.append(region[0][0])
        s1.append(region[0][1])
        e1.append(region[0][2])
        c2.append(region[1][0])
        s2.append(region[1][1])
        e2.append(region[1][2])        
       
        mat_bal=c.matrix(balance=True).fetch(region[0], region[1])
        mat_raw=c.matrix(balance=False).fetch(region[0], region[1])

        n_valid= np.count_nonzero(~np.isnan(mat_bal))
        balance_sum = np.nansum(mat_bal)
        n.append(n_valid)
        cmed.append(np.nanmedian(mat_raw))
        bmed.append(np.nanmedian(mat_bal))
        csum.append(np.nansum(mat_raw))
        bsum.append(balance_sum)
        baverage.append(balance_sum/n_valid)
        
    df = pd.DataFrame()
    df['chrom1'] = c1
    df['start1'] = s1
    df['end1'] = e1
    df['chrom2'] = c2
    df['start2'] = s2
    df['end2'] = e2
    df['n_valid'] = n
    df['count.sum'] = csum
    df['count.median'] = cmed
    df['balanced.sum'] = bsum
    df['balanced.median'] = bmed
    df['balanced.average'] = baverage
    
    df.to_csv(outfile,sep='\t', index=False, na_rep='nan')
    
    
def cis_diag_regions_stats (c, regions, outfile):
    c1, s1, e1, n, cmed, bmed, csum, bsum, baverage = ([] for i in range(9))
    for region in regions:
        c1.append(region[0])
        s1.append(region[1])
        e1.append(region[2])        
       
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
        
    df = pd.DataFrame()
    df['chrom'] = c1
    df['start'] = s1
    df['end'] = e1
    df['n_valid'] = n
    df['count.sum'] = csum
    df['count.median'] = cmed
    df['balanced.sum'] = bsum
    df['balanced.median'] = bmed
    df['balanced.average'] = baverage
    
    df.to_csv(outfile,sep='\t', index=False, na_rep='nan')

for stage in stages:

    coolfile='/data05/combined_spermatocytes/' + stage + '/fragment_db/' + stage + '_combined_10kb.cool'

    c = cooler.Cooler(coolfile)  
    chrs=c.chromnames
    chrs.remove('chrM')
    chrs.remove('chrX')
    chrs.remove('chrY')

    chroms=pd.DataFrame({'chrom': chrs, 'start': [0]*len(c.chromsizes[chrs]), 'end': list(c.chromsizes[chrs])})

    cen_end= pd.DataFrame({'chrom': chrs, 'start': [0]*len(chrs), 'end': [dist_to_end]*len(chrs)})
    tel_end= pd.DataFrame({'chrom': chrs, 'start': list(c.chromsizes[chrs]-dist_to_end), 'end': list(c.chromsizes[chrs])})

    cen_regions= list(cen_end.itertuples(index=False, name=None))
    tel_regions= list(tel_end.itertuples(index=False, name=None))

    trans_cen_cen = list(itertools.combinations(cen_regions,2))
    trans_tel_tel = list(itertools.combinations(tel_regions,2))

    cen_tel = list(itertools.product(cen_regions, tel_regions))
    trans_cen_tel = [ x for x in cen_tel if x[0][0]!= x[1][0]] 
    cis_cen_tel = [ x for x in cen_tel if x[0][0]== x[1][0]]

    cis_diag_regions_stats (c, cen_regions, stage+'_cis_cen_'+ str(dist_to_end)+ '.tsv')
    cis_diag_regions_stats (c, tel_regions, stage+'_cis_tel_'+ str(dist_to_end)+ '.tsv')
    cis_offdiag_regions_stats (c, cis_cen_tel, stage+'_cis_cen_tel_'+ str(dist_to_end)+ '.tsv') 
    trans_region_stats (c, trans_cen_cen, stage+'_trans_cen_cen_'+ str(dist_to_end)+ '.tsv')
    trans_region_stats (c, trans_tel_tel, stage+'_trans_tel_tel_'+ str(dist_to_end)+ '.tsv')   
    trans_region_stats (c, trans_cen_tel, stage+'_trans_cen_tel_'+ str(dist_to_end)+ '.tsv')   
   
