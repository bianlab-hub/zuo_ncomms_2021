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

#calculate expected interaction at each distance on all autosomes for A/B compartment domains larger than 1Mb 

stages=["sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII"]


def cal_region_expected(regions, c, outfile, ignorediag=2):

    with multiprocess.Pool(20) as pool:
        expected = cooltools.expected.diagsum(c,regions,transforms={'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']},map=pool.map, ignore_diags=ignorediag)

    expected_df = pd.concat([exp.reset_index().assign(chrom=reg[0], start=reg[1], end=reg[2]) for reg, exp in expected.items()])

    expected_df = expected_df.groupby(('chrom','diag')).aggregate({'n_valid':'sum','count.sum':'sum','balanced.sum':'sum'}).reset_index()

    expected_df['balanced.avg'] = expected_df['balanced.sum'] / expected_df['n_valid']

    expected_df.to_csv(outfile,sep='\t', index=False, na_rep='nan')
    
def combine_expected(exp_tsv, sumfile, combine_by='diag'):
    expected_df=pd.read_csv(exp_tsv, sep='\t')
    expected_df = expected_df.groupby(combine_by).aggregate({'n_valid':'sum','count.sum':'sum','balanced.sum':'sum'}).reset_index()
    expected_df['balanced.avg'] = expected_df['balanced.sum'] / expected_df['n_valid']
    expected_df.to_csv(sumfile,sep='\t', index=False, na_rep='nan')
    
    
    
for stage in stages:
    coolfile='/data05/combined_spermatocytes/' + stage + '/fragment_db/' + stage + '_combined_10kb.cool'
    clr = cooler.Cooler(coolfile)
    
    
    a_bed= '/data05/combined_spermatocytes/A_B_compartment_expected/' +stage +'_50kb_A_compartment_domains_larger_than_1Mb.bed'
    b_bed= '/data05/combined_spermatocytes/A_B_compartment_expected/' +stage +'_50kb_B_compartment_domains_larger_than_1Mb.bed'

    a_domain=pd.read_csv(a_bed, sep='\t', header=None, usecols=[0,1,2])
    a_domain=a_domain[a_domain[0]!='chrX']

    b_domain=pd.read_csv(b_bed, sep='\t', header=None, usecols=[0,1,2])
    b_domain=b_domain[b_domain[0]!='chrX']

    a_regions= list(a_domain.itertuples(index=False, name=None))
    b_regions= list(b_domain.itertuples(index=False, name=None))

    
    a_exp_file= stage + '_A_compartment_expected_10kb.tsv'
    b_exp_file= stage + '_B_compartment_expected_10kb.tsv'
    
    cal_region_expected(a_regions, clr, a_exp_file, 0)
    cal_region_expected(b_regions, clr, b_exp_file, 0)

    a_exp_combined= stage + '_combined_autosomes_A_compartment_expected_10kb.tsv'
    b_exp_combined= stage + '_combined_autosomes_B_compartment_expected_10kb.tsv'
    
    combine_expected (a_exp_file, a_exp_combined, 'diag')
    combine_expected (b_exp_file, b_exp_combined, 'diag')
    
