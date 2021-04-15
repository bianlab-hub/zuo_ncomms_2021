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


stages=["sertoli", "PreL", "L", "Z", "P", "D"]



def cal_each_region_expected(regions, c, outfile, ignorediag=2):

    with multiprocess.Pool(20) as pool:
        expected = cooltools.expected.diagsum(c,regions,transforms={'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']},map=pool.map, ignore_diags=ignorediag)

    expected_df = pd.concat([exp.reset_index().assign(chrom=reg[0], start=reg[1], end=reg[2]) for reg, exp in expected.items()])

    expected_df.to_csv(outfile,sep='\t', index=False, na_rep='nan')
    

    
    
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

    
    a_exp_file= stage + '_A_compartment_each_region_expected_10kb.tsv'
    b_exp_file= stage + '_B_compartment_each_region_expected_10kb.tsv'
    
    cal_each_region_expected(a_regions, clr, a_exp_file)
    cal_each_region_expected(b_regions, clr, b_exp_file)
    