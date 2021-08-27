import numpy as np
import pandas as pd
import pybedtools
import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt
import itertools
import cooler
import cooltools
import cooltools.expected
import cooltools.snipping
import multiprocess


halfwindow = 500000
binsize = 50000
hw = int(halfwindow/binsize)
    
coolfile='/data03/haplo_hic/nsmb_mousehic/GSE122622_pachynema_interhomolog_50kb.cool'
c = cooler.Cooler(coolfile)

chroms=c.chromnames
chroms.remove('M')
chroms.remove('X')
chroms.remove('Y')

chr=[]
start=[]
end=[]
score=[]

for chrom in chroms:
    mat = c.matrix().fetch(chrom)
    avg = np.nanmean(mat)
    l =mat.shape [0]
    chr_1= [chrom] * l
    start_1 = [i * binsize for i in list(range(l))]
    end_1 = [i * binsize for i in list(range(1, l+1))]
    
    score_1 = []
    for i in list(range(l)):
        square = mat[max(0, i-hw) : min(i+hw, l), max(0, i-hw) : min(i+hw, l)]
        ratio = np.nanmean(square) / avg
        score_1.append(ratio)
        
    chr= chr + chr_1
    start = start + start_1
    end = end + end_1
    score = score + score_1
    
interhomolog = pd.DataFrame( {'chr': chr, 'start': start,'end': end,'alignmentscore': score})

interhomolog.to_csv('Z_interhomolog_alignment_ratio_halfwindow_500kb.tsv',sep='\t', index=True, na_rep='nan')



