library(dplyr)
library(ggplot2)

####################################################################################################
Z_500kb_1Mb_co <- read.table('Z_compartment_co_dsb_summary.tsv', sep='\t', header=TRUE)

table(Z_500kb_1Mb_co$A_region_groups)

#A1  A2  A3   B 
#76  73  74 377
###################################################################################################

# co_davis_dsb_ratio
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=co_davis_dsb_ratio, fill=A_region_groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("red", "orange", "orange", "dodgerblue"))+
  theme_bw() + 
  ylim(0,1)
  #+geom_jitter(shape=16, position=position_jitter(0.2))

ggsave('Fig4g_hinch_co_davis_dsb_ratio_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)
wilcox.test(Z_500kb_1Mb_co$co_davis_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$co_davis_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A2'], alternative = 'greater')
#p-value = 0.02754

#co_smagulova_dsb_ratio
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=co_smagulova_dsb_ratio, fill=A_region_groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("red", "orange", "orange", "dodgerblue"))+
  theme_bw() + 
  ylim(0,1)
#+geom_jitter(shape=16, position=position_jitter(0.2))

ggsave('Fig4h_hinch_co_smagulova_dsb_ratio_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)
wilcox.test(Z_500kb_1Mb_co$co_smagulova_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$co_smagulova_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A2'], alternative = 'greater')
#p-value = 0.01854


#hinch_co_density
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=hinch_co_density, fill=A_region_groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("red", "orange", "orange", "dodgerblue"))+  
  theme_bw()+  
  ylim(0, 15)
  #geom_jitter(shape=16, position=position_jitter(0.2))
ggsave('Fig4f_hinch_co_density_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)

wilcox.test(Z_500kb_1Mb_co$hinch_co_density[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$hinch_co_density[Z_500kb_1Mb_co$A_region_groups=='A3'], alternative = 'greater')
#p-value = 0.1054

#davis_dsb_density
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=davis_dsb_density, fill=A_region_groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("red", "orange", "orange", "dodgerblue"))+  
  theme_bw()+  
  ylim(0, 15) 
  #geom_jitter(shape=16, position=position_jitter(0.2))
ggsave('Fig4d_davis_dsb_density_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)

wilcox.test(Z_500kb_1Mb_co$davis_dsb_density[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$davis_dsb_density[Z_500kb_1Mb_co$A_region_groups=='A3'], alternative = 'greater')

#p-value = 0.1794


#smagulova_dsb_density
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=smagulova_dsb_density, fill=A_region_groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("red", "orange", "orange", "dodgerblue"))+  
  theme_bw()+  
  ylim(0, 15) 
  # geom_jitter(shape=16, position=position_jitter(0.2))
ggsave('Fig4e_smagulova_dsb_density_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)

wilcox.test(Z_500kb_1Mb_co$smagulova_dsb_density[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$smagulova_dsb_density[Z_500kb_1Mb_co$A_region_groups=='A3'], alternative = 'greater')
#p-value = 0.06014