library(dplyr)
library(ggplot2)

#prep data

Z_500kb_1Mb <- read.table('Z_compartment_region_500kb_1Mb_interaction_ratio.tsv', sep='\t', header=TRUE) 
colnames(Z_500kb_1Mb)[5] <- 'int_freq_500kb'
colnames(Z_500kb_1Mb)[7] <- 'int_freq_1Mb'
colnames(Z_500kb_1Mb)[8] <- 'ratio_500kb_1Mb'

hinch_co_center <- read.table('hinch_co_center.bed', sep='\t', header =FALSE)

davis_dsb <- read.table('Davis_DMC1_hotspot.txt', sep='\t', header =TRUE)
davis_dsb <- davis_dsb %>% mutate(chr= paste('chr', chr, sep='')) %>% filter(allele == 'CAST')
davis_dsb_bg <- davis_dsb %>% mutate(start = dmc1_center, end = dmc1_center+1) %>% select(chr, start, end, dmc1_heat)
write.table(davis_dsb_bg, 'davis_dsb.bedGraph', sep='\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
#table(davis_dsb$allele)

smagulova_dsb <- read.table('GSM1954839_B6fXCASTm_hotspots.tab.bed', sep='\t', header =FALSE)
colnames(smagulova_dsb) <- c('chr', 'start', 'end', 'dmc1')
smagulova_dsb$center <- round((smagulova_dsb$start + smagulova_dsb$end)/2)
smagulova_dsb_bg <- smagulova_dsb %>% mutate(centerplusone = center+1) %>% select(chr, center, centerplusone, dmc1)
write.table(smagulova_dsb_bg, 'smagulova_dsb.bedGraph', sep='\t', row.names = FALSE, col.names = FALSE, quote=FALSE)

############################################################################################################################
#Count the number of DSB and CO occurrence for each region

#Z stage
c1 <- apply(Z_500kb_1Mb, 1, function(x) sum((as.character(hinch_co_center$V1)==as.character(x[2])) &  
                                              (hinch_co_center$V2 > as.numeric(x[3])) & 
                                              (hinch_co_center$V3 < as.numeric(x[4]))))
Z_500kb_1Mb$hinch_co_count <- c1

c2 <- apply(Z_500kb_1Mb, 1, function(x) sum((as.character(davis_dsb$chr)==as.character(x[2])) &  
                                              (davis_dsb$dmc1_center > as.numeric(x[3])) & 
                                              (davis_dsb$dmc1_center < as.numeric(x[4]))))
Z_500kb_1Mb$davis_dsb_count <- c2
d2 <- apply(Z_500kb_1Mb, 1, function(x) sum(davis_dsb$dmc1_heat[(as.character(davis_dsb$chr)==as.character(x[2])) &  
                                                                  (davis_dsb$dmc1_center > as.numeric(x[3])) & 
                                                                  (davis_dsb$dmc1_center < as.numeric(x[4]))]))
Z_500kb_1Mb$davis_DMC1_sum <- d2

c3 <- apply(Z_500kb_1Mb, 1, function(x) sum((as.character(smagulova_dsb$chr)==as.character(x[2])) &  
                                              (smagulova_dsb$start > as.numeric(x[3])) & 
                                              (smagulova_dsb$end < as.numeric(x[4]))))
Z_500kb_1Mb$smagulova_dsb_count <- c3
d3 <- apply(Z_500kb_1Mb, 1, function(x) sum(smagulova_dsb$dmc1[(as.character(smagulova_dsb$chr)==as.character(x[2])) &  
                                                                 (smagulova_dsb$start > as.numeric(x[3])) & 
                                                                 (smagulova_dsb$end < as.numeric(x[4]))]))
Z_500kb_1Mb$smagulova_DMC1_sum <- d3

Z_500kb_1Mb$region_length <- Z_500kb_1Mb$end - Z_500kb_1Mb$start
write.table(Z_500kb_1Mb, 'Z_500kb_1Mb_dsb_co.tsv',sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE )

######################################################################################################################
# Correlation between loopsize and CO probability
chromsizes <- read.table('mm10.chrom.sizes')
colnames(chromsizes) <- c('chrom', 'chrlength')

Z_500kb_1Mb <- read.table('Z_500kb_1Mb_dsb_co.tsv', header=TRUE)
L_500kb_1Mb <- read.table('L_500kb_1Mb_dsb_co.tsv', header=TRUE)

df_co_dsb <- function(df) {
  #normalize co, dsb etc by length of each region
  df <- df %>% mutate(hinch_co_density = hinch_co_count/region_length* 1000000, 
                      davis_dsb_density= davis_dsb_count/region_length* 1000000,
                      smagulova_dsb_density= smagulova_dsb_count/region_length* 1000000,
                      davis_dmc1_coverage_density= davis_DMC1_sum/region_length* 1000000,
                      smagulova_dmc1_coverage_density= smagulova_DMC1_sum/region_length* 1000000,
                      co_davis_dsb_ratio = hinch_co_count/davis_dsb_count,
                      co_smagulova_dsb_ratio = hinch_co_count/smagulova_dsb_count,
                      co_davis_dmc1_ratio = hinch_co_count/davis_DMC1_sum,  
                      co_smagulova_dmc1_ratio = hinch_co_count/smagulova_DMC1_sum)
  
  df <- merge(df, chromsizes, by='chrom')
  df <- df %>% mutate(chrposition = (start+end)/2/chrlength)
  df$near_end <- ifelse(df$chrposition>0.15 & df$chrposition<0.85, FALSE, TRUE)
  
  #divide A compartment regions into 3 groups based on the values of 500kb interactions
  df$A_region_groups <- as.character(df$compartment)
  breaks <- as.numeric(quantile(df$int_freq_500kb[df$compartment=='A'], probs = seq(0, 1, 0.33)))
  df$A_region_groups[df$compartment=='A'] <- 'A2'
  df$A_region_groups[df$compartment=='A' & df$int_freq_500kb> breaks[3] ] <- 'A1'
  df$A_region_groups[df$compartment=='A' & df$int_freq_500kb< breaks[2] ] <- 'A3'
  df$A_region_groups = factor(df$A_region_groups, levels = c('A1', 'A2', 'A3', 'B'))
  return(df)
  
}

Z_500kb_1Mb_co <- df_co_dsb(Z_500kb_1Mb)

write.table(Z_500kb_1Mb_co, 'Z_compartment_co_dsb_summary.tsv', sep='\t', row.names = FALSE, col.names = TRUE,quote=FALSE)


ggplot(Z_500kb_1Mb_co[Z_500kb_1Mb_co$compartment=='A' ,], 
       aes(x= int_freq_500kb, y=co_smagulova_dsb_ratio, color=near_end)) +
  geom_point() +
  stat_smooth(method = "lm", col = "#C42126", se = FALSE,size = 1)+
  theme_bw()

# co_davis_dsb_ratio
p1 = ggplot(Z_500kb_1Mb_co, 
       aes(x= A_region_groups, y=co_davis_dsb_ratio, fill=A_region_groups)) +
  geom_boxplot() +
  theme_bw() + 
  ylim(0,2)+
  geom_jitter(shape=16, position=position_jitter(0.2))

ggsave('hinch_co_davis_dsb_ratio_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)
wilcox.test(Z_500kb_1Mb_co$co_davis_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$co_davis_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A2'])
#p-value = 0.02754

#co_smagulova_dsb_ratio
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=co_smagulova_dsb_ratio, fill=A_region_groups)) +
  geom_boxplot() +
  theme_bw() + 
  ylim(0,2)+
  geom_jitter(shape=16, position=position_jitter(0.2))

ggsave('hinch_co_smagulova_dsb_ratio_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)
wilcox.test(Z_500kb_1Mb_co$co_smagulova_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$co_smagulova_dsb_ratio[Z_500kb_1Mb_co$A_region_groups=='A2'])
#p-value = 0.01854

# co_davis_dmc1_ratio
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=co_davis_dmc1_ratio, fill=A_region_groups)) +
  geom_boxplot() +
  theme_bw() + 
  ylim(0,2)+
  geom_jitter(shape=16, position=position_jitter(0.2))

ggsave('hinch_co_davis_dmc1_ratio_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)
wilcox.test(Z_500kb_1Mb_co$co_davis_dmc1_ratio[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$co_davis_dmc1_ratio[Z_500kb_1Mb_co$A_region_groups=='A2'])
#p-value = 0.03704

#co_smagulova_dmc1_ratio
p1 = ggplot(Z_500kb_1Mb_co, 
            aes(x= A_region_groups, y=co_smagulova_dmc1_ratio, fill=A_region_groups)) +
  geom_boxplot() +
  theme_bw() + 
  ylim(0,2)+
  geom_jitter(shape=16, position=position_jitter(0.2))

ggsave('hinch_co_smagulova_dmc1_ratio_vs_compartment_groups.pdf', p1, width=5, height=5, dpi=300)
wilcox.test(Z_500kb_1Mb_co$co_smagulova_dmc1_ratio[Z_500kb_1Mb_co$A_region_groups=='A1'],
            Z_500kb_1Mb_co$co_smagulova_dmc1_ratio[Z_500kb_1Mb_co$A_region_groups=='A2'])
#p-value = 0.03502