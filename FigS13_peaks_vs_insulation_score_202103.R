library(ggplot2)
##################################################
#prep data

ctcf_rec8 <- read.table('PD_CTCF_REC8_overlapping_10kb.bed')
rec8_only <- read.table('PD_REC8_only_10kb.bed')
ctcf_only <- read.table('PD_CTCF_only_10kb.bed')

#read insulation profiles
#sertoli_insulation <- read.table('sertoli--is500001--nt0.1--ids200001--ss1--immean.insulation.bedGraph', skip = 1, sep='\t')
#preL_insulation <- read.table('preL--is500001--nt0.1--ids200001--ss1--immean.insulation.bedGraph', skip = 1, sep='\t')

sertoli_insulation <- read.table('Sertoli_10kb_insulation_500kb_window.tsv', header=TRUE, sep='\t')
preL_insulation <- read.table('PreL_10kb_insulation_500kb_window.tsv', header=TRUE, sep='\t')

ins_data <- function(insulation_filename, stage) {
  name <- paste(stage, 'insulation', sep='_')
  insulation <- read.table(insulation_filename, header=TRUE, sep='\t')
  insulation$peak_name <- paste(insulation$chrom, insulation$start, insulation$end, sep='_')
  insulation[name] <- insulation$log2_insulation_score_500000
  return(insulation)
}


#merge data

ctcf_only$peak_name <- paste(ctcf_only$V1, ctcf_only$V2, ctcf_only$V3, sep='_')
ctcf_only$peak_type <- 'CTCF_only'

rec8_only$peak_name <- paste(rec8_only$V1, rec8_only$V2, rec8_only$V3, sep='_')
rec8_only$peak_type <- 'REC8_only'

ctcf_rec8$peak_name <- paste(ctcf_rec8$V1, ctcf_rec8$V2, ctcf_rec8$V3, sep='_')
ctcf_rec8$peak_type <- 'both_CTCF_REC8'

all_peaks<- rbind(ctcf_only, rec8_only, ctcf_rec8)


sertoli_insulation <- ins_data('Sertoli_10kb_insulation_500kb_window.tsv', 'sertoli')
spg_insulation <- ins_data('spg_10kb_insulation_500kb_window.tsv', 'spg')
preL_insulation <- ins_data('PreL_10kb_insulation_500kb_window.tsv', 'preL')
L_insulation <- ins_data('L_10kb_insulation_500kb_window.tsv', 'L')
Z_insulation <- ins_data('Z_10kb_insulation_500kb_window.tsv', 'Z')
P_insulation <- ins_data('P_10kb_insulation_500kb_window.tsv', 'P')
D_insulation <- ins_data('D_10kb_insulation_500kb_window.tsv', 'D')


peaks_insulation <- merge(sertoli_insulation, all_peaks, by='peak_name', all.x = TRUE)
peaks_insulation$peak_type [is.na(peaks_insulation$peak_type)] <- 'Not Peaks'
peaks_insulation <- merge(peaks_insulation, spg_insulation, by='peak_name')
peaks_insulation <- merge(peaks_insulation, preL_insulation, by='peak_name')
peaks_insulation <- merge(peaks_insulation, L_insulation, by='peak_name')
peaks_insulation <- merge(peaks_insulation, Z_insulation, by='peak_name')
peaks_insulation <- merge(peaks_insulation, P_insulation, by='peak_name')
peaks_insulation <- merge(peaks_insulation, D_insulation, by='peak_name')

peaks_insulation <- peaks_insulation[, c('peak_name', 'peak_type', 'sertoli_insulation','spg_insulation','preL_insulation',
                                         'L_insulation', 'Z_insulation','P_insulation','D_insulation')]
write.table(peaks_insulation, 'peaks_vs_insulation_202103.tsv', sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

#################################################################################################
peaks_insulation <- read.table('peaks_vs_insulation_202103.tsv', sep='\t', header=TRUE, row.names = NULL)
peaks_insulation$insulation_change <- peaks_insulation$preL_insulation - peaks_insulation$sertoli_insulation
peaks_insulation$peak_type <- factor(peaks_insulation$peak_type,
                       levels = c('both_CTCF_REC8','CTCF_only', 'REC8_only', 'Not Peaks'),ordered = TRUE)

table(peaks_insulation$peak_type)
#both_CTCF_REC8      CTCF_only      REC8_only      Not Peaks 
#5745                12653           9421         244747

p1= ggplot(peaks_insulation, aes(x=peak_type, y=sertoli_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")

ggsave('peaks_insulation_score_sertoli_202103.pdf', p1, width=5, height=8, dpi=300)


wilcox.test(peaks_insulation$sertoli_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$sertoli_insulation[peaks_insulation$peak_type=='CTCF_only'])
# p-value = 1.989e-06

p2= ggplot(peaks_insulation, aes(x=peak_type, y=preL_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")
ggsave('peaks_insulation_score_preL_202103.pdf', p2, width=5, height=8, dpi=300)

wilcox.test(peaks_insulation$preL_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$preL_insulation[peaks_insulation$peak_type=='CTCF_only'])
# p-value < 2.2e-16

p= ggplot(peaks_insulation, aes(x=peak_type, y=insulation_change, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")
ggsave('peaks_insulation_change_sertoli_vs_preL_202103.pdf', p, width=5, height=8, dpi=300)


wilcox.test(peaks_insulation$insulation_change[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$insulation_change[peaks_insulation$peak_type=='CTCF_only'])
# p-value < 2.2e-16

#L stage
p3= ggplot(peaks_insulation, aes(x=peak_type, y=L_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")
ggsave('peaks_insulation_score_L_202103.pdf', p3, width=5, height=8, dpi=300)

wilcox.test(peaks_insulation$L_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$L_insulation[peaks_insulation$peak_type=='CTCF_only'])
#p-value < 2.2e-16


#Z stage
p4= ggplot(peaks_insulation, aes(x=peak_type, y=Z_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")
ggsave('peaks_insulation_score_Z_202103.pdf', p4, width=5, height=8, dpi=300)


wilcox.test(peaks_insulation$Z_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$Z_insulation[peaks_insulation$peak_type=='CTCF_only'])
#p-value < 2.2e-16
#p-value < 2.2e-16

#P stage
p5= ggplot(peaks_insulation, aes(x=peak_type, y=P_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")
ggsave('peaks_insulation_score_P_202103.pdf', p5, width=5, height=8, dpi=300)


wilcox.test(peaks_insulation$P_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$P_insulation[peaks_insulation$peak_type=='CTCF_only'])

#p-value < 2.2e-16


#D stage
p6= ggplot(peaks_insulation, aes(x=peak_type, y=D_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")
ggsave('peaks_insulation_score_D_202103.pdf', p6, width=5, height=8, dpi=300)

wilcox.test(peaks_insulation$D_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$D_insulation[peaks_insulation$peak_type=='CTCF_only'])
#p-value < 2.2e-16

#spg
p7= ggplot(peaks_insulation, aes(x=peak_type, y=spg_insulation, fill=peak_type)) + 
  geom_boxplot() +
  ylim(-1, 1)+
  theme_bw()+
  scale_fill_manual(values=c('firebrick1', 'darkorange1', 'steelblue3', 'grey'))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")

ggsave('peaks_insulation_score_spg_202103.pdf', p7, width=5, height=8, dpi=300)


wilcox.test(peaks_insulation$spg_insulation[peaks_insulation$peak_type=='both_CTCF_REC8'],
            peaks_insulation$spg_insulation[peaks_insulation$peak_type=='CTCF_only'])
# p-value = 1.989e-06
