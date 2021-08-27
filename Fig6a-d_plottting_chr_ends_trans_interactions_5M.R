#########################################################################################
# data preparation
#########################################################################################

autosomes=paste('chr', c(1:19), sep='')
stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")
trans_int=c('trans_cen_cen', 'trans_cen_tel', 'trans_tel_tel')

for (trans_int_type in trans_int){
  tmpfile<-read.csv(paste('L', trans_int_type, '5000000.tsv', sep='_'), sep='\t', header = TRUE)
  chrom1<-tmpfile$chrom1
  chrom2<-tmpfile$chrom2
  interactions_summary<-data.frame(chrom1, chrom2)
  for (stage in stages) {
    filename=paste(stage, trans_int_type, '5000000.tsv', sep='_')
    df=read.csv(filename, sep='\t', header = TRUE)
    interactions_summary[stage] <- df$balanced.average
  }
  save_file=paste(trans_int_type, 'autosome_summary_2.tsv', sep='_')
  write.table(interactions_summary, save_file, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}


df1=read.csv('trans_cen_cen_autosome_summary_2.tsv', sep='\t', header = TRUE)
df2=read.csv('trans_tel_tel_autosome_summary_2.tsv', sep='\t', header = TRUE)
df3=read.csv('trans_cen_tel_autosome_summary_2.tsv', sep='\t', header = TRUE)


chrom1<-as.vector(df1$chrom1)
chrom2<-as.vector(df1$chrom2)
trans_polarity<-data.frame(chrom1, chrom2)
for (stage in stages) {
  cen_cen <- df1[stage]
  tel_tel <- df2[stage]
  cen_tel <-c()
  for ( n in 1:nrow(cen_cen)) {
    c1 <- chrom1[n]
    c2 <- chrom2[n]
    cen_tel_tmp <- unlist(df3[stage])
    cen_tel<- c(cen_tel, mean(cen_tel_tmp[(df3$chrom1 == c1 & df3$chrom2 == c2) | (df3$chrom1 == c2 & df3$chrom2 == c1)]))
  }
  polarity=(cen_cen+tel_tel )/2/cen_tel 
  trans_polarity[stage] <- polarity
}
save_file=paste('trans_polarity_autosome_summary_2.tsv', sep='_')
write.table(trans_polarity, save_file, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


##############################################################################
# Plot interactions through stages for each interaction type
##############################################################################
trans_int=c('trans_cen_cen', 'trans_cen_tel', 'trans_tel_tel')
int_type=  trans_int

for (int in int_type) {
  filename<-paste(int, 'autosome_summary_2.tsv', sep='_')
  data<- read.csv(filename, sep='\t', header = TRUE)
  
  sertoli <- data$sertoli
  spg <- data$spg
  preL <- data$PreL
  L <- data$L
  Z <- data$Z
  
  P <- data$P
  D <- data$D
  
  pdf(file= paste(int, '_interaction_frequency_by_chrpair.pdf', sep=''), width=12, height=8)
  boxplot(sertoli, spg, preL, L, Z, P, D, names=c("sertoli", "spg", "PreL", "L", "Z", "P", "D"),
          col="gold", outline=FALSE, 
          xlab = "stages", ylab = "average interaction frequency", main = paste(int, '_interaction_frequency', sep=''))
  dev.off()  
  
}

#############################################################
#trans plots

int='trans_cen_cen'
filename<-paste(int, 'autosome_summary_2.tsv', sep='_')
data<- read.csv(filename, sep='\t', header = TRUE)

sertoli <- data$sertoli
spg <- data$spg
preL <- data$PreL
L <- data$L
Z <- data$Z
P <- data$P
D <- data$D

pdf(file= paste(int, '_interaction_frequency_by_chrpair_v3.pdf', sep=''), width=12, height=8)
boxplot(sertoli, spg, preL, L, Z, P, D, names=c("sertoli", "spg", "PreL", "L", "Z", "P", "D"),
        col="red", outline=FALSE, ylim=c(1e-6, 1e-5),
        xlab = "stages", ylab = "average interaction frequency", main = paste(int, '_interaction_frequency', sep=''))
dev.off() 


int='trans_tel_tel'
filename<-paste(int, 'autosome_summary_2.tsv', sep='_')
data<- read.csv(filename, sep='\t', header = TRUE)

sertoli <- data$sertoli
spg <- data$spg
preL <- data$PreL
L <- data$L
Z <- data$Z
P <- data$P
D <- data$D

pdf(file= paste(int, '_interaction_frequency_by_chrpair_v3.pdf', sep=''), width=12, height=8)
boxplot(sertoli, spg, preL, L, Z, P, D, names=c("sertoli", "spg", "PreL", "L", "Z", "P", "D"),
        col="gold", outline=FALSE, ylim=c(1e-6, 8e-6),
        xlab = "stages", ylab = "average interaction frequency", main = paste(int, '_interaction_frequency', sep=''))
dev.off() 

int='trans_cen_tel'
filename<-paste(int, 'autosome_summary_2.tsv', sep='_')
data<- read.csv(filename, sep='\t', header = TRUE)

sertoli <- data$sertoli
spg <- data$spg
preL <- data$PreL
L <- data$L
Z <- data$Z
P <- data$P
D <- data$D

pdf(file= paste(int, '_interaction_frequency_by_chrpair_v3.pdf', sep=''), width=12, height=6)
boxplot(sertoli, spg, preL, L, Z, P, D, names=c("sertoli", "spg", "PreL", "L", "Z", "P", "D"),
        col="dodgerblue", outline=FALSE, ylim=c(1e-6, 4e-6),
        xlab = "stages", ylab = "average interaction frequency", main = paste(int, '_interaction_frequency', sep=''))
dev.off() 

###################################################################################

###################################################################################
# plot polarity for cis and trans

data<- read.csv('trans_polarity_autosome_summary_2.tsv', sep='\t', header = TRUE)

sertoli <- data$sertoli
spg <- data$spg
preL <- data$PreL
L <- data$L
Z <- data$Z
P <- data$P
D <- data$D

pdf(file= paste('trans interaction polarity by chrpair.pdf', sep=''), width=12, height=6)
boxplot(sertoli, spg, preL, L, Z, P, D, names=c("sertoli", "spg", "PreL", "L", "Z", "P", "D"),
        col="cyan", outline=FALSE,
        xlab = "stages", ylab = "(CEN_CEN + TEL_TEL)/(2*CEL_TEL)", main = paste('trans interaction polarity', sep=''))
dev.off() 


#####################################################################################################
#trans plots wt vs sun-1

cis_int=c('cis_cen', 'cis_cen_tel', 'cis_tel')
trans_int=c('trans_cen_cen', 'trans_tel_tel','trans_cen_tel')

chrs=c()
wt=c()
sun1=c()
int_type=c()

for (int in trans_int) {
  filename<-paste(int, 'autosome_summary_2.tsv', sep='_')
  data<- read.csv(filename, sep='\t', header = TRUE)
  
  chrs=c(chrs, as.vector(data$autosomes))
  wt=c(wt, data$Z)
  sun1=c(sun1, data$Sun1)
  int_type= c(int_type, rep(int, nrow(data)))
}

df<- data.frame(chrs, wt, sun1, int_type)
colnames(df) <- c('chrs', 'WT_Z', 'Sun1_Z_like', "interaction_category")
  
  pdf(file= paste(int, '_interaction_frequency_by_chrpair.pdf', sep=''), width=12, height=8)
  boxplot(sertoli, spg, preL, L, Z, P, D, names=c("sertoli", "spg", "PreL", "L", "Z", "P", "D"),
          col="gold", outline=FALSE, 
          xlab = "stages", ylab = "average interaction frequency", main = paste(int, '_interaction_frequency', sep=''))
  dev.off()  

library(reshape2) 
  
df.melt<- melt(df, measure.vars = c("WT_Z", "Sun1_Z_like"))

colnames(df.melt)[3]<-'genotype'
df.melt$interaction_category <- ordered(df.melt$interaction_category, levels=c('trans_cen_cen', 'trans_tel_tel','trans_cen_tel'))


library(ggplot2)
library(extrafont)
loadfonts(device = "win")


pdf('wt_vs_sun1_chromosome_ends_trans_interaction_by_chrpair.pdf', height=8, width=6)
ggplot(df.melt, aes(x=interaction_category, y=value)) + 
  geom_boxplot(aes(fill = genotype), position="dodge", outlier.shape = NA)+
  theme_bw()+
  labs(x='', y='Average interaction frequency per 10kb bin')+
  theme(legend.position = 'right')
dev.off()

#legend.title = element_text(family='Arial Black', size=12),
#legend.text = element_text(family='Arial', size=10),
#axis.text.x = element_text(family='Arial', size=10)


