setwd("~/Documents/project/meiosis/scaling/")
# replicates scaling
stages=c('Sertoli1', 'Sertoli2')
stages=c('PreL1', 'PreL2')
stages=c('L1', 'L2')
stages=c('Z1', 'Z2')
stages=c('P1', 'P2')
stages=c('D1', 'D2')
stages=c('spg1', 'spg2')
stages=c('MII1', 'MII2')

bins=unique(round(1.12**c(1:87)))
distances=bins*10000

scaling_summary<-data.frame(distances)
head(scaling_summary)
for (stage in stages) {
  filename=paste(stage, '_autosomes_interaction_versus_distance_keep_diag.tsv', sep='')
  expected_df=read.csv(filename, sep='\t', header = TRUE)
  frequency=expected_df$frequency[expected_df$diag %in% bins]
  scaling_summary[stage] <- frequency
}
head(scaling_summary)
colnames(scaling_summary)[c(2:3)]<-c('rep1','rep2')
write.table(scaling_summary, 'autosomes_scaling_summary_keep_diag1_MII.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


bins=unique(round(1.12**c(1:86)))
distances=bins*10000

scaling_summary_x<-data.frame(distances)

for (stage in stages) {
  filename=paste(stage, '_x_interaction_versus_distance_keep_diag.tsv', sep='')
  expected_df=read.csv(filename, sep='\t', header = TRUE)
  frequency=expected_df$frequency[expected_df$diag %in% bins]
  scaling_summary_x[stage] <- frequency
}
head(scaling_summary_x)
colnames(scaling_summary_x)[c(2:3)]<-c('rep1','rep2')
write.table(scaling_summary_x, 'x_scaling_summary_keep_MII.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

################################
library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

#########################################################################
#Making autosomal scaling plot
#########################################################################

df=read.csv('autosomes_scaling_summary_keep_diag1_Sertoli.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_spg.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_PreL.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_MII.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_L.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_Z.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_P.tsv', sep='\t', header = TRUE)
df=read.csv('autosomes_scaling_summary_keep_diag1_D.tsv', sep='\t', header = TRUE)

df['Psline1']=(df['distances']^ -0.6)*15

cols <- c( "rep1" = "red","rep2"="darkorchid1", "Psline1"="grey60", "Psline2"="grey40")

for (n in c(2:3)) {
  df[,n] <- df[,n] * 0.01 / df[df$distances==50000, n] 
}

p=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= rep1, colour='rep1'), linetype=1, size=1.2)+
  geom_line(aes(y= rep2, colour='rep2'), linetype=1, size=1.2)+
  geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('rep1', 'rep2', 'Psline1'),
                      labels= c('rep1', 'rep2', 'P~s^(-0.6)')) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e4, 1e8)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-5, 1e-1)) +
  theme_bw()+
  labs(x='Separation (bp)', y='Contact Probability', title='Autosomes scaling plot')+
  theme(legend.position = 'right',
        legend.text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5), 
        axis.text.x.bottom = element_text(size=15),
        axis.text.y.left = element_text(size=15),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15))

pdf(file='autosomes_scaling_D.pdf', width = 10, height = 10)
p
dev.off()

#########################################################################
#Making X chr scaling plot
#########################################################################
df=read.csv('x_scaling_summary_keep_Sertoli.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_spg.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_MII.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_PreL.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_L.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_Z.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_P.tsv', sep='\t', header = TRUE)
df=read.csv('x_scaling_summary_keep_D.tsv', sep='\t', header = TRUE)

# normalize at 50kb
df['Psline1']=(df['distances']^ -0.6)*15

cols <- c( "rep1" = "red","rep2"="darkorchid1", "Psline1"="grey60", "Psline2"="grey40")

for (n in c(2:3)) {
  df[,n] <- df[,n] * 0.01 / df[df$distances==50000, n] 
}

p=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= rep1, colour='rep1'), linetype=1, size=1.2)+
  geom_line(aes(y= rep2, colour='rep2'), linetype=1, size=1.2)+
  geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('rep1', 'rep2', 'Psline1'),
                      labels= c('rep1', 'rep2', 'P~s^(-0.6)')) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e4, 1e8)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-5, 1e-1)) +
  theme_bw()+
  labs(x='Separation (bp)', y='Contact Probability', title='X chromosome scaling plot')+
  theme(legend.position = 'right',
        legend.text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5), 
        axis.text.x.bottom = element_text(size=15),
        axis.text.y.left = element_text(size=15),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15))

pdf(file='X_chromosome_scaling_D.pdf', width = 10, height = 10)
p
dev.off()
