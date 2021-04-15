#######################################################################
bins=unique(round(1.12**c(1:86)))
distances=bins*10000

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

scaling_summary_x<-data.frame(distances)

for (stage in stages) {
  filename=paste(stage, '_x_interaction_versus_distance_keep_diag.tsv', sep='')
  expected_df=read.csv(filename, sep='\t', header = TRUE)
  frequency=c()
  for (n in 1:length(bins)) {
    start <- ifelse(n>1, bins[n-1], 0)
    end <- bins[n]
    total_int=sum(expected_df$balanced.sum[expected_df$diag > start & expected_df$diag <= end])
    total_bins=sum(expected_df$n_valid[expected_df$diag > start & expected_df$diag <= end])
    ave_int = total_int / total_bins
    frequency <- c(frequency, ave_int)
  }
  
  scaling_summary[stage] <- frequency
  
  
}

write.table(scaling_summary_x, 'x_scaling_summary_keep_diag_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
#######################################################################

library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

#########################################################################
#Making X chr scaling plot
#########################################################################
df=read.csv('x_scaling_summary_keep_diag.tsv', sep='\t', header = TRUE)

for (n in c(2:12)) {
  df[,n] <- df[,n] * 0.01 / df[df$distances==10000, n] 
}

df['Psline1']=(df['distances']^ -0.6)*15
#df['Psline2']=(df['distances']^ -1.5)*10000

cols <- c("sertoli" = "red", "spg" = "orange", "PreL" = "darkorchid1", "L" = "lightblue1", "Z"= "cyan", "P"="dodgerblue", "D"="navy", "Sun1"="green", "Psline1"="grey60", "Psline2"="grey40")

# sertoli, spg, PreL, L
p=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= sertoli, colour='sertoli'), linetype=1, size=1.2)+
  geom_line(aes(y= spg, colour='spg'), linetype=1, size=1.2)+
  geom_line(aes(y= PreL, colour='PreL'), linetype=1, size=1.2)+
  geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
  geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('sertoli', 'spg', 'PreL', 'L', 'Psline1'),
                      labels= c('sertoli', 'spg', 'PreL', 'L', 'P~s^(-0.6)')) +
  #geom_vline(xintercept = 1e5, linetype=2, size=1)+
  #geom_vline(xintercept = 2e6, linetype=2, size=1)+
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

p

# L, Z, P, D

q=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
  geom_line(aes(y= Z, colour='Z'), linetype=1, size=1.2)+
  geom_line(aes(y= P, colour='P'), linetype=1, size=1.2)+
  geom_line(aes(y= D, colour='D'), linetype=1, size=1.2)+
  geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('L', 'Z', 'P', 'D', 'Psline1'),
                      labels= c('L', 'Z', 'P', 'D', 'P~s^(-0.6)')) +
  #geom_vline(xintercept = 1e5, linetype=2, size=1)+
  #geom_vline(xintercept = 2e6, linetype=2, size=1)+
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

q

pdf(file='X_scaling1.pdf', width = 10, height = 10)
p
dev.off()

pdf(file='X_scaling2.pdf', width = 10, height = 10)
q
dev.off()

s=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= Z, colour='Z'), linetype=1, size=1.2)+
  geom_line(aes(y= Sun1, colour='Sun1'), linetype=1, size=1.2)+
  geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('Z', 'Sun1', 'Psline1'),
                      labels= c('Z', 'Sun1', 'P~s^(-0.6)')) +
  #geom_vline(xintercept = 1e5, linetype=2, size=1)+
  #geom_vline(xintercept = 2e6, linetype=2, size=1)+
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


pdf(file='X_scaling1_202005.pdf', width = 10, height = 10)
p
dev.off()

pdf(file='X_scaling2_202005.pdf', width = 10, height = 10)
q
dev.off()

pdf(file='X_scaling_Z_vs_Sun1_202005.pdf', width = 10, height = 10)
s
dev.off()