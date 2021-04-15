library(pspline)
library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

#################################################################
df=read.csv('autosomes_scaling_summary_keep_diag_202005.tsv', sep='\t', header = TRUE)

#calculate log10 of distances and interaction frequencies
lgdf<- log10(df)
colnames(lgdf)[1] <- 'log_dist'
lgdf$distances <- df$distances

#calculate 1st order derivative at each distance
derivatives<- lgdf
stages<- c("sertoli", "spg", "PreL", "L",  "Z",  "Sun1", "P", "D", "MII", "RS", "sperm")

for (stage in stages) {
  x<- lgdf$log_dist
  y<- lgdf[,stage]
  d<- predict(sm.spline(x, y), x, 1)
  derivatives[,stage] <- d
  
}

write.table(derivatives, 'autosomes_scaling_derivatives_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#####################################################################################
#Smoothen derivatives

df_smooth=read.csv('autosomes_scaling_derivatives_202005.tsv', sep='\t', header = TRUE)

for (n in c(2:12)) {
  df_smooth[,n] <- smooth(df_smooth[,n]) 
}

write.table(df_smooth, 'autosomes_scaling_derivatives_turkey_smooth_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#loess smoothen
df=read.csv('autosomes_scaling_derivatives_202005.tsv', sep='\t', header = TRUE)

sertoli=loess.smooth(df_smooth$log_dist, df_smooth$sertoli, span = 1/3, degree = 1,
             family ="gaussian", evaluation = 77)

df_loess <- data.frame(x=sertoli$x)
df_loess$distances <- round(10^df_loess$x)
stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

for (stage in stages) {
  x<- df$log_dist
  y<- df[,stage]
  loessfit <- loess.smooth(x, y, span = 1/3, degree = 1, family ="gaussian", evaluation = 77)
  df_loess[,stage] <- loessfit$y
}

write.table(df_loess, 'autosomes_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)



#####################################################################################
#plot derivatives

df=read.csv('autosomes_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', header = TRUE)

prefix='autosomes_scaling_derivatives_loess_smooth'

cols <- c("sertoli" = "red", "spg" = "orange", "PreL" = "darkorchid1", "L" = "lightblue1", "Z"= "cyan", "P"="dodgerblue", "D"="navy", "Sun1"="green", "Psline1"="grey60", "Psline2"="grey40")

# sertoli, spg, PreL, L
p=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= sertoli, colour='sertoli'), linetype=1, size=1.2)+
  geom_line(aes(y= spg, colour='spg'), linetype=1, size=1.2)+
  geom_line(aes(y= PreL, colour='PreL'), linetype=1, size=1.2)+
  geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('sertoli', 'spg', 'PreL', 'L'),
                      labels= c('sertoli', 'spg', 'PreL', 'L')) +
  #geom_vline(xintercept = 1e5, linetype=2, size=1)+
  #geom_vline(xintercept = 2e6, linetype=2, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e4, 1e8)) +
  scale_y_continuous(limits = c(-2.5, -0.2)) +
  theme_bw()+
  labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling plot')+
  theme(legend.position = 'right',
        legend.text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5), 
        axis.text.x.bottom = element_text(size=15),
        axis.text.y.left = element_text(size=15),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15))

# L, Z, P, D

q=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
  geom_line(aes(y= Z, colour='Z'), linetype=1, size=1.2)+
  geom_line(aes(y= P, colour='P'), linetype=1, size=1.2)+
  geom_line(aes(y= D, colour='D'), linetype=1, size=1.2)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('L', 'Z', 'P', 'D'),
                      labels= c('L', 'Z', 'P', 'D')) +
  #geom_vline(xintercept = 1e5, linetype=2, size=1)+
  #geom_vline(xintercept = 2e6, linetype=2, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e4, 1e8)) +
  scale_y_continuous(limits = c(-2.5, -0.2)) +
  theme_bw()+
  labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling plot')+
  theme(legend.position = 'right',
        legend.text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5), 
        axis.text.x.bottom = element_text(size=15),
        axis.text.y.left = element_text(size=15),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15))

s=ggplot(data=df, aes(x = distances))+
  geom_line(aes(y= Z, colour='Z'), linetype=1, size=1.2)+
  geom_line(aes(y= Sun1, colour='Sun1'), linetype=1, size=1.2)+
  #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
  scale_colour_manual("", 
                      values = cols, 
                      breaks= c('Z', 'Sun1'),
                      labels= c('Z', 'Sun1')) +
  #geom_vline(xintercept = 1e5, linetype=2, size=1)+
  #geom_vline(xintercept = 2e6, linetype=2, size=1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e4, 1e8)) +
  scale_y_continuous(limits = c(-2.5, -0.2)) +
  theme_bw()+
  labs(x='Separation (bp)', y='Scaling slope', title='Autosomes scaling plot')+
  theme(legend.position = 'right',
        legend.text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5), 
        axis.text.x.bottom = element_text(size=15),
        axis.text.y.left = element_text(size=15),
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15))

pdf(paste(prefix, '_1.pdf', sep=''), width = 20, height = 10)
p
dev.off()

pdf(paste(prefix, '_2.pdf', sep=''), width = 20, height = 10)
q
dev.off()

pdf(paste(prefix, '_WT_Sun1_Z.pdf', sep=''), width = 20, height = 10)
s
dev.off()


#########################################################################

lgdf2<- lgdf[lgdf$distances>500000 & lgdf$distances<3000000 , ]
lm(lgdf2$sertoli ~ lgdf2$log_dist)$coeff[[2]]
# [1] -1.176305
