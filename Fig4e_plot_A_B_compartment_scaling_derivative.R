library(pspline)
library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

#################################################################

#calculate log10 of distances and interaction frequencies

calc_derivatives <- function(df) {
  
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
  
  return(derivatives)
}

a=read.csv('Autosomes_A_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)
b=read.csv('Autosomes_B_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)

a<- a[c(1:50), ]
b<- b[c(1:50), ]


write.table(calc_derivatives(a), 'autosomes_A_compartment_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(calc_derivatives(b), 'autosomes_B_compartment_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#####################################################################################
#Smoothen derivatives
df_smooth=read.csv('autosomes_A_compartment_scaling_derivatives.tsv', sep='\t', header = TRUE)
for (n in c(2:12)) {
  df_smooth[,n] <- smooth(df_smooth[,n]) 
}
write.table(df_smooth, 'autosomes_A_compartment_scaling_derivatives_turkey_smooth_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

df_smooth=read.csv('autosomes_B_compartment_scaling_derivatives.tsv', sep='\t', header = TRUE)
for (n in c(2:12)) {
  df_smooth[,n] <- smooth(df_smooth[,n]) 
}
write.table(df_smooth, 'autosomes_B_compartment_scaling_derivatives_turkey_smooth_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


#######################################################################################################################
#loess smoothen derivatives

loess_smooth <- function(df) {
  sertoli=loess.smooth(df$log_dist, df$sertoli, span = 1/3, degree = 1,
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
  return (df_loess)
}

a_derivatives=read.csv('autosomes_A_compartment_scaling_derivatives.tsv', sep='\t', header = TRUE)
b_derivatives=read.csv('autosomes_B_compartment_scaling_derivatives.tsv', sep='\t', header = TRUE)

write.table(loess_smooth(a_derivatives), 'autosomes_A_compartment_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(loess_smooth(b_derivatives), 'autosomes_B_compartment_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)



#####################################################################################
#plot derivatives
# A vs B for each stage

a=read.csv('autosomes_A_compartment_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', header = TRUE)
b=read.csv('autosomes_B_compartment_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', header = TRUE)

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

cols <- c("A" = "red", "B" = "dodgerblue")
for (n in c(1:length(stages))) {
  stage=stages[n]
  df<- data.frame(a$distances)
  colnames(df) <-'distances'
  df['A']<- a[stage]
  df['B']<- b[stage]
  
  p=ggplot(data=df, aes(x = distances))+
    geom_line(aes(y= A, colour='A'), linetype=1, size=1.2)+
    geom_line(aes(y= B, colour='B'), linetype=1, size=1.2)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A', 'B'),
                        labels= c(paste(stage, ", A Compartment", sep=""), paste(stage, ", B Compartment", sep=""))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e4, 5e6)) +
    scale_y_continuous(limits = c(-1.5, 0.5)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling slope')+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_A_B_scaling_slope.pdf', sep=''), width = 10, height = 5)
  print(p)
  dev.off()
  
}


########################################################################################################

# Make plots across stages
plot_scaling <- function(df, plot_prefix){
  cols <- c("sertoli" = "red", "spg" = "orange", "PreL" = "darkorchid1", "L" = "lightblue1", "Z"= "cyan", "P"="dodgerblue", "D"="navy", "Sun1"="green", "Psline1"="grey60", "Psline2"="grey40")
  
  # sertoli, spg, PreL, L
  p=ggplot(data=df, aes(x = distances))+
    geom_line(aes(y= sertoli, colour='sertoli'), linetype=1, size=1.2)+
    geom_line(aes(y= spg, colour='spg'), linetype=1, size=1.2)+
    geom_line(aes(y= PreL, colour='PreL'), linetype=1, size=1.2)+
    geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('sertoli', 'spg', 'PreL', 'L', 'Psline1'),
                        labels= c('sertoli', 'spg', 'PreL', 'L', 'P~s^(-0.6)')) +

    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A', 'B'),
                        labels= c(paste(stage, ", A Compartment", sep=""), paste(stage, ", B Compartment", sep=""))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e4, 5e6)) +
    scale_y_continuous(limits = c(-1.5, 0.5)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling slope')+
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
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('L', 'Z', 'P', 'D', 'Psline1'),
                        labels= c('L', 'Z', 'P', 'D', 'P~s^(-0.6)')) +
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A', 'B'),
                        labels= c(paste(stage, ", A Compartment", sep=""), paste(stage, ", B Compartment", sep=""))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e4, 5e6)) +
    scale_y_continuous(limits = c(-1.5, 0.5)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling slope')+
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
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('Z', 'Sun1', 'Psline1'),
                        labels= c('Z', 'Sun1', 'P~s^(-0.6)')) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A', 'B'),
                        labels= c(paste(stage, ", A Compartment", sep=""), paste(stage, ", B Compartment", sep=""))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e4, 5e6)) +
    scale_y_continuous(limits = c(-1.5, 0.5)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling slope')+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(paste(plot_prefix, '_scaling_slope_SSPL.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
  pdf(paste(plot_prefix, '_scaling_slope_LZPD.pdf', sep=''), width = 10, height = 10)
  print(q)
  dev.off()
  
  pdf(paste(plot_prefix, '_scaling_slope_Sun1_WT_Z.pdf', sep=''), width = 10, height = 10)
  print(s)
  dev.off()
  
}

a=read.csv('autosomes_A_compartment_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', header = TRUE)
b=read.csv('autosomes_B_compartment_scaling_derivatives_loess_smooth_202005.tsv', sep='\t', header = TRUE)
plot_scaling(a, 'Autosomes_A_compartment')
plot_scaling(b, 'Autosomes_B_compartment')
