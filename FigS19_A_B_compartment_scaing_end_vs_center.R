library(pspline)
library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

#########################################################################
#read data

a_left=read.csv('Autosomes_A_compartment_left_end_expected_vs_distances.tsv', sep='\t', header = TRUE)
a_right=read.csv('Autosomes_A_compartment_right_end_expected_vs_distances.tsv', sep='\t', header = TRUE)
a_middle=read.csv('Autosomes_A_compartment_middle_expected_vs_distances.tsv', sep='\t', header = TRUE)

b_left=read.csv('Autosomes_B_compartment_left_end_expected_vs_distances.tsv', sep='\t', header = TRUE)
b_right=read.csv('Autosomes_B_compartment_right_end_expected_vs_distances.tsv', sep='\t', header = TRUE)
b_middle=read.csv('Autosomes_B_compartment_middle_expected_vs_distances.tsv', sep='\t', header = TRUE)

#########################################################################
#Plot end vs middle scaling plot for each stage

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

for (stage in stages) {
  
  df<- data.frame(a$distances)
  colnames(df) <-'distances'
  df['A_left']<- a_left[stage]
  df['A_right']<- a_right[stage]
  df['A_middle']<- a_middle[stage]
  
  cols <- c("A_left" = "orange", "A_middle" = "dodgerblue", "A_right" = "red")
  
  #df['Psline1']=(df['distances']^ -0.6)*15
  #df['Psline2']=(df['distances']^ -1.5)*10000
  
  # sertoli, spg, PreL, L
  p=ggplot(data=df, aes(x = distances))+
    geom_line(aes(y= A_left, colour='A_left'), linetype=1, size=1.2)+
    geom_line(aes(y= A_right, colour='A_right'), linetype=1, size=1.2)+
    geom_line(aes(y= A_middle, colour='A_middle'), linetype=1, size=1.2)+
    #geom_line(aes(y= A, colour='A'), linetype=1, size=1.2)+
    #geom_line(aes(y= B, colour='B'), linetype=1, size=1.2)+
    #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A_left', 'A_middle', 'A_right'),
                        labels= c(paste(stage, ", A_left", sep=""), 
                                  paste(stage, ", A_middle", sep=""),
                                  paste(stage, ", A_right", sep=""))) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e5, 5e6)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e-4, 5e-3)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Average interaction frequency per 10kb bin', title=paste(stage, 'Interaction frequency vs distances'))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_A_chromosome_ends_vs_center_scaling_50k-5Mb.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
}


#######################################################################################
# scaling plots from 50kb to 5Mb, loess smooth

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

for (stage in stages) {
  
  df<- data.frame(a$distances)
  colnames(df) <-'distances'
  df['A_left']<- a_left[stage]
  df['A_right']<- a_right[stage]
  df['A_middle']<- a_middle[stage]
  
  cols <- c("A_left" = "orange", "A_middle" = "dodgerblue", "A_right" = "red")
  
  #df['Psline1']=(df['distances']^ -0.6)*15
  #df['Psline2']=(df['distances']^ -1.5)*10000
  
  # sertoli, spg, PreL, L
  p=ggplot(data=df, aes(x = distances))+
    geom_smooth(aes(y= A_left, colour='A_left'), linetype=1, size=1.2, method='loess', se=F)+
    geom_smooth(aes(y= A_right, colour='A_right'), linetype=1, size=1.2, method='loess', se=F)+
    geom_smooth(aes(y= A_middle, colour='A_middle'), linetype=1, size=1.2, method='loess', se=F)+
    #geom_line(aes(y= A, colour='A'), linetype=1, size=1.2)+
    #geom_line(aes(y= B, colour='B'), linetype=1, size=1.2)+
    #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A_left', 'A_middle', 'A_right'),
                        labels= c(paste(stage, ", A_left", sep=""), 
                                  paste(stage, ", A_middle", sep=""),
                                  paste(stage, ", A_right", sep=""))) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e5, 5e6)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e-4, 5e-3)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Average interaction frequency per 10kb bin', title=paste(stage, 'Interaction frequency vs distances'))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_A_chromosome_ends_vs_center_scaling_50k-5Mb_loess.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
}

######################################################################################
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

write.table(calc_derivatives(a_left[c(1:49), ]), 'autosomes_A_compartment_left_end_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(calc_derivatives(a_right[c(1:49), ]), 'autosomes_A_compartment_right_end_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(calc_derivatives(a_middle[c(1:49), ]), 'autosomes_A_compartment_middle_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(calc_derivatives(b_left[c(1:45), ]), 'autosomes_B_compartment_left_end_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(calc_derivatives(b_right[c(1:45), ]), 'autosomes_B_compartment_right_end_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(calc_derivatives(b_middle[c(1:45), ]), 'autosomes_B_compartment_middle_scaling_derivatives.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

######################################################################################################
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

a_left_derivatives=read.csv('autosomes_A_compartment_left_end_scaling_derivatives.tsv', sep='\t', header = TRUE)
a_right_derivatives=read.csv('autosomes_A_compartment_right_end_scaling_derivatives.tsv', sep='\t', header = TRUE)
a_middle_derivatives=read.csv('autosomes_A_compartment_middle_scaling_derivatives.tsv', sep='\t', header = TRUE)
b_left_derivatives=read.csv('autosomes_B_compartment_left_end_scaling_derivatives.tsv', sep='\t', header = TRUE)
b_right_derivatives=read.csv('autosomes_B_compartment_right_end_scaling_derivatives.tsv', sep='\t', header = TRUE)
b_middle_derivatives=read.csv('autosomes_B_compartment_middle_scaling_derivatives.tsv', sep='\t', header = TRUE)


write.table(loess_smooth(a_left_derivatives), 'autosomes_A_compartment_left_end_scaling_derivatives_loess_smooth.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(loess_smooth(a_right_derivatives), 'autosomes_A_compartment_right_end_scaling_derivatives_loess_smooth.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(loess_smooth(a_middle_derivatives), 'autosomes_A_compartment_middle_scaling_derivatives_loess_smooth.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(loess_smooth(b_left_derivatives), 'autosomes_B_compartment_left_end_scaling_derivatives_loess_smooth.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(loess_smooth(b_right_derivatives), 'autosomes_B_compartment_right_end_scaling_derivatives_loess_smooth.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(loess_smooth(b_middle_derivatives), 'autosomes_B_compartment_middle_scaling_derivatives_loess_smooth.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

###################################################################################################################
#read loess smoothed file
a_left_derivatives=read.csv('autosomes_A_compartment_left_end_scaling_derivatives_loess_smooth.tsv', sep='\t', header = TRUE)
a_right_derivatives=read.csv('autosomes_A_compartment_right_end_scaling_derivatives_loess_smooth.tsv', sep='\t', header = TRUE)
a_middle_derivatives=read.csv('autosomes_A_compartment_middle_scaling_derivatives_loess_smooth.tsv', sep='\t', header = TRUE)
b_left_derivatives=read.csv('autosomes_B_compartment_left_end_scaling_derivatives_loess_smooth.tsv', sep='\t', header = TRUE)
b_right_derivatives=read.csv('autosomes_B_compartment_right_end_scaling_derivatives_loess_smooth.tsv', sep='\t', header = TRUE)
b_middle_derivatives=read.csv('autosomes_B_compartment_middle_scaling_derivatives_loess_smooth.tsv', sep='\t', header = TRUE)






stages=c("L", "Z", "Sun1")

for (stage in stages) {
  
  df<- data.frame(a_left_derivatives$distances)
  colnames(df) <-'distances'
  df['A_left']<- a_left_derivatives[stage]
  df['A_right']<- a_right_derivatives[stage]
  df['A_middle']<- a_middle_derivatives[stage]
  
  cols <- c("A_left" = "orange", "A_middle" = "dodgerblue", "A_right" = "red")
  
  #df['Psline1']=(df['distances']^ -0.6)*15
  #df['Psline2']=(df['distances']^ -1.5)*10000
  
  # sertoli, spg, PreL, L
  p=ggplot(data=df, aes(x = distances))+
    geom_smooth(aes(y= A_left, colour='A_left'), linetype=1, size=1.2, method='loess', se=F)+
    geom_smooth(aes(y= A_right, colour='A_right'), linetype=1, size=1.2, method='loess', se=F)+
    geom_smooth(aes(y= A_middle, colour='A_middle'), linetype=1, size=1.2, method='loess', se=F)+
    #geom_line(aes(y= A, colour='A'), linetype=1, size=1.2)+
    #geom_line(aes(y= B, colour='B'), linetype=1, size=1.2)+
    #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A_left', 'A_middle', 'A_right'),
                        labels= c(paste(stage, ", A_left", sep=""), 
                                  paste(stage, ", A_middle", sep=""),
                                  paste(stage, ", A_right", sep=""))) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e4, 5e6)) +
    scale_y_continuous(limits = c(-1.2, 0.2)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling slope')+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_A_chromosome_ends_vs_center_scaling_slope_50k-3Mb_loess.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
}



for (stage in stages) {
  
  df<- data.frame(b_left_derivatives$distances)
  colnames(df) <-'distances'
  df['B_left']<- b_left_derivatives[stage]
  df['B_right']<- b_right_derivatives[stage]
  df['B_middle']<- b_middle_derivatives[stage]
  
  cols <- c("B_left" = "orange", "B_middle" = "dodgerblue", "B_right" = "red")
  
  #df['Psline1']=(df['distances']^ -0.6)*15
  #df['Psline2']=(df['distances']^ -1.5)*10000
  
  # sertoli, spg, PreL, L
  p=ggplot(data=df, aes(x = distances))+
    geom_smooth(aes(y= B_left, colour='B_left'), linetype=1, size=1.2, method='loess', se=F)+
    geom_smooth(aes(y= B_right, colour='B_right'), linetype=1, size=1.2, method='loess', se=F)+
    geom_smooth(aes(y= B_middle, colour='B_middle'), linetype=1, size=1.2, method='loess', se=F)+
    #geom_line(aes(y= A, colour='A'), linetype=1, size=1.2)+
    #geom_line(aes(y= B, colour='B'), linetype=1, size=1.2)+
    #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('B_left', 'B_middle', 'B_right'),
                        labels= c(paste(stage, ", B_left", sep=""), 
                                  paste(stage, ", B_middle", sep=""),
                                  paste(stage, ", B_right", sep=""))) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e4, 5e6)) +
    scale_y_continuous(limits = c(-1.2, 0.2)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title='Autosomes scaling slope')+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_B_chromosome_ends_vs_center_scaling_slope_50k-3Mb_loess.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
}
