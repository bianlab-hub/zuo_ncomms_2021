#######################################################################
bins=unique(round(1.12**c(1:82)))
distances=bins*10000

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII")

scaling_summary<-data.frame(distances)

for (stage in stages) {
  filename=paste(stage, '_combined_autosomes_A_compartment_expected_10kb.tsv', sep='')
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

write.table(scaling_summary, 'Autosomes_A_compartment_expected_vs_distances_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
########################################################################
bins=unique(round(1.12**c(1:82)))
distances=bins*10000

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII")

scaling_summary<-data.frame(distances)

for (stage in stages) {
  filename=paste(stage, '_combined_autosomes_B_compartment_expected_10kb.tsv', sep='')
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

write.table(scaling_summary, 'Autosomes_B_compartment_expected_vs_distances_202005.tsv', sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

############################################################################
#######################################################################


library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

#########################################################################
#Making autosomal scaling plot
#########################################################################
#Plot interaction frequency with loess scaling
a=read.csv('Autosomes_A_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)
b=read.csv('Autosomes_B_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII")

for (stage in stages) {
  
  df<- data.frame(a$distances)
  colnames(df) <-'distances'
  df['A']<- a[stage]
  df['B']<- b[stage]
  
  cols <- c("A" = "red", "B" = "dodgerblue")
  
  #df['Psline1']=(df['distances']^ -0.6)*15
  #df['Psline2']=(df['distances']^ -1.5)*10000
  
  #cols <- c("sertoli" = "red", "spg" = "orange", "PreL" = "darkorchid1", "L" = "lightblue1", "Z"= "cyan", "P"="dodgerblue", "D"="navy", "Sun1"="green", "Psline1"="grey60", "Psline2"="grey40")
  
  # sertoli, spg, PreL, L
  p=ggplot(data=df, aes(x = distances))+
    geom_smooth(aes(y= A, colour='A'), linetype=1, size=1.2, method = "loess")+
    geom_smooth(aes(y= B, colour='B'), linetype=1, size=1.2, method = "loess")+
    #geom_line(aes(y= A, colour='A'), linetype=1, size=1.2)+
    #geom_line(aes(y= B, colour='B'), linetype=1, size=1.2)+
    #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('A', 'B'),
                        labels= c(paste(stage, ", A Compartment", sep=""), paste(stage, ", B Compartment", sep=""))) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e5, 1e7)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e-4, 1e-2)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Average interaction frequency per 10kb bin', title=paste(stage, 'Interaction frequency vs distances'))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_A_B_scaling_loess.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
}


#######################################################################################
# scaling plots from 50kb to 5Mb

a=read.csv('Autosomes_A_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)
b=read.csv('Autosomes_B_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII")

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
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e-5, 5e-3)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Average interaction frequency per 10kb bin', title=paste(stage, 'Interaction frequency vs distances'))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(file=paste(stage, '_autosomes_A_B_scaling_50k-5Mb.pdf', sep=''), width = 10, height = 5)
  print(p)
  dev.off()
  
}

######################################################################################

# Make plots across stages
plot_scaling <- function(df, plot_prefix){
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
                  limits = c(1e4, 5e6)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e-5, 1e-1)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Contact Probability', title='Autosomes scaling plot')+
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
                  limits = c(1e4, 5e6)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e-5, 1e-1)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Contact Probability', title='Autosomes scaling plot')+
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
                  limits = c(1e4, 5e6)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(5e-5, 1e-1)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Contact Probability', title='Autosomes scaling plot')+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  pdf(paste(plot_prefix, '_scaling_SSPL.pdf', sep=''), width = 10, height = 10)
  print(p)
  dev.off()
  
  pdf(paste(plot_prefix, '_scaling_LZPD.pdf', sep=''), width = 10, height = 10)
  print(q)
  dev.off()
  
  pdf(paste(plot_prefix, '_scaling_Sun1_WT_Z.pdf', sep=''), width = 10, height = 10)
  print(s)
  dev.off()
  
}


a=read.csv('Autosomes_A_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)
b=read.csv('Autosomes_B_compartment_expected_vs_distances_202005.tsv', sep='\t', header = TRUE)


plot_scaling(a, 'Autosomes_A_compartment')
plot_scaling(b, 'Autosomes_B_compartment')
