chroms=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX')

stages=c('Sertoli','MII','spg','PreL','L','Z','P','D')

chroms=c('chr1')#87
chroms=c('chr2','chrX')#86
chroms=c('chr3','chr4')#85
chroms=c('chr5','chr6','chr7')#84
chroms=c('chr8','chr9','chr10','chr11','chr14')#83
chroms=c('chr12','chr13')#82
chroms=c('chr15','chr16')#81
chroms=c('chr17','chr18')#80
chroms=c('chr19')#76

for (chrom in chroms){
  bins=unique(round(1.12**c(1:76)))
  distances=bins*10000
  scaling_summary<-data.frame(distances)
  for (stage in stages){
    filename=paste(stage,'_',chrom,'_interaction_versus_distance_keep_diag.tsv', sep='')
    expected_df=read.csv(filename, sep='\t', header = TRUE)
    frequency=expected_df$frequency[expected_df$diag %in% bins]
    scaling_summary[stage] <- frequency
  }
  filename1=paste('autosomes_scaling_summary_keep_diag1_',chrom,'.tsv', sep='')
  write.table(scaling_summary, filename1, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}

## plot
setwd("~/Documents/project/meiosis/scaling/")

library(ggplot2)
library(scales)
#install.packages("extrafont")
library(extrafont)
font_import()

cols <- c("Sertoli" = "red", "spg" = "orange", "PreL" = "darkorchid1", "L" = "lightblue1", "Z"= "cyan", "P"="dodgerblue", "D"="navy", "MII"="green", "Psline1"="grey60", "Psline2"="grey40")

for (i in paste("chr", as.character(19),sep = "")){
  filename=paste('autosomes_scaling_summary_keep_diag1_',i,'.tsv',sep = "")
  df=read.csv(filename, sep='\t', header = TRUE)
  df= df[-1,]
  for (n in c(2:9)) {
    df[,n] <- df[,n] * 0.01 / df[df$distances==50000, n] 
  }
  df['Psline1']=(df['distances']^ -0.6)*15
  p=ggplot(data=df, aes(x = distances))+
    geom_line(aes(y= Sertoli, colour='Sertoli'), linetype=1, size=1.2)+
    geom_line(aes(y= spg, colour='spg'), linetype=1, size=1.2)+
    geom_line(aes(y= PreL, colour='PreL'), linetype=1, size=1.2)+
    geom_line(aes(y= MII, colour='MII'), linetype=1, size=1.2)+
    geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
    geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('Sertoli', 'spg', 'PreL', 'L','MII', 'Psline1'),
                        labels= c('Sertoli', 'spg', 'PreL', 'L','MII', 'P~s^(-0.6)')) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e4, 1e8)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e-5, 1e-1)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Contact Probability', title=paste(i,'_scaling plot',sep = ''))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  
  q=ggplot(data=df, aes(x = distances))+
    geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
    geom_line(aes(y= Z, colour='Z'), linetype=1, size=1.2)+
    geom_line(aes(y= P, colour='P'), linetype=1, size=1.2)+
    geom_line(aes(y= D, colour='D'), linetype=1, size=1.2)+
    geom_line(aes(y= Psline1, colour='Psline1'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('L', 'Z', 'P', 'D', 'Psline1'),
                        labels= c('L', 'Z', 'P', 'D', 'P~s^(-0.6)')) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e4, 1e8)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e-5, 1e-1)) +
    theme_bw()+
    labs(x='Separation (bp)', y='Contact Probability', title=paste(i,'_scaling plot',sep = ''))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
  pdf(file=paste(i,'_sacling_1.pdf',sep = ""), width = 10, height = 10)
  p
  dev.off()

  pdf(file=paste(i,'_sacling_2.pdf',sep = ""), width = 10, height = 10)
  q
  dev.off()
}

