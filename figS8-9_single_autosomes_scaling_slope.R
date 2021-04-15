for (i in paste("chr", as.character(c(1:19)),sep = "")){
#for (i in paste("chr", 'X',sep = "")){
  filename=paste('autosomes_scaling_summary_keep_diag1_',i,'.tsv',sep = "")
  df=read.csv(filename, sep='\t', header = TRUE)
  df= df[-1,]
  #calculate log10 of distances and interaction frequencies
  lgdf<- log10(df)
  colnames(lgdf)[1] <- 'log_dist'
  lgdf$distances <- df$distances
  
  #calculate 1st order derivative at each distance
  lgdf<- lgdf[1:length(df$distances)-1,]
  derivatives<- lgdf
  stages<- c('Sertoli','MII','spg','PreL','L','Z','P','D')
  
  for (stage in stages) {
    x<- lgdf$log_dist
    y<- lgdf[,stage]
    d<- predict(sm.spline(x, y), x, 1)
    derivatives[,stage] <- d
  }
  
  write.table(derivatives, paste(i,'_autosomes_scaling_derivatives.tsv',sep = ""), sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  df_smooth=derivatives
  
  for (n in c(2:9)) {
    df_smooth[,n] <- smooth(df_smooth[,n]) 
  }
  
  write.table(df_smooth, paste(i,'_autosomes_scaling_slopes.tsv',sep = ""), sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #loess smoothen
  df=derivatives
  
  Sertoli=loess.smooth(df_smooth$log_dist, df_smooth$Sertoli, span = 1/3, degree = 1,
                       family ="gaussian", evaluation = 77)
  
  df_loess <- data.frame(x=Sertoli$x)
  df_loess$distances <- round(10^df_loess$x)
  
  for (stage in stages) {
    x<- df$log_dist
    y<- df[,stage]
    loessfit <- loess.smooth(x, y, span = 1/3, degree = 1, family ="gaussian", evaluation = 77)
    df_loess[,stage] <- loessfit$y
  }
  
  write.table(df_loess, paste(i,'_autosomes_scaling_slopes_smoothened.tsv'), sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}


for (i in paste("chr", as.character(19),sep = "")){
#for (i in paste("chr",'X',sep = "")){
  filename=paste(i,' _autosomes_scaling_slopes_smoothened.tsv',sep = "")
  df=read.csv(filename, sep='\t', header = TRUE)
  p=ggplot(data=df, aes(x = distances))+
    geom_line(aes(y= Sertoli, colour='Sertoli'), linetype=1, size=1.2)+
    geom_line(aes(y= spg, colour='spg'), linetype=1, size=1.2)+
    geom_line(aes(y= PreL, colour='PreL'), linetype=1, size=1.2)+
    geom_line(aes(y= L, colour='L'), linetype=1, size=1.2)+
    geom_line(aes(y= MII, colour='MII'), linetype=1, size=1.2)+
    #geom_line(aes(y= Psline2, colour='Psline2'), linetype=5, size=1.5)+
    scale_colour_manual("", 
                        values = cols, 
                        breaks= c('Sertoli', 'spg', 'PreL', 'L','MII'),
                        labels= c('Sertoli', 'spg', 'PreL', 'L','MII')) +
    #geom_vline(xintercept = 1e5, linetype=2, size=1)+
    #geom_vline(xintercept = 2e6, linetype=2, size=1)+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1e4, 1e8)) +
    scale_y_continuous(limits = c(-2, -0.2)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title=paste(i,'_Autosomes scaling plot',sep = ""))+
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
    scale_y_continuous(limits = c(-2, -0.2)) +
    theme_bw()+
    labs(x='Separation (bp)', y='scaling slope', title=paste(i,'_Autosomes scaling plot',sep = ""))+
    theme(legend.position = 'right',
          legend.text = element_text(size=20),
          plot.title = element_text(size = 20, hjust = 0.5), 
          axis.text.x.bottom = element_text(size=15),
          axis.text.y.left = element_text(size=15),
          axis.title.x = element_text(size=15), 
          axis.title.y = element_text(size=15))
    pdf(file=paste(i,'_sacling_slop_smoothed_1.pdf',sep = ""), width = 10, height = 10)
    p
    dev.off()

    pdf(file=paste(i,'_sacling_slop_smoothed_2.pdf',sep = ""), width = 10, height = 10)
    q
    dev.off()
}


