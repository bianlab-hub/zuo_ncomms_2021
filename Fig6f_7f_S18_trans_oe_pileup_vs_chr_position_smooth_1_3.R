##############################################################################################################
#Trans OE
#Plot loess smoothed diagnol
#slide a 5x5 window along diagnonal of the 500x500 trans OE average heatmaps to get diagnol values
diag_sliding_window <- function(matrix, halfwindow) {
  s <- dim(matrix)[1]
  line <- c()
  
  for (i in c(1:s)) {
    left <- max(0, i-halfwindow)
    right <- min (i+halfwindow, s)
    line <- c(line, mean(matrix[c(left:right), c(left:right)]))
    
  }
  return(line)
}

#read a trans obs matrix and get smoothed diagnol

trans_oe_diag <- function(stage, halfwindow=5) {
  
  #read trans_obs matrix for a given stage
  trans_oe <- read.table(paste(stage, '_average_trans_oe_500x500.matrix', sep=''), sep='\t', header=FALSE)
  trans_oe <- as.matrix(trans_oe)
  
  #loess smooth a series

  diag <- diag_sliding_window(trans_oe, halfwindow)
  smoothed <- loess.smooth(c(20:500)/500, diag[20:500], span = 1/3, degree = 1, family ="gaussian", evaluation = 100)
  return (smoothed)
  
}

preL <- trans_oe_diag('PreL')
L <- trans_oe_diag('L')
Z <- trans_oe_diag('Z')
P <- trans_oe_diag('P')
D <- trans_oe_diag('D')
sun1 <- trans_obs_diag('Sun1')

pdf('L_Z_trans_oe_pileup_alignment.pdf', width=10, height=6)
plot(c(1:500)/500, c(1:500)/500, col='white', ylab='trans obs/exp',xlab='position along chromosome',
     main='inter-chromosomal alignment', type='l',cex.lab=0.8,cex.main=1.2,cex.axis=0.8,
     xlim=c(0,1), ylim <- c(0.5,3))
lines(L$x, L$y,type='l', lwd=3, col='orange', xlim=c(0,1))
#abline(h=median(L$y), lty=2, lwd=2, col='orange')
lines(Z$x, Z$y,type='l', lwd=3, col='red', xlim=c(0,1))
#abline(h=median(Z$y),  lty=2, lwd=2, col='red')
#lines(sun1$x, sun1$y,type='l', lwd=3, col='black', xlim=c(0,1))
legend('top', legend=c('L', 'Z') ,col=c('red', 'orange') , lwd=2 , bty='n', cex=1.5)
dev.off()

pdf('WT_vs_Sun1_Z_trans_oe_pileup_alignment.pdf', width=10, height=6)
plot(c(1:500)/500, c(1:500)/500, col='white', ylab='trans obs/exp',xlab='position along chromosome',
     main='inter-chromosomal alignment', type='l',cex.lab=0.8,cex.main=1.2,cex.axis=0.8,
     xlim=c(0,1), ylim <- c(0.5,3))
lines(Z$x, Z$y,type='l', lwd=3, col='red', xlim=c(0,1))
lines(sun1$x, sun1$y,type='l', lwd=3, col='dodgerblue', xlim=c(0,1))
legend('top', legend=c('WT', 'Sun1 Mutant') ,col=c('red', 'dodgerblue') , lwd=2 , bty='n', cex=1.5)
dev.off()


stages= c('sertoli', 'spg', 'PreL', 'L', 'Z', 'P', 'D', 'Sun1', 'MII')


for (stage in stages) {
  plotdata <- preL <- trans_oe_diag(stage)

  pdf(paste(stage, '_trans_oe_pileup_alignment.pdf', sep=''), width=10, height=6)
  plot(c(1:500)/500, c(1:500)/500, col='white', ylab='trans obs/exp',xlab='position along chromosome',
       main='inter-chromosomal alignment', type='l',cex.lab=0.8,cex.main=1.2,cex.axis=0.8,
       xlim=c(0,1), ylim <- c(0.5,3))
  lines(plotdata$x, plotdata$y,type='l', lwd=3, col='red', xlim=c(0,1))
  legend('top', legend=stage ,col='red', lwd=2 , bty='n', cex=1.5)
  dev.off()
  
}
