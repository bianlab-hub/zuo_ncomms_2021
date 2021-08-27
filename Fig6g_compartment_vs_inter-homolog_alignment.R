# get A/B compartment annotation for preL

compartment <- read.csv('PreL_combined_50kb.cis.vecs.tsv', sep='\t', header = TRUE)

compartment <- compartment[!is.na(compartment$E1), ]

compartment$comp <- 'A'
compartment$comp[compartment$E1<0] <- 'B'

e1 <- compartment[, c(1,2,3,5)]

#
alignment <- read.csv('Z_interhomolog_alignment_ratio_halfwindow_500kb.tsv', sep='\t', header = TRUE)
alignment$chr <- paste('chr',alignment$chr, sep='' )
alignment <- alignment[!is.na(alignment$alignmentscore), c(2,3,4,5)]

compartment$binname <- paste(compartment$chrom, compartment$start, compartment$end, sep='_')

#summarize alignment scores
alignment$binname <- paste(alignment$chr, alignment$start, alignment$end, sep='_')

alignment <- merge(alignment, compartment[, c(8,9)], by='binname')
alignment <- alignment[, c(2:6)]

chromsizes <- read.table('mm10.chrom.sizes')
chrsizes <- chromsizes$V2
names(chrsizes) <- chromsizes$V1

alignment$chrlen <- chrsizes[alignment$chr]

alignment$position_index <- (alignment$start + alignment$end)/2 / alignment$chrlen

alignment <- alignment[order(alignment$chr, alignment$start), ]

alignment$position_index_bin <- round(alignment$position_index / 0.01) * 0.01


align_a <- alignment[alignment$comp == 'A', ]
align_b <- alignment[alignment$comp == 'B', ]


align_position <- data.frame(position_index=c(0:100)*0.01)
align_position$A_compartment <- NA
align_position$B_compartment <- NA

for ( n in c(0:100)) {
  
  align_position$A_compartment[n] <- median(align_a$alignmentscore[align_a$position_index_bin == n*0.01 ], na.rm=TRUE)
  align_position$B_compartment[n] <- median(align_b$alignmentscore[align_b$position_index_bin == n*0.01 ], na.rm=TRUE)
  
}


smoothen<-function(list, windowsize=10){
  halfwidth=round((windowsize-1)/2)
  newlist=c()
  for (i in c(1:length(list))) {
    x=median(list[max(0, i-halfwidth) : min(i+halfwidth, length(list))], na.rm = TRUE)
    newlist<-c(newlist, x)
  }
  return(newlist)
}


loess_smooth <- function(df) {
  sm <-loess.smooth(df[,1], df[,2], span = 1/3, degree = 1, family ="gaussian", evaluation = 100)
  return(sm)
}


pdf('interhomolog alignment A vs B compartment position index loess smooth.pdf', width=10, height=8)
plot(align_position$position_index,align_position$A_compartment, col='white', ylab='alignmentscore',xlab='position along chromosome',
     main='interhomolog alignment vs position along chromosome', type='l',cex.lab=0.8,cex.main=1.2,cex.axis=0.8, ylim=c(0,30))
lines(loess_smooth(align_position[, c(1,2)])$x,loess_smooth(align_position[, c(1,2)])$y, type='l', lwd=3, col='red')
lines(loess_smooth(align_position[, c(1,3)])$x,loess_smooth(align_position[, c(1,3)])$y, type='l', lwd=3, col='dodgerblue')
legend('topright', legend=c('A compartment', 'B compartment') ,col=c('red', 'dodgerblue') , lwd=2 , bty='n', cex=1.5)
dev.off()

