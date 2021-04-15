library(rtracklayer)
library(GenomicRanges)

stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

for (stage in stages) {
  
  a_domain <-read.table(paste(stage,'_50kb_A_compartment_domains.bed', sep=''), sep='\t', header = FALSE, quote = NULL)
  a_domain <- a_domain[(a_domain$V3 - a_domain$V2)>= 1000000, ]
  write.table(a_domain, paste(stage,'_50kb_A_compartment_domains_larger_than_1Mb.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  
  
  b_domain <-read.table(paste(stage,'_50kb_B_compartment_domains.bed', sep=''), sep='\t', header = FALSE, quote = NULL)
  b_domain <- b_domain[(b_domain$V3 - b_domain$V2)>= 1000000, ]
  write.table(b_domain, paste(stage,'_50kb_B_compartment_domains_larger_than_1Mb.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  
}


###############################################################################################################

#select domains >1Mb at at the right end of chromosome (position index >0.8)
#select domains >1Mb at in the middle of chromosome (position index >0.4 and < 0.6)

chromsizes <- read.table('mm10.chrom.sizes')
chrsizes <- chromsizes$V2
names(chrsizes) <- chromsizes$V1


stages=c("sertoli", "spg", "PreL", "L", "Z", "Sun1", "P", "D", "MII", "RS", "sperm")

for (stage in stages) {
  
  a_domain <-read.table(paste(stage,'_50kb_A_compartment_domains.bed', sep=''), sep='\t', header = FALSE, quote = NULL)
  a_domain_big <- a_domain[(a_domain$V3 - a_domain$V2)>= 1000000, ]
  
  a_domain_big$position_index <- (a_domain_big$V2 + a_domain_big$V3)/2 / chrsizes[as.character(a_domain_big$V1)]
  a_domain_middle <- a_domain_big[a_domain_big$position_index >0.4 & a_domain_big$position_index <0.6, ]
  a_domain_right_end <- a_domain_big[a_domain_big$position_index >0.8, ]
  a_domain_left_end <- a_domain_big[a_domain_big$position_index <0.2, ]
  
  write.table(a_domain_middle, paste(stage,'_50kb_A_compartment_domains_larger_than_1Mb_middle.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  write.table(a_domain_right_end, paste(stage,'_50kb_A_compartment_domains_larger_than_1Mb_right_end.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  write.table(a_domain_left_end, paste(stage,'_50kb_A_compartment_domains_larger_than_1Mb_left_end.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  
  
  b_domain <-read.table(paste(stage,'_50kb_B_compartment_domains.bed', sep=''), sep='\t', header = FALSE, quote = NULL)
  b_domain_big <- b_domain[(b_domain$V3 - b_domain$V2)>= 1000000, ]
  
  b_domain_big$position_index <- (b_domain_big$V2 + b_domain_big$V3)/2 / chrsizes[as.character(b_domain_big$V1)]
  b_domain_middle <- b_domain_big[b_domain_big$position_index >0.4 & b_domain_big$position_index <0.6, ]
  b_domain_right_end <- b_domain_big[b_domain_big$position_index >0.8, ]
  b_domain_left_end <- b_domain_big[b_domain_big$position_index <0.2, ]
  
  write.table(b_domain_middle, paste(stage,'_50kb_B_compartment_domains_larger_than_1Mb_middle.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  write.table(b_domain_right_end, paste(stage,'_50kb_B_compartment_domains_larger_than_1Mb_right_end.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  write.table(b_domain_left_end, paste(stage,'_50kb_B_compartment_domains_larger_than_1Mb_left_end.bed', sep=''), sep='\t', quote = FALSE, row.names = FALSE, col.names =FALSE)
  
}


