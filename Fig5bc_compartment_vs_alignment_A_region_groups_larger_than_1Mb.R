library(dplyr)

# get A/B compartment annotation

compartment <- read.csv('Z_combined_50kb.cis.vecs.tsv', sep='\t', header = TRUE)

compartment <- compartment[!is.na(compartment$E1), ]

compartment$comp <- 'A'
compartment$comp[compartment$E1<0] <- 'B'
compartment$binname <- paste(compartment$chrom, compartment$start, compartment$end, sep='_')

#get A1 A2 A3 B groups

df <- read.table('Z_compartment_region_500kb_interaction.tsv', header=TRUE)

df$A_region_groups <- as.character(df$compartment)
breaks <- as.numeric(quantile(df$balanced.average[df$compartment=='A'], probs = seq(0, 1, 0.33)))
df$A_region_groups[df$compartment=='A'] <- 'A2'
df$A_region_groups[df$compartment=='A' & df$balanced.average> breaks[3] ] <- 'A1'
df$A_region_groups[df$compartment=='A' & df$balanced.average< breaks[2] ] <- 'A3'
df$A_region_groups = factor(df$A_region_groups, levels = c('A1', 'A2', 'A3', 'B'))

A1 <- df %>% filter(A_region_groups=='A1') %>% select(chrom, start, end)
A2 <- df %>% filter(A_region_groups=='A2') %>% select(chrom, start, end)
A3 <- df %>% filter(A_region_groups=='A3') %>% select(chrom, start, end)
B <- df %>% filter(A_region_groups=='B') %>% select(chrom, start, end)


table(df$A_region_groups)

get_rownumbers_domains <- function(df, domains) {
  rownumbers <- c()
  for ( i in 1: nrow(domains)) {
    chri = as.character(domains$chrom[i])
    starti = domains$start [i]
    endi = domains$end [i]
    rownumbers <- c(rownumbers, which(df$chr== chri & df$start >= starti & df$end <= endi))
    
  }
  
  return(unique(rownumbers))
}
#########################################################################################################
windowsize='100kb'
#
alignment <- read.csv(paste('Z_interhomolog_alignment_ratio_halfwindow_',windowsize, '.tsv', sep=''), sep='\t', header = TRUE)
alignment$chr <- paste('chr',alignment$chr, sep='' )
alignment <- alignment[!is.na(alignment$alignmentscore), c(2,3,4,5)]

#summarize alignment scores
alignment$binname <- paste(alignment$chr, alignment$start, alignment$end, sep='_')

alignment <- merge(alignment, compartment[, c(8,9)], by='binname')
alignment <- alignment[, c(2:6)]

chromsizes <- read.table('mm10.chrom.sizes')
chrsizes <- chromsizes$V2
names(chrsizes) <- chromsizes$V1

alignment$chrlen <- chrsizes[alignment$chr]

alignment$dist_to_left <- alignment$start
alignment$dist_to_right <- alignment$chrlen - alignment$end

alignment$position_index <- (alignment$start + alignment$end)/2 / alignment$chrlen

alignment <- alignment[order(alignment$chr, alignment$start), ]

alignment$position_index_bin <- round(alignment$position_index / 0.01) * 0.01


#Fig4b
# Make boxplots A vs B

pdf(file=paste('alignment_score_',windowsize,'_A_B_202104.pdf', sep=''))
scores <- boxplot(alignment$alignmentscore[alignment$comp=='A'],alignment$alignmentscore[alignment$comp=='B'],
                  col=c('red','dodgerblue'),
                  outline=F, names=c('A','B'),ylab='alignment score',
                  main="zygotene stage homolog alignment score")
text(x=1:2, y=scores$stats[4,]+1, labels=paste("n=", scores$n))
text(x=1:2, y=scores$stats[3,]-1, labels=round(scores$stats[3,],digits=3))
dev.off()

wilcox.test(alignment$alignmentscore[alignment$comp=='A'], 
            alignment$alignmentscore[alignment$comp=='B'])$p.value

#Fig4c
# Make boxplots A compartment groups

pdf(file=paste('alignment_score_',windowsize,'_A_groups_202104.pdf', sep=''))
scores <- boxplot(alignment$alignmentscore[get_rownumbers_domains(alignment, A1)],
                  alignment$alignmentscore[get_rownumbers_domains(alignment, A2)],
                  alignment$alignmentscore[get_rownumbers_domains(alignment, A3)],
                  alignment$alignmentscore[get_rownumbers_domains(alignment, B)],
                  col=c('red', 'orange','orange','dodgerblue'),
                  outline=F, names=c('A1','A2','A3','B'),ylab='alignment score',
                  main="zygotene stage homolog alignment score")
text(x=1:4, y=scores$stats[4,]+1, labels=paste("n=", scores$n))
text(x=1:4, y=scores$stats[3,]-1, labels=round(scores$stats[3,],digits=3))
dev.off()

wilcox.test(alignment$alignmentscore[get_rownumbers_domains(alignment, A1)], 
            alignment$alignmentscore[get_rownumbers_domains(alignment, A2)])$p.value

#[1] 2.540057e-16

wilcox.test(alignment$alignmentscore[get_rownumbers_domains(alignment, A2)], 
            alignment$alignmentscore[get_rownumbers_domains(alignment, A3)])$p.value
#[1] 6.122008e-07

wilcox.test(alignment$alignmentscore[get_rownumbers_domains(alignment, A3)], 
            alignment$alignmentscore[get_rownumbers_domains(alignment, B)])$p.value

#[1] 1.136943e-254

