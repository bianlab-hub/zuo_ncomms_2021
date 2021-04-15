library(dplyr)
library(ggplot2)

#TAD boundaries are calculated using the diamond-insulation utility in cooltools and saved as tsv files

#select TAD boundaries using a score > 0.1 cutoff

filenames <- list.files(pattern = '*.tsv')

for (filename in filenames) {
  
  prefix <- strsplit(filename, '_insulation_500kb_window.tsv')[[1]]
  insulation <- read.table(filename, sep='\t', header=TRUE)
  boundaries <- insulation %>% filter(boundary_strength_500000 >0) %>% select (chrom, start, end, boundary_strength_500000)
  write.table(boundaries, paste(prefix, '_TAD_boundary.bedGraph', sep=''), sep='\t', row.names=FALSE, col.names=FALSE,quote=FALSE)
  write.table(boundaries[c(1,2,3)], paste(prefix, '_TAD_boundary.bed', sep=''), sep='\t', row.names=FALSE, col.names=FALSE,quote=FALSE)
  
  boundaries_0.1 <- insulation %>% filter(boundary_strength_500000 >0.1) %>% select (chrom, start, end, boundary_strength_500000)
  write.table(boundaries_0.1, paste(prefix, '_TAD_boundary_0.1.bedGraph', sep=''), sep='\t', row.names=FALSE, col.names=FALSE,quote=FALSE)
  write.table(boundaries_0.1[c(1,2,3)], paste(prefix, '_TAD_boundary_0.1.bed', sep=''), sep='\t', row.names=FALSE, col.names=FALSE,quote=FALSE)
  
  
}


stages <- c('Sertoli', 'spg', 'PreL', 'L', 'Z', 'P', 'D', 'MII')

boundary_strength <- data.frame(stage=c(), strength =c())
for (stage in stages) {
  df <- read.table(paste(stage, '_10kb_TAD_boundary_0.1.bedGraph', sep=''), sep='\t')
  df$stage <- stage 
  df$strength <- df$V4 
  boundary_strength <- rbind(boundary_strength , df[c('stage', 'strength')])
}

#Make boxplots

boundary_strength$stage <- factor(boundary_strength$stage,
                                     levels = c('Sertoli', 'spg', 'PreL', 'L', 'Z', 'P', 'D', 'MII'),ordered = TRUE)

boundary_sum <- boundary_strength %>%
  group_by(stage) %>%
  summarise(median = median(strength), n = n())

p1= ggplot(boundary_strength, aes(x=stage, y=strength, fill=stage)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 2.5)+
  theme_bw()+
  scale_fill_manual(values=c(rep('firebrick1', 2), rep('steelblue3', 6)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave('cooltools_TAD_boundary_strength_202104.pdf', p1, dpi=300, width = 5, height=5)
