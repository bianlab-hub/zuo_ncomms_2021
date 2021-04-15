library(dplyr)
library(pspline)
library(ggplot2)
library(scales)

filename <- 'Z_A_compartment_each_region_expected_10kb.tsv'

# A/ B compartment regions larger than 1 Mb

average_interaction_distance <- function(filename, dist) {
  expected <- read.table(filename, sep='\t',header=TRUE)
  expected <- expected %>% mutate (region =paste(chrom, start, end, sep='_'), distance = diag*10000, balanced.average = balanced.sum/ n_valid) %>% 
    select (chrom, start, end, region, distance, balanced.average)
  balanced_average <- expected %>% filter(distance == dist) %>% select(chrom, start, end, balanced.average)
  return(balanced_average)
}

L_A_500kb <- average_interaction_distance('L_A_compartment_each_region_expected_10kb.tsv', 500000)
L_A_500kb$compartment <- 'A'
L_B_500kb <- average_interaction_distance('L_B_compartment_each_region_expected_10kb.tsv', 500000)
L_B_500kb$compartment <- 'B'
L_500kb <- rbind(L_A_500kb, L_B_500kb)
write.table(L_500kb, 'L_compartment_region_500kb_interaction.bedGraph', sep='\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(L_500kb, 'L_compartment_region_500kb_interaction.tsv', sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)

Z_A_500kb <- average_interaction_distance('Z_A_compartment_each_region_expected_10kb.tsv', 500000)
Z_A_500kb$compartment <- 'A'
Z_B_500kb <- average_interaction_distance('Z_B_compartment_each_region_expected_10kb.tsv', 500000)
Z_B_500kb$compartment <- 'B'
Z_500kb <- rbind(Z_A_500kb, Z_B_500kb)
write.table(Z_500kb, 'Z_compartment_region_500kb_interaction.bedGraph', sep='\t', row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(Z_500kb, 'Z_compartment_region_500kb_interaction.tsv', sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)


################################################################################

Z_500kb <- read.table('Z_compartment_region_500kb_interaction.bedGraph', sep='\t', header=TRUE)

Z_500kb$balanced.average.normalized <- Z_500kb$balanced.average/ median(Z_500kb$balanced.average)

p= ggplot(Z_500kb, aes(x=compartment, y=balanced.average.normalized, fill=compartment)) + 
  geom_boxplot() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")

