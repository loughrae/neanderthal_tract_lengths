library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpubr)
library(D3GB)

#### Read in three recombination files, haplotype file and gene_map file ###
map <- fread('processed_map.txt')
haps <- fread('count_gene_overlaps.bed')
colnames(haps) <- c('Chromosome', 'Start', 'End', 'ID', 'bp', 'cM', 'row', 'genes', 'distinct_genes')
haps$pop <- substr(haps$ID, start = 10, stop = 12)
haps$grp <- ifelse(haps$pop %in% c('CHB', 'CHS', 'JPT'), 'EAS', 'EUR')

win <- fread('windows_rec-all.bed')

hr <- fread('haps_rec-all.bed') %>% 
  as.data.frame() %>% 
  separate(V10, into = c(NA, NA, 'pop'), sep = '_') %>% 
  separate(pop, into = c('pop', NA), sep = '\\.')

gr <- fread('genes_rec-all.bed') %>% 
  as.data.frame() %>% 
  separate(V9, into = c(NA, NA, 'pop'), sep = c('_')) %>% 
  separate(pop, into = c('pop', NA), sep = '\\.')
win$par <- ifelse(win$V2 >= 60001 & win$V3 <= 2699520 & win$V1 == 'chrX', 'PAR1', ifelse(win$V2 >= 154931044 & win$V3 <= 155260560 & win$V1 == 'chrX', 'PAR2', 'Neither'))

p1 <- haps %>% 
  mutate(chr = gsub('chr', '', Chromosome)) %>% 
  count(grp, chr) %>% 
  left_join(GRCh37, by = c('chr')) %>% 
  mutate(chr = ifelse(chr == 'X', 23, chr)) %>% 
  mutate(chr = as.numeric(chr)) %>% 
  mutate(tracts_per_mb = n*1000000 / end) %>% 
  ggplot(aes(x = chr, y = tracts_per_mb, color = chr == 23)) + geom_point() + facet_wrap(~grp) + ylab('Nea tracts per MB') + xlab('Chromosome') + theme(legend.position = 'none') + ggtitle('A')

p2 <- map %>% 
  mutate(autosome = ifelse(Chromosome == 'chrX', 'X', 'Autosome')) %>% 
  ggplot(aes(x = autosome, fill = n_uniq >= 1)) + geom_bar(position = 'fill') + facet_wrap(~pop, nrow = 2) + ylab('Proportion of genes introgressed') + xlab('Chromosome') + theme(legend.position = 'none') + ggtitle('B')

map %>% group_by(Chromosome == 'chrX') %>% summarize(sum(n_uniq >= 1)/n())

p3 <- haps %>% 
  mutate(autosome = ifelse(Chromosome == 'chrX', 'X', 'Autosome')) %>% 
  ggplot(aes(x = autosome, y = bp, fill = autosome)) + geom_boxplot() + stat_compare_means() + ylab('Tract Length (log10)') + scale_y_continuous(trans = 'log10') + xlab('Chromosome') + ggtitle('C') + theme(legend.position = 'none')

s1 <- haps %>% ggplot(aes(x = genes > 0, y = bp, fill = genes > 0)) + geom_boxplot() + stat_compare_means() + scale_y_continuous(trans = 'log10') + ylab('Tract Length') + xlab('Tract Overlaps Genes') + ggtitle('A') + theme(legend.position = 'none')
s2 <- win %>% filter(V5 > 0) %>% ggplot(aes(x = V4 > 0, y = V7, fill = V4 > 0)) + geom_boxplot() + stat_compare_means() + scale_y_continuous(trans = 'log10') + ylab('Median length of tracts overlapping window') + xlab('Window Genic?') + ggtitle('B') + theme(legend.position = 'none')
s3 <- win %>% ggplot(aes(x = V4 > 0, fill = V5 > 0)) + geom_bar(position = 'fill') + xlab('Window Genic?') + ylab('Proportion of windows with a tract') + theme(legend.position = 'none') + ggtitle('C')

p4 = hr %>% 
  filter(V8 != -999) %>% 
  ggplot(aes(x = V8, y = V5)) + geom_point() + geom_smooth() + stat_cor(aes(label = ..r.label..)) + xlab('Recombination rate per base per generation') + ylab('Nea tract length') + ggtitle('D')

p5 = map %>% 
  left_join(gr, by = c('pop', 'Gene' = 'V4')) %>% 
  filter(n_uniq >= 3, V7 != -999) %>% 
  ggplot(aes(x = V7, y = bp_pres)) + geom_point() + geom_smooth() + stat_cor(aes(label = ..r.label..), hjust = -1) + xlab('Recombination rate per base per generation') + ylab('Mean overlapping tract length') + ggtitle('E')

p6 <- win %>% 
  filter(V9 != 'fuck', V5 > 0) %>% 
  ggplot(aes(x = as.numeric(V9), y = V7)) + geom_point() + stat_cor(aes(label = ..r.label..), hjust = -2) + xlab('Recombination rate per base per generation') + ylab('Median length of Nea tracts in window') + ggtitle('F')

s4 <- win %>% 
  filter(V9 != 'fuck') %>% 
  ggplot(aes(y = as.numeric(V9), x = V5 > 0, fill = V5 > 0)) + geom_boxplot() + stat_compare_means() + ylab('Recombination rate per base per generation, log10') + xlab('Tract present in window') + ggtitle('D') + scale_y_continuous(trans = 'log10') + theme(legend.position = 'none')


fig2 <- (p1 | p2 | p3) / (p4 | p5 | p6)
ggsave(fig2, file = 'fig3emilia.png', width = 20, height = 15)

supp2 <- (s1 | s2) / (s3 | s4)
ggsave(supp2, file = 'supp3emilia.png', width = 20, height = 15)



