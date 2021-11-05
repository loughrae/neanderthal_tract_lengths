library(tidyverse)
library(data.table)

map <- fread('processed_map.txt')
haps <- fread('count_gene_overlaps.bed') 
colnames(haps) <- c('Chromosome', 'Start', 'End', 'ID', 'bp', 'cM', 'row', 'genes', 'distinct_genes')
haps$pop <- substr(haps$ID, start = 10, stop = 12)
haps$grp <- ifelse(haps$pop %in% c('CHB', 'CHS', 'JPT'), 'EAS', 'EUR')

#### Figure 1 #### 
map %>% 
  filter(GeneName %in% c('OAS1', 'OAS2', 'OAS3')) %>% 
  select(pop_gene, AI, bp_pres, cM_pres, Start, End, pop, GeneName, tracts)  %>% 
  separate_rows(tracts, sep = ',') %>% 
  mutate(tracts = as.numeric(tracts)) %>% 
  left_join(haps, by = c('pop', 'tracts' = 'row'))  %>% 
  mutate(tracts = as.factor(tracts)) %>% 
  ggplot() + geom_rect(aes(xmin = Start.x, xmax = End.x, fill = GeneName), ymin = -Inf, ymax = Inf, alpha = 0.08) + geom_segment(aes(x = Start.y, xend = End.y, y = ID, yend = ID)) + ylab('') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab('Position on chromosome 12') + facet_wrap(~pop, nrow = 2, scales = 'free') + guides(fill = guide_legend(override.aes = list(alpha = 1) ) ) + labs(fill = 'Gene')
ggsave(file = 'fig1_oas.png', width = 20, height = 15)


