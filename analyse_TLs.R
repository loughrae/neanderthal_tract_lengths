library(tidyverse)
library(data.table)
library(ggpubr)

ids <- fread('ids/ids_auto.txt', header = T) %>% as.data.frame() %>% mutate(match_chr = ifelse(chr == '.chr-X.ids', 'chrX', 'chrOther'))
ai_genes <- read.table('AI_genes.txt') %>% pull(V1) 


maps <- fread('mapped_tracts-all.bed', header = F, stringsAsFactors = F) %>% as.data.frame() 
colnames(maps) <- c('Chromosome', 'Start', 'End', 'Gene', 'GeneName', 'Biotype', 'bp_sum', 'cM_sum', 'n_tracts', 'n_uniq', 'haploids', 'distinct_haploids', 'tracts', 'ID')

mop <- maps %>% mutate(ID = gsub("mapped_tracts_", "", ID)) %>% mutate(pop = gsub(".bed", "", ID)) %>% select(-ID) %>% mutate(grp = ifelse(pop %in% c('CHB', 'CHS', 'JPT'), 'EAS', 'EUR')) %>%
	 mutate(match_chr = ifelse(Chromosome == 'chrX', 'chrX', 'chrOther')) %>% left_join(ids, by = c('match_chr', 'pop')) %>%
	 mutate(bp_all = bp_sum / num, bp_pres = bp_sum / n_uniq, cM_pres = cM_sum / n_uniq) %>%
	 mutate(pop_gene = paste(pop, Gene, sep = '_')) %>%
	 mutate(AI = ifelse(pop_gene %in% ai_genes, 'AI', 'Not AI')) 

mop %>% count(pop, n_uniq >= 3) %>% ggplot(aes(x = pop, fill = `n_uniq >= 3`, y = n)) + geom_col() 
ggsave('included_genes_by_pop.png', width = 12, height = 10)
	
mop %>% fwrite(file = 'processed_map.txt', col.names = T, row.names = F, quote = F, sep = '\t')

print(table(mop$AI, mop$Biotype))




