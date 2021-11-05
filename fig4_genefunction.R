set.seed(123)

library(gprofiler2)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(data.table)
library(ggrepel)
library(patchwork)
library(ggpubr)

gr <- fread('genes_rec-all.bed') %>% rename(recrate = V8) %>% filter(recrate != -999) %>% mutate(pop = substr(V9, start = 11, stop = 13)) %>% mutate(pop_gene = paste(pop, V4, sep = '_'))
map <- fread('processed_map.txt') %>% mutate(len = End - Start) %>% left_join(gr, by = c('pop_gene', 'pop')) %>% filter(n_uniq >= 3)

#### Create set of control non-AI genes with similar recombination and gene length distributions to the AI genes ####
lens = map %>% 
  filter(AI == 'AI') %>%
  mutate(len_quants = cut(len, 10)) %>% 
  pull(len_quants) %>%
  unique() %>%
  gsub('\\(', '', .) %>% 
  gsub('\\]', '', .)
  
recs = map %>% 
  filter(AI == 'AI') %>%
  filter(!is.na(recrate)) %>% #there's an AI gene on chrX and our recombination rate dataset doesn't include chrX
  mutate(rec_quants = cut(recrate, 10)) %>% 
  pull(rec_quants) %>%
  unique() %>%
  gsub('\\(', '', .) %>% 
  gsub('\\]', '', .)

quants = expand.grid(len = lens, rec = recs) %>% #make combinations of all bins (4x8 = 32 rows)
  separate(len, into = c('min_len', 'max_len'), sep = ',') %>% 
  separate(rec, into = c('min_rec', 'max_rec'), sep = ',') %>% 
  mutate_all(function(x) as.numeric(as.character(x)))

lis <- vector(mode = 'list')
for (row in 1:nrow(quants)) {
  qq = quants[row,]
  n_ai <- map %>% filter(AI == 'AI') %>% filter(len > qq[,1] & len <= qq[,2] & recrate > qq[,3] & recrate <= qq[,4]) %>% nrow()
  print(n_ai)
  cont = map %>% filter(AI == 'Not AI') %>% filter(len > qq[,1] & len <= qq[,2] & recrate > qq[,3] & recrate <= qq[,4])
  if (nrow(cont) > 5*n_ai) {
    cont <- sample_n(cont, 5*n_ai)
  }
  lis[[row]] <- cont
  print(row)
}

controlled <- lis %>% bind_rows() %>% rbind(map[map$AI == 'AI'])

#verify that the recombination rate and gene length distributions look similar
s1 <- controlled %>% ggplot(aes(x = recrate, fill = AI)) + geom_density(alpha = 0.3) + facet_wrap(~AI, scales = 'free_y') + xlab('Recombination rate per base per generation') + theme(legend.position = 'none') + ggtitle('A')
s2 <- controlled %>% ggplot(aes(x = len, fill = AI)) + geom_density(alpha = 0.3) + facet_wrap(~AI, scales = 'free_y') + xlab('Gene Length') + theme(legend.position = 'none') + ggtitle('B')

#compare AI genes to control non-AI genes by mean tract length (bp)
p1 <- controlled %>% ggplot(aes(x = AI, y = bp_pres, fill = AI)) + geom_boxplot() + scale_y_continuous(trans = 'log10') + ylab('Mean tract length (bp, log10)') + xlab('') + stat_compare_means() + ggtitle('A') + theme(legend.position = 'none')
#same but cM instead of bp
s3 <- controlled %>% ggplot(aes(x = AI, y = cM_pres)) + geom_boxplot() + scale_y_continuous(trans = 'log10') + ylab('Mean tract length (cM, log10)') + xlab('') + stat_compare_means() + ggtitle('C')

map <- fread('processed_map.txt') #note refreshed map

#compare AI genes to all other genes by mean tract length
p2 <- map %>% filter(n_uniq >= 3) %>% ggplot(aes(x = AI, y = bp_pres, fill = AI)) + geom_boxplot() + stat_compare_means() + ggtitle('B') + ylab('Mean tract length (bp, log10)') + scale_y_continuous(trans = 'log10') + theme(legend.position = 'none') 
#compare pseudogenes to coding genes by mean tract length
p3 <- map %>% 
 filter(n_uniq >= 3) %>% 
 filter(Biotype %in% c('pseudogene', 'protein_coding')) %>% 
 ggplot(aes(x = pop, fill = Biotype, y = bp_pres)) + geom_boxplot() + scale_y_continuous(trans = 'log10') + xlab('Population') + ylab('Mean tract length (bp, log10)')

#### GO Enrichment Analysis ####
#gprofiler for GO, ranked:
gost(query = c(map %>% filter(grp == 'EAS', n_uniq >= 3) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = T, significant = TRUE, sources = 'GO:BP')[[1]]
gost(query = c(map %>% filter(grp == 'EUR', n_uniq >= 3) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = T, significant = TRUE, sources = 'GO:BP')[[1]]

#gprofiler for GO, overrepresentation:
gost(query = c(map %>% filter(grp == 'EAS', n_uniq >= 1) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = F, significant = TRUE, sources = 'GO:BP')[[1]]
gost(query = c(map %>% filter(grp == 'EUR', n_uniq >= 1) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = F, significant = TRUE, sources = 'GO:BP')[[1]]

#gprofiler for GO, underrepresentation:
gost(query = c(map %>% filter(grp == 'EAS', n_uniq >= 1) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = F, significant = TRUE, sources = 'GO:BP', measure_underrepresentation = T)[[1]]
gost(query = c(map %>% filter(grp == 'EUR', n_uniq >= 1) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = F, significant = TRUE, sources = 'GO:BP', measure_underrepresentation = T)[[1]]


#gprofiler for KEGG, ranked:
gost(query = c(map %>% filter(grp == 'EAS', n_uniq >= 3) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = TRUE, significant = TRUE, sources = 'KEGG')
gost(query = c(map %>% filter(grp == 'EUR', n_uniq >= 3) %>% group_by(Gene) %>% summarize(median_len = median(bp_pres)) %>% arrange(desc(median_len)) %>% pull(Gene)), organism = 'hsapiens', ordered_query = TRUE, significant = TRUE, sources = 'KEGG')


#prep named gene vectors for FGSEA
ea = map %>% filter(grp == 'EAS', n_uniq >= 3) %>% group_by(Gene) %>% summarize(median(bp_pres)) %>% deframe()
eu = map %>% filter(grp == 'EUR', n_uniq >= 3) %>% group_by(Gene) %>% summarize(median(bp_pres)) %>% deframe()
#retrieve GO terms
gogo <- msigdbr(species = "human", category = "C5", subcategory = 'GO:BP')
gogo <- split(x = gogo$ensembl_gene, f = gogo$gs_name)
#fgsea for GO, ranked:
fgsea(pathways = gogo, stats = ea, scoreType = 'pos') %>% arrange(padj) %>% head(20)
fgsea(pathways = gogo, stats = eu, scoreType = 'pos') %>% arrange(padj) %>% head(20)
#retrieve KEGG terms from MSIGDB:
kegg <- msigdbr(species = "human", category = "C2", subcategory = 'CP:KEGG') 
kegg <- split(x = kegg$ensembl_gene, f = kegg$gs_name)
#fgsea for KEGG, ranked:
fgsea(pathways = kegg, stats = ea, scoreType = 'pos') %>% arrange(padj) %>% head(20)
fgsea(pathways = kegg, stats = eu, scoreType = 'pos') %>% arrange(padj) %>% head(20)

#### Broad GO Terms ####
map %>% filter(pop == 'CHB') %>% pull(Gene) %>% write.table('allgenes_forGOmapper.txt', quote = F, sep = '\t', col.names = F, row.names = F)

broad <- fread('~/Downloads/Broad_GO_terms.txt', header = T) %>%
  as.data.frame() %>%
  separate_rows(ANNOTATED_GENES, sep = ',') %>%
  mutate(ANNOTATED_GENES = trimws(ANNOTATED_GENES)) #nrow = 89,381

## GO Boxplots ##
s4 <- map %>% filter(n_uniq >= 3, grp == 'EAS') %>% left_join(broad, by = c('Gene' = 'ANNOTATED_GENES')) %>% filter(!is.na(TERM)) %>%  ggplot(aes(x = reorder(TERM, bp_pres, FUN = median), y = bp_pres, fill = TERM)) + geom_boxplot() + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10') + xlab('GO Biological Process') + ylab('Mean tract length') + coord_flip() + labs(subtitle = 'EAS') + ggtitle('D')
s5 <- map %>% filter(n_uniq >= 3, grp == 'EUR') %>% left_join(broad, by = c('Gene' = 'ANNOTATED_GENES')) %>% filter(!is.na(TERM)) %>%  ggplot(aes(x = reorder(TERM, bp_pres, FUN = median), y = bp_pres, fill = TERM)) + geom_boxplot() + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10') + xlab('GO Biological Process') + ylab('Mean tract length') + coord_flip() + labs(subtitle = 'EUR') + ggtitle('E')
## GO EAS v EUR scatter plot
p4 <- map %>% filter(n_uniq >= 3) %>% left_join(broad, by = c('Gene' = 'ANNOTATED_GENES')) %>% filter(!is.na(TERM)) %>% group_by(TERM, grp) %>% summarize(me = median(bp_pres)) %>% pivot_wider(names_from = grp, values_from = me) %>% ggplot(aes(x = EAS, y = EUR)) + geom_point(aes(color = TERM)) + stat_cor(hjust = -0.5, vjust = 3) + geom_text_repel(aes(label = TERM, color = TERM)) + theme(legend.position = 'none') + xlab('Mean tract length in EAS') + ylab('Mean tract length in EUR') + ggtitle('D')


#### Pseudogenes ####
map %>% 
  filter(n_uniq >= 3) %>% 
  filter(Biotype %in% c('pseudogene', 'protein_coding')) %>% 
  group_by(Biotype) %>% 
  summarize(median(bp_pres))


##### Monoallelic Expression ####
mae <- read.csv('~/Downloads/mae_genes.csv') %>% distinct(Ensemble.ID, .keep_all = T)
mae$intro <- ifelse(mae$Ensemble.ID %in% map[map$n_uniq >= 1,]$Gene, 'Int', 'Not')
prop.table(table(mae$MAE.1_BAE.0, mae$intro),1)
chisq.test(table(mae$MAE.1_BAE.0, mae$intro))

#graph mean tract length of monoallelically- vs biallelically-expressed genes
s6 <- map %>% 
  filter(n_uniq >= 3) %>% 
  left_join(mae, by = c('Gene' = 'Ensemble.ID')) %>% 
  filter(!is.na(MAE.1_BAE.0)) %>% 
  ggplot(aes(x = as.factor(MAE.1_BAE.0), y = bp_pres)) + geom_boxplot() + stat_compare_means() + scale_y_continuous(trans = 'log10') + ggtitle('F') + xlab('Biallelic (0) vs Monoallelic (1) Expression') + ylab('Mean tract length (bp, log10)')

fig4 <- (p1 | p2 ) / (p3 | p4)
ggsave(fig4, file = 'fig4emil.png', width = 20, height = 15)
supp4 <- (s1 | s2 | s3) / (s4 | s5 | s6)
ggsave(supp4, file = 'supp4emil.png', width = 20, height = 15)
