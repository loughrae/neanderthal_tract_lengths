set.seed(123)
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggVennDiagram)
library(patchwork)

map <- fread('processed_map.txt')
haps <- fread('count_gene_overlaps.bed')
colnames(haps) <- c('Chromosome', 'Start', 'End', 'ID', 'bp', 'cM', 'row', 'genes', 'distinct_genes')
haps$pop <- substr(haps$ID, start = 10, stop = 12)
haps$grp <- ifelse(haps$pop %in% c('CHB', 'CHS', 'JPT'), 'EAS', 'EUR')

s1 <- haps %>% group_by(ID, grp) %>% summarize(tot = sum(bp)) %>% ggplot(aes(x = grp, y = tot, fill = grp)) + geom_boxplot() + stat_compare_means() + ylab('Total tract length per individual') + xlab('') + ggtitle('A') + theme(legend.position = 'none')

p1 <- haps %>% count(ID, grp) %>% ggplot(aes(x = grp, y = n, fill = grp)) + geom_boxplot() + stat_compare_means() + ylab('Number of tracts per individual') + xlab('') + ggtitle('A') + theme(legend.position = 'none')
p2 <- haps %>% ggplot(aes(x = grp, y = bp, fill = grp)) + geom_violin() + ylab('Tract length (log10)') + xlab('') + ggtitle('B') + theme(legend.position = 'none') + scale_y_continuous(trans = 'log10') + stat_compare_means()

p3 <- map %>% 
  filter(n_uniq >= 3) %>% 
  group_by(grp, Gene) %>% 
  summarize(mtl = median(bp_pres)) %>%
  pivot_wider(names_from = grp, values_from = mtl) %>% 
  filter(!is.na(EUR), !is.na(EAS)) %>% 
  ggplot(aes(x = EUR, y = EAS)) + geom_point() + stat_cor() + xlab('Mean tract length in EUR (bp)') + ylab('Mean tract length in EAS (bp)') + ggtitle('C')

s2 <- map %>% 
  filter(n_uniq >= 1) %>% #changing to 3 makes it r = 0.24
  group_by(grp, Gene) %>% 
  summarize(nu = median(n_uniq)) %>%
  pivot_wider(names_from = grp, values_from = nu) %>% 
  filter(!is.na(EUR), !is.na(EAS)) %>% 
  ggplot(aes(x = EUR, y = EAS)) + geom_point() + stat_cor() + xlab('Median number of individuals with tract in gene in EUR populations') + ylab('Median number of individuals with tract in gene in EAS populations') + ggtitle('B')

goh <- fread('map_genes_onto_haps.bed') %>% 
  separate(V4, into = c(NA, 'pop', 'indiv'), remove = F, sep = '_') %>% 
  separate(pop, into = c('pop', NA), sep = '\\.') %>% 
  mutate(grp = ifelse(pop %in% c("CHB", "CHS", "JPT"), 'EAS', 'EUR')) %>%
  rename(Gene = V10, ind_pop = V4) %>%
  mutate(indiv = as.numeric(indiv)) %>%
  mutate(person = floor(indiv/2)) %>% 
  mutate(which = indiv %% 2) %>%
  filter(Gene != '.') %>%
  separate_rows(Gene, sep = ',') 

#### Overall Sharing of Introgressed Genes Between Populations ####
eur <- goh %>% filter(grp == 'EUR') 
eas <- goh %>% filter(grp == 'EAS') 
p4 <- ggVennDiagram(list(EUR = unique(eur$Gene), EAS = unique(eas$Gene))) + theme(legend.position = 'none')  + ggtitle('D') + scale_fill_gradient(low = "white", high = "white") 

#### Sharing Between Equally-Sized Populations ####
holder <- vector(mode = 'list', length = 2000) 
coun = 1
for (siz in c(1,5,10,25,50,100,150,200,500)) {
  for (iter in 1:100) {
    sam_eur <- sample(unique(eur$ind_pop), size = siz)
    sam_eas <- sample(unique(eas$ind_pop), size = siz)
    genes_eur <- eur %>% filter(ind_pop %in% sam_eur) %>% pull(Gene) %>% unique()
    genes_eas <- eas %>% filter(ind_pop %in% sam_eas) %>% pull(Gene) %>% unique()
    holder[[coun]] <- data.frame(samplesize = siz, it = iter, n_eur = NROW(genes_eur), n_eas = NROW(genes_eas), inter = NROW(intersect(genes_eur, genes_eas)))
    print(coun)
    coun = coun + 1
  }
}

sharing <- bind_rows(holder) %>% 
  mutate(perc = (100*inter) / (n_eas + n_eur - inter)) %>% 
  mutate(samplesize = as.factor(samplesize))

p5 <- sharing %>% ggplot(aes(x = samplesize, y = n_eur/n_eas, fill = samplesize)) + geom_boxplot() + xlab('Sample Size') + ylab('EUR/EAS ratio of introgressed genes') + theme(legend.position = 'none') + geom_hline(yintercept = 1) + ggtitle('E')
s3 <- sharing %>% ggplot(aes(x = samplesize, y = perc, fill = samplesize)) + geom_boxplot() + theme(legend.position = 'none') + ylab('% of introgressed genes shared') + xlab('Sample Size') + ggtitle('C')

#### Homozygosity for Introgressed Genes ####
left <- goh %>% filter(which == 0) 
right <- goh %>% filter(which == 1) 
people <- left %>% left_join(right, by = c('grp', 'pop', 'person', 'Gene')) %>%
  distinct(pop, person, Gene, .keep_all = T) %>% 
  group_by(grp, pop, person) %>% 
  summarize(genes_with_match = sum(!is.na(which.y)), n_genes = n(), perc_homo = 100*genes_with_match/n_genes) 

people %>% 
  group_by(grp) %>% 
  summarize(median(perc_homo))

p6 <- people %>% ggplot(aes(x = grp, y = perc_homo, fill = grp)) + geom_violin() + stat_compare_means() + theme(legend.position = 'none') + ylab('Percentage of Genes Homozygous for Introgression') + xlab('') + ggtitle('F')


#### Sharing Between Random Individuals in a Population: EAS v EUR ####
pair_share <- function(df, set) {
holder2 <- vector(mode = 'list', length = 150)
for (i in 1:100) {
  sam1 <- sample(unique(df$ind_pop), size = 1)
  res1 <- df %>% filter(ind_pop %in% sam1) %>% pull(Gene) %>% unique()
  sam2 <- sample(df %>% filter(ind_pop != sam1) %>% pull(ind_pop) %>% unique(), size = 1)
  res2 <- df %>% filter(ind_pop %in% sam2) %>% pull(Gene) %>% unique()
  holder2[[i]] <- data.frame(it = i, n_one = NROW(res1), n_two = NROW(res2), inter = NROW(intersect(res1, res2)), set = set)
}
holder2 %>% bind_rows() %>% mutate(perc_inter = (100*inter) / (n_one + n_two - inter)) %>% return() 
}

eas_pairs <- pair_share(df = goh %>% filter(grp == 'EAS'), set = 'EAS')
eur_pairs <- pair_share(df = goh %>% filter(grp == 'EUR'), set = 'EUR')
s4 <- rbind(eas_pairs, eur_pairs) %>% ggplot(aes(x = set, y = perc_inter, fill = set)) + geom_boxplot() + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('% of introgressed genes shared between two indivdiuals') + ggtitle('D')

#### Sharing Between Random Individuals in a Population: CHB v CEU ####
chb_pairs <- pair_share(df = goh %>% filter(pop == 'CHB'), set = 'CHB')
ceu_pairs <- pair_share(df = goh %>% filter(pop == 'CEU'), set = 'CEU')
s5 <- rbind(chb_pairs, ceu_pairs) %>% ggplot(aes(x = set, y = perc_inter, fill = set)) + geom_boxplot() + stat_compare_means() + theme(legend.position = 'none') + xlab('') + ylab('% of introgressed genes shared between two indivdiuals') + ggtitle('E')

fig1 <- (p1 | p2 | p3) / (p4 | p5 | p6)
ggsave(fig1, file = 'fig2emilia.png', width = 20, height = 15)
supp1 <- (s1 | s2) / (s3 | s4 | s5)
ggsave(supp1, file = 'supp1emilia.png', width = 20, height = 15)



