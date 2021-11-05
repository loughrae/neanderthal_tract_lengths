library(tidyverse)
library(data.table)
library(ggpubr)
library(ggrepel)
library(patchwork)

#read in all raw output from grid search
daim <- fread('all_DAIM_results_sep.txt')
#and reformat
grid = daim %>% 
  as.data.frame() %>% 
  separate(V15, into = c('numeric.AF', 'numeric.sel', 'numeric.t1', 'numeric.t2'), sep = ' ') %>% 
  separate(V7, into = c('logit.AF', 'logit.sel', 'logit.t1', 'logit.t2'), sep = ' ') %>% 
  mutate(numeric.AF = as.numeric(numeric.AF), sel = as.numeric(numeric.sel), t1 = as.numeric(numeric.t1), t2 = as.numeric(numeric.t2)) %>% 
  mutate(numeric.cM = V13*100) %>% 
  mutate(numeric.AF = numeric.AF, digits = 2)


gaf <- fread('gene_archaic_af-all.bed') #maximum Neanderthal allele frequency in each gene-pop
map <- fread('processed_map.txt', header = T) #mean tract lengths and other info for each gene
broad <- fread('~/Downloads/Broad_GO_terms.txt', header = T) %>% #broad Gene Ontology Biological Processes 
  as.data.frame() %>%
  separate_rows(ANNOTATED_GENES, sep = ',') %>%
  mutate(ANNOTATED_GENES = trimws(ANNOTATED_GENES))

#table containing all Neanderthal tracts
haps <- fread('count_gene_overlaps.bed')
colnames(haps) <- c('Chromosome', 'Start', 'End', 'ID', 'bp', 'cM', 'row', 'genes', 'distinct_genes')
haps$pop <- substr(haps$ID, start = 10, stop = 12)
haps$grp <- ifelse(haps$pop %in% c('CHB', 'CHS', 'JPT'), 'EAS', 'EUR')

#filter gene-pops to be used in the grid search
df <- gaf %>%
  mutate(pop = substr(V13, start = 17, stop = 19)) %>% #extract population from filename
  mutate(pop_gene = paste(pop, V4, sep = '_')) %>%
  left_join(map, by = c('pop_gene')) %>% 
  filter(n_uniq >= 3) %>% #has an introgressed tract in at least 3 haploid individuals in its population
  filter(V10 > 0) %>% #has a listed Neanderthal allele frequency
  filter(n_uniq == n_tracts) %>% #tracts aren't broken
  mutate(empirical_af = as.numeric(V7)) #convert AF to numeric and rename
 

 #remove gene-populations that aren't within the range of the grid search's predictions (mean tract length or allele frequency) for any of the three admixture times
  #this is a first pass; we will later focus on admixture time = 1931 generations ago and remove additional gene-pops that don't fit within that range
df_glob <- df %>%
  filter(empirical_af >= min(grid$numeric.AF), empirical_af <= max(grid$numeric.AF)) %>% 
  filter(cM_pres >= min(grid$numeric.cM), cM_pres <= max(grid$numeric.cM)) 

admixtimes <- c(1621, 1931, 2241) #generations ago
lis <- vector(mode = 'list', length = 3*NROW(unique(df_glob$pop_gene))) #three admixture times
names(lis) <- paste(admixtimes, rep(unique(df_glob$pop_gene), 3), sep = '-') #combinations of the three admixture times and pop-genes to name the list
counter = 0
st = Sys.time()

for (admix_time in admixtimes) {
  gr = grid %>% filter(numeric.t2 == admix_time)
  for (pg in unique(df_glob$pop_gene)) {
    #cut the dataframe of pop genes to just one and pull out the allele frequency, mean tract length (in centiMorgans), AI status and gene name
    pdf <- df_glob %>% filter(pop_gene == pg)
    emp_af = pdf %>% pull(empirical_af)
    emp_cM = pdf %>% pull(cM_pres)
    ai_status = pdf %>% pull(AI)
    gene_name = pdf %>% pull(GeneName)
    ind = paste(admix_time, pg, sep = '-') #name of list element
    #matching grid search output to gene-populations by minimising the sum of absolute differences between each predicted mean tract length and AF and the empirical versions
    #normalised by the empirical versions
    lis[[ind]] <- gr %>%
      mutate(diff = (abs(numeric.cM - emp_cM)/emp_cM) + (abs(numeric.AF - emp_af)/emp_af)) %>%
      arrange(diff) %>%
      head(1) %>% #get the closest match by arranging the df in ascending order of the distance and taking the first row (containing the selection parameters)
      mutate(eAF = emp_af, eCM = emp_cM, pop_gene = pg, row = row_number(), AI = ai_status, genename = gene_name)
    
    counter = counter + 1
    if (counter %% 100 == 0) { print(counter) }
  }
}

tim = Sys.time() - st
print(tim)
bind_rows(lis) %>% fwrite('all_inferences.txt', sep = '\t', quote = F, col.names = T, row.names = F)

#### Graphs from Raw Grid Search Output ####
e1 <- grid %>% ggplot(aes(x = as.numeric(numeric.sel), y = as.numeric(numeric.t1), colour = numeric.cM)) + geom_point()  + facet_wrap(~numeric.t2) + labs(color = 'MTL (cM)') + ylab('Selection Time') + xlab('Selection Coefficient') + ggtitle('A')
e2 <- grid %>% ggplot(aes(x = as.numeric(numeric.sel), y = as.numeric(numeric.t1), colour = numeric.AF)) + geom_point()  + facet_wrap(~numeric.t2) + labs(color = 'AF') + ylab('Selection Time') + xlab('Selection Coefficient') + ggtitle('B')

#renaming and rejigging columns in grid to make the stained-graph density plots
modgrid <- grid %>% 
  rename(Logit = V5, Numeric = numeric.cM) %>%
  mutate(Logit = 100*Logit) %>%
  filter(numeric.t2 == 1931) %>%
  pivot_longer(cols = c(Numeric, Logit), names_to = 'Method',  values_to = 'cM_pres') %>%
  dplyr::select(Method, cM_pres)

#make stained-glass density plots by taking the AI and non-AI empirical mean tract lengths from map and grafting on the predicted mean tract lengths from the logistic and numerical versions of DAIM
#storing in long format under the name Method
e3 <- map %>%
  filter(n_uniq >= 3) %>%
  dplyr::select(AI, cM_pres) %>% 
  rename(Method = AI) %>%
  rbind(modgrid) %>%
  mutate(Method = factor(Method, levels = c('AI', 'Not AI', 'Logit', 'Numeric'))) %>%
  ggplot(aes(x = cM_pres, fill = Method)) + geom_density(alpha = 0.3) + xlim(0,0.5) + xlab('Mean tract length (cM)') + ylab('Density') + ggtitle('C')

#make exclusions barplot: take the list of pop-genes with all previous filters (e.g. on broken tracts) and for each admixture time see how many would fall outside the range
#by having their AF or mean tract length either too low or too high
lo = vector(mode = 'list', length = length(admixtimes))
for (admix_time in admixtimes) {
  gr <- grid %>% filter(numeric.t2 == admix_time) 
  lo[[admix_time]] <- df %>% 
    mutate(t2 = admix_time) %>%
    mutate(valid = ifelse(cM_pres > max(gr$numeric.cM) | empirical_af > max(gr$numeric.AF), 'Too High', ifelse(cM_pres < min(gr$numeric.cM) | empirical_af < min(gr$numeric.AF), 'Too Low', 'Valid')))
}
e4 <- lo %>% bind_rows() %>% ggplot(aes(x = AI, fill = valid)) + geom_bar(position = 'fill') + facet_wrap(~t2) + xlab('') + ylab('Proportion of Gene-Populations') + ggtitle('D') + labs(fill = '')


patch5 <- (e1 | e2) / (e3 | e4)
ggsave(file = 'fig5emilia.png', width = 17, height = 15)

#### Graphs on Inferred Selection Parameters #### 
lims <- grid %>% filter(numeric.t2 == 1931) %>% summarize(minlen = min(numeric.cM), maxlen = max(numeric.cM), minaf = min(numeric.AF), maxaf = max(numeric.AF))

#filter all inferences to keep those valid for admixture 1931 generations ago
fer = fread('all_inferences.txt')
intermed <- fer %>%
  filter(numeric.t2 == 1931) %>%
  filter(eCM >= lims[,1] & eCM <= lims[,2] & eAF >= lims[,3] & eAF <= lims[,4]) %>%
  mutate(pop = substr(pop_gene, start = 1, stop = 3)) %>%
  mutate(grp = ifelse(pop %in% c('CHB', 'CHS', 'JPT'), 'EAS', 'EUR')) %>%
  left_join(map %>% select(pop_gene, Gene), by = c('pop_gene')) %>%
  mutate(immune = ifelse(Gene %in% broad[broad$TERM == 'immune system process',]$ANNOTATED_GENES, 'Immune', 'Other'))

#violin plot of selection coefficients for AI genes vs other genes
f1 <- intermed %>% 
  ggplot(aes(x = AI, y = numeric.sel, fill = AI)) + geom_violin() + ylab('Selection') + xlab('') + theme(legend.position = 'none') + stat_compare_means() + ggtitle('A')
#violin plot of selection times for AI genes vs other genes
f2 <- intermed %>% 
  ggplot(aes(x = AI, y = numeric.t1, fill = AI)) + geom_violin() + ylab('Selection Time') + xlab('') + theme(legend.position = 'none') + stat_compare_means(hjust = -1) + ggtitle('B')
#scatterplot showing selection and selection coefficient for each gene (AI and non-AI separately)
f3 <- intermed %>% 
  ggplot(aes(x = numeric.sel, y = numeric.t1)) + geom_jitter(aes(colour = pop)) + xlab('Selection') + ylab('Selection time (generations ago)') + geom_text_repel(aes(label = genename)) + facet_wrap(~AI) + ggtitle('C')

#find candidates for adaptive introgression i.e. genes with selection coefficient >= 0.004 that don't have a known Eurasian AI gene on the same chromosome in the same population
candidates <- map %>% left_join(intermed, by = c('pop_gene')) %>% group_by(pop.x, Chromosome) %>% mutate(sai = sum(AI.x == 'AI')) %>% filter(numeric.sel >= 0.004) %>% dplyr::select(-haploids, -distinct_haploids)  %>% filter(AI.x == 'Not AI') %>% ungroup() %>% filter(sai < 1)
top5 <- candidates %>% filter(Biotype == 'protein_coding') %>% arrange(desc(sel)) %>% head(5) %>% pull(pop_gene)
f4 <- map %>% 
  filter(pop_gene %in% top5) %>% 
  select(pop_gene, AI, bp_pres, cM_pres, Start, End, pop, GeneName, tracts) %>% 
  separate_rows(tracts, sep = ',') %>% 
  mutate(tracts = as.numeric(tracts)) %>% 
  left_join(haps, by = c('pop', 'tracts' = 'row')) %>% 
  mutate(tracts = as.factor(tracts)) %>% 
  ggplot() + 
  geom_rect(aes(xmin = Start.x, xmax = End.x, fill = GeneName), ymin = -Inf, ymax = Inf, alpha = 0.08) + 
  geom_segment(aes(x = Start.y, xend = End.y, y = tracts, yend = tracts)) + ylab('') + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab('Position on chromosome') + 
  facet_wrap(~pop_gene, scales = 'free', nrow = 1) + 
  guides(fill = guide_legend(override.aes = list(alpha = 1))) + labs(fill = 'Gene') + ggtitle('D')

#investigate the selection time for immune genes vs other genes
f5 <- intermed %>% ggplot(aes(x = immune, y = t1, fill = immune)) + geom_violin() + stat_compare_means() + theme(legend.position = 'none') + ylab('Selection time') + xlab('') + ggtitle('F')
f6 <- intermed %>% ggplot(aes(x = immune, fill = t1 > 1431)) + geom_bar(position = 'fill') + facet_wrap(~AI) + xlab('') + ylab('Proportion under early selection') + labs(fill = 'Selected within 500 generations') + theme(legend.position = 'bottom')
patch6 = (((f1 | f2 ) / f3 ) / f4) / ( f5 | f6)
ggsave(file = 'fig6emilia.png', width = 25, height = 30)

#Fig S6
map %>% 
  filter(pop_gene %in% candidates$pop_gene) %>% 
  select(pop_gene, AI, bp_pres, cM_pres, Start, End, pop, GeneName, tracts) %>% 
  separate_rows(tracts, sep = ',') %>% 
  mutate(tracts = as.numeric(tracts)) %>% 
  left_join(haps, by = c('pop', 'tracts' = 'row')) %>% 
  mutate(tracts = as.factor(tracts)) %>% 
  ggplot() + 
  geom_rect(aes(xmin = Start.x, xmax = End.x, fill = GeneName), ymin = -Inf, ymax = Inf, alpha = 0.08) + 
  geom_segment(aes(x = Start.y, xend = End.y, y = tracts, yend = tracts)) + ylab('') + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab('Position on chromosome') + 
  facet_wrap(~pop_gene, scales = 'free') + 
  guides(fill = guide_legend(override.aes = list(alpha = 1))) + labs(fill = 'Gene') + ggtitle('D')
ggsave('all_candidates.png', width = 30, height = 30)

