library(biomaRt)
library(tidyverse) 

ens <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'ensembl_gene_id', 'external_gene_name', 'gene_biotype'), mart = ens)

genes %>% filter(chromosome_name %in% c(1:22, 'X')) %>% write.table('genes_hg19.txt', col.names = F, row.names = F, sep = '\t', quote = F)  
