module load R/4.1.0
module load bedtools

Rscript download_genes.R
#can use existing output to avoiding querying Ensembl every time

echo "===Making gene file into bed file===" 
awk '{print "chr"$1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' genes_hg19.txt | bedtools slop -b 1000 -g human.hg19.genome.1 | sort -k 1,1 -k2,2n > genes_hg19_extended_sorted.bed

echo "### Get AI genes from AI regions ###"
Rscript process_AI_regions.R
sort -k 1,1 -k2,2n AI_regions.bed > AI_regions_sorted.bed
#map human genes onto AI regions with a left outer join. awk removes AI regions without an overlapping gene and then prints the pasted-together population and gene ID.
bedtools intersect -a AI_regions_sorted.bed -b genes_hg19_extended_sorted.bed -loj | awk '$5 != "."'  | awk '{print $4"_"$8}' > AI_genes.txt


echo ==="Make haplotype .bed files and map the haplotypes onto genes==="
for pop in CHB CHS JPT IBS GBR FIN CEU TSI
do
awk '{print $0 "\t" FILENAME}' allhapls_${pop}.txt |  awk '{print "chr"$1 "\t" $3 "\t" $4 "\t" $10"_"$2 "\t" $5 "\t" $6}' | sort -k 1,1 -k2,2n | awk '{print $0 "\t" NR}'  > allhapls_${pop}_sorted.bed
bedtools map -a genes_hg19_extended_sorted.bed -b allhapls_${pop}_sorted.bed -c 5,6,4,4,4,4,7 -o sum,sum,count,count_distinct,collapse,distinct,collapse -null -999 > mapped_tracts_${pop}.bed
done

echo "===Concatenate mapped files and add filename===" 
awk '{print $0 "\t" FILENAME}'  mapped_tracts_*.bed > mapped_tracts-all.bed


echo "### Process mapped genes to get mean tract lengths ###"
Rscript analyse_TLs.R

echo "===Concatenate haplotype files and count number of genes overlapped by each haplotype==="
cat allhapls_*_sorted.bed | sort -k 1,1 -k2,2n > all_hapls-all_sorted.bed
bedtools map -a all_hapls-all_sorted.bed -b genes_hg19_extended_sorted.bed -c 4,4 -o count,count_distinct > count_gene_overlaps.bed

echo "### Map genes onto haplotypes ###"
bedtools map -a all_hapls-all_sorted.bed -b genes_hg19_extended_sorted.bed -c 4 -o count,count_distinct,collapse > map_genes_onto_haps.bed

echo "===Make 10 kb windows across genome and calculate how many genes and Nea haplotypes overlap==="
bedtools makewindows -g human.hg19.genome.1 -w 10000 | grep -v "random" | grep -v "hap" | grep -v "chrUn" | grep -v "chrM" | sort -k 1,1 -k2,2n > genome_windows.bed
bedtools map -a genome_windows.bed -b genes_hg19_extended_sorted.bed -c 4 -o count > windows_genes.bed
bedtools map -a windows_genes.bed -b count_gene_overlaps.bed -c 4,5,5 -o count,mean,median -null -999 > windows_haplotypes-all.bed

echo "### Calculate recombination rates for genes, tracts and genomic windows ###"

for pop in CHB CHS JPT IBS GBR FIN CEU TSI
do
bedtools map -a windows_genes.bed -b allhapls_${pop}_sorted.bed -c 4,5,5 -o count,mean,median -null -999 > windows_haplotypes_${pop}.bed
awk '{print $0 "\t" FILENAME}' ../recombination_maps/${pop}_recombination_map_hg19_chr*.bed | grep -v "#" | sort -k 1,1 -k2,2n > recombination_map_${pop}.bed
bedtools map -a genes_hg19_extended_sorted.bed -b recombination_map_${pop}.bed -c 4,4 -o mean,median -null -999 > genes_rec_${pop}.bed
bedtools map -a allhapls_${pop}_sorted.bed -b recombination_map_${pop}.bed -c 4,4 -o mean,median -null -999 > haps_rec_${pop}.bed
bedtools map -a windows_haplotypes_${pop}.bed -b recombination_map_${pop}.bed -c 4,4 -o mean,median -null "fuck" > windows_rec_${pop}.bed
done

awk '{print $0 "\t" FILENAME}' windows_rec_*.bed > windows_rec-all.bed
awk '{print $0 "\t" FILENAME}' genes_rec_*.bed > genes_rec-all.bed
awk '{print $0 "\t" FILENAME}' haps_rec_*.bed > haps_rec-all.bed

echo "### Getting Neanderthal allele frequencies for each gene ###"
for pop in CHB CHS JPT FIN GBR CEU IBS TSI
do
grep -v "Chrom" ../archaic_genomes/ArchaicUniquelysharedsites/${pop}_0.01_Der_all.tsv | awk '{print "chr"$1 "\t" $2-1 "\t"  $2 "\t" FILENAME "\t" $3 "\t" $4}' | awk '$6 != 0' | sort -k 1,1 -k2,2n > archa$
bedtools map -a genes_hg19_extended_sorted.bed -b archaic_af_${pop}.bed -c 5,5,5,5,5,3 -o max,min,mean,median,count,collapse > gene_archaic_af_${pop}.bed
bedtools map -a allhapls_${pop}_sorted.bed -b archaic_af_${pop}.bed -c 5,5,5,5,5,3 -o max,min,mean,median,count,collapse > haps_archaic_af_${pop}.bed
done

awk '{print $0 "\t" FILENAME}' gene_archaic_af_*.bed > gene_archaic_af-all.bed
awk '{print $0 "\t" FILENAME}' haps_archaic_af_*.bed > haps_archaic_af-all.bed

echo "### Combining grid search results ###"
awk '{print $0 "\t" FILENAME}' grid/logit_results_nov_*.txt | grep -v "t1"  > all_logit_results_nov.txt
awk '{print $0 "\t" FILENAME}' grid/numeric_results_sep_*.txt | grep -v "t1"  > all_numeric_results_nov.txt
paste all_logit_results_nov.txt all_numeric_results_nov.txt > all_DAIM_results_nov.txt

echo '### Make figures ###'
Rscript fig1_oas.R
Rscript fig2_populations.R
Rscript fig3_location.R
Rscript fig4_genefunction.R
Rscript fig5and6_gridsearch.R

