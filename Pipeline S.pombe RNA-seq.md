# Pipeline: S.pombe RNA-seq 

## 1. Sample information

| Source_name | Strain  | Genotype                                                   | Assay Type        | Rep  |
| ----------- | ------- | ---------------------------------------------------------- | ----------------- | ---- |
| WT1         | YYX43sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216               | RNA-seq/pA enrich | Rep2 |
| WT1         | YYX43sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216               | RNA-seq/pA enrich | Rep3 |
| delCph1     | YYX46sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph1Δ::KANMX | RNA-seq/pA enrich | Rep2 |
| delCph1     | YYX46sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph1Δ::KANMX | RNA-seq/pA enrich | Rep3 |
| delCph2     | YYX47sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph2Δ::KANMX | RNA-seq/pA enrich | Rep2 |
| delCph2     | YYX47sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph2Δ::KANMX | RNA-seq/pA enrich | Rep3 |
| WT2         | YYX41sp | wt (parental strain of set2 deletion)                      | RNA-seq/pA enrich | Rep2 |
| WT2         | YYX41sp | wt (parental strain of set2 deletion)                      | RNA-seq/pA enrich | Rep3 |
| delSet2     | YYX42sp | set2:null:KANMX                                            | RNA-seq/pA enrich | Rep2 |
| delSet2     | YYX42sp | set2:null:KANMX                                            | RNA-seq/pA enrich | Rep3 |

## 2. Pipeline

```bash
##== linux command ==##
ml STAR  #STAR/2.7.10a
ml samtools #samtools/1.16.1 
ml deeptools #deeptools/3.5.1
#the following analysis requires the indexes of genome fasta and gtf files, so I have rebuiled the reference genome to perform the meta analysis in R
index=/data/qiz/2_genome/pombe/STAR_index/
dir=/data/yixuanp/0_data/rnaseq/pombe_230911/
fastq_dir=/data/yixuanp/0_data/rnaseq/pombe_230911/fastq/
STAR_dir=${dir}STAR/
bam_dir=${dir}bamfiles/
bigwig_dir=${dir}bigwig/
FILES=${fastq_dir}*_1.clean.fq.gz
mkdir -p $STAR_dir $bigwig_dir $bam_dir
for fq in $FILES
do
    echo $fq
    base=$(basename $fq "_1.clean.fq.gz")
    echo ${fastq_dir}${base}_2.clean.fq.gz
    R1=${fastq_dir}${base}_1.clean.fq.gz
    R2=${fastq_dir}${base}_2.clean.fq.gz
    echo $base $R1 $R2
STAR --runThreadN 6 --runMode alignReads --genomeDir $index --readFilesIn $R1, $R2 --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outFilterMismatchNmax 2 --outFilterMultimapScoreRange 0 --alignIntronMax 500000 --outSAMunmapped None --outSAMtype BAM Unsorted --outFileNamePrefix ${STAR_dir}$base"_"
    samtools sort -o ${STAR_dir}${base}.sorted.bam ${STAR_dir}${base}_Aligned.out.bam
##MarkDuplicates define duplicates and give a marker
     java -jar /data/qiz/pipline/picard.jar MarkDuplicates INPUT=${STAR_dir}${base}.sorted.bam OUTPUT=${STAR_dir}${base}.sorted.marked.bam METRICS_FILE=${STAR_dir}${base}.sorted.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT
##index the bam file
     samtools index ${STAR_dir}${base}.sorted.marked.bam
bamCoverage -p 20 --bam ${STAR_dir}${base}.sorted.marked.bam --skipNonCoveredRegions --normalizeUsing RPKM -o ${bigwig_dir}${base}.bw
done
echo "done!"
```

## 3. Calculate counts per genes: featureCounts

```bash
##== linux command ==##
module load subread #subread/2.0.1
files=/data/yixuanp/0_data/rnaseq/pombe_230911/STAR/*.sorted.marked.bam
featureCounts -a /data/yixuanp/1_reference/pombe/Schizosaccharomyces_pombe.ASM294v2.57.gtf -t exon -g gene_id -T 10 -o counts.txt $files
```

## 4. Visualization 

### PCA analysis using DE_seq

```R
##== R command ==##
counts <- read.table("counts.txt",row.names = 1,header = T)
count<-counts
counts<-counts[,-c(1,2,3,4,5,6,9,12,15,18)]
colnames(counts) <- gsub("X.data.qiz.4_otherdata.PYX.20231007_DE_Pombe.bamfiles.", "", colnames(counts))
colnames(counts) <- gsub(".sorted.marked.bam", "", colnames(counts))
```

```R
##== R command ==##
index1<-data.frame(
  name= c("dcph1_rep2","dcph1_rep3","dcph2_rep2","dcph2_rep3","wt1_rep2","wt1_rep3"),
  group=c("dcph1","dcph1", "dcph2", "dcph2","wt1","wt1")
)
index1$group<-factor(index1$group)
```

```R
##== R command ==##
library(DESeq2)
library(DEGreport)
library(ggplot2)
library(ggplotify)
```

```R
##== R command ==##
dds <- DESeqDataSetFromMatrix(countData=counts2, 
                              colData=index2, 
                              design=~group)
keep <- rowSums(counts(dds)) >= 10
dds<- dds[keep,]
vsdata <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsdata, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca<-ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=2,shape=16) +
  geom_text_repel(aes(label=name),size=3, max.overlaps = 30)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_cartesian(xlim=c(-30,30),ylim=c(-10,10))
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
pdf("PCA_all_rep2_rep3.pdf", width = 4.7, height =4)
pca
dev.off()
```

### Differential expression analysis

```R
##== R command ==##
cph1_rep2<-data.frame(dcph1_rep2=counts$dcph1_rep2,row.names = row.names(counts))
cph1_rep3<-data.frame(dchp1_rep3=counts$dcph1_rep3,row.names = row.names(counts))
cph2_rep2<-data.frame(dchp2_rep2=counts$dcph2_rep2,row.names = row.names(counts))
cph2_rep3<-data.frame(dchp2_rep3=counts$dcph2_rep3,row.names = row.names(counts))
cphwt_rep2<-data.frame(chpwt_rep2=counts$wt1_rep2,row.names = row.names(counts))
cphwt_rep3<-data.frame(chpwt_rep3=counts$wt1_rep3,row.names = row.names(counts))
set2_rep2<-data.frame(set2_rep2=counts$dset2_rep2,row.names = row.names(counts))
set2_rep3<-data.frame(set2_rep3=counts$dset2_rep3,row.names = row.names(counts))
setwt_rep2<-data.frame(set2wt_rep2=counts$wt2_rep2,row.names = row.names(counts))
setwt_rep3<-data.frame(set2wt_rep3=counts$wt2_rep3,row.names = row.names(counts))
write.csv(cph1_rep2,"cph1_rep2.csv")
write.csv(cph1_rep3,"cph1_rep3.csv")
write.csv(cph2_rep2,"cph2_rep2.csv")
write.csv(cph2_rep3,"cph2_rep3.csv")
write.csv(cphwt_rep2,"cphwt_rep2.csv")
write.csv(cphwt_rep3,"cphwt_rep3.csv")
write.csv(set2d_rep2,"set2_rep2.csv")
write.csv(set2d_rep3,"set2_rep3.csv")
write.csv(setwt_rep2,"setwt_rep2.csv")
write.csv(setwt_rep3,"setwt_rep3.csv")
```

```R
##== R command ==##
genes <- c("cph1", "cph2", "set2")

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through each time point
for (i in 1:length(genes)) {
  time <- genes[i]
  
  # Construct file names for _1 and _2 files
  if (time == "set2") {
    file_1 <- paste0("setwt_rep2.csv")
    file_2 <- paste0("setwt_rep3.csv")
  } else {
    file_1 <- paste0("cphwt_rep2.csv")
    file_2 <- paste0("cphwt_rep3.csv")
  }
  file_other_1 <- paste0(time, "_rep2.csv")
  file_other_2 <- paste0(time, "_rep3.csv")
  
  # Read data for _1 and _2 files
  data_1 <- read.csv(file_1,header = T)
  data_2 <- read.csv(file_2, header = T)
  data_other_1 <- read.csv(file_other_1)
  data_other_2 <- read.csv(file_other_2)
  
  # Use the first column of file_1 as row names
  rownames(data_1) <- data_1[, 1]
  
  # Create a data frame for the current time point
  current_df <- data.frame(
    Gene = rownames(data_1),
    wt_1 = data_1[, 2],
    wt_2 = data_2[, 2],
    treat_1 = data_other_1[, 2],
    treat_2= data_other_2[, 2]
  )
  col_name_other_1 <- gsub(".csv", "", file_other_1)
  col_name_other_2 <- gsub(".csv", "", file_other_2)
  colnames(current_df)[colnames(current_df) == "treat_1"] <- col_name_other_1
  colnames(current_df)[colnames(current_df) == "treat_2"] <- col_name_other_2
  
  # Store the data frame in the list
  data_frames[[time]] <- current_df
}

# Write each data frame to a separate CSV file
for (i in 1:length(data_frames)) {
  time <- genes[i]
  output_file <- paste0("Raw_counts_", time, ".csv")
  write.csv(data_frames[[time]], output_file, row.names = FALSE)
}
```

```R
##== R command ==##
library(edgeR)
# List of time points
genes <- c("cph1", "cph2", "set2")

# Loop through each time point
for (gene in genes) {
  # Load raw counts data
  file_name <- paste0("Raw_counts_", gene, ".csv")
  counts <- read.csv(file_name, header = TRUE, row.names = 1)
  # Filtering based on read counts
  countData <- counts[rowSums(counts) > 20, ]
  
  # Correlation heatmap
  cor_matrix <- cor(countData)
  pdf(paste0("correlation_heatmap_", gene, ".pdf"), height = 8, width = 8)
  heatmap(cor_matrix)
  dev.off()
  
  # Create DGEList and normalize data
  y <- DGEList(countData, group = c(1, 1, 2, 2))
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  
  # Exact test for differential expression
  et <- exactTest(y)
  # Filtering based on p-value and logFC
  sig.up <- et$table$PValue < 0.05 & et$table$logFC > 1
  sig.dn <- et$table$PValue < 0.05 & et$table$logFC < -1
  
  # Adjust p-values
  et$table$padj <- p.adjust(et$table$PValue, method = "BH")
  test <- as.data.frame(et)
  write.csv(test, file = paste0("fold_", gene, ".csv"), row.names = T)
  # Get lists of enriched and depleted genes
  Del.enriched <- et$table$logFC > 1 & et$table$padj < 0.01
  Del.depleted <- et$table$logFC < -1 & et$table$padj < 0.01
  
  # Store the results in a data frame
  DE.exp.data <- countData[Del.enriched | Del.depleted, ]
  DE.exp.data <- as.matrix(DE.exp.data)
  
  # Create DE_gene data frame
  Del.genes <- rownames(test[which(test$logFC > 1 & test$padj < 0.01), ])
  WT.genes <- rownames(test[which(test$logFC < -1 & test$padj < 0.01), ])
  
  up <- data.frame(ENSEMBL = Del.genes, group = paste("up_", gene, sep = ""))
  dn <- data.frame(ENSEMBL = WT.genes, group = paste("dn_", gene, sep = ""))
  DE_gene <- rbind(up, dn)
  
  # Save the DE_gene data frame
  write.csv(DE_gene, file = paste0("DE_gene_", gene, ".csv"), row.names = FALSE)
```

```R
##== R command ==##
library(VennDiagram)

# List of time points and corresponding files
time_points <- c("cph1", "cph2", "set2")

# Create lists to store up and dn genes for each time point
up_gene_lists <- list()
dn_gene_lists <- list()
venn_colors <- c("#4DBBD5", "#A6A6A6", "#F39B7F")
# Loop through each time point
for (time in time_points) {
  # Read DE gene data
  de_gene <- read.csv(paste0("DE_gene_", time, ".csv"), header = TRUE)
  
  # Define DE gene categories
  up_genes <- de_gene[which(de_gene$group == paste0("up_", time)), 1]
  dn_genes <- de_gene[which(de_gene$group == paste0("dn_", time)), 1]
  
  # Store up and dn genes in lists
  up_gene_lists[[time]] <- up_genes
  dn_gene_lists[[time]] <- dn_genes
}

# Generate Venn diagrams for up genes
venn_file_up <- "venn_diagram_up.png"
venn.diagram(
  x = up_gene_lists,
  filename = venn_file_up,
  output = TRUE,
  imagetype = "png",
  main = "Venn Diagram of Upregulated Gene Lists",
  col = "transparent",
  fill = venn_colors,
  alpha = 0.75,
  col.label = venn_colors,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = venn_colors,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.cex = 1.2,
  cex = 1.5,
  margin = 0.05
)

# Generate Venn diagrams for dn genes
venn_file_dn <- "venn_diagram_dn.png"
venn.diagram(
  x = dn_gene_lists,
  filename = venn_file_dn,
  output = TRUE,
  imagetype = "png",
  main = "Venn Diagram of Downregulated Gene Lists",
  col = "transparent",
  fill = venn_colors,
  alpha = 0.75,
  col.label = venn_colors,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = venn_colors,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.cex = 1.2,
  cex = 1.5,
  margin = 0.05
)
```

```R
##== R command ==##
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
# List of time points and corresponding files
time_points <- c("cph1", "cph2", "set2")
file_names <- c("fold_cph1.csv", "fold_cph2.csv", "fold_set2.csv")
time_colors <- list(
  "cph1" = c( "#4DBBD5", "#A6A6A6", "#F39B7F"),
  "cph2" = c("#4DBBD5", "#A6A6A6", "#F39B7F"),
  "set2" = c("#4DBBD5", "#A6A6A6", "#F39B7F")
)
for (i in seq_along(time_points)) {
  time <- time_points[i]
  file_name <- file_names[i]
  
  # Read fold data
  fold_data <- read.csv(file_name, header = TRUE)
  colnames(fold_data)[1]<-"ENSEMBL"
#  fold_data$padj <- p.adjust(fold_data$PValue, method = "BH")
  
  fold_dn <- fold_data[which(fold_data$logFC < -1 & fold_data$padj < 0.01),]
  fold_up <- fold_data[which(fold_data$logFC > 1 & fold_data$padj < 0.01),]
  fold_unchange <- fold_data[which(fold_data$logFC >= -1 & fold_data$logFC <= 1),]
up_ten <- fold_up[order(-fold_up$logFC), ][1:5, ]
dn_ten <- fold_dn[order(fold_dn$logFC), ][1:5, ]
valcano_p<-ggplot() + 
  geom_vline(xintercept = c(-1,1),color="black",size=0.25,lty=5)+
  geom_point(data=fold_unchange,aes(y= -log10(PValue),x=logFC),size=1,shape=16,col=time_colors[[time]][2],alpha=0.5)+
  geom_point(data=fold_up,aes(y= -log10(PValue),x=logFC),size=1,shape=16,col=time_colors[[time]][3],alpha=0.5)+
  geom_point(data=fold_dn,aes(y= -log10(PValue),x=logFC),size=1,shape=16,col=time_colors[[time]][1],alpha=0.5)+
  geom_text_repel(data=up_ten,aes(y= -log10(PValue),x=logFC,label=ENSEMBL),size=2)+
  geom_text_repel(data=dn_ten,aes(y= -log10(PValue),x=logFC,label=ENSEMBL),size=2)+
  coord_cartesian(xlim=c(-10,10))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(),axis.text.x = element_text(size = 13),axis.text.y = element_text(size = 13))+
  labs(title =paste(time, " vs WT") , x = paste("logFC",sep = ""), y = "-log10(PValue)")
  ggsave(paste("DE_gene_Valcano_label_", time, "_1008.pdf", sep = ""), plot = valcano_p, height = 3, width = 3.5)
#  pdf(paste("Length_express_", time, "_merge_0816.pdf", sep = ""), height = 3.5, width = 3)
 
}
```

## 5. Packages 

```R
library(DESeq2)
library(DEGreport)
library(ggplot2)
library(ggplotify)
library(edgeR)
library(VennDiagram)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
sessionInfo()
## R version 4.2.2 (2022-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 22621)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
## [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
## [3] LC_MONETARY=Chinese (Simplified)_China.utf8
## [4] LC_NUMERIC=C                               
## [5] LC_TIME=Chinese (Simplified)_China.utf8    
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] RColorBrewer_1.1-3          ggrepel_0.9.3              
##  [3] ggpubr_0.6.0                VennDiagram_1.7.3          
##  [5] futile.logger_1.4.3         edgeR_3.38.4               
##  [7] limma_3.52.4                ggplotify_0.1.1            
##  [9] ggplot2_3.4.2               DEGreport_1.32.0           
## [11] DESeq2_1.36.0               SummarizedExperiment_1.26.1
## [13] Biobase_2.56.0              MatrixGenerics_1.8.1       
## [15] matrixStats_1.0.0           GenomicRanges_1.48.0       
## [17] GenomeInfoDb_1.32.4         IRanges_2.30.1             
## [19] S4Vectors_0.34.0            BiocGenerics_0.42.0        
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_2.1-0            ggsignif_0.6.4             
##   [3] rjson_0.2.21                circlize_0.4.15            
##   [5] XVector_0.36.0              GlobalOptions_0.1.2        
##   [7] ggdendro_0.1.23             clue_0.3-64                
##   [9] rstudioapi_0.15.0           bit64_4.0.5                
##  [11] AnnotationDbi_1.58.0        fansi_1.0.4                
##  [13] codetools_0.2-19            splines_4.2.2              
##  [15] logging_0.10-108            mnormt_2.1.1               
##  [17] doParallel_1.0.17           cachem_1.0.8               
##  [19] geneplotter_1.74.0          knitr_1.43                 
##  [21] Nozzle.R1_1.1-1.1           broom_1.0.5                
##  [23] annotate_1.74.0             cluster_2.1.4              
##  [25] png_0.1-8                   compiler_4.2.2             
##  [27] httr_1.4.6                  backports_1.4.1            
##  [29] Matrix_1.6-0                fastmap_1.1.1              
##  [31] cli_3.6.1                   formatR_1.14               
##  [33] lasso2_1.2-22               htmltools_0.5.5            
##  [35] tools_4.2.2                 gtable_0.3.3               
##  [37] glue_1.6.2                  GenomeInfoDbData_1.2.8     
##  [39] dplyr_1.1.2                 Rcpp_1.0.11                
##  [41] carData_3.0-5               vctrs_0.6.3                
##  [43] Biostrings_2.64.1           nlme_3.1-162               
##  [45] iterators_1.0.14            psych_2.3.6                
##  [47] xfun_0.39                   stringr_1.5.0              
##  [49] lifecycle_1.0.3             rstatix_0.7.2              
##  [51] XML_3.99-0.14               zlibbioc_1.42.0            
##  [53] MASS_7.3-60                 scales_1.2.1               
##  [55] parallel_4.2.2              lambda.r_1.2.4             
##  [57] ComplexHeatmap_2.12.1       yaml_2.3.7                 
##  [59] memoise_2.0.1               yulab.utils_0.0.6          
##  [61] reshape_0.8.9               stringi_1.7.12             
##  [63] RSQLite_2.3.1               genefilter_1.78.0          
##  [65] foreach_1.5.2               BiocParallel_1.30.4        
##  [67] shape_1.4.6                 rlang_1.1.1                
##  [69] pkgconfig_2.0.3             bitops_1.0-7               
##  [71] evaluate_0.21               lattice_0.21-8             
##  [73] purrr_1.0.1                 cowplot_1.1.1              
##  [75] bit_4.0.5                   tidyselect_1.2.0           
##  [77] plyr_1.8.8                  magrittr_2.0.3             
##  [79] R6_2.5.1                    generics_0.1.3             
##  [81] DelayedArray_0.22.0         DBI_1.1.3                  
##  [83] pillar_1.9.0                withr_2.5.0                
##  [85] abind_1.4-5                 survival_3.5-5             
##  [87] KEGGREST_1.36.3             RCurl_1.98-1.12            
##  [89] tibble_3.2.1                car_3.1-2                  
##  [91] crayon_1.5.2                futile.options_1.0.1       
##  [93] utf8_1.2.3                  rmarkdown_2.23             
##  [95] GetoptLong_1.0.5            locfit_1.5-9.8             
##  [97] blob_1.2.4                  ConsensusClusterPlus_1.60.0
##  [99] digest_0.6.33               xtable_1.8-4               
## [101] tidyr_1.3.0                 gridGraphics_0.5-1         
## [103] munsell_0.5.0
```



