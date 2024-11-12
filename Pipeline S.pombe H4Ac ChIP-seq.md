# Pipeline: S.pombe H4Ac ChIP-seq 

## 1. Sample information

| Source_name | Strain  | Genotype                                                   | Assay Type    | Rep  |
| ----------- | ------- | ---------------------------------------------------------- | ------------- | ---- |
| 5% Input    | YYX43sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216               | ChIP-seq/H4Ac | Rep1 |
| 5% Input    | YYX43sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216               | ChIP-seq/H4Ac | Rep2 |
| WT          | YYX43sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216               | ChIP-seq/H4Ac | Rep1 |
| WT          | YYX43sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216               | ChIP-seq/H4Ac | Rep2 |
| cph1Δ       | YYX46sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph1Δ::KANMX | ChIP-seq/H4Ac | Rep1 |
| cph1Δ       | YYX46sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph1Δ::KANMX | ChIP-seq/H4Ac | Rep2 |
| cph2Δ       | YYX47sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph2Δ::KANMX | ChIP-seq/H4Ac | Rep1 |
| cph2Δ       | YYX47sp | dy1779 h- his3 d1 ura4-d18 leu1-32 ade6-m216  cph2Δ::KANMX | ChIP-seq/H4Ac | Rep2 |

## 2. Pipeline

```bash
##== linux command ==##
ml trim_galore ## trim_galore/0.6.7 
ml samtools ## samtools/1.16.1 
ml deeptools ## deeptools/3.5.1
ml bowtie2 ## bowtie2/2.2.5
dir=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/
fasta=/data/yixuanp/0_data/chipseq/pombe_h4ac_230911/fastq/
FILES=${fasta}*_1.clean.fq.gz
index=/data/qiz/2_genome/pombe/bowtie2_index/
trim_fq_dir=${dir}trim_galore/
bowtie2_dir=${dir}bowtie2/
peak_dir=${dir}peak/
bw_dir=${dir}bigwig/
mkdir -p $trim_fq_dir $bowtie2_dir $peak_dir $bw_dir
echo "trim_galore=====================================>"
tg="trim_galore"
for f in $FILES
do
echo $f 
    base=$(basename $f "_1.clean.fq.gz")
    echo ${fasta}${base}_2.clean.fq.gz
    f1=$f
    f2=${fasta}${base}_2.clean.fq.gz
$tg --paired --quality 20 --length 25 -e 0.1 --stringency 4 --dont_gzip -o $trim_fq_dir $f1 $f2
echo "trimgalore done!"
R1=$trim_fq_dir${base}_1.clean_val_1.fq
R2=$trim_fq_dir${base}_2.clean_val_2.fq
##perform align
echo "align================================================================>"
bowtie2 -p 20 -x $index -1 $R1 -2 $R2 | samtools view -S -b -o ${bowtie2_dir}${base}.raw.bam
echo "sort bam================================================================>"
samtools sort -o ${bowtie2_dir}${base}.sorted.bam ${bowtie2_dir}${base}.raw.bam
echo "mark and remove duplicate bam================================================================>"
java -jar /data/qiz/pipline/picard.jar MarkDuplicates INPUT=${bowtie2_dir}${base}.sorted.bam OUTPUT=${bowtie2_dir}${base}.sorted.marked.bam METRICS_FILE=${bowtie2_dir}${base}.sorted.marked.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT
samtools markdup -r ${bowtie2_dir}${base}.sorted.last.bam ${bowtie2_dir}${base}.sorted.markdup.bam 
samtools index ${bowtie2_dir}${base}.sorted.marked.bam
echo "bamcoverage================================================================>"
bamCoverage -p 20 --bam ${bowtie2_dir}${base}.sorted.marked.bam --skipNonCoveredRegions --normalizeUsing RPKM -o $bw_dir${base}.bw
done
echo "complete!"
```

## 3. Visualization 

### 3.1 Genome browser tracks (IGV visualization)

Open bigwig files (product from last step) in IGV and take representative views of browser tracks. 

- [ ] > *Where: /Users/yixuan/Desktop/FangCloud/LiLabData/Yixuan_data/11_Manuscript/Rpd3S/PDU/Figure 7B-C Pombe_ChIP-seq*/bigwig

### 3.2 Heatmap over transcription units

#### `BamCompare` from deepTools was applied to calculate the log2 fold change of mutant versus wt. 

```bash
##== linux command ==##
module load deeptools ## deeptools/3.5.1
inp_wt_rep1=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/wt_inp_rep1.sorted.marked.bam
inp_wt_rep2=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/wt_inp_rep2.sorted.marked.bam
delcph1_rep1=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/delcph1_rep1.sorted.marked.bam
delcph1_rep2=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/delcph1_rep2.sorted.marked.bam
delcph2_rep1=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/delcph2_rep1.sorted.marked.bam
delcph2_rep2=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/delcph2_rep2.sorted.marked.bam
wt_rep1=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/wt_rep1.sorted.marked.bam
wt_rep2=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bowtie2/wt_rep2.sorted.marked.bam
bamCompare -b1 $delcph1_rep1 -b2 $wt_rep1 -o h4ac_delcph1_rep1_vs_wt.bw
bamCompare -b1 $delcph2_rep1 -b2 $wt_rep1 -o h4ac_delcph2_rep1_vs_wt.bw
bamCompare -b1 $delcph1_rep2 -b2 $wt_rep2 -o h4ac_delcph1_rep2_vs_wt.bw
bamCompare -b1 $delcph2_rep2 -b2 $wt_rep2 -o h4ac_delcph2_rep2_vs_wt.bw
echo "done"
```



#### Then `computeMatrix`, `plotProfile` and `plotHeatmap` functions from deepTools to generate the metagene plot and heatmap. 

```bash
##== linux command ==##
module load deeptools ## deeptools/3.5.1
bw=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bamcompare/*.bw
bed=/data/qiz/2_genome/pombe/pombe.bed
computeMatrix scale-regions -S $bw -R $bed --beforeRegionStartLength 1000 --regionBodyLength 3000 --afterRegionStartLength 1000 --skipZeros -p 10 -o h4ac_genobed_log.gz
gz=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/heatmap_log2/h4ac_genobed_log.gz
tss=/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/heatmap_log2/h4ac_TSS_log.gz
plotProfile -m $gz --perGroup -out profile_h4ac_genobed_log.pdf
echo "done"
#plotHeatmap -m $gz --colorMap Greens --heatmapHeight 12 --heatmapWidth 3 -out h4ac_genobed_log.pdf
#plotProfile -m $tss --perGroup -out profile_h4ac_tss_log.pdf
#plotHeatmap -m $tss --colorMap Greens --heatmapHeight 12 --heatmapWidth 3 -out h4ac_tss_log.pdf
```

### 3.3 Assess replicate reproducibility

To study the reproducibility between replicates and across conditions, the genome is split into 1000 bp bins, and a Pearson correlation of the log2-transformed values of read counts in each bin is calculated between replicate datasets. 

```bash
##== linux command ==##
##/data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bam/
module load deeptools ## deeptools/3.5.1
multiBamSummary bins --bamfiles *.bam --binSize 1000 --numberOfProcessors 10 --outRawCounts repbin1000.txt -o repbin1000.npz
plotCorrelation -in /data/yixuanp/0_data/chipseq/pombe_h4ac_index2/bam/repbin1000.npz --corMethod pearson --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --removeOutliers -o heatmap_pearsonCorr_readCounts_bin1000.pdf --outFileCorMatrix pearsonCorr_readCounts_bin1000.tab
```


