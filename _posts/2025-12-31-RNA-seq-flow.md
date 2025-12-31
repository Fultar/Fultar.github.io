---
title: "RNA-seq上游分析流程"
date: 2025-12-31
categories: 
 - R
 - RNA-seq
tags: 
 - R
 - RNA-seq
---

RNA的上游分析流程比SNP calling稍微简单点，我这里用自己的鲈鱼样本数据做个参考流程。

# 一、QC质控

一般给的测序数据都是`fasta`格式，那么我们利用fastQC软件对获得的fastq序列文件进行质量分析，生成html格式的结果报告，记得把每份样品的报告都浏览一遍：
```sh
nohup fastqc *.raw.fastq.gz> fastqc.log 2>fastqc.error.log &
```

# 二、序列比对

- 先对鲈鱼的参考基因组文件构建索引：

```sh
# gff转格式为gtf
gffread -T GCA_022435785.1_ASM2243578v1_genomic.cds.gff  -o GCA_022435785.1_ASM2243578v1_genomic.cds.gtf

# 建立索引
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir GCA_022435785.1_star_index --genomeFastaFiles /home/yzhou/MS_2025/genome/GCA_022435785.1_ASM2243578v1_genomic.fna --sjdbGTFfile GCA_022435785.1_ASM2243578v1_genomic.cds.gtf --sjdbOverhang 149


--runThreadN：线程数。

--runMode genomeGenerate：构建基因组索引。

--genomeDir：索引目录。（index_dir一定要是存在的文件夹，需提前建好）

--genomeFastaFiles：基因组文件。

--sjdbGTFfile：基因组注释文件。可选项，STAR可以从中提取剪切位点，可以提高比对准确性，推荐！

--sjdbOverhang：reads长度减1。
```



- 索引建立好后开始正式的比对，我这里使用的是`star`比对软件，`star`可以直接输出排序文件`sort.bam`
  
```sh
单样本测试：
STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 6 --genomeDir /home/yzhou/MS_2025/genome/GCA_022435785.1_star_index --alignIntronMin 20 --alignIntronMax 50000 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA --outFilterMismatchNmax 2 --outSJfilterReads Unique --outSAMmultNmax 1 --outFileNamePrefix /home/yzhou/MS_2025/MS_RNA_analysis/MS25_RNA_ --outSAMmapqUnique 60 --readFilesCommand gunzip -c --readFilesIn /home/yzhou/MS_2025/rawdata/MS24_RNA_601-700/MS24_P2_RNA601.R1.raw.fastq.gz /home/yzhou/MS_2025/rawdata/MS24_RNA_601-700/MS24_P2_RNA601.R2.raw.fastq.gz
```

批处理的bash脚本：
```sh
#!/bin/bash

GENOME_DIR="/home/yzhou/MS_2025/genome/GCA_022435785.1_star_index"
OUTDIR="/home/yzhou/MS_2025/MS_RNA_analysis"

# 因为我的样本数字后缀是600-700，分为4个批次跑，减轻服务器压力

# 循环 601 到 625 
for i in $(seq 601 625)
do
    sample="MS24_P2_RNA${i}"
    R1=${sample}.R1.raw.fastq.gz
    R2=${sample}.R2.raw.fastq.gz
    
    if [[ -f $R1 && -f $R2 ]]; then
        echo "Processing sample: $sample"

        STAR \
            --twopassMode Basic \
            --quantMode TranscriptomeSAM GeneCounts \
            --runThreadN 6 \
            --genomeDir $GENOME_DIR \
            --alignIntronMin 20 \
            --alignIntronMax 50000 \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbOverhang 149 \
            --outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA \
            --outFilterMismatchNmax 2 \
            --outSJfilterReads Unique \
            --outSAMmultNmax 1 \
            --outFileNamePrefix $OUTDIR/${sample}_ \
            --outSAMmapqUnique 60 \
            --readFilesCommand gunzip -c \
            --readFilesIn $R1 $R2
    else
        echo "Warning: Missing files for $sample"
    fi
done



nohup bash star_aline_MS601-625.sh >aline_MS601-625.log&
nohup bash star_aline_MS626-650.sh >aline_MS626-650.log&
nohup bash star_aline_MS651-675.sh >aline_MS651-675.log&
nohup bash star_aline_MS676-700.sh >aline_MS676-700.log&


# 修改star比对输出的结果文件名：
for file in `ls *_Aligned.sortedByCoord.out.bam`; do x=${file/Aligned.sortedByCoord.out.bam/}sorted.bam; mv $file $x; done

```

这里讲下`star`软件的参数：
```sh
--twopassMode Basic：使用two-pass模式进行reads比对。简单来说就是先按索引进行第一次比对，而后把第一次比对发现的新剪切位点信息加入到索引中进行第二次比对。

--quantMode TranscriptomeSAM GeneCounts：生成Aligned.toTranscriptome.out.bam后缀文件；reads 映射到 转录本坐标系 的 BAM 文件，主要用在 RSEM 等转录本定量，不是用来做基因组可视化的。

--runThreadN：线程数。

--genomeDir：索引目录。

--alignIntronMin：最短的内含子长度。（根据GTF文件计算）

--alignIntronMax：最长的内含子长度。（根据GTF文件计算）

--outSAMtype BAM SortedByCoordinate：输出BAM文件并进行排序。

--sjdbOverhang：reads长度减1。

--outSAMattrRGline：ID代表样本ID，SM代表样本名称，PL为测序平台。在使用GATK进行SNP Calling时同一SM的样本可以合并在一起。

--outFilterMismatchNmax：比对时允许的最大错配数。

--outSJfilterReads Unique：对于跨越剪切位点的reads(junction reads)，只考虑跨越唯一剪切位点的reads。

--outSAMmultNmax：每条reads输出比对结果的数量。

--outFileNamePrefix：输出文件前缀。

--outSAMmapqUnique 60：将uniquely mapping reads的MAPQ值调整为60，满足下游使用GATK进行分析的需要。

--readFilesCommand：对FASTQ文件进行操作。

--readFilesIn：输入FASTQ文件的路径。

注意：
1、star支持fasta和fastq格式的文件
2、如果read文件是压缩形式的，如果是*.gz, --readFilesCommand zcat或者 --readFilesCommand gunzip -c，如果是bzip2格式的， --readFilesCommand bunzip2 -c
3、对于readfilesin参数，支持多个测序文件一起比对
如果是单端测序的多测序文件一起比对，则每个文件需要用","隔开，如--readFilesIn sample1.fq,sample2.fq,sample3.fq     注意中间没有空格！！！
如果是双端测序的多测序文件比对，格式为S1read1,S2read1,S3read1 S1read2,S2read2,S3read2
即read1和read2内部用逗号隔开，read1和read2之间用空格隔开！
```

# 三、reads定量

需要用到上一步比对完成后的`sort.bam`后缀文件，以及定量分析软件`featureCounts`

```shell
# 检查bam是否已排序
samtools view -H /home/yzhou/MS_2025/MS_RNA_analysis/BAM_Sorted/rawdata_bam/MS24_P2_RNA601_sorted.bam | grep '^@HD'

# 定量计算
nohup featureCounts -T 8 -p -t exon -g gene_id -a /home/yzhou/MS_2025/genome/GCA_022435785.1_ASM2243578v1_genomic.cds.gtf -o MS24_P2_RNA_counts.txt /home/yzhou/MS_2025/MS_RNA_analysis/BAM_Sorted/rawdata_bam/*.bam >gene_counts.log 2>&1 &

```

# 四、表达量计算

表达量计算使用R语言，通过定量得到的counts数计算每个基因的TPM

```R
library(data.table)

rm(list=ls())
gc()

# 读取文件
file_counts <- fread("MS24_P2_RNA_counts.txt", header=TRUE)

# 提取基因长度（单位：bp）
gene_length <- file_counts$Length

# 提取 count 矩阵（第7列开始是样本）
counts <- as.matrix(file_counts[, 7:ncol(file_counts)])
rownames(counts) <- file_counts$Geneid

# 将长度转换为 kilobases
length_kb <- gene_length / 1000

# 计算 RPK = counts / length_kb
rpk <- sweep(counts, 1, length_kb, "/")

# 计算每个样本的 scaling factor（所有基因 RPK 之和 / 1e6）
scaling_factor <- colSums(rpk) / 1e6

# 计算 TPM = RPK / scaling_factor
tpm <- sweep(rpk, 2, scaling_factor, "/")

# 保存结果
write.table(tpm, file = "MS24_P2_RNA_TPM_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)

# TPM 标准化
tpm_log2 <- log2(tpm + 1)

write.table(tpm_log2, file = "MS24_P2_RNA_TPM_log2_matrix.txt", sep = "\t", quote = FALSE, col.names = NA)

```

输出的表达量TPM矩阵可用于RNA-seq下游分析，比如寻找差异基因、火山图、KEGG和GO分析等。


