---
title: "使用rMATS进行可变剪切分析"
date: 2026-04-23
categories: 
 - Linux
 - RNA-seq
tags: 
 - RNA-seq
 - Linux
 - Alternative splicing
 - rMATS
 - sashimiplot
---

# 一、rMATS安装

推荐使用`conda`安装：
```bash
conda install bioconda::rmats
```

也可以从官网下载压缩包，解压后自行安装：
```bash
# 网站
https://github.com/Xinglab/rmats-turbo

tar -xvf rmats_turbo_v4_3_0.tar.gz

cd ./rmats_turbo_v4_3_0

./build_rmats

# 然后即可运行
rmats.py {arguments}
```

如果报错：
```shell
/home/yzhou/miniforge3/rMATS/rMATS_C/rMATSexe: error while loading shared libraries: libgsl.so.25: cannot open shared object file: No such file or directory
Traceback (most recent call last):
  File "/home/yzhou/miniforge3/rMATS/rMATS_P/FDR.py", line 53, in <module>
    ifile=open(sys.argv[1]);title=ifile.readline();
          ^^^^^^^^^^^^^^^^^

```

说明可能缺少GSL文件，那么我们来安装下：
```bash
conda install -c conda-forge gsl
```
安装好后就能正常运行了

# 二、rMATS分析可变剪切

## 2.1 所需的分析文件

1. 待分析样本RNA测序后的`bam`格式文件，且构建好索引；
   

2. 待分析物种的`gft`格式的注释文件；
如果没有gft文件，可以用`gffread`软件将`gff`文件转换为`gtf`格式
```bash
gffread -T genomic.gff  -o genomic.gtf
```

3. 样本分组的信息
`rMATS`支持两组间可变剪切比对，所以需要对照组和实验组样本的名单。

例如，在`b1.txt`中写入对照组样本的路径，样本路径间用`,`分割，不需要空格：
```bash
/home/xiang/alternative/01.gonad_transcriptome/10dph/36-10D-F1_S46_R.bam,/home/xiang/alternative/01.gonad_transcriptome/10dph/36-10D-F2_S47_R.bam,/home/xiang/alternative/01.gonad_transcriptome/10dph/36-10D-F3_S48_R.bam
```

也可以使用以下代码快速构建样本名单：
```bash
# grep用于排除.bai索引文件；head表示前40个bam文件作为对照组，tail表示后20个bam文件作为实验组
ls ../BAM_Sorted/rawdata_bam/*.bam | grep -v '\.bai$' | head -n 40 |paste -sd ',' > b1.txt
ls ../BAM_Sorted/rawdata_bam/*.bam | grep -v '\.bai$' | tail -n 20 |paste -sd ',' > b2.txt
```



## 2.2 可变剪切分析

准备好以上文件后，使用`rMATS`进行可变剪切分析：
```bash
rmats.py --b1 /home/xiang/alternative/01.gonad_transcriptome/rmats/b1.txt --b2 /home/xiang/alternative/01.gonad_transcriptome/rmats/b2.txt --gtf /home/xiang/alternative/genome/Oreochromis_niloticus.O_niloticus_UMD_NMBU.100.gtf -t paired --readLength 151 --variable-read-length --nthread 20 --novelSS --od /home/xiang/alternative/01.gonad_transcriptome/rmats/demo --tmp /home/xiang/alternative/tmp

#--b1 为组别1（对照组）的bam文件的路径，若有生物学重复则bam文件路径用逗号隔开；
#--b2 为组别2（实验组）的bam文件的路径，若有生物学重复则bam文件路径用逗号隔开；
# 为单比较组时，仅给b1即可

#--gtf 为已知的基因及转录本的gtf文件；
#--od 即为输出路径；
#-t 测序类型为单端或者双端;
#--readLength 测序读长，可通过测序报告获得，一般为151；
#--libType 文库类型，可选择是否为链特异性；
#--tmp 缓存目录，每次重新运行时需要删除上一轮的内容，否则有概率报错；
```

## 2.3 结果文件各列含义

生成的结果中，一般看`.MATS.JC.txt`后缀的结果文件。

**ID**：rMATS 事件的ID；

**GeneID**：Gene ID；

**geneSymbol**：Gene 名称；

**chr**：染色体；

**strand**：基因的正负链情况；

**riExonStart**:RI事件的起始位置(被滞留的外显子的起始前一碱基位置)；

**riExonEnd**：RI事件的结束位置(被滞留的外显子的结束后一碱基位置)；

**shortES**：A5SS和A3SS中特有，指AS事件后，被裁剪的短外显子的起始位置；

**shortEE**：A5SS和A3SS中特有，指AS事件后，被裁剪的短外显子的结束位置；

**flankingES**：A5SS和A3SS中特有，指AS事件后，与短外显子的剪切端相连的前一外显子的起始位置；

**flankingEE**：A5SS和A3SS中特有，指AS事件后，与短外显子的剪切端相连的前一外显子的结束位置；

**upstreamES/EE**：发生ES事件上游exon的起始/终止位置；

**downstreamES/EE**：发生ES事件下游exon的起始/终止位置；

**1st/2ndExonstart/end**：MXE事件中，互斥的第一个/第二个外显子的起始/终止位置；

**IJC_SAMPLE_1**：sample 1中包含剪切区域的reads数，生物学重复以逗号分隔；

**SJC_SAMPLE_1**：sample 1中不包含剪切区域(skipping junction counts)的reads数，生物学重复以逗号分隔；

**IJC_SAMPLE_2**：sample 2中包含剪切区域的reads数，生物学重复以逗号分隔；

**SJC_SAMPLE_2**：sample 2中不包含剪切区域的reads数(skipping junction counts)，生物学重复以逗号分隔；

**IncFormLen**：包含区域的长度，用于校正；

**SkipFormLen**：跳过区域的长度，用于校正；

**PValue**：两个比较组可变剪切差异的显著性（仅在使用statistical model时存在）；

**FDR**：由 p-value计算的错误发现率（仅在使用statistical model时存在）；

**IncLevel1**：由校正后reads数得到的sample 1的区域等级，生物学重复以逗号分隔；

**IncLevel2**：由校正后reads数得到的sample 2的区域等级，生物学重复以逗号分隔；

**IncLevelDifference**：average(IncLevel1) - average(IncLevel2)。

# 三、可视化可变剪切结果

使用与`rMATS`配套的`rmats2sashimiplot`进行可视化。

推荐使用`conda`安装：
```bash
conda install bioconda::rmats2sashimiplot
```

也可在官网自行下载、解压、安装：
```bash
https://github.com/Xinglab/rmats2sashimiplot
```

`rmats2sashimiplot`可视化所需文件与`rMATS`相同，此外还需要`rMATS`的结果文件，结果文件要用筛选出的关键基因的剪切事件，不然会把所有事件都绘图输出

```bash
rmats2sashimiplot --b1 /home/xiang/alternative/01.gonad_transcriptome/rmats/5dph_control.txt --b2 /home/xiang/alternative/01.gonad_transcriptome/rmats/5dph_treat.txt --event-type RI -e /home/xiang/alternative/01.gonad_transcriptome/rmats/sashimiplot/RI.MATS.JC.TOP20.txt --l1 control --l2 case --exon_s 1 --intron_s 6 -o RI_plot.test

# --b1 对照组样本名单
# --b2 实验组样本名单
# --event-type 可变剪切类型
# -e rMATS输出的可变剪切结果
# --l1/l2 定义b1/b2组名称
```

生成的`pdf`文件：
![alt text](/pictures/sashimiplot_sample)

纵轴`RPKM`为表达量高低；每个色块间的连接线表示检测到多少reads没有这段区域，该区域可能发生了剪切事件