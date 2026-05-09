---
title: "HiBLUP计算个体的GBLUP育种值"
date: 2026-05-08
categories: 
 - Linux
 - Breeding Value
tags: 
 - Blup
 - GBLUP
 - Linux
---

HiBLUP，中文名“天权”，是华中农业大学赵书红教授团队开发的一款针对农业动物遗传育种的全基因组选择软件。

相比于之前的一系列育种值计算软件，Hiblup最大的优势是计算速度快，此外可选择的计算模型多，不仅能计算BLUP和GBLUP，还可以计算PBLUP、SSBLUP，还能够根据单性状和多性状计算育种值。总之，比`rrBLUP `R包的可拓展性高多了。

官网地址：`https://www.hiblup.com/`

下载方式在官网也有，支持Linux、Windows、MacOS等版本，我用的是Linux版本，数据处理也都是在Linux中处理的。

## 一、GBLUP训练集制作

Hiblup计算群体中个体的GBLUP育种值，需要准备训练集和验证集。

下面先介绍怎么从原始的`vcf`制作符合Hiblup标准格式的训练集数据。


### 1.1 训练集格式转换与snp位点提取
先从`MS353_SNP_delJAK_beagle.vcf`MS353样本集的`2354733`个SNP位点中提取与5k芯片相同位点的SNPs。

```sh
# 把下机的基因型数据转为vcf格式
plink --file All.GT_nonoChr --recode-vcf --out Huzhou2000

# 填充缺失基因型
java -Xmx100g -jar /home/yzhou/biosoft/beagle_v5.5/beagle.27Feb25.75f.jar window=50 overlap=5 gt=Huzhou2000.vcf out=Huzhou2000.beagle

# 提取5k芯片的染色体和位点坐标
grep -v "^#" Huzhou2000.beagle.vcf | cut -f1,2 > 5k_snp_pos.txt

# 从353样本集中提取5k芯片位点
vcftools --vcf /home/yzhou/MS_2025/MS24_DNA_vcf/vcf_combine/MS353_SNP_delJAK_beagle.vcf --positions 5k_snp_pos.txt --recode --out MS353_5k_chip
# After filtering, kept 6094 out of a possible 2354724 Sites
# 基本全提出来了，还是可以的
```
### 1.2 改为Hiblup需要的标准格式
格式转换，转换后需要给文件加表头
```sh
plink --vcf ../rawdata/MS353_5k_chip.recode.vcf --make-bed --out MS353_5k_chip --allow-extra-chr

plink --bfile MS353_5k_chip --recode --out MS353_5k_chip --allow-extra-chr

# 修改ped，只要前4列
awk '{print $1,$3,$4}' ../MS353_5k_chip.ped >MS353_5k_chip.blup.ped

# bim文件snp名重命名
awk 'BEGIN{OFS="\t"} {$2=$1"_"$4; print}' ../MS353_5k_chip.bim > MS353_5k_chip.blup.bim

# 除了plink结果文件，还要表型文件，需要有表头，还可以加上各种协变量
# 表型文件格式，缺失数据可以用“NA”表示
id	sex	season	day	bornweight	location	dam	T1	T2	T3
IND1001	Male	Winter	92	1.2	l32	IND0921	4.76582022911475	-3.6176788560136	24.4309950025418
IND1002	Male	Spring	88	2.7	l36	IND0921	12.4097715906115	10.6740955422296	20.7120836359468
IND1003	Male	Spring	91	1	l17	IND0968	4.85449880306195	0.602879070640462	17.4291590558733
IND1004	Male	Autumn	93	1	l37	IND0968	33.2216999796794	18.1082437974961	22.665404399452
IND1005	Male	Winter	93	2.7	l19	IND0983	13.9741699934432	12.8330653509005	2.73368458760639

# 默认情况下，HIBLUP 会取第二列的性状进行分析，用户可以通过使用 --pheno-pos n 指定不同的列，例如，--pheno-pos 8 表示第八列进行分析

```

### 1.3 训练集SNP效应值估计
对MS353训练集中的每个SNP打分，评估它们对表型（耐热性）的效应值：
```sh
hiblup --single-trait --pheno phenotype_353.txt --pheno-pos 3 --bfile MS353_5k_chip.blup --pedigree MS353_5k_chip.blup.ped --add --snp-effect --thread 32 --out MS353_5k_chip.blup

# 之后会输出一系列的文件：
# 随机效应的方差组分和遗传率（estimated variance components）储存在 *.var 
# 协变量和固定效应的估计系数（Coefficients of all covariates）储存在 *.beta ；截距（Mu）是总体平均值
# 环境随机效应、遗传随机效应（Random effects of all individuals）储存在 *.rand 
# snp效应值（SNP effects）储存在 *.snpeff ，可用于后续的验证集计算
```


## 二、验证集

处理完训练集后，另取一份其他群体样本作为验证集。

### 2.1 湖州2000个体验证集
完成训练集的SNP效应值估计后，开始把huzhou 2k个体的基因型数据作为验证集，计算2k个体的GEBV值。

首先当然要把2k vcf转为 plink 格式：
```sh
# 提取MS353中的SNP pos号码，保证353与2k个体的snp完全一致，否则会报错，无法估计育种值
grep -v "^#" ../rawdata/MS353_5k_chip.recode.vcf | cut -f1,2 > MS353_5k_snp_pos.txt

vcftools --vcf ../rawdata/Huzhou2000.beagle.vcf --positions MS353_5k_snp_pos.txt --recode --out Huzhou2000.beagle.check

# 后续分析中出现报错：2k个体中部分snp缺少基因型，在vcf中显示“.”，在bim文件中基因型错误显示为[0/T]。
# 最简单的办法是删除这些缺失型的位点：
bcftools view -e 'REF="." || ALT="."' Huzhou2000.beagle.check.recode.vcf -Oz -o Huzhou2000.beagle.check.vcf

plink --vcf Huzhou2000.beagle.check.vcf --make-bed --out Huzhou2000.beagle.check --allow-extra-chr

plink --bfile Huzhou2000.beagle.check --recode --out Huzhou2000.beagle.check --allow-extra-chr

# 验证集不用修改ped文件格式
```

之后，用训练集的snp值预测2k个体的GEBV：
```sh
hiblup --pred --bfile Huzhou2000.beagle.check --score ../../train_data_MS353/hiblup_output/MS353_5k_chip.blup.snpeff --threads 10 --out Huzhou2000.beagle.check

# GEBV值储存在 *.bv 文件中
```

### 2.2 结果检验

计算完成后，我们可以查看`*.bv`后缀的文件，第二列为估计的GBLUP：
```bash
id	add_a1
N1509	62.3118
N1649	61.9679
N1136	60.6082
N401	59.2897
N181	49.6425
……
```

再与`rrBLUP`R包估计的GBLUP值作比较，好像`rrBLUP`多了一段截距：
```bash
ID	GEBV_LOE(Ref353)
N1649	162.1159151
N738	153.6514605
N1509	152.8801641
N1705	152.3269435
N181	148.850447
N401	147.9798587
N508	147.1178807
N1972	146.7340417
……
```

计算两份数据的相关系数`r = 0.99942 `。

不用关注两者数值的大小，主要看样本的排名。

可见两种软件计算得到的个体的GBLUP值完全相同。