---
title: "注释GWAS筛选到的显著性SNPs"
date: 2026-07-06
categories: 
 - Linux
 - GWAS
tags: 
 - GWAS
 - SNPs
 - Linux
 - bedtools
---

在GWAS分析中，我们通过`GEMMA`、`GCTA`等线性模型软件，可以筛选到许多显著性SNPs。但这些SNPs结果中，只有染色体和位置信息，并没有SNPs位于哪个基因的内含子、外显子或基因间区的结构信息，这需要我们自己来注释。

我在这里分享一个根据`GEMMA`和`GCTA`的结果文件，使用`bedtools`快速注释SNPs结构信息的脚本。

## 一、SNPs输入文件格式整理

我们用`GEMMA`得到的显著性SNPs结果文件格式如下：
```bash
chr	rs	ps	n_miss	allele1	allele0	af	beta	se	logl_H1	l_remle	p_wald	logP
4	.	2050970	0	C	T	0.37	-0.4857806	0.07621047	-37.37568	0.6038327	7.451639e-09	8.12774819297407
8	.	8945253	0	A	T	0.37	-0.3254752	0.05141566	-37.0684	1e+05	9.081277e-09	8.04185307713328
4	.	2050891	0	T	G	0.401	-0.4159068	0.06812961	-38.82142	0.4660491	2.489329e-08	7.60391770144599
8	.	26818072	0	G	C	0.089	-0.5015317	0.08473756	-42.19864	1e-05	5.657212e-08	7.24739754607751
8	.	26815040	0	T	C	0.094	-0.5184381	0.08795359	-42.30195	1e-05	6.290239e-08	7.20133285305882
4	.	24294676	0	A	C	0.182	0.4618176	0.08082621	-39.78295	1.068385	1.381443e-07	6.85966702989055
4	.	24294677	0	A	C	0.182	0.4618176	0.08082621	-39.78295	1.068385	1.381443e-07	6.85966702989055
```

我们需要从中提取三列信息（染色体、SNP位置）转为`snp.bed`格式：
这三列依次是SNPs的染色体（chr）、位置起始（pos-1）和位置结束（pos）
```bash
# snp.bed
4	2050970	2050971
8	8945253	8945254
4	2050891	2050892
8	26818072	26818073
8	26815040	26815041
```

## 二、从GFF3提取外显子/内含子区间

需要准备好物种的`gtf`注释文件。

如果没有`gtf`文件，可以从`gff`格式转为`gtf`：
```bash
gffread -T filename.gff  -o filename.gtf
```

运行以下脚本从`gtf`中拆分基因结构：
```bash

GTF="filename.gtf"   # 替换为你的GTF文件名


# ── 提取基因区间──
grep -v "^#" $GTF | awk '$3=="transcript"' | \
  awk '{
    match($0, /gene_id "([^"]+)"/, arr)
    gene_id = arr[1]
    print $1"\t"$4-1"\t"$5"\t"gene_id
  }' > fish_genes.bed

# ── 提取外显子区间 ──────────────────────────────────
grep -v "^#" $GTF | awk '$3=="exon"' | \
  awk '{
    match($0, /gene_id "([^"]+)"/, g)
    match($0, /transcript_id "([^"]+)"/, t)
    print $1"\t"$4-1"\t"$5"\t"g[1]"|"t[1]
  }' > fish_exons.bed

```


## 三、为SNPs注释基因结构

输入文件为上一步拆分的`gtf`各结构信息文件，以及`snp.bed`文件。

使用的软件是`bedtools`，如果没安装需要安装下。

```bash
# ── 1. SNP落在外显子 ────────────────────────────────
bedtools intersect \
  -a sigSNPS.bed \
  -b fish_exons.bed \
  -wa -wb \
  > fish_sigSNPs_exonic.txt
# 最后一列为：基因名|转录本名

# ── 2. SNP落在基因区间内（外显子+内含子）────────────
bedtools intersect \
  -a sigSNPS.bed \
  -b fish_genes.bed \
  -wa -wb \
  > fish_sigSNPs_in_genes.txt

# ── 3. SNP位于基因间区（不与任何基因重叠）────────────
bedtools intersect \
  -a sigSNPS.bed \
  -b fish_genes.bed \
  -v \
  > fish_sigSNPs_intergenic.txt

# ── 4. SNP位于内含子（基因内 但 不在外显子）──────────
bedtools intersect \
  -a fish_sigSNPs_in_genes.txt \
  -b fish_exons.bed \
  -v \
  > fish_sigSNPs_intronic.txt

echo "BEDTools注释完成"
wc -l fish_sigSNPs_exonic.txt fish_sigSNPs_intronic.txt fish_sigSNPs_intergenic.txt
```

## 四、python脚本合并、整理注释结果：

使用以下脚本把上一步注释得到的文件，整理到一个文件中，方便整合和查看。

```py
#!/usr/bin/env python3
# merge_annotation.py

def parse_exonic(file):
    """同一SNP多个转录本合并到一行，区间取第一次出现的"""
    results = {}
    with open(file) as f:
        for line in f:
            p = line.strip().split('\t')
            key = (p[0], p[1], p[2])          # chr, start, end
            interval = (p[4], p[5])           # 区间 start, end
            gene_trans = p[6]                 # ENSONIG...|ENSONIT...
            parts = gene_trans.split('|')
            gene_id  = parts[0]
            trans_id = parts[1] if len(parts) > 1 else '.'

            if key not in results:
                results[key] = {
                    'interval': interval,
                    'gene_id':  gene_id,
                    'transcripts': [trans_id]
                }
            else:
                results[key]['transcripts'].append(trans_id)
    return results

def parse_intronic(file):
    results = {}
    with open(file) as f:
        for line in f:
            p = line.strip().split('\t')
            key = (p[0], p[1], p[2])
            results[key] = {
                'interval': (p[4], p[5]),
                'gene_id':  p[6]
            }
    return results

def parse_intergenic(file):
    results = {}
    with open(file) as f:
        for line in f:
            p = line.strip().split('\t')
            key = (p[0], p[1], p[2])
            results[key] = {}
    return results

# ── 读取三个文件 ──────────────────────────────────────
exonic     = parse_exonic("fish_sigSNPs_exonic.txt")
intronic   = parse_intronic("fish_sigSNPs_intronic.txt")
intergenic = parse_intergenic("fish_sigSNPs_intergenic.txt")

# ── 写出结果 ──────────────────────────────────────────
with open("snp_annotation_result.txt", "w") as out:
    out.write("chr\tpos\tregion_start\tregion_end\tannotation\tname\n")

    for key, v in sorted(exonic.items(), key=lambda x: int(x[0][2])):
        chr_, start, end = key
        pos = end                              # BED end = 原始pos
        region_start, region_end = v['interval']
        gene_id = v['gene_id']
        trans    = ','.join(sorted(set(v['transcripts'])))
        name = f"{gene_id}|{trans}"
        out.write(f"{chr_}\t{pos}\t{region_start}\t{region_end}\texonic\t{name}\n")

    for key, v in sorted(intronic.items(), key=lambda x: int(x[0][2])):
        chr_, start, end = key
        pos = end
        region_start, region_end = v['interval']
        out.write(f"{chr_}\t{pos}\t{region_start}\t{region_end}\tintronic\t{v['gene_id']}\n")

    for key in sorted(intergenic.keys(), key=lambda x: int(x[2])):
        chr_, start, end = key
        pos = end
        out.write(f"{chr_}\t{pos}\t.\t.\tintergenic\t.\n")

# ── 统计摘要 ──────────────────────────────────────────
print(f"exonic:     {len(exonic)} SNPs")
print(f"intronic:   {len(intronic)} SNPs")
print(f"intergenic: {len(intergenic)} SNPs")
print(f"total:      {len(exonic)+len(intronic)+len(intergenic)} SNPs")
print("输出：snp_annotation_result.txt")

```

最后合并得到的SNPs注释文件内容大致如下：
```bash
chr	pos	region_start	region_end	annotation	name
4	2050859	2050572	2050882	exonic	EVM0007878.1.path1|EVM0007878.1.mrna1
4	2050971	2050926	2051219	exonic	EVM0007878.1.path1|EVM0007878.1.mrna1
4	2051080	2050926	2051219	exonic	EVM0007878.1.path1|EVM0007878.1.mrna1
8	26832215	26832060	26832229	exonic	EVM0026242.3.path1|EVM0026242.3.mrna1
4	2050892	2050572	2051219	intronic	EVM0007878.1.path1
8	8945254	8942619	8961605	intronic	EVM0002321.1.path1
8	27328234	27309330	27333683	intronic	EVM0020052.1.path1
4	24294677	.	.	intergenic	.
4	24294678	.	.	intergenic	.
8	26815041	.	.	intergenic	.
8	26816739	.	.	intergenic	.
8	26818073	.	.	intergenic	.
8	26821686	.	.	intergenic	.
```





