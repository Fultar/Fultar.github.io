---
title: "GWAS qq图和曼哈顿图绘制"
date: 2025-12-10
categories: 
 - R
 - GWAS
tags: 
 - R
 - GWAS
 - CMplot
 - qqPlot
 - manhattanPlot
---

# 一、数据准备

在完成GWAS的LMM线性模型分析后会得到各位点的数据，包括位点的等位基因和频率、效应值（beta）、显著性P值等。
```shell
# 通过GEMMA LMM模型分析得到的结果文件格式：
$head gwas_results.assoc.txt
chr	rs	ps	n_miss	allele1	allele0	af	beta	se	logl_H1	l_remle	l_mle	p_wald	p_lrt	p_score
1	.	3529	0	C	A	0.427	1.100208e+00	1.872614e+00	-3.510330e+02	1.000000e+05	1.000000e+05	5.583062e-01	5.465892e-01	5.484665e-01
1	.	3565	0	T	A	0.495	6.400726e-01	2.872317e+00	-3.511885e+02	1.000000e+05	1.000000e+05	8.241580e-01	8.189848e-01	8.195225e-01
1	.	3571	0	A	T	0.490	1.450959e-01	2.836265e+00	-3.512133e+02	1.000000e+05	1.000000e+05	9.593122e-01	9.580955e-01	9.582110e-01
1	.	38714	0	T	C	0.083	2.197605e+00	2.231193e+00	-3.507057e+02	1.000000e+05	1.000000e+05	3.272621e-01	3.129926e-01	3.169427e-01
1	.	82709	0	A	C	0.422	1.206520e+00	1.412927e+00	-3.508316e+02	1.000000e+05	1.000000e+05	3.953942e-01	3.814045e-01	3.846531e-01
```
# 二、绘图代码

这里我推荐用R包‘CMplot’来可视化位点结果。
```R
#install.packages('CMplot')
library(CMplot)
library(data.table) # 加快读取文件的速度
library(tidyverse)


rm(list = ls())
gc()

LMM_result = fread("gwas_results.assoc.txt")

# 计算P值，方便比较
LMM_result$logP <- -log10(LMM_result$`p_wald`)

# cmplot只需要4列：snpid、染色体、pos、p值
LMM_df = LMM_result %>% select(rs, chr, ps, p_wald)
str(LMM_df)
table(LMM_df$chr)

# 如果想用bonferroni检验可用这个公式计算
pvalue <- c(0.05/nrow(LMM_df))

# 筛选大于阈值的点
sig_snps <- subset(LMM_result, p_wald < 1e-5)           # 具体小于多少p值可以自己决定
sig_snps <- sig_snps[order(sig_snps$p_wald, decreasing = F)]        # 排序
write.table(sig_snps, 'GEMMA结果-显著性snps位点信息.txt',
            quote = F, row.names = F, sep = '\t')

#qq图
CMplot(LMM_df, plot.type = "q",cex=c(0.5,0.5,0.5), col = c('#377EB8'),
       file="jpg", lab.font=2, lab.cex=1.5, main.cex = 1.5,
       main = "Q-Q Plot", dpi = 600,
       chr.labels.angle=90, file.name="GWAS_qqplot")


#manhatten图
# threshold = 1e-5 这段阈值可以自己划定
CMplot(LMM_df, plot.type = "m",threshold = 1e-5,col = c('#377EB8','#f8c120'),
       threshold.col=c('black'),
       threshold.lty = c(1,2),threshold.lwd = c(1,1), amplify = F,
       signal.cex = c(0.8,0.8), signal.pch = c(10,10),cex=c(0.5,0.5,0.5),
       file="jpg", lab.font=2, lab.cex=1.5, main.cex = 1.5,
       main = "Manhatten Plot", dpi = 600,
       chr.labels.angle=0, file.name="GWAS_manhattan")

#SNP-Density Plotting 染色体密度图
CMplot(LMM_df,plot.type = "d",bin.size = 1e6, 
        col = c("#3D9970","#f8c120","#F012BE"), 
        file.name="snp_density")

#带染色体密度的环形曼哈顿图
CMplot(LMM_df,plot.type="c",r=0.4,col=c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
       chr.labels=paste("Chr",c(1:23),sep=""),      # 染色体标签
       threshold=c(1e-6,1e-7),      # 阈值线
       cir.chr.h=1.5,
       amplify=F,
       threshold.lty=c(1,2),
       threshold.col=c("red","blue"),
       signal.line=1,
       signal.col=c("red","green"),
       chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,
       outward=FALSE,
       file="jpg", file.name="cir_manhattan", 
       dpi=300,file.output=TRUE,verbose=TRUE)

```