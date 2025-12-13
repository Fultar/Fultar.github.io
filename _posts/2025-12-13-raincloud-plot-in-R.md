---
title: "R语言云雨图绘制"
date: 2025-12-13
categories: 
 - R
 - scientific plot
tags: 
 - R
 - ggplot2
---
云雨图（Raincloud plot）常用于组间比较，能够同时展示数据的表型、基因型、表达量分布、组间差异等信息，是一种高信息密度和高效率的可视化方法。

在R语言中，我们可以通过半小提琴图+箱型图+散点图的方法绘制云雨图。

# 一、数据准备

这里我用鲈鱼的积温数据表型作为横轴，样本的基因型作为纵轴，探究表型与基因型之间的分布关系。
```R
library(ggplot2)
library(data.table)

rm(list = ls())
gc()

file_pheno <- fread('样本积温.txt')

> file_pheno[1:3,]
        V1    V2    V3     V4    V5
MS24_P2_601 221.0    34      Y     1
MS24_P2_602 225.3    34      Y     1
MS24_P2_603 193.4    34      Y     1
# V5列的高温下存活时间是我们关心的积温表型

# 接着需要读取基因型数据，我采用从vcf中读入并选取需要的SNP位点的基因型
file_vcf <- fread('MS24_P2.gwas.vcf')
snp_13046948 <- subset(file_vcf, POS == c('13046948'))
snp_13047310 <- subset(file_vcf, POS == c('13047310'))
snp_13047491 <- subset(file_vcf, POS == c('13047491'))

# 合并这3个位点的数据, 并转置
gene_vcf <- t(rbind(snp_13046948, snp_13047310, snp_13047491))
colnames(gene_vcf) <- c('snp_13046948', 'snp_13047310', 'snp_13047491')
> gene_vcf[10:20, ]
            snp_13046948 snp_13047310 snp_13047491
MS24_P2_601 "1|0"        "1|0"        "1|0"       
MS24_P2_602 "1|0"        "1|0"        "1|0"       
MS24_P2_603 "0|1"        "0|1"        "0|1"

# 把0|1都改为1|0
gene_vcf[gene_vcf == "0|1"] <- "1|0"

# 合并表型和基因型:表型只要V5列，vcf只要10列以后的基因型
gene_plot <- cbind(file_pheno$V5, gene_vcf[10:nrow(gene_vcf), ])

## 最终的绘图数据，第一列是表型，第二列以后的是基因型：
> gene_plot[1:3,]
                snp_13046948 snp_13047310 snp_13047491
MS24_P2_601 "1" "1|1"        "1|0"        "0|0"       
MS24_P2_602 "1" "1|0"        "0|0"        "1|1"       
MS24_P2_603 "1" "1|0"        "1|0"        "1|0" 
```

# 二、云雨图绘制
第一部分整理好了表型与基因型数据，这里再计算下基因型散点的x轴坐标，使其分布在x轴上，另外可以选择性地加上t检验组间差异
```R
# 绘图要用到包
library(gghalves)
library(ggpubr)
library(dplyr)

# 制作基因型散点的坐标
gene_plot$snp_13046948_xpos <- as.numeric(factor(gene_plot$snp_13046948))

gene_plot$snp_13047310_xpos <- as.numeric(factor(gene_plot$snp_13047310))

gene_plot$snp_13047491_xpos <- as.numeric(factor(gene_plot$snp_13047491))

# t检验显著性
comparisons <- list(
  c("1|1", "1|0"),
  c("1|0", "0|0"),
  c("0|0", "1|1"))

```

下面使用ggplot包开始绘图：
```R
ggplot(gene_plot, aes(x = snp_13047310, y = V1) ) + 
  # 箱型图
  geom_boxplot(aes(fill = snp_13047310, color = snp_13047310),
               outlier.shape = NA,
               width = 0.1, 
               alpha = 0.7) +
  # 半小提琴图
  geom_half_violin(aes(fill = snp_13047310, color = snp_13047310),
                   position = position_nudge(x = 0.1, y = 0),   
                   # 设置小提琴图的一侧（'R'表示右侧）
                   side = 'R', adjust = 1.2, trim = F, 
                    alpha = 0.8) +  
  # 散点图, 注意不同组间减去的值不一样
  geom_point(aes(x = snp_13046948_xpos - 0.1, y = V2), 
             position = position_jitter(width = 0.03),   
             size = 1.5,  color = '#0B2631' ,
             shape = 4) +  
  
  geom_point(aes(x = snp_13047310_xpos - 0.2, y = V2), 
             position = position_jitter(width = 0.03),   
             size = 1.5,   color = '#18445D' ,
             shape = 20) +  
  geom_point(aes(x = snp_13047491_xpos - 0.3, y = V2), 
             position = position_jitter(width = 0.03),   
             size = 1.5,   color = '#1F7799' ,
             shape = 17) +  
  
  stat_compare_means(comparisons = comparisons, method = "t.test") +
  labs(x = "Genotype", y = "Phenotype", title = "GeneName") +
  scale_x_discrete(labels = c("A2/A2", "A1/A2", "A1/A1")) +
  scale_color_manual(values=c("#001F3F",'#178236',"#85144B"))+  # 手动调整点的颜色
  scale_fill_manual(values=c("#001F3F", '#31C950',"#85144B"))+   # 手动调整填充颜色
  coord_flip()+  # 翻转坐标轴
  theme_bw() +
  theme(## 调整画布大小，避免文字显示不全
    plot.margin=unit(rep(0.5,4),'cm'),
    ## 隐藏图例
    legend.position="none",
    ## 调整标题
    plot.title = element_text(hjust = 0.5,face = "bold",size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),# 去除网格线
    ## 调整边框线条
    axis.line.y = element_line(size = 0.5, color = "black"),
    axis.line.x = element_line(size = 0.5, color = "black"),
    ## 调整坐标轴字体
    axis.title.x = element_text(size = 14,face = "bold",color = "black"), # 调整坐标轴字体
    axis.title.y = element_text(size = 14,face = "bold",color = "black"),
    axis.text = element_text(size = 10,color = "black"),
    axis.text.x = element_text(angle = 0,hjust = 0.5 ), # 调整x轴刻度
    ## 移除外边框
    panel.border = element_rect(size = 1.5, color = "black"),
    ## 刻度线调整
    axis.ticks.y = element_line(color = "black",size = 1),
    axis.ticks.x = element_line(color = "black",size = 1),
    axis.ticks.length.x = unit(0.2,'cm'), 
    axis.ticks.length.y = unit(0.2,'cm')
  )

ggsave('filename.png', height = 8, width = 10, dpi =600)

```
最后会输出类似的图片：![alt text](/pictures/raincloud_test.png)


参考网站：
https://blog.csdn.net/letsdrinkwithme/article/details/135442326
