---
title: "R语言绘制PCA主成分分析图"
date: 2025-12-05
categories: 
 - R
 - GWAS
tags: 
 - R
 - PCA
 - ggplot
 - GWAS
 - principal components analysis
 - linux
 - plink
 - bash
 - shell
---

# 一、文件准备

PCA首先需要群体的数据，可以是表型或基因型矩阵。

在GWAS分析中，vcf格式的基因型文件，需要经过plink过滤后，最终转为plink格式，才能进行PCA群体分层。


# 二、PCA分析

可以选用Plink和GCTA做PCA。二者不同的是，GCTA会计算样本的亲缘关系矩阵，并将其引入PCA的计算中，导致两款软件最后的计算结果略有不同。

plink计算pca：
```bash
plink  \
    --bfile samplename \    # plink格式的样本
    --pca n \               # n 为需要计算的pca数量，一般计算前3个主成分就够了
    --out pca_output        # 输出文件名
```

GCTA计算pca，需要先计算亲缘关系矩阵，再计算pca：
```bash
gcta \
--bfile samplename \        # plink格式样本
--make-grm \                # 输出kinship矩阵
--make-grm-alg 0 \          # 可以选0或1，0,为Yang的方法，1为Van的方法
--out kinship_output
 
gcta \
--grm kinship_output \      # 选择上一步生成的kinship矩阵
--pca n \                   # n为主成分个数
--out pca_output
```

以上两款软件最后的pca分析都会生成两个结果文件：`.eigenval`和`.eigenvec`。

`.eigenval`记录的是特征向量，用于绘图；`.eigenvec`记录的是特征方差，即每个主成分PC占解释率的百分比。

# 三、PCA绘图

使用R语言的ggplot包绘制散点图：

```R
library(ggplot2)
library(ggsci)

rm(list = ls())
gc()

# PCA 图
## 1. 读取数据文件
file_pca_val <- read.table("pca.eigenval",header = F)
file_pca_vec <- read.table("pca.eigenvec",header = F)

# > head(file_pca_vec)
#          V1          V2        V3         V4         V5 
# MS24_P2_601 MS24_P2_601 0.1118100  0.2149710 0.13412900    
# MS24_P2_602 MS24_P2_602 0.0831065  0.1877040 0.15384700      
# MS24_P2_603 MS24_P2_603 0.1211640  0.2145500 0.14093700      
# MS24_P2_604 MS24_P2_604 0.0992557  0.1947370 0.16821500      
# MS24_P2_605 MS24_P2_605 0.0713111 -0.0752131 0.00596257     
# MS24_P2_606 MS24_P2_606 0.0622080 -0.0869438 0.02502900     

### 表型类群、家系文件
file_group <- read.table("样本分类.txt",header = T)

# > head(file_group)
#         sample family
# MS24_P2_RNA601      Y 
# MS24_P2_RNA602      Y 
# MS24_P2_RNA603      Y 
# MS24_P2_RNA604      Y 
# MS24_P2_RNA605     AC 
# MS24_P2_RNA606     AC 

#### 为pca添加家系表型
file_pca_vec$family <- file_group$family

## 2. 计算每个主成分的方差解释比例
pve <- file_pca_val$V1 / sum(file_pca_val$V1) * 100

## 3. ggplot绘图

ggplot(file_pca_vec, 
       aes(x = V3, y = V4, color = family)) +
  geom_point(size = 3) +
  # 添加坐标原点
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  # 添加散点的置信区间椭圆
  stat_ellipse(aes(x = V3, y = V4), linetype = 2, size = 0.5, level = 0.95) + 
  # 修改坐标轴标签，添加主成分百分比
  labs(
    x = paste0("PC1 (", round(pve[1], 2), "%)"),
    y = paste0("PC2 (", round(pve[2], 2), "%)")
  ) +
  theme_bw()+
  # 设置配色
  scale_color_d3('category20')

ggsave('主成分分析图.png', width = 12, height = 8, dpi = 300)

```

可以看到通过前两个主成分已经很好地把各家系样本区分开来，在后续的GWAS分析中，我们可以把PCA作为协变量加入模型。
![alt text](../pictures/pca_plot_1205.png)



部分信息可能参考以下网站：
https://www.cnblogs.com/wzbzk/p/18086804
https://cloud.tencent.com/developer/article/1940981