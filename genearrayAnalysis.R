
###
# 基因表达矩阵-找差异基因
# edger
library(edgeR)
library(openxlsx)

# 1. 载入数据 读取read count数
data <-  read.delim("circleRNA_expression.txt", row.names=1, stringsAsFactors=FALSE) 
# openxlsx::read.xlsx()

head(data)
dim(data)

#2. 构建分组变量
#分为 Control组和Case组  分别为4个和3个重复
targets <- data.frame(Lane = c(1:6,8), Treatment = c(rep("Control",4),rep("Case",3)),
                      Label = c(paste("Con", 1:4, sep=""), paste("Case", 1:3, sep="")));


#3. 创建基因表达列表 进行标准化因子计算
y <- DGEList(counts=data[,1:7], group=targets$Treatment, genes=data.frame(Length=data[,8]));
colnames(y) <- targets$Label;
dim(y);

#过滤表达量偏低的基因
#基因在至少3个样本中得count per million（cpm）要大于1
keep <- rowSums(cpm(y)>1) >= 3;
y <- y[keep,];
dim(y)
# [1] 16494 7
#重新计算库大小
y$samples$lib.size <- colSums(y$counts);

#3. 进行标准化因子计算 默认使用TMM方法
y <- calcNormFactors(y);
y


#这里主要是通过图形的方式来展示实验组与对照组样本是否能明显的分开
#以及同一组内样本是否能聚的比较近 即重复样本是否具有一致性
plotMDS(y);

#4. 估计离散度
y <- estimateCommonDisp(y, verbose=TRUE)
# Disp = 0.02002 , BCV = 0.1415 
y <- estimateTagwiseDisp(y);

plotBCV(y);

#5. 通过检验来鉴别差异表达基因
et <- exactTest(y);
top <- topTags(et);
top

#6. 定义差异表达基因与基本统计
summary(de <- decideTestsDGE(et)); # 默认选取FDR = 0.05为阈值

#输出
# [,1] 
# -1  2094   #显著下调
# 0  12060   #没有显著差异
# 1   2340   #显著上调

#图形展示检验结果
detags <- rownames(y)[as.logical(de)];
plotSmear(et, de.tags=detags);
abline(h=c(-1, 1), col="blue");




###
# 差异基因的GSEA分析




###
# GO/KEGG富集分析





###
# 特殊基因集的其他分析