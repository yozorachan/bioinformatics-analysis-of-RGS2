setwd('E:/R/draw')

rm(list = ls())
library(tidyverse)
library(survminer)
library(survival)
library(ggplot2)
library(grid)
library(fBasics)
library(tidyr)
library(agricolae)
library(Rmisc)
library(gdata)
library(rlist)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(estimate)
library(reshape2)
library(plyr)
library(ImageGP)
library(ggunchained)
library(RColorBrewer)
library(preprocessCore)
library(ggsci)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#读取数据（下载于https://xenabrowser.net/datapages/）
ALLdata <- data.table::fread("TCGA-LIHC.htseq_counts.tsv",data.table = F)
ALLdata[1:5,1:5]

#读取基因注释文件（下载于https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/）
gtf <- rtracklayer::import('Homo_sapiens.GRCh38.109.chr.gtf.gz')
#转换为数据框
gtf <- as.data.frame(gtf)
#保存
save(gtf,file = "gtf.Rdata")

#去掉基因名小数点及后面的数字方便下一步转换
colnames(ALLdata)[1]<-'gene_id'
ALLdata1<-separate(ALLdata,gene_id,into = c("gene_id"),sep="\\.") 
ALLdata1[1:5,1:5]

#提取编码蛋白，与表达谱进行合并以转换基因名
mRNA<-dplyr::filter(gtf,type=="gene",gene_biotype=="protein_coding")%>%#选择编码蛋白
  select(gene_name,gene_id,gene_biotype)%>%#选择有用的三列
  inner_join(ALLdata1,by ="gene_id")%>%#与表达谱合并
  select(-gene_id,-gene_biotype)%>%
  distinct(gene_name,.keep_all = T)
mRNA[1:5,1:5]

#若row.names有遗漏值
sum(is.na(mRNA$gene_name))
which(is.na(mRNA$gene_name))
mRNA<-mRNA[-79,]
sum(is.na(mRNA$gene_name))

#将gene_name这一列作为行名
row.names(mRNA)<-mRNA[,1]
mRNA<-mRNA[,-1]
mRNA[1:5,1:5]

#保存
save(mRNA,file = 'pancancer_mRNA.Rdata')

#重新开始
rm(list = ls())
load('pancancer_mRNA.Rdata')

#将表达谱倒置，方便后续与临床资料的合并
mRNA<-as.data.frame(t(mRNA))
mRNA[1:5,1:5]

#设计Normal和Cancer的分组，根据TCGA数据ID的特性可区分Normal和Cancer
group_list=ifelse(as.numeric(substr(rownames(mRNA),14,15)) < 10,'tumor','normal')
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=rownames(mRNA)
head(design)

#给design数据数据变成数据框，添加一列ID
class(design)
design<-as.data.frame(design)
design$ID<-row.names(design)
design[1:3,1:3]

#将design分组中tumor这一列的1换成Tumor，0换成Normal
design$tumor[design$tumor==1]  <-  "Tumor"#1换成Tumor
design$tumor[design$tumor==0]  <-  "Normal"#0换成Normal
design<-design[,-1]#去掉normal这一列
colnames(design)[1]<-'Type'#cancer这一列列名改为Type
head(design)

#合并数据
mRNA$ID<-row.names(mRNA)
mRNA=inner_join(mRNA,design,by ="ID",copy=T)
mRNA[1:5,1:5]
#调整列的顺序便于观察
mRNA<-select(mRNA,ID,Type,everything())
mRNA[1:5,1:5]

#读取临床信息
clin <- data.table::fread("TCGA-LIHC.GDC_phenotype.tsv",data.table = F)
clin[1:5,1:5]
#筛选临床资料
clin1 <-clin %>% 
  select(submitter_id.samples,disease_code,gender.demographic,tumor_stage.diagnoses,vital_status.demographic,days_to_last_follow_up.diagnoses) %>% #选取这些列
  rename(ID=submitter_id.samples,Cancer=disease_code,Gender=gender.demographic,Stage=tumor_stage.diagnoses,status=vital_status.demographic,time=days_to_last_follow_up.diagnoses)
clin1[1:5,1:5]
clin1$status[clin1$status=='Alive'] <- TRUE  #Alive换成1
clin1$status[clin1$status=='Dead'] <- FALSE   #Dead换成0
clin1$status<-as.logical(clin1$status)
clin1[1:5,1:5]
save(clin1,file = 'pancancer_clin.Rdata')

load('pancancer_clin.Rdata')

#合并
drawdata<-dplyr::inner_join(mRNA,clin1,by ="ID",copy=T)
drawdata[1:5,1:5]
#调整列名
drawdata<-select(drawdata,ID,Cancer,Gender,Stage,status,time,everything())
drawdata[1:5,1:10]
#保存数据
save(drawdata,file = "pancancer_drawdata.Rdata")

rm(list = ls())
load('pancancer_drawdata.Rdata')
#去除正常样本
OSdata<-filter(drawdata,drawdata$Type == 'Tumor')
#建立分组信息
OSdata <- OSdata[order(OSdata$RGS2), ]
OSdata[,"RGS2"]
drawdata<-OSdata[-c(76:185),]
save(drawdata,file = "alldata.Rdata")

group<-c(rep("RGS_low",times=75),rep("RGS_high",times=75))
groups<-data.frame(drawdata$ID,group)
save(groups,file = "group.Rdata")

rm(list = ls())
#生存曲线分析由（http://gepia.cancer-pku.cn/）获得
load('alldata.Rdata')
load('group.Rdata')

#RGS2与临床分期之间的关系
OSdata=drawdata %>% mutate(Stage=case_when(Stage == "stage iiia" ~ "stage iii",
                                           Stage == "stage iiib" ~ "stage iii",
                                           Stage == "stage iiic" ~ "stage iii",
                                           Stage == "stage iii" ~ "stage iii",
                                           Stage == "stage ii" ~ "stage ii",
                                           Stage == "stage i" ~ "stage i"))
OSdata<-filter(OSdata,OSdata$Stage!='NA')
OSdata$Stage
save(OSdata,file = "stage_draw.Rdata")

load("stage_draw.Rdata")
OSdata$Stage=factor(OSdata$Stage,levels =c("stage i","stage ii","stage iii"))
ggboxplot(OSdata, x = "Stage", y = "RGS2",
          fill = 'Stage') +
  stat_compare_means(label.y = 4.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "stage i")
#ESIMATE肿瘤纯度分析
rm(list = ls())
load('alldata.Rdata')
load('group.Rdata')

mRNA <- select(drawdata,-c(Cancer,Gender,Stage,time,status,Type))
rownames(mRNA)<-mRNA[,1]
mRNA<-mRNA[,-1]
mRNA<-t(mRNA)

#保存为Tab键分割的格式，供estimate包使用
sp_writeTable(mRNA, file="mRNA.tsv")
#输入txt格式的表达矩阵，输出ESIMATE计算结果
filterCommonGenes(input.f= "mRNA.tsv", 
                  output.f="mRNA.gct", id="GeneSymbol")
estimateScore(input.ds = "mRNA.gct",
              output.ds = "mRNA_score.gct", 
              platform="illumina")
ESTI_score <- read.table("mRNA_score.gct",skip = 2,header = T,row.names = 1)
ESTI_score <- as.data.frame(t(ESTI_score[2:ncol(ESTI_score)]))
head(ESTI_score)
#融合ESIMATE计算结果与分组信息
groups$drawdata.ID <- gsub("-",".",groups$drawdata.ID)
ESTI_score$drawdata.ID<-row.names(ESTI_score)
melt<-merge(ESTI_score,groups,by='drawdata.ID')
melt1<-melt(melt)
colnames(melt1)=c("ID","group","status","score")  #设置行名
save(melt1,file = "ESTI_score.Rdata")

load('ESTI_score.Rdata')
#计算误差线
ESTI_Data_summary <- summarySE(melt1, measurevar="score", groupvars=c("group","status"))
head(ESTI_Data_summary)

ESTI_split_violin <- ggplot(melt1,aes(x= status,y= score,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = ESTI_Data_summary,aes(x= status, y= score),pch=19,
             position=position_dodge(0.4),size= 1)+ #绘制均值为点图
  geom_errorbar(data = ESTI_Data_summary,aes(ymin = score-ci, ymax= score+ci), 
                width= 0.05, 
                position= position_dodge(0.4), 
                color="black",
                alpha = 0.8,
                linewidth= 0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  labs(y=("TME score"),x=NULL) + 
  theme_bw() +
  scale_x_discrete(labels=c("Stromal","Immune","ESTIMATE")) +
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox",
                     label.y = max(melt1$score),
                     hide.ns = T)
ESTI_split_violin

#CIBERSORT免疫浸润分析
rm(list = ls())
load('alldata.Rdata')
load('group.Rdata')
source('Cibersort.R')

LM22<- "LM22.txt"

#构建基因表达矩阵
drawdata <- select(drawdata,-c(Cancer,Gender,Stage,time,status,Type))
rownames(drawdata)<-drawdata[,1]
drawdata<-drawdata[,-1]
drawdata<-t(drawdata)
write.table (drawdata, file ="tumor_mRNA.txt", sep ="\t", row.names =TRUE, col.names =TRUE, quote =TRUE)
mRNA<- "tumor_mRNA.txt"

#perm表示置换次数, QN如果是芯片设置为T，如果是测序就设置为F
cibersort <- CIBERSORT(LM22,mRNA, perm = 50, QN = F)
write.table(cibersort, "CIBERSORT.txt")
#提取cibersort前22列数据，23-25列为 P-value， P-value，RMSE数据
cibersort_data <- as.data.frame(cibersort[,1:22])
#将行名转化为列名
cibersort_data<-rownames_to_column(cibersort_data,var="Sample")
#读取样本分组文件,包括两列信息，第一列为样本名，第二列为分组信息。
colnames(groups)[1]<-"Sample"
#连接分组信息和cibersort结果文件
cibersort1<-left_join(cibersort_data,groups,by="Sample")
save(cibersort1,file = "cibersort1.Rdata")

#长宽数据转换
cibersort2<- melt(cibersort1,id.vars=c("Sample","group"))
#设置行名
colnames(cibersort2)<-c("Sample","Group","celltype","composition")
save(cibersort2,file="cibersort2.Rdata")

#绘制差异分析箱线图
boxplot_cibersort<- ggplot(cibersort2, aes(x = celltype, y = composition))+ 
  labs(y="Cell composition",x= "")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5)+ 
  scale_fill_npg()+
  #修改主题
  theme_bw() + 
  theme(axis.title = element_text(size = 12,color ="black"), 
        axis.text = element_text(size= 12,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)
  ) +
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)#隐藏不显著的
boxplot_cibersort

#TIDE评分
rm(list = ls())
load('alldata.Rdata')
load('group.Rdata')

mRNA<-read.table ("tumor_mRNA.txt", sep ="\t")
# 均值标准化：表达量减去每个基因所在样本的均值（即按行计算均值，再用每个表达量-均值）
TIDE <- t(apply(mRNA, 1, function(x){x-(mean(x))}))
TIDE[1:5,1:5]

write.table(TIDE, file = 'mRNA_TIDE.txt', sep = "\t", quote = F, row.names = T)

#TIDE评分计算直接使用官网（http://tide.dfci.harvard.edu/login/）
result <- read.csv('TIDE_result.csv')
colnames(result)
#将TIDE计算结果与分组信息融合
groups$drawdata.ID <- gsub("-",".",groups$drawdata.ID)
colnames(groups)[1]<-"Patient"
result<-merge(result,groups,by='Patient')

# 小提琴图展示结果：  
# 1.TIDE小提琴图：  
result$group <- factor(result$group, levels = c("RGS_high", "RGS_low"))
my_comparisons <- list( c("RGS_high", "RGS_low")) # 添加比较分组  
p1 <- ggviolin(result, x = 'group', y = 'TIDE', fill = 'group',  
               palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')  
p1

# 2.Dysfunction小提琴图：
# dysfunction score的计算原理：免疫失调作用的基因拥有更高的权重，再乘以表达量  
p2 <- ggviolin(result, x = 'group', y = 'Dysfunction', fill = 'group',  
               palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')  
p2

# Exclusion小提琴图：  
# exclusion score是由免疫排斥的基因拥有更高的权重，再乘以表达量得到
p3 <- ggviolin(result, x = 'group', y = 'Exclusion', fill = 'group',  
               palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')  
p3

# MSI小提琴图：  
colnames(result)[6]  
colnames(result)[6] <- c('MSI') # 简化一下列名  
p4 <- ggviolin(result, x = 'group', y = 'MSI', fill = 'group',  
               palette = c("#2E9FDF", "#E7B800"),  
               add = 'boxplot', add.params = list(fill = "white")) +  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')  
p4
p <- p1 + p2 + p3 + p4  
p
save(p,file = "TIDE_paint.Rdata")

#CD8 T cell分析
rm(list = ls())
load('alldata.Rdata')
load('group.Rdata')
#构建标志基因列表
marker<- list('CD8A',c('TOX','PDCD1','LAG3','ENTPD1'),'TCF7')
names(marker)<-c('CD8 T cell','up-regulated gene','down-regulated gene')
save(marker,file = "marker.Rdata")
#构建基因表达矩阵
colnames(groups)[1]<-'ID'
data<-merge(drawdata,groups,by='ID')
mRNA <- select(data,-c(Cancer,Gender,Stage,time,status,Type,group))
rownames(mRNA)<-mRNA[,1]
mRNA<-mRNA[,-1]
mRNA<-t(mRNA)
save(mRNA,file = 'mRNA_CD8T.Rdata')
#建立分组信息
group <-data %>% 
  select(ID,group)      #选取这些列
row.names(group) <- group[,1]
group <- group[,-1]       #使mygroup最终为character
save(group, file = 'group_CD8T.Rdata')

#进行gsva分析
re <- gsva(mRNA, marker, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)
draw_boxplot(re,group,color = c("#e5171a","#1d4a9b"),ylab="Value", xlab = 'CD8 T')

#差异基因筛选&富集分析
#差异基因筛选
rm(list = ls())
load('alldata.Rdata')
load('group.Rdata')

#构建基因表达矩阵
colnames(groups)[1]<-'ID'
data<-merge(drawdata,groups,by='ID')
mRNA <-dplyr::select(data,-c(Cancer,Gender,Stage,time,status,Type,group))
rownames(mRNA)<-mRNA[,1]
mRNA<-mRNA[,-1]
mRNA<-t(mRNA)
mRNA<-2^mRNA-1
mRNA<-round(mRNA, digits = 0)
mRNA <- mRNA[rowMeans(mRNA)>1,]
save(mRNA,file = 'mRNA_kegg.Rdata')

load('mRNA_kegg.Rdata')
#指定分组因子顺序
colData <- data.frame(row.names=colnames(mRNA), condition = factor(data[,'group']))
#构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = mRNA, colData = colData, design= ~condition)
head(dds)
#第二步，计算差异倍数并获得 p 值
#备注：parallel = TRUE 可以多线程运行，在数据量较大时建议开启
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = TRUE)
#注意，RGS_high 在前，RGS_low 在后，意为 RGS_high 相较于 RGS_low 中哪些基因上调/下调
res <- results(dds1, contrast = c('condition', 'RGS_high', 'RGS_low'))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'DESeq2_gene.txt', col.names = NA, sep = '\t', quote = FALSE)

##筛选差异表达基因
#首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & padj<0.05 标识 up，代表显著上调的基因
#log2FC≤-1 & padj<0.05 标识 down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'

res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'DESeq2_select_gene.txt', sep = '\t', col.names = NA, quote = FALSE)

res1_select<- read.table('DESeq2_select_gene.txt',sep = '\t',header = TRUE,row.names = 1)
#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'RGS2', color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-12, 12) + ylim(0, 35)  #定义刻度边界
p

#对upgene和downgene分别富集
res1_up <- subset(res1_select, sig %in% 'up')
res1_down <- subset(res1_select, sig %in% 'down')
write.table(res1_up, file = 'DESeq2_up_gene.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'DESeq2_down_gene.txt', sep = '\t', col.names = NA, quote = FALSE)
#富集分析
upgene_list <- rownames(res1_up)
mRNA_up <- mRNA[upgene_list,]

#symbol to ID
upgene.id <- bitr(
  upgene_list, fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db) 
#enrichKEGG通路富集
upkegg <- enrichKEGG(
  gene = upgene.id$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
head(upkegg)
dotplot(upkegg)

#富集分析
downgene_list <- rownames(res1_down)
mRNA_down <- mRNA[downgene_list,]

#symbol to ID
downgene.id <- bitr(
  downgene_list, fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db) 
#enrichKEGG通路富集
downkegg <- enrichKEGG(
  gene = downgene.id$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
head(downkegg)
dotplot(downkegg)
#Enrichment Map 根据通路之间是否有基因交叠来确定通路间是否存在互作边
edo_up <- pairwise_termsim(upkegg)
emapplot(edo_up, layout="kk") 
#对通路关系网络进行聚类展示
emapplot_cluster(edo_up, node_scale=1.5, layout="kk") 

#Enrichment Map 根据通路之间是否有基因交叠来确定通路间是否存在互作边
edo_down <- pairwise_termsim(downkegg)
emapplot(edo_down, layout="kk") 
#对通路关系网络进行聚类展示
emapplot_cluster(edo_down, node_scale=1.5, layout="kk")

#构建差异基因矩阵
selectgene_list <- rownames(res1_select)
mRNA_select <- mRNA[selectgene_list,]
mRNA_select <- t(mRNA_select)

#表达相关性分析
corr.result<-cor(mRNA_select,method = 'pearson') #计算相关性系数
ggcorrplot(corr= corr.result,
           hc.order = TRUE, 
           tl.cex = 0.1, 
           outline = FALSE,
           )
