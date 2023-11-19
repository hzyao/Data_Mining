

############################### GSEA 富集分析 ##################################

# 加载GSEA分析所需要的包
library(clusterProfiler) # GSEA富集/基因集读取


# 导入差异分析后的数据，以便后续使用logFC进行基因排序
load("./gsea_data/DEG_limma.Rdata")
head(DEG_limma)
#             logFC  AveExpr         t      P.Value    adj.P.Val        B
# ELOVL6  1.0386031 3.382675  80.14137 1.250078e-08 0.0002957671 9.555464
# PAICS   1.0957881 4.859392  74.97349 1.715935e-08 0.0002957671 9.432400
# DSCC1   1.0066575 3.714065  64.41613 3.529092e-08 0.0003055031 9.109090
# SLBP    0.7234739 5.742791  64.32998 3.551600e-08 0.0003055031 9.105962
# CDKN2B -2.0770246 4.617310 -59.21432 5.264576e-08 0.0003055031 8.902549
# CDKN3   0.8263725 4.852420  58.56735 5.546570e-08 0.0003055031 8.874130

# 加载基因集，基因集介绍往下滑！
geneSet_go <- read.gmt("./gsea_data/c5.go.bp.v2023.2.Hs.symbols.gmt")
head(geneSet_go)
#                                               term     gene
# 1 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS AASDHPPT
# 2 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS  ALDH1L1
# 3 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS  ALDH1L2
# 4 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS   MTHFD1
# 5 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS  MTHFD1L
# 6 GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS  MTHFD2L

# 这里展示的是GO的基因集
# 可以看到共两列，第一列为term，可以理解为通路名称，第二列为基因



# 我们再展示一个基因集，它来自c6，
geneSet_onco <- read.gmt("./gsea_data/c6.all.v2023.2.Hs.symbols.gmt")
head(geneSet_onco)
#           term    gene
# 1 AKT_UP.V1_DN   ACKR3
# 2 AKT_UP.V1_DN  ADGRL1
# 3 AKT_UP.V1_DN  ADHFE1
# 4 AKT_UP.V1_DN    ALPL
# 5 AKT_UP.V1_DN   AMPD2
# 6 AKT_UP.V1_DN ANGPTL4



# 接下来我们进行基因排序
geneList <- DEG_limma$logFC                 # 获取GeneList
names(geneList) <- rownames(DEG_limma)      # 对GeneList命名
geneList <- sort(geneList, decreasing = T)  # 从高到低排序

# 排序后的基因列表
head(geneList)
#     ETV4    MYEOV     ETV5    DUSP6      FGG    FOSL1 
# 2.893571 2.789140 2.617158 2.500177 2.442258 2.355249


# 开始GSEA富集分析
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet_onco, # 基因集
                        pvalueCutoff = 0.05,      # P值阈值
                        minGSSize = 20,           # 最小基因数量
                        maxGSSize = 1000,         # 最大基因数量
                        eps = 0,                  # P值边界
                        pAdjustMethod = "BH")     # 校正P值的计算方法

result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)
# [1] 114  11

# 可以看到，最终富集到114条通路，结果有11列，咱们来看看每列都是什么吧！



# 接下来，我们进行富集结果可视化

# 加载必要包
library(enrichplot) # 富集结果可视化

# 特定通路作图——单个通路
gseaplot2(GSEA_enrichment, "KRAS.600.LUNG.BREAST_UP.V1_DN", color = "red3", pvalue_table = T) 


# 特定通路绘图——多个通路
gseaplot2(GSEA_enrichment, c("VEGF_A_UP.V1_UP", "VEGF_A_UP.V1_DN"), color = c("red3", "blue4"), pvalue_table = T) 


# 展示富集到的通路，我们这里选择展示前15个
dotplot(GSEA_enrichment, showCategory = 15, color = "p.adjust") 


# 将通路分为激活和抑制两个部分进行展示
library(ggplot2)     # 画图图
dotplot(GSEA_enrichment, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10,color = "black"), 
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))



library(tidyverse)
rank_Signal2Noise <- read_tsv("./gsea_data/ranked_Signal2Noise_gene_list_lapa_versus_ctrl_1698913651364.tsv")
head(rank_Signal2Noise)
rank_Signal2Noise <- rank_Signal2Noise[,-4]
saveRDS(rank_Signal2Noise, "./gsea_data/rank_Signal2Noise.rds")


rank_Signal2Noise <- readRDS("./gsea_data/rank_Signal2Noise.rds")
head(rank_Signal2Noise)
# # A tibble: 6 × 3
#   NAME    TITLE                                                                            SCORE
#   <chr>   <chr>                                                                            <dbl>
# 1 SLC15A2 solute carrier family 15 member 2 [Source:HGNC Symbol;Acc:HGNC:10921]             4.46
# 2 MUC19   mucin 19, oligomeric [Source:HGNC Symbol;Acc:HGNC:14362]                          4.32
# 3 ELSPBP1 epididymal sperm binding protein 1 [Source:HGNC Symbol;Acc:HGNC:14417]            4.29
# 4 AGXT    alanine--glyoxylate aminotransferase [Source:HGNC Symbol;Acc:HGNC:341]            4.28
# 5 CYP4F8  cytochrome P450 family 4 subfamily F member 8 [Source:HGNC Symbol;Acc:HGNC:2648]  3.91
# 6 MMP13   matrix metallopeptidase 13 [Source:HGNC Symbol;Acc:HGNC:7159]                     3.80



