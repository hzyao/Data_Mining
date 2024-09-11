

############################## PPI 网络 ##################################

# 数据准备，这里我用了我自己随便搞的 gene list，标注上下调基因用到了之前在差异分析中得到的结果，有需要的小伙伴们可以按提示获取（见上文）！
# 大家也可以直接使用自己的 gene list 去 string 进行分析！

# 这是我的 gene list
gene_list <- read.table("./data/gene_list.txt")
gene_list
#           V1
# 1      ADH1B
# 2     ADIPOQ
# 3    ANGPTL7
# 4       APOB
# 5       AQP7
# 6     ATP1A2
# 7  C14orf180
# 8        CA4
# 9    CD300LG
# 10    CHRDL1
# 11     CIDEA
# 12     CIDEC
# 13   COL11A1
# 14   COL10A1
# 15      DLK1
# 16     FABP4
# 17     GLYAT
# 18      GPD1
# 19      IBSP
# 20      MMP1
# 21     MMP11
# 22     MMP13
# 23      PCK1
# 24     PLIN1
# 25     PLIN4
# 26      RBP4


# 这个用来后续标注上下调基因
load("./data/DEG_limma_voom.Rdata")
head(DEG_limma_voom)
#             logFC    AveExpr         t       P.Value     adj.P.Val        B
# FIGF    -5.984847 -0.7193930 -51.67041 1.843289e-309 4.938355e-305 698.5939
# CA4     -6.844833 -2.5701167 -44.96985 3.340380e-261 4.474605e-257 587.5876
# PAMR1   -3.989305  2.3605059 -44.85958 2.161519e-260 1.665003e-256 585.9261
# LYVE1   -4.786578  1.3531474 -44.85132 2.485914e-260 1.665003e-256 585.7724
# CD300LG -6.537456 -0.0898487 -43.57667 6.384798e-251 3.421102e-247 564.1320
# SDPR    -4.600471  2.7186631 -43.38389 1.712581e-249 7.646961e-246 560.8745



##################### STRING 处理结束，开始画图图啦！###########################

# 首先加载我们需要的包包，没有安装的小伙伴们记得安装一下哟！
library(igraph)

# 读取从 STRING 得到的蛋白互作数据
links <- read.delim("./data/string_interactions.tsv")

# 构建无向图，仅选择节点列和边权重列
network <- graph_from_data_frame(d = links[, c(1, 2, 13)], directed = FALSE)

# 计算节点的度数（连接数量）
deg <- degree(network, mode = "all")

# 提取唯一的基因列表
genes <- unique(links$X.node1) %>% as.data.frame()
colnames(genes) <- "gene"

# 筛选出与基因列表匹配的差异表达基因
DEG_filtered <- DEG_limma_voom %>%
  filter(rownames(DEG_limma_voom) %in% genes$gene) %>%
  mutate(gene = rownames(.)) %>%
  select(gene, logFC)

# 结合基因列表与 DEG 数据，基于 logFC 值分组（上下调）
nodes <- genes %>%
  inner_join(DEG_filtered, by = "gene") %>%
  mutate(group = ifelse(logFC > 0, "up", "down"))

# 使用新的节点数据构建网络图
network <- graph_from_data_frame(d = links, vertices = nodes, directed = FALSE)

# 为节点分配颜色，基于分组信息
my_color <- c("#66C2A5", "#FC8D62", "#8DA0CB")[as.numeric(as.factor(V(network)$group))]

# 保存网络图为 PDF 文件
pdf(file = "./figure/ppi_network.pdf", width = 8, height = 8)

# 设置图形参数并绘制网络图
par(bg = "white", mar = c(0, 0, 0, 0))

plot(network, 
     vertex.size = deg,  # 根据度数调整节点大小
     layout = layout.circle,  # 使用圆形布局
     vertex.color = my_color,  # 节点颜色
     vertex.label.cex = 0.9,  # 标签字体大小
     vertex.label.color = "black",  # 标签颜色
     vertex.frame.color = "transparent",  # 节点边框颜色
     edge.width = E(network)$combined_score * 3,  # 边宽度
     edge.curved = 0.3)  # 边的弯曲程度

# 添加图例，显示基因分组信息
legend(x = 1, y = 1, 
       legend = unique(V(network)$group), 
       col = unique(my_color), 
       bty = "n",  # 无边框
       pch = 20,  # 点形状
       pt.cex = 3,  # 点大小
       cex = 0.8,  # 文本大小
       text.col = "black",  # 文本颜色
       horiz = FALSE,  # 垂直排列图例
       y.intersp = 1.5)  # 图例项间距

# 关闭设备，完成 PDF 文件保存
dev.off()

