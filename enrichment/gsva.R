

############################### GSVA 富集分析 ##################################


# 加载包，没有安装的记得安装一下哟！
library(clusterProfiler)
library(GSVA)

# 加载表达矩阵，这里的数据是log2(fpkm+1)
# 不同数据类型用到参数会有所不同，大家要注意哟！咱们后文会一一介绍！
gsva_data <- readRDS("./gsva_data/exp_dlbc_fpkm.rds")
head(gsva_data)[1:5, 1:5]
#          TCGA-RQ-A6JB-01A TCGA-FF-8046-01A TCGA-FF-A7CW-01A TCGA-RQ-A68N-01A TCGA-FM-8000-01A
# TSPAN6           1.026771       0.76890645       0.35573821        0.7864696        1.0575194
# TNMD             0.000000       0.07154948       0.02051499        0.0000000        0.1362624
# DPM1             5.134474       5.01494933       5.28879599        5.1281748        4.6654726
# SCYL3            1.534159       1.18431712       1.91343387        1.2004038        1.5156969
# C1orf112         1.894523       1.72454346       2.24504933        1.2145814        1.7478664

# 加载基因集
geneSet <- read.gmt("./gsva_data/c6.all.v2023.2.Hs.symbols.gmt")
head(geneSet)
#           term    gene
# 1 AKT_UP.V1_DN   ACKR3
# 2 AKT_UP.V1_DN  ADGRL1
# 3 AKT_UP.V1_DN  ADHFE1
# 4 AKT_UP.V1_DN    ALPL
# 5 AKT_UP.V1_DN   AMPD2
# 6 AKT_UP.V1_DN ANGPTL4

# 可以看到共两列，第一列为term，可以理解为通路名称，第二列为基因

# 我们需要简单处理一下，按照通路名对基因进行分组
geneSet <- split(geneSet$gene, geneSet$term)


# 接下来就可以直接开启GSVA分析啦！
gsva_result <- gsva(as.matrix(gsva_data), geneSet, method = "gsva", min.sz = 1,
                    max.sz = Inf, kcdf = "Gaussian", parallel.sz = 1L)
# Estimating GSVA scores for 189 gene sets.
# Estimating ECDFs with Gaussian kernels
# |==============================================================================================================================| 100%
# 
# Warning messages:
# 1: useNames = NA is deprecated. Instead, specify either useNames = TRUE or useNames = FALSE. 
# 2: In .filterFeatures(expr, method) :
#   3900 genes with constant expression values throuhgout the samples.
# 3: In .filterFeatures(expr, method) :
#   Since argument method!="ssgsea", genes with constant expression values are discarded.


# 我们看一下结果长啥样！
head(gsva_result)[1:5, 1:5]
#                      TCGA-RQ-A6JB-01A TCGA-FF-8046-01A TCGA-FF-A7CW-01A TCGA-RQ-A68N-01A TCGA-FM-8000-01A
# AKT_UP.V1_DN               0.16805899       0.03322358      -0.27828409       0.03211608      -0.01535361
# AKT_UP.V1_UP               0.06711940      -0.23781825      -0.14112397       0.17413258      -0.28218023
# AKT_UP_MTOR_DN.V1_DN       0.05688383       0.11518378      -0.08668666       0.16164634      -0.05672208
# AKT_UP_MTOR_DN.V1_UP       0.03973035      -0.26006328      -0.14883800       0.16617699      -0.18649161
# ALK_DN.V1_DN               0.02174707       0.12221061      -0.13024908       0.03811105       0.06346741

# 可以看到，行为通路，列为样本，其中的值就是样本在对应通路中的GSVA得分，也就是通路活性估计值。




