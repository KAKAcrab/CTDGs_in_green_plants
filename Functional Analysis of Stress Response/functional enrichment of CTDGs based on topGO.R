library(AnnotationHub)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(graph)
library(GO.db)
library(SparseM)
library(ggplot2)
library(DOSE)
library(topGO)
library(grid)
library(Rgraphviz)
library(gridExtra)
library(openxlsx)
setwd("~/Desktop/project/CTDG_of_207_spe/draft_fig/")
## 使用topGO对差异表达基因进行功能富集，适用于cluster profile无法得到显著结果的差异表达基因
## 需要准备差异表达基因列表和包括所有基因的 gene2GO 文件

#########################################################################################################################
#########################################################################################################################

## 加载成簇差异表达基因列表，由DESeq2生成；加载该物种的所有基因，作为gene2go的mapping文件
diff <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/1_salt_DE-CTDGs_genes_list",header = F)
diff <- as.character(diff[,1]) #第一列转换为字符类型
head(diff)
all_genes <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/all_genes_list",header = F)
all_genes <- as.character(all_genes[,1])
head(all_genes)
##gene2go数据导入
gene2go <- readMappings("~/R  project/207spe_ctdg/1.expression/16.Aty/gene2go",sep = '\t')
genenames <- names(gene2go)
genelist <- factor(as.integer(genenames %in% diff))
# view(genelist)
names(genelist) <- genenames
head(genelist)
# view(genelist)

#########################################################################################################################
###################################使用不同ontology进行富集分析,MF#######################################################

GOdata_MF <- new("topGOdata", ontology="MF", allGenes = genelist,
                 annot = annFUN.gene2GO, gene2GO = gene2go)
##统计检验
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata_MF, test.stat)
test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata_MF, test.stat)
test.stat <- new("weightCount", testStatistic = GOFisherTest, name="Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata_MF, test.stat)
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata_MF, test.stat)

elim.ks <- runTest(GOdata_MF, algorithm = "elim", statistic = "ks")
allRes_MF <- GenTable(GOdata_MF, classic=elim.ks, KS=resultKS, weight = resultWeight,
                      orderBy = "weight", ranksOf = "classic", topNodes =15) ##只保留前15个
allRes_MF <- as.data.frame(allRes_MF)
##添加物种分组信息，保存文件
allRes_MF$Group <- "salt"
allRes_MF$ontology <- "MF"

#########################################################################################################################
###################################使用不同ontology进行富集分析,BP#######################################################

GOdata_BP <- new("topGOdata", ontology="BP", allGenes = genelist,
                 annot = annFUN.gene2GO, gene2GO = gene2go)
##统计检验
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata_BP, test.stat)
test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata_BP, test.stat)
test.stat <- new("weightCount", testStatistic = GOFisherTest, name="Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata_BP, test.stat)
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata_BP, test.stat)

elim.ks <- runTest(GOdata_BP, algorithm = "elim", statistic = "ks")
allRes_BP <- GenTable(GOdata_BP, classic=elim.ks, KS=resultKS, weight = resultWeight,
                      orderBy = "weight", ranksOf = "classic", topNodes =15) ##只保留前15个
allRes_BP <- as.data.frame(allRes_BP)
##添加物种分组信息，保存文件
allRes_BP$Group <- "salt"
allRes_BP$ontology <- "BP"

#########################################################################################################################
###################################使用不同ontology进行富集分析,CC#######################################################

GOdata_CC <- new("topGOdata", ontology="CC", allGenes = genelist,
                 annot = annFUN.gene2GO, gene2GO = gene2go)
GOdata_CC <- na.omit(GOdata_CC)
##统计检验
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata_CC, test.stat)
test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata_CC, test.stat)
test.stat <- new("weightCount", testStatistic = GOFisherTest, name="Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata_CC, test.stat)
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata_CC, test.stat)

elim.ks <- runTest(GOdata_CC, algorithm = "elim", statistic = "ks")
allRes_CC <- GenTable(GOdata_CC, classic=elim.ks, KS=resultKS, weight = resultWeight,
                      orderBy = "weight", ranksOf = "classic", topNodes =15) ##只保留前15个
allRes_CC <- as.data.frame(allRes_CC)
##添加物种分组信息，保存文件
allRes_CC$Group <- "salt"
allRes_CC$ontology <- "CC"

#########################################################################################################################
#########################################################################################################################

allRes <- rbind(allRes_BP,allRes_MF,allRes_CC)
# print(allRes)
write.table(allRes, file = "~/R  project/207spe_ctdg/1.expression/16.Aty/DE-CTDGs_salt_topGO.tsv",
            row.name = FALSE, col.names=TRUE,quote = F,sep = '\t')

#########################################################################################################################
#########################################将物种合并绘图##################################################################

desired_order <-c("cold","drought","heat","light","salt")
df_1 <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/DE-CTDGs_cold_topGO.tsv",sep = '\t')
df_2 <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/DE-CTDGs_drought_topGO.tsv",sep = '\t')
df_3 <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/DE-CTDGs_heat_topGO.tsv",sep = '\t')
df_4 <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/DE-CTDGs_light_topGO.tsv",sep = '\t')
df_5 <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/DE-CTDGs_salt_topGO.tsv",sep = '\t')
df_combined <- rbind(df_1, df_2, df_3, df_4, df_5)

df_combined$Term <- factor(df_combined$Term, levels = unique(df_combined$Term))
df_combined$Group <- factor(df_combined$Group, levels = desired_order)
df_combined$ontology <- factor(df_combined$ontology,levels = c("MF","CC","BP")) ##根据ontology排列term，按照 ontology 列的顺序对 Term 进行因子化
df_combined <- df_combined %>%
  arrange(ontology) %>%
  mutate(Term = factor(Term, levels = unique(Term)))
# head(df_combined)
# view(df_combined)
# write.xlsx(df_combined,"~/Desktop/project/CTDG_of_207_spe/supplymentary tables/Table S5.xlsx")
# write.xlsx(df_combined,"~/R  project/207spe_ctdg/1.expression/16.Aty/GO_enrichment_of_DE-CTDGs.xlsx")
# showSigOfNodes(GOdata, score(resultWeight), firstSigNodes = 10, useInfo = "all")
##使用ggplot绘图
# 只保留ontology为“MF”的行
df_combined_1 <- df_combined %>% filter(ontology == "MF")
p1 <- ggplot(df_combined_1,aes(Group,Term))+
  geom_point(aes(color=KS, size=Significant))+theme_bw()+
  theme()+
  scale_color_gradient(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))+
  theme(legend.direction = "vertical",
        legend.position = "right",
        axis.text.x = element_text(size = 8,colour = "black",face = "bold"),
        axis.text.y = element_text(size = 8,color = "black")
        
        # legend.title = element_blank()
  )+
  # facet_wrap(~Group, nrow = 4,scales = "free_y")+
  scale_y_discrete(position = "left");p1

ggsave("DE-CTDGs_enrichment_of_MF_8.26*9.pdf", 
       plot = p1, width = 8.26, height = 9, units = "in")



