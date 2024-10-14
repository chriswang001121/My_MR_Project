# ----install package----

# install.packages("data.table")

# install.packages("ggplot2")

# install.packages("ggrepel")

# install.packages("grid")

# install.packages("readr")

# install.packages("forestploter")

# install.packages("colorspace")

# install.packages("stringi")

# install.packages("ggplot2")

# install.packages("circlize")

# install.packages("RColorBrewer")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

# BiocManager::install("org.Hs.eg.db")

# BiocManager::install("DOSE")

# BiocManager::install("clusterProfiler")

# BiocManager::install("enrichplot")

# BiocManager::install("ComplexHeatmap")




# ----library package----

library(data.table)
library(tidyverse)
library(foreach)
library(ggplot2)                         
library(ggrepel)
library(grid)
library(readr)
library(forestploter)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

# ----filter by IVW & pleiotropy & direction----

MR_Result_OR <- fread(file="MR_Result_OR.csv",sep=",",data.table = F) # input MR result
MR_Result_OR_adjust <- MR_Result_OR %>%
  mutate(adjust.p = pval*nsnp)

MR_pleiotropy <- fread(file="MR_pleiotropy.csv",sep=",",data.table = F) # pleiotropy

AllGene <- unique(MR_Result_OR$exposure)

Positive_Result_adjust.p <- data.frame()
foreach(i=AllGene, .errorhandling = "pass") %do%{
  SingleResult <- MR_Result_OR_adjust %>%
    dplyr::filter(exposure == i)
  SingleResult_ivw <- SingleResult %>%
    dplyr::filter(method == 'Inverse variance weighted')
  pleiotropy <- MR_pleiotropy %>%
    dplyr::filter(exposure == i)
  if(SingleResult_ivw$adjust.p > 0.05 | pleiotropy$pval < 0.05){
    SingleResult <- data.frame()
  }
  if(sum(SingleResult$or>1) == nrow(SingleResult) | sum(SingleResult$or<1)==nrow(SingleResult)){
    Positive_Result_adjust.p <- rbind(Positive_Result_adjust.p,SingleResult)
  }
}
write.csv(Positive_Result_adjust.p, file="Positive_Result_adjust.p.csv", row.names=F)

Positive_Result_pval <- data.frame()
foreach(i=AllGene, .errorhandling = "pass") %do%{
  SingleResult <- MR_Result_OR_adjust %>%
    dplyr::filter(exposure == i)
  SingleResult_ivw <- SingleResult %>%
    dplyr::filter(method == 'Inverse variance weighted')
  pleiotropy <- MR_pleiotropy %>%
    dplyr::filter(exposure == i)
  if(SingleResult_ivw$pval > 0.05 | pleiotropy$pval < 0.05){
    SingleResult <- data.frame()
  }
  if(sum(SingleResult$or>1) == nrow(SingleResult) | sum(SingleResult$or<1)==nrow(SingleResult)){
    Positive_Result_pval <- rbind(Positive_Result_pval ,SingleResult)
  }
}
write.csv(Positive_Result_pval, file="Positive_Result_pval.csv", row.names=F)


pleiotropy_direct <- data.frame()
foreach(i=AllGene, .errorhandling = "pass") %do%{
  SingleResult <- MR_Result_OR_adjust %>%
    dplyr::filter(exposure == i)
  SingleResult_ivw <- SingleResult %>%
    dplyr::filter(method == 'Inverse variance weighted')
  pleiotropy <- MR_pleiotropy %>%
    dplyr::filter(exposure == i)
  if(pleiotropy$pval < 0.05){
    SingleResult <- data.frame()
  }
  if(sum(SingleResult$or>1) == nrow(SingleResult) | sum(SingleResult$or<1)==nrow(SingleResult)){
    pleiotropy_direct <- rbind(pleiotropy_direct,SingleResult)
  }
}
write.csv(pleiotropy_direct, file="pleiotropy_direct.csv", row.names=F)

# ----volcano----

data <- pleiotropy_direct %>%
  dplyr::filter(method == 'Inverse variance weighted')

data$risk <- "not"
data$risk[data$pval < 0.05 & data$b > 0] <- "high"
data$risk[data$pval < 0.05 & data$b < 0] <- "low"
table(data$risk)

p <- ggplot(data,aes(b,-1*log10(pval))) +
  geom_point(aes(color=risk),size=2) +
  scale_color_manual(values = c("#E64B35B2","#00A087B2","#3C5488B2")) +
  geom_hline(yintercept = -1*log10(0.05),linetype=4,size=0.3) +
  geom_vline(xintercept = c(-0.3,0.3),linetype=4,size=0.3) +
  xlim(-1,1) +
  ylim(0,5) +
  theme_classic() +
  theme(title=element_text(size = 18),text = element_text(size=18)) +
  labs(x="Beta",y="-log10(pvalue)")

data_label <- data %>%
  filter(adjust.p < 0.05) # 添加标签

p1 <- p + geom_label_repel(data=data_label,
                     aes(label=exposure))

ggsave(p1,file = 'volcano.tiff',width = 40,height = 30,units = 'cm')


# ----forest----

data <- Positive_Result_adjust.p %>% # select method
  dplyr::filter(method == 'Inverse variance weighted' | method == 'Weighted median' | method == 'MR Egger')

data$' ' <- paste(rep(" ", 10), collapse = " ")

data$'OR(95% CI)' <- ifelse(is.na(data$or), 
                            "", 
                            sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))
        
data$pval <- ifelse(data$pval<0.001, 
                   "<0.001", 
                   sprintf("%.3f", data$pval))

data$exposure <- ifelse(is.na(data$exposure), 
                        "", 
                        data$exposure)

data$nsnp <- ifelse(is.na(data$nsnp), 
                    "", 
                    data$nsnp)

data[duplicated(data$exposure),]$exposure <- ""

head(data)

# the data is too much, keep first 15 line

data <- data[1:15,]

## 1. Simple forest plot

p <- forest(data[, c("exposure","nsnp","method","pval"," ","OR(95% CI)")],
            est = data$or,
            lower = data$or_lci95,
            upper = data$or_uci95,
            # sizes = 0.5,
            ci_column = 5,  # 可信区间所在的列
            ref_line = 1,  # 参考线条的位置，虚线
            arrow_lab = c("low", "high"),
            xlim = c(0, 3), # X轴的范围
            ticks_at = c(0.5, 1, 2, 3),
            footnote = "exposure on outcome")# 图形的参数

# Print plot
plot(p)

## 1. Simple forest plot

## 2. Change theme

tm <- forest_theme(base_size = 18, # 图形整体的大小
                   # 可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                   ci_pch = 16, 
                   ci_lty = 1, 
                   ci_lwd = 1.5, 
                   ci_col = "black", 
                   ci_alpha = 0.8,
                   ci_fill = "#E64B35",
                   ci_Theight = 0.2, # Set a T end at the end of CI
                   # 参考线条的形状、宽度、颜色
                   refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey20"),
                   # x轴刻度字体的大小
                   xaxis_cex=0.8,
                   # Vertical line width/type/color
                   vertline_lwd = 10, # Line width for extra vertical line.
                   vertline_lty = "dashed",
                   vertline_col = "orange",
                   # Change summary color for filling and borders
                   summary_fill = "#762a83",
                   summary_col = "#762a83",
                   #脚注大小、颜色
                   footnote_gp = gpar(cex = 0.6, fontface = "italic", col = "blue"))

pt <- forest(data[, c("exposure","nsnp","method","pval"," ","OR(95% CI)")],
            est = data$or,
            lower = data$or_lci95,
            upper = data$or_uci95,
            # sizes = 0.5,
            is_summary = c(rep(FALSE, nrow(data)-1), TRUE), # 最后一行特殊标记
            ci_column = 5,  # 可信区间所在的列
            ref_line = 1,  # 参考线条的位置，虚线
            arrow_lab = c("low", "high"),
            xlim = c(0, 3), # X轴的范围
            ticks_at = c(0.5, 1, 2, 3),
            footnote = "exposure on outcome", # 图形的参数
            theme = tm)
# Print plot
plot(pt)

## 3. edit text and plot

tm <- forest_theme(base_size = 18, # 图形整体的大小
                   # 可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                   ci_pch = 16, 
                   ci_lty = 1, 
                   ci_lwd = 1.5, 
                   ci_col = "black", 
                   ci_alpha = 0.8,
                   ci_fill = "black",
                   ci_Theight = 0.2, # Set a T end at the end of CI
                   # 参考线条的形状、宽度、颜色
                   refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey20"),
                   # x轴刻度字体的大小
                   xaxis_gp= gpar(cex = 0.8),
                   # Vertical line width/type/color
                   vertline_lwd = 10, # Line width for extra vertical line.
                   vertline_lty = "dashed",
                   vertline_col = "orange",
                   # Change summary color for filling and borders
                   summary_fill = "#762a83",
                   summary_col = "#762a83",
                   #脚注大小、颜色
                   footnote_gp = gpar(cex = 0.6, fontface = "italic", col = "blue"))

pt <- forest(data[, c("exposure","nsnp","method","pval"," ","OR(95% CI)")],
             est = data$or,
             lower = data$or_lci95,
             upper = data$or_uci95,
             # sizes = 0.5,
             # is_summary = c(rep(FALSE, nrow(data)-1), TRUE), # 最后一行特殊标记
             ci_column = 5,  # 可信区间所在的列
             ref_line = 1,  # 参考线条的位置，虚线
             arrow_lab = c("low", "high"),
             xlim = c(0, 3), # X轴的范围
             ticks_at = c(0.5, 1, 2, 3),
             footnote = "exposure on outcome", # 图形的参数
             theme = tm)

# edit box color
boxcolor <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148") 
boxcolor <- boxcolor[as.numeric(as.factor(data$method))]  

foreach(i = 1:nrow(data), .errorhandling = 'pass') %do% {
  pt <- edit_plot(pt, col=5,row = i, which = "ci", gp = gpar(fill = boxcolor[i],fontsize=25))
  ptc <- pt
}

# bold pval < 0.05
pos_bold_pval = which(as.numeric(gsub('<',"",data$pval))<0.05) # find row with p < 0.05
if(length(pos_bold_pval)>0){
  foreach(i = pos_bold_pval, .errorhandling = 'pass') %do% {
    ptc <- edit_plot(ptc, col=4,row = i, which = "text", gp = gpar(fontface="bold"))
    ptcp <- ptc
  }
}

# add border
lineVec <- which(data$exposure != '')
ptcpb <- add_border(ptcp, part = "header", row =1,where = "top",gp = gpar(lwd =2))|>
  add_border(part = "header", row = lineVec, gp = gpar(lwd =1))

# text align
ptcpbt <- edit_plot(ptcpb, col=1:ncol(data),row = 1:nrow(data), # all text size 
                   which = "text", gp = gpar(fontsize=12)) |>
  edit_plot(col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"),part="header",
                  x = unit(0.5, "npc")) |>
  edit_plot(col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),part="body",
                  x = unit(0.5, "npc"))

ggsave(ptcpbt,file = 'forest.tiff',width = 40,height = 30,units = 'cm')

# ----GO analyse and visualize----

Positive_Result_adjust.p <- fread(file = "Positive_Result_adjust.p.csv",data.table = FALSE)

data_GO <- Positive_Result_adjust.p %>%
  dplyr::filter(method == 'Inverse variance weighted')

Gene_GO <- unique(data_GO$exposure)

enrich_GO <- enrichGO(gene = Gene_GO,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "all",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable=T)
result_GO <- as.data.frame(enrich_GO)

result_GO_adjust.p <- result_GO %>%
  dplyr::filter(p.adjust < 0.05)

write.table(result_GO_adjust.p, file="result_GO_adjust.p.txt", sep="\t", quote=F, row.names = F)

# barplot
GO_barplot <- barplot(enrich_GO,showCategory = 10,label_format=100, 
                      split="ONTOLOGY",drop=TRUE,color="p.adjust") +
  facet_grid(ONTOLOGY~., scale='free')
ggsave(GO_barplot,file = 'GO_barplot.tiff',width = 40,height = 30,units = 'cm')

# dotplot
GO_dotplot <- dotplot(enrich_GO,showCategory = 10,label_format=100, 
                      split="ONTOLOGY",orderBy="GeneRatio",color="p.adjust") +
  facet_grid(ONTOLOGY~., scale='free')
ggsave(GO_dotplot,file = 'GO_dotplot.tiff',width = 40,height = 30,units = 'cm')

# circle----
ontology.col=c("#E64B35B2","#00A087B2","#3C5488B2")

data_circle <- result_GO_adjust.p %>%
  arrange(p.adjust)
  
datasig=data_circle[data_circle$p.adjust<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf(file="GO.circlize.pdf", width=10, height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制富集显著性pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()




# ----KEGG analyse and visualize----
Positive_Result_adjust.p <- fread(file = "Positive_Result_adjust.p.csv",data.table = FALSE)

data_KEGG <- Positive_Result_adjust.p %>%
  dplyr::filter(method == 'Inverse variance weighted')

Gene_KEGG <- unique(data_GO$exposure)

# symbol to id
Gene_KEGG_id <- mapIds(x = org.Hs.eg.db,
                       keys =  Gene_KEGG ,
                       keytype = "SYMBOL",
                       column = "ENTREZID") |>
  na.omit()

enrich_KEGG <- enrichKEGG(gene = Gene_KEGG_id,
                          organism = "hsa", pvalueCutoff=1, qvalueCutoff=1)

enrich_KEGG_symbol <- setReadable(enrich_KEGG, 
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID") # id to symbol
result_KEGG <- enrich_KEGG_symbol@result

result_KEGG_adjust.p <- result_KEGG %>%
  dplyr::filter(pvalue < 0.05)

write.table(result_KEGG_adjust.p, file="result_KEGG_adjust.p.txt", sep="\t", quote=F, row.names = F)

# barplot
KEGG_barplot <- barplot(enrich_KEGG_symbol, drop=TRUE, 
             showCategory=10, label_format=100, color='pvalue')
ggsave(KEGG_barplot,file = 'KEGG_barplot.tiff',width = 40,height = 30,units = 'cm')

# dotplot
KEGG_dotplot <- dotplot(enrich_KEGG_symbol, showCategory=10, 
                        orderBy="GeneRatio", label_format=100, color='pvalue')
ggsave(KEGG_dotplot,file = 'KEGG_dotplot.tiff',width = 40,height = 30,units = 'cm')



