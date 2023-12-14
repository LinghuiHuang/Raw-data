############################## WGCNA
##### cluster gene
library(WGCNA)
################# GSE57065
load("raw_datas/GEO/gse57065.exp.RData")
load("raw_datas/GEO/gse57065.cli.RData")
all(colnames(gse57065.exp)==rownames(gse57065.cli))
gse57065.cli_dt=data.frame(gse57065.cli,NETosis=gse57065.pcd.ssgsea[rownames(gse57065.cli),'Netotic.cell.death'])
table(gse57065.cli_dt$disease)
# healthy  sepsis 
# 25      82
dim(gse57065.cli_dt) ## 107*9
####################
dim(gse57065.exp) ## 20549*107
mads=apply(gse57065.exp[,rownames(gse57065.cli_dt)],1,mad)
pcg.selected=names(mads[order(mads,decreasing = T)[1:(length(mads)*0.7)]])
length(pcg.selected) ## 14384
###############
edata.wgcna=gse57065.exp[pcg.selected,rownames(gse57065.cli_dt)]
edata.wgcna=t(edata.wgcna)
dim(edata.wgcna) ## 107 * 14384
range(edata.wgcna)
## 1.334591 14.234079

########
mg_wgcna_get_power=function(exp,RsquaredCut=0.85,height=0,net_type='unsigned',blockSize = 7000){
  library(WGCNA)
  logs=c()
  #print('=====')
  filter.exp=exp
  sampleTree = hclust(dist(filter.exp), method = "average")
  #pdf(file = paste0(MG_GENE_FOLDER,'/SampleCluster.pdf'),width = 12,height = 6)
  #plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  #dev.off()
  #height=140000
  #print(height)
  #print(min(sampleTree$height))
  if(height>min(sampleTree$height)){
    cte=cutree(tree = sampleTree,h = height)
    #plot(sampleTree)
    #abline(h=height)
    s.inds=match(names(which(cte==names(which.max(table(cte))))),row.names(filter.exp))
    #sum(sampleTree$height<height)
    logs=c(logs,paste0('remove samples:',(nrow(filter.exp)-length(s.inds)),' row count')) 
    filter.exp=filter.exp[s.inds,]
    #print('==0')
  }
  #print('===1')
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  #print('===2')
  #print(head(filter.exp))
  sft = pickSoftThreshold(filter.exp, powerVector=powers, 
                          RsquaredCut=RsquaredCut,
                          networkType=net_type, verbose=5,blockSize = blockSize)
  #print('===3')
  cutPower=sft$powerEstimate
  if(is.na(cutPower)){
    cutPower=0
  }
  print(paste0('power succed:',cutPower,',starting plot'))
  logs=c(logs,paste0('power succed:',cutPower,',starting plot')) 
  
  #library(customLayout)
  #lay1 <- lay_new(
  #  matrix(1:1, nc = 1),
  #  widths = c(5),
  #  heights = c(4))
  #lay2 <- lay_new(
  #  matrix(1:2, nc = 2),
  #  widths = c(2.5,2.5),
  #  heights = c(4))
  
  #cl = lay_bind_row(lay1, lay2,heights = c(0.8,1))
  #lay_set(cl)
  
  layout(matrix(c(1,1,2,3),2,2,byrow=T),widths=c(1,1),heights=c(0.8,1))
  
  mai=par('mai')
  
  #pdf(file = paste0(MG_GENE_FOLDER,'/powers.pdf'),width = 9,height = 8)
  mai1=mai
  mai1[1]=0
  par(mai=mai1)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  if(height>min(sampleTree$height)){
    abline(h=height,col='red')
  }
  par(mai=mai)
  
  #lay_show(cl)
  #par(mfrow = c(1,2))
  
  cex1 = rep(0.9,length(sft$fitIndices[,1]))
  col1=rep('red',length(sft$fitIndices[,1]))
  if(cutPower>0){
    cex1[which(sft$fitIndices[,1]==cutPower)]=1.5
    col1[which(sft$fitIndices[,1]==cutPower)]='blue'
  }
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col=col1)
  # abline(h=0.85,col="red")
  abline(h=RsquaredCut ,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
       cex=cex1, col=col1)
  arg=paste0('height=',height,',net_type=',net_type,',blockSize = ',blockSize)
  return(list(cutPower=cutPower,log=logs,Exp=filter.exp,arg=arg))  
}

pdf('analysis/3_WGCNA/Fig3ABC.pdf',width = 8,height = 8)
power=mg_wgcna_get_power(edata.wgcna,RsquaredCut=0.9)
dev.off()

power$cutPower ## 7
net=mg_WGCNA_getModule(edata.wgcna,power = power$cutPower
                       , deepSplit = 1, mergeCutHeight = 0.25
                       , minModuleSize = 60)
length(table(net$Modules[,2])) ## 19

pdf('analysis/3_WGCNA/Fig3D.pdf',height = 6,width = 8)
plotDendroAndColors(net$Tree, net$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

fig3e=mg_barplot_point(labels = names(table(net$Modules[,2]))
                       ,values = as.numeric(table(net$Modules[,2]))
                       ,point_sizes = 2
                       ,point_cols = names(table(net$Modules[,2]))
                       ,xlab = 'Number of Genes',legend.pos = NULL)
fig3e=fig3e+theme(axis.text.y = element_text(family = "Times", face = "plain"),
                  axis.title.y = element_text(family = "Times", face = "plain"),
                  panel.background = element_rect(fill = "white",
                                                  colour = "black"),
                  plot.background = element_rect(fill = "white",
                                                 colour = "black")
)
fig3e 
savePDF('analysis/3_WGCNA/Fig3E.pdf',fig3e,height = 5,width = 4)

#### 模块特征向量聚类
# Calculate eigengenes
MEs = net$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('analysis/3_WGCNA/Fig3F.pdf',height = 6,width = 10,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()

############### GSE57065
dim(gse57065.cli_dt)
datTraits0 = data.frame(type=gse57065.cli_dt$disease,NETosis=gse57065.cli_dt$NETosis,gse57065.cli_dt[,c("age",'gender')],stringsAsFactors = F)
datTraits0$type=factor(datTraits0$type,levels = c('healthy','sepsis'))
datTraits0$gender=factor(datTraits0$gender,levels = c("Female","Male"))
str(datTraits0)
######
datTraits=datTraits0
datTraits[,c('type',"gender")]=sapply(datTraits[, c('type',"gender")], function(x)as.numeric(as.factor(x)))
head(datTraits)
head(datTraits0)
dim(datTraits)
####### Calculate module eigengenes
MEs<-net$MEs
dim(MEs)
## Define numbers of genes and samples
nGenes = ncol(edata.wgcna)
nSamples = nrow(edata.wgcna)
## 计算module和clinical data的相关性.
modTraitCor = WGCNA::cor(MEs[,rownames(MEDiss)[METree$order]]
                         , datTraits
                         , use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, nSamples)
## 对correlation进行可视化:
textMatrix = paste(signif(modTraitCor, 2), "(", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

dev.off()
pdf('analysis/3_WGCNA/Fig3G.pdf',width = 6,height = 6)
labeledHeatmap(Matrix = data.frame(modTraitCor), 
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = rownames(modTraitCor), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix), setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1,1),xLabelsAngle =45,xLabelsAdj =1,
               main = paste("Module-trait relationships"))
dev.off()

#基因的模块成员度（module membership）计算 
#即各基因表达值与相应模块特征基因的相关性，其衡量了基因在全局网络中的位置
## Notes: signedKME函数的列名是根据输入数据的列名称从第3个字母开始截取，所以如果前面的模块的名称前面的ME如果已经去除，记得在前面加上两个字母
##计算每个基因和module的相关性
dim(edata.wgcna) ## 107*14384
dim(net$MEs) ## 107 * 19
geneModuleMembership <- as.data.frame(signedKME(edata.wgcna, data.frame(net$MEs), outputColumnName = ""))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

#计算基因和clinical data的相关性
all(rownames(edata.wgcna)==rownames(datTraits))
geneTraitSignificance <- as.data.frame(cor(edata.wgcna, datTraits, use = 'pairwise.complete.obs'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

###模块内分析：鉴定具有高GS和MM的基因
### Intramodular analysis: identifying genes with high GS and MM
############
get_module_hub_genes=function(net=NULL,module=NULL,trait=NULL,output=NULL,MM=0.6,GS=0.5,pval=0.05){
  modNames = substring(names(net$MEs), 3)
  moduleColors = unname(net$Modules[,2])
  #########
  module = module
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  pdf(output,height = 6,width = 6)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, trait]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",trait),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                     , col = module,lwd=2)
  abline(v = MM, col = "red", lwd = 2, lty = 2)
  abline(h = GS, col = "red", lwd = 2, lty = 2)
  dev.off()
  print(colnames(geneModuleMembership)[column])
  inds=(abs(geneModuleMembership[moduleGenes, column])>MM & MMPvalue[moduleGenes, column]<pval) & (abs(geneTraitSignificance[moduleGenes, trait])>GS & GSPvalue[moduleGenes, trait]<pval)
  hub.genes=(rownames(geneModuleMembership)[moduleGenes])[inds]
  print(length(hub.genes))
  return(hub.genes)
}

table(net$Modules[,2])

rownames(modTraitCor)
#### module-GENE
table(net$Modules[,2])
blue.module.genes=rownames(net$Modules)[net$Modules[,2]=='blue']
length(blue.module.genes) ## 1828

df=net$Modules[,2][net$Modules[,2] %in% c('blue')]
df=data.frame(df)
colnames(df)='Module'
df=df %>% tibble::rownames_to_column('Gene')
dim(df) ## 1828*2
head(df)
write.table(df,file = 'results/Files/gse65682.wgcna.blue.module.genes.txt',sep = '\t',quote = F,row.names = F,col.names = T)
write.table(df,file = 'analysis/3_WGCNA/gse65682.wgcna.blue.module.genes.txt',sep = '\t',quote = F,row.names = F,col.names = T)

########### module-blue
blue.hub.genes=get_module_hub_genes(net = net,module = "blue",trait = "NETosis"
                                         ,output = "analysis/3_WGCNA/wgcna.blue.scatterplot.pdf"
                                         ,MM=0.4,GS=0.4,pval=0.05)
length(blue.hub.genes) ## 735 相关
#########################
modNames = substring(names(MEs), 3)
moduleColors = unname(net$Modules[,2])
modNames
moduleColors

module = "blue"
column = match(module, modNames)
column
moduleGenes = moduleColors==module
trait="NETosis"
table(moduleColors)[module]
colnames(geneModuleMembership)[column]

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, trait]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ",trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
table(moduleColors)[module]

######################### blue module
dev.off()
pdf('analysis/3_WGCNA/gse57065.module.degs.venn.pdf',height = 7,width = 7)
mg_venn_plot(list(Up=gse57065.degs.up,Blue=blue.hub.genes,Down=gse57065.degs.dwn),fill=c(mycolors[4],'blue',mycolors[5]))
dev.off()

###################
wgcna.degs.genes=Reduce(intersect,list(Blue=blue.hub.genes,DEGs=c(gse57065.degs.up,gse57065.degs.dwn)))
length(wgcna.degs.genes) ## 231
write.table(data.frame(wgcna.degs.genes),'analysis/4_PPI/wgcna.degs.genes.txt',sep = "\t",row.names = F,col.names = F,quote = F)
