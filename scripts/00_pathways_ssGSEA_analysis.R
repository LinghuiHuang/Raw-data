
################################################
########### metabolism pathway
load("raw_datas/KEGG_Pathways/metabolism.paths.genelist.RData")
length(paths.genelist) ## 84
metabolism.genelist=paths.genelist
#### hallmark基因集
h <- clusterProfiler::read.gmt("raw_datas/h.all.v2023.1.Hs.symbols.gmt")
h.all.genelist<-lapply(unique(h$term),function(x){h$gene[h$term==x]})
names(h.all.genelist) <- unique(h$term)
################################# 程序性死亡基因集
pcd.genesets=readxl::read_excel('raw_datas/PCD.geneSets.PMID36341760.xlsx')
pcd.genesets=data.frame(pcd.genesets)
pcd.genesets=pcd.genesets[,-1]
head(pcd.genesets)

pcd.list=list()
pcd.genesets.df=c()
for(i in colnames(pcd.genesets)){
  pcd.list[[i]]=as.character(na.omit(pcd.genesets[,i]))
  pcd.genesets.df=rbind(pcd.genesets.df,data.frame(PCD=i,Symbol=pcd.genesets[,i],check.names = F,stringsAsFactors = F))
}

#######################

TILs.genesets=TILs.genesets[,-1]
colnames(TILs.genesets)=NULL
head(TILs.genesets)
################
TILs.list=apply(TILs.genesets, 1, function(x){(x[which(x!="")])})
length(TILs.list)

###############
load("raw_datas/GEO/gse57065.exp.RData")
########## 代谢通路富集评分
gse57065.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse57065.exp)
                                                           , genelist = metabolism.genelist)
save(gse57065.metabolism.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse57065.metabolism.ssgsea.RData')

########## hallmark富集评分
gse57065.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse57065.exp)
                                                      , genelist = h.all.genelist)
save(gse57065.h.all.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse57065.h.all.ssgsea.RData')
########## PCD程序性死亡富集评分
gse57065.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse57065.exp
                                                    ,genelist = pcd.list)
save(gse57065.pcd.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse57065.pcd.ssgsea.RData')
########## TILs富集评分
gse57065.TILs.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse57065.exp
                                                     ,genelist = TILs.list)
save(gse57065.TILs.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse57065.TILs.ssgsea.RData')

##################
library(IOBR)
#### ESTIMATE
gse57065.exp.estimate<-deconvo_estimate(eset=gse57065.exp)
save(gse57065.exp.estimate,file='raw_datas/Pathways_ssGSEA/gse57065.exp.estimate.RData')
### CIBERSORT
gse57065.exp.cibersort<-deconvo_cibersort(eset=gse57065.exp,arrays=T)
save(gse57065.exp.cibersort,file='raw_datas/Pathways_ssGSEA/gse57065.exp.cibersort.RData')
#################################
###############
load("raw_datas/GEO/gse65682.exp.RData")
########## 代谢通路富集评分
gse65682.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse65682.exp)
                                                           , genelist = metabolism.genelist)
save(gse65682.metabolism.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse65682.metabolism.ssgsea.RData')

########## hallmark富集评分
gse65682.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse65682.exp)
                                                      , genelist = h.all.genelist)
save(gse65682.h.all.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse65682.h.all.ssgsea.RData')
########## PCD程序性死亡富集评分
gse65682.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse65682.exp
                                                    ,genelist = pcd.list)
save(gse65682.pcd.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse65682.pcd.ssgsea.RData')
########## TILs富集评分
gse65682.TILs.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse65682.exp
                                                     ,genelist = TILs.list)
save(gse65682.TILs.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse65682.TILs.ssgsea.RData')

##################
library(IOBR)
#### ESTIMATE
gse65682.exp.estimate<-deconvo_estimate(eset=gse65682.exp)
save(gse65682.exp.estimate,file='raw_datas/Pathways_ssGSEA/gse65682.exp.estimate.RData')
### CIBERSORT
gse65682.exp.cibersort<-deconvo_cibersort(eset=gse65682.exp,arrays=T)
save(gse65682.exp.cibersort,file='raw_datas/Pathways_ssGSEA/gse65682.exp.cibersort.RData')
#################################
###############
load("raw_datas/GEO/gse145227.exp.RData")
########## 代谢通路富集评分
gse145227.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse145227.exp)
                                                            , genelist = metabolism.genelist)
save(gse145227.metabolism.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse145227.metabolism.ssgsea.RData')

########## hallmark富集评分
gse145227.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse145227.exp)
                                                       , genelist = h.all.genelist)
save(gse145227.h.all.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse145227.h.all.ssgsea.RData')
########## PCD程序性死亡富集评分
gse145227.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse145227.exp
                                                     ,genelist = pcd.list)
save(gse145227.pcd.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse145227.pcd.ssgsea.RData')
########## TILs富集评分
gse145227.TILs.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse145227.exp
                                                      ,genelist = TILs.list)
save(gse145227.TILs.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse145227.TILs.ssgsea.RData')

##################
library(IOBR)
#### ESTIMATE
gse145227.exp.estimate<-deconvo_estimate(eset=gse145227.exp)
save(gse145227.exp.estimate,file='raw_datas/Pathways_ssGSEA/gse145227.exp.estimate.RData')
### CIBERSORT
# dim(gse145227.exp)
# gse145227.exp.cibersort<-deconvo_cibersort(eset=gse145227.exp,arrays=T)
# save(gse145227.exp.cibersort,file='raw_datas/Pathways_ssGSEA/gse145227.exp.cibersort.RData')
#################################
###############
load("raw_datas/GEO/gse54514.exp.RData")
########## 代谢通路富集评分
gse54514.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse54514.exp)
                                                           , genelist = metabolism.genelist)
save(gse54514.metabolism.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse54514.metabolism.ssgsea.RData')

########## hallmark富集评分
gse54514.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse54514.exp)
                                                      , genelist = h.all.genelist)
save(gse54514.h.all.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse54514.h.all.ssgsea.RData')
########## PCD程序性死亡富集评分
gse54514.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse54514.exp
                                                    ,genelist = pcd.list)
save(gse54514.pcd.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse54514.pcd.ssgsea.RData')
########## TILs富集评分
gse54514.TILs.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse54514.exp
                                                     ,genelist = TILs.list)
save(gse54514.TILs.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse54514.TILs.ssgsea.RData')

##################
library(IOBR)
#### ESTIMATE
gse54514.exp.estimate<-deconvo_estimate(eset=gse54514.exp)
save(gse54514.exp.estimate,file='raw_datas/Pathways_ssGSEA/gse54514.exp.estimate.RData')
### CIBERSORT
gse54514.exp.cibersort<-deconvo_cibersort(eset=gse54514.exp,arrays=T)
save(gse54514.exp.cibersort,file='raw_datas/Pathways_ssGSEA/gse54514.exp.cibersort.RData')
#################################
#################################
###############
load("raw_datas/GEO/gse95233.exp.RData")
########## 代谢通路富集评分
gse95233.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse95233.exp)
                                                           , genelist = metabolism.genelist)
save(gse95233.metabolism.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse95233.metabolism.ssgsea.RData')

########## hallmark富集评分
gse95233.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse95233.exp)
                                                      , genelist = h.all.genelist)
save(gse95233.h.all.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse95233.h.all.ssgsea.RData')
########## PCD程序性死亡富集评分
gse95233.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse95233.exp
                                                    ,genelist = pcd.list)
save(gse95233.pcd.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse95233.pcd.ssgsea.RData')
########## TILs富集评分
gse95233.TILs.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse95233.exp
                                                     ,genelist = TILs.list)
save(gse95233.TILs.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse95233.TILs.ssgsea.RData')

##################
library(IOBR)
#### ESTIMATE
gse95233.exp.estimate<-deconvo_estimate(eset=gse95233.exp)
save(gse95233.exp.estimate,file='raw_datas/Pathways_ssGSEA/gse95233.exp.estimate.RData')
### CIBERSORT
gse95233.exp.cibersort<-deconvo_cibersort(eset=gse95233.exp,arrays=T)
save(gse95233.exp.cibersort,file='raw_datas/Pathways_ssGSEA/gse95233.exp.cibersort.RData')
#################################
############### GS185263
load("raw_datas/GEO/gse185263.exp.RData")
dim(gse185263.exp)
########## 代谢通路富集评分
gse185263.metabolism.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse185263.exp)
                                                            , genelist = metabolism.genelist)
save(gse185263.metabolism.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse185263.metabolism.ssgsea.RData')

########## hallmark富集评分
gse185263.h.all.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = as.matrix(gse185263.exp)
                                                       , genelist = h.all.genelist)
save(gse185263.h.all.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse185263.h.all.ssgsea.RData')
########## PCD程序性死亡富集评分
gse185263.pcd.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse185263.exp
                                                     ,genelist = pcd.list)
save(gse185263.pcd.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse185263.pcd.ssgsea.RData')
########## TILs富集评分
gse185263.TILs.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = gse185263.exp
                                                      ,genelist = TILs.list)
save(gse185263.TILs.ssgsea,file = 'raw_datas/Pathways_ssGSEA/gse185263.TILs.ssgsea.RData')

