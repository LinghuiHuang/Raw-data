
########################################
## GSE57065 dataset
# GSE57065=getGEOExpDataByCel('GSE57065')
# save(GSE57065,file='raw_datas/GEO/GSE57065.RData')
load('raw_datas/GEO/GSE57065.RData')
########## 基因表达谱整理
# gse57065.exp<-exp_probe2symbol_v2(GSE57065$Exp$GPL570_Data,GPL='GPL570')
# save(gse57065.exp,file='raw_datas/GEO/gse57065.exp.RData')
load('raw_datas/GEO/gse57065.exp.RData')
gse57065.exp=gse57065.exp[rownames(gse57065.exp) %in% gencode.pcg$gene_name,]
dim(gse57065.exp) ## 17283*107

########## 临床信息整理
gse57065.cli<-GSE57065$Sample
rownames(gse57065.cli)=gse57065.cli$Acc
###########
gse57065.cli$patient=gsub("(.*),.*","\\1",gse57065.cli$Source)
gse57065.cli$disease=gsub("(.*?) .*","\\1",gse57065.cli$Source)
gse57065.cli$disease[which(gse57065.cli$disease=='patient')]='sepsis'
gse57065.cli$disease[which(gse57065.cli$disease=='Temoin')]='healthy'
##################
table(gse57065.cli$disease,gse57065.cli$sapsii)
#           High Low NA
# healthy    0   0 25
# sepsis    40  42  0
gse57065.cli=gse57065.cli[,c('Acc','disease','patient','age','gender','sapsii','collection time')]
gse57065.cli$gender[which(gse57065.cli$gender=="FEMALE")]='Female'
gse57065.cli$gender[which(gse57065.cli$gender=="MALE")]='Male'
gse57065.cli$sapsii[which(gse57065.cli$sapsii=='SAPSII-High')]='High'
gse57065.cli$sapsii[which(gse57065.cli$sapsii=='SAPSII-Low')]='Low'
gse57065.cli$`collection time`[which(gse57065.cli$`collection time`=='NULL')]='healthy'
gse57065.cli[which(gse57065.cli=="NA",arr.ind = T)]=NA
gse57065.cli$Age1=ifelse(gse57065.cli$age>65,">65","<=65")
############
table(gse57065.cli$disease)
head(gse57065.cli)
save(gse57065.cli,file = 'raw_datas/GEO/gse57065.cli.RData')

########################################
## GSE54514 dataset
# GSE54514=getGEOExpData('GSE54514')
# save(GSE54514,file='raw_datas/GEO/GSE54514.RData')
load('raw_datas/GEO/GSE54514.RData')
########## 基因表达谱整理
# gse54514.exp<-exp_probe2symbol_v2(GSE54514$Exp$GPL6947_24840_Data_col1,GPL='GPL6947')
# save(gse54514.exp,file='raw_datas/GEO/gse54514.exp.RData')
load('raw_datas/GEO/gse54514.exp.RData')
gse54514.exp=gse54514.exp[rownames(gse54514.exp) %in% gencode.pcg$gene_name,]
dim(gse54514.exp) ## 13114*163

########## 临床信息整理
gse54514.cli<-GSE54514$Sample
rownames(gse54514.cli)=gse54514.cli$Acc
###########
gse54514.cli=gse54514.cli[,c('Acc','age (years)','gender','disease status','group_day','group_id','severity (apacheii)','neutrophil proportion','site of infection')]
gse54514.cli=reshape::rename(gse54514.cli,c('age (years)'='age','disease status'='disease','group_id'='patient','group_day'='collection time'))
head(gse54514.cli)
table(gse54514.cli$`collection time`)
gse54514.cli$`collection time`=gsub(".*_(.*)","\\1",gse54514.cli$`collection time`)
head(gse54514.cli)

gse54514.cli$OS.time=gse54514.cli$`collection time`
gse54514.cli$OS.time=gsub("D","",gse54514.cli$OS.time)
mode(gse54514.cli$OS.time)="numeric"

gse54514.cli$OS=gse54514.cli$disease
gse54514.cli$OS[which(gse54514.cli$OS=='healthy')]=NA
gse54514.cli$OS[which(gse54514.cli$OS=='sepsis nonsurvivor')]=1
gse54514.cli$OS[which(gse54514.cli$OS=='sepsis survivor')]=0
mode(gse54514.cli$OS)="numeric"
head(gse54514.cli)
save(gse54514.cli,file = 'raw_datas/GEO/gse54514.cli.RData')

gse54514.exp=gse54514.exp[,rownames(gse54514.cli)]
dim(gse54514.exp) ## 13114*163

########################################
## GSE65682 dataset
# GSE65682=getGEOExpDataByCel('GSE65682')
# save(GSE65682,file='raw_datas/GEO/GSE65682.RData')
load('raw_datas/GEO/GSE65682.RData')
########## 基因表达谱整理
# gse65682.exp<-exp_probe2symbol_v2(GSE65682$Exp$GPL13667_Data,GPL='GPL13667')
# save(gse65682.exp,file='raw_datas/GEO/gse65682.exp.RData')
load('raw_datas/GEO/gse65682.exp.RData')
gse65682.exp=gse65682.exp[rownames(gse65682.exp) %in% gencode.pcg$gene_name,]
dim(gse65682.exp) ## 18028*802

########## 临床信息整理
gse65682.cli<-GSE65682$Sample
rownames(gse65682.cli)=gse65682.cli$Acc
###########
table(gse65682.cli$endotype_cohort)
# discovery         NA validation 
# 263        323        216 
table(gse65682.cli$endotype_cohort,gse65682.cli$`pneumonia diagnoses`)
#             cap hap  NA no-cap
# discovery   64  41 158      0
# NA           2   7 281     33
# validation  42  36 138      0
gse65682.cli$disease=gsub(" \\[.*\\]","",gse65682.cli$Title)
table(gse65682.cli$disease)
# healthy subject intensive-care unit patient 
# 42                         760 
gse65682.cli$disease[which(gse65682.cli$disease=='intensive-care unit patient')]='sepsis'
gse65682.cli$disease[which(gse65682.cli$disease=='healthy subject')]='healthy'

# 4个脓毒症亚型，即Mars 1~ 4型，其中Mars 1型病情最重，发生休克的比例最高，病死率也显著高于其他3个亚型
gse65682.cli=gse65682.cli[,c('Acc','disease','age','gender','pneumonia diagnoses','thrombocytopenia',
                             'endotype_cohort','endotype_class','mortality_event_28days','time_to_event_28days'
                             ,'icu_acquired_infection','icu_acquired_infection_paired','diabetes_mellitus','abdominal_sepsis_and_controls')]
head(gse65682.cli)
table(gse65682.cli$gender)
gse65682.cli$gender[which(gse65682.cli$gender=='male')]='Male'
gse65682.cli$gender[which(gse65682.cli$gender=='female')]='Female'
gse65682.cli[which(gse65682.cli=='NA',arr.ind = T)]=NA

gse65682.cli$survival=gse65682.cli$mortality_event_28days
gse65682.cli$survival[which(gse65682.cli$survival==0)]='Survivor'
gse65682.cli$survival[which(gse65682.cli$survival==1)]='Non survivor'

gse65682.cli$OS.time=gse65682.cli$time_to_event_28days
gse65682.cli$OS=gse65682.cli$mortality_event_28days

gse65682.cli$Age1=ifelse(gse65682.cli$age>65,">65","<=65")
head(gse65682.cli)
save(gse65682.cli,file = 'raw_datas/GEO/gse65682.cli.RData')
writexl::write_xlsx(gse65682.cli,path = 'raw_datas/GEO/gse65682.cli.xlsx')

gse65682.exp=gse65682.exp[,rownames(gse65682.cli)]
dim(gse65682.exp) ## 18028*802
########################################
## GSE145227 dataset
# GSE145227=getGEOExpData('GSE145227')
# save(GSE145227,file='raw_datas/GEO/GSE145227.RData')
load('raw_datas/GEO/GSE145227.RData')
########## 基因表达谱整理
# dt <- GSE145227$Anno$GPL23178
# dt$genes <- stringr::str_split_fixed(dt$V7, ' // ', 3)[, 2]
# gse145227.exp<-exp_probe2symbol_v2(GSE145227$Exp$GPL23178_95722_Data_col1,anno = dt[,c(1,12)])
# save(gse145227.exp,file='raw_datas/GEO/gse145227.exp.RData')
# load('raw_datas/GEO/gse145227.exp.RData')
# gse145227.exp=gse145227.exp[rownames(gse145227.exp) %in% gencode.pcg$gene_name,]
# dim(gse145227.exp) ## 17595*22
# save(gse145227.exp,file='raw_datas/GEO/gse145227.exp.RData')
load('raw_datas/GEO/gse145227.exp.RData')

########## 临床信息整理
gse145227.cli<-GSE145227$Sample
rownames(gse145227.cli)=gse145227.cli$Acc
###########
gse145227.cli=gse145227.cli[,c('Acc','Source')]
gse145227.cli$Source=gsub("[0-9]+","",gse145227.cli$Source)
head(gse145227.cli)
save(gse145227.cli,file = 'raw_datas/GEO/gse145227.cli.RData')

gse145227.exp=gse145227.exp[,rownames(gse145227.cli)]
dim(gse145227.exp) ## 17283*107

########################################
## GSE95233 dataset
# GSE95233=getGEOExpDataByCel('GSE95233')
# save(GSE95233,file='raw_datas/GEO/GSE95233.RData')
load('raw_datas/GEO/GSE95233.RData')
########## 基因表达谱整理
# gse95233.exp<-exp_probe2symbol_v2(GSE95233$Exp$GPL570_Data,GPL='GPL570')
# save(gse95233.exp,file='raw_datas/GEO/gse95233.exp.RData')
load('raw_datas/GEO/gse95233.exp.RData')
dim(gse95233.exp)
gse95233.exp=gse95233.exp[rownames(gse95233.exp) %in% gencode.pcg$gene_name,]
dim(gse95233.exp) ## 13114*163

########## 临床信息整理
gse95233.cli<-GSE95233$Sample
rownames(gse95233.cli)=gse95233.cli$Acc
###########
gse95233.cli=gse95233.cli[,c('Acc','Source','age','gender','time point','survival')]
gse95233.cli$Source=gsub("(.*?) .*","\\1",gse95233.cli$Source)
gse95233.cli$Source[which(gse95233.cli$Source=='Patient')]='sepsis'
gse95233.cli$Source[which(gse95233.cli$Source=='Control')]='healthy'
gse95233.cli=reshape::rename(gse95233.cli,c('Source'='disease'))
table(gse95233.cli$disease)
head(gse95233.cli)

gse95233.cli$OS.time=gse95233.cli$`time point`
table(gse95233.cli$OS.time)
gse95233.cli$OS.time[which(gse95233.cli$OS.time=='D01')]=1
gse95233.cli$OS.time[which(gse95233.cli$OS.time=='D02')]=2
gse95233.cli$OS.time[which(gse95233.cli$OS.time=='D03')]=3
gse95233.cli$OS.time[which(gse95233.cli$OS.time=='NA')]=NA
gse95233.cli$OS.time=as.numeric(gse95233.cli$OS.time)

gse95233.cli$OS=gse95233.cli$survival
table(gse95233.cli$OS)
gse95233.cli$OS[which(gse95233.cli$OS=="Survivor")]=0
gse95233.cli$OS[which(gse95233.cli$OS=="Non Survivor")]=1
gse95233.cli$OS[which(gse95233.cli$OS=="NA")]=NA
gse95233.cli$OS=as.numeric(gse95233.cli$OS)
gse95233.cli$Age1=ifelse(gse95233.cli$age > 65,">65","<=65")
gse95233.cli$gender[which(gse95233.cli$gender=="F")]="Female"
gse95233.cli$gender[which(gse95233.cli$gender=="M")]="Male"
gse95233.cli$gender[which(gse95233.cli$gender=="H")]=NA
save(gse95233.cli,file = 'raw_datas/GEO/gse95233.cli.RData')

gse95233.exp=gse95233.exp[,rownames(gse95233.cli)]
dim(gse95233.exp) ## 17283*124

## GSE185263 dataset
# library(GEOquery)
# gse=getGEO('GSE185263')
# gse185263.cli=pData(gse[[1]])
# gse185263.cli=gse185263.cli[,c(1:2,46:53)]
# colnames(gse185263.cli)=gsub(":.*","",colnames(gse185263.cli))
# gse185263.cli=reshape::rename(gse185263.cli,c("disease state"='disease'))
# save(gse185263.cli,file = 'raw_datas/GEO/gse185263.cli.online.RData')
load("raw_datas/GEO/gse185263.cli.online.RData")
gse185263.cli=reshape::rename(gse185263.cli,c("in hospital mortality"="survival"))
table(gse185263.cli$survival)
gse185263.cli$survival[which(gse185263.cli$survival=='Survived')]="Survivor"
gse185263.cli$survival[which(gse185263.cli$survival=='Died')]="Non Survivor"
gse185263.cli$survival[which(gse185263.cli$survival=='NA')]=NA

gse185263.cli$OS=gse185263.cli$survival
table(gse185263.cli$OS)
gse185263.cli$OS[which(gse185263.cli$OS=="Non Survivor")]=1
gse185263.cli$OS[which(gse185263.cli$OS=="Survivor")]=0

gse185263.cli$OS.time=1
gse185263.cli$OS=as.numeric(gse185263.cli$OS)
gse185263.cli$Age1=ifelse(gse185263.cli$age>65,">65","<=65")
head(gse185263.cli)
gse185263.cli=reshape::rename(gse185263.cli,c('Sex'='gender',"sofa 24h post admisssion"="Sofa score"))
save(gse185263.cli,file = 'raw_datas/GEO/gse185263.cli.RData')

########## 基因表达谱整理
# gse185263.exp=read.csv('raw_datas/GEO/GSE185263_raw_counts.csv')
# rownames(gse185263.exp)=gse185263.exp$ensembl_gene_id
# gse185263.exp=gse185263.exp[,-1]
# gse185263.exp=mg_count2TPMs(gse185263.exp)
# gse185263.exp=log2(gse185263.exp+1)
# gse185263.exp=exp_ensg2symbol(gse185263.exp)
# colnames(gse185263.exp)=gse185263.cli$geo_accession[match(colnames(gse185263.exp),gse185263.cli$title)]
# dim(gse185263.exp)
# gse185263.exp[1:4,1:5]
# save(gse185263.exp,file = 'raw_datas/GEO/gse185263.exp.RData')
load('raw_datas/GEO/gse185263.exp.RData')
dim(gse185263.exp)
gse185263.exp=gse185263.exp[rownames(gse185263.exp) %in% gencode.pcg$gene_name,]
dim(gse185263.exp) ## 19040*392

## GSE110487 dataset
# library(GEOquery)
# gse=getGEO('GSE110487')
# all(colnames(pData(gse[[1]]))==colnames(pData(gse[[2]])))
# gse110487.cli=rbind(pData(gse[[1]]),pData(gse[[2]]))
# gse110487.cli=gse110487.cli[,c(1:2,42:44)]
# colnames(gse110487.cli)=gsub(":.*","",colnames(gse110487.cli))
# gse110487.cli=reshape::rename(gse110487.cli,c("clinical classification"='Response','title'='ID'))
# save(gse110487.cli,file = 'raw_datas/GEO/gse110487.cli.RData')
load("raw_datas/GEO/gse110487.cli.RData")

########## 基因表达谱整理
# gse110487.exp=readxl::read_xlsx('raw_datas/GEO/GSE110487_rawcounts.xlsx')
# gse110487.exp=data.frame(gse110487.exp)
# rownames(gse110487.exp)=gse110487.exp$Geneid
# gse110487.exp=gse110487.exp[,-1]
# gse110487.exp=mg_count2TPMs(gse110487.exp)
# gse110487.exp=log2(gse110487.exp+1)
# gse110487.exp=exp_ensg2symbol(gse110487.exp)
# colnames(gse110487.exp)=gse110487.cli$geo_accession[match(colnames(gse110487.exp),gse110487.cli$ID)]
# dim(gse110487.exp)
# gse110487.exp[1:4,1:5]
# save(gse110487.exp,file = 'raw_datas/GEO/gse110487.exp.RData')
load('raw_datas/GEO/gse110487.exp.RData')
dim(gse110487.exp) ## 34250*62
gse110487.exp=gse110487.exp[rownames(gse110487.exp) %in% gencode.pcg$gene_name,]
dim(gse110487.exp) ## 19040*392
