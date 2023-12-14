mg_lasso_cox_use=function(dat,time,event,nfolds=3,lambda.min=T,show_text=T,figLabels=c('A','B')){
  library("glmnet") 
  library('survival')
  t.inds=which(!is.na(time)&!is.na(event)&time>0)
  dat=dat[t.inds,]
  time=as.numeric(time[t.inds])
  event=as.numeric(event[t.inds])
  y=Surv(time,event)
  # set.seed(123456)
  set.seed(56565656)
  fit1_cv = cv.glmnet(as.matrix(dat), y, family = "cox", nfolds=nfolds
                      ,nlambda=100, alpha=1
  )
  fit<-glmnet(dat, y, family = "cox")
  if(lambda.min){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=fit1_cv$lambda.1se
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  genes=row.names(coefficients)[Active.Index]
  Active.coefficients<-coefficients[Active.Index]  
  g=mg_plot_lasso(fit,fit1_cv,lambda = lambda,show_text=show_text,figLabels=figLabels)
  return(list(Mode1=fit,Model2=fit1_cv,Genes=genes,Coef=Active.coefficients,lambda=lambda,plot=g))
}

mg_plot_lasso_use <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
  if(is.null(lambda)){
    lmda=cv_fit$lambda.min
  }else{
    lmda=lambda
  }
  fit.coef=fit$beta[(apply(fit$beta,1,function(x){
    return(sum(x!=0))
  })>0),]
  
  fit.coef=as.matrix(fit.coef)
  colnames(fit.coef)=fit$lambda
  #fit$lambda==cv_fit$lambda
  library(ggplot2)
  dat=data.table::melt(t(as.matrix(fit.coef)))
  dat_z=dat[which(dat$value==0),]
  dat=dat[which(dat$value!=0),]
  dat.sv=rbind()
  for (u in unique(dat_z[,2])) {
    t.z=dat_z[which(dat_z[,2]==u),1]
    t.zx=max(t.z)
    dat.sv=rbind(dat.sv,c(t.zx,u,0))
    t.zn=min(t.z)
    if(t.zx!=t.zn){
      dat.sv=rbind(dat.sv,c(t.zn,u,0))
    }
  }
  colnames(dat.sv)=colnames(dat_z)
  #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
  dat=crbind2DataFrame(rbind(dat,dat.sv))
  mn=min(-log(dat$Var1))
  mx=max(-log(dat$Var1))
  if(show_text){
    mx=(mx-mn)*0.1+mx
  }
  p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
  p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
  if(show_text){
    fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
    for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
    p=p+ggrepel::geom_label_repel(
      aes(label = Var2,color=Var2),
      data = for_label,hjust = 0
    )
  }
  p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
  p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
  tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                 ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
  p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
    geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
    geom_point(aes(colour=col))
  p1=p1+theme_bw()+theme(legend.position = "none")
  # gal=p+p1
  gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                        ,align = "hv"
                        ,labels = figLabels)
  return(gal)
}

mg_surv_pROC.temp=function(time,status,score,mks=c(1,3,5)){
  # mks=mg_predict_time_ymd(time,mks)
  library(pROC)
  library(survival)
  dat=data.frame(time,status,score)
  dat=dat[!is.na(time)&!is.na(status)&!is.na(score),]
  colnames(dat)=c('time','status','score')
  plots=list()
  AUCs=c()
  cols=mg_colors
  print(mks)
  for (mk in mks) {
    dat1=dat
    dat1$year=rep(NA,nrow(dat))
    dat1$year[dat1$time >=mk & dat1$status == 0] =0
    dat1$year[dat1$time >=mk & dat1$status == 1] =0
    dat1$year[dat1$time < mk & dat1$status == 0] =NA
    dat1$year[dat1$time < mk & dat1$status == 1] =1
    dat1=dat1[which(!is.na(dat1$year)),]
    #print(head(dat1))
    #print(length(plots)>0)
    pt=mg_surv_pROC_smooth_trycatch(dat1$year,dat1$score,cols[length(plots)+1],add = length(plots)>0)
    #print(pt)
    AUCs=c(AUCs,round(pt$auc,2))
    #print(AUCs)
    plots=c(plots,list(pt))
    
  }
  if(max(mks)<20){
    lb=paste0(mks,'-Years')
  }else if(max(mks)<365){
    lb=paste0(round(mks/12,0),'-Years')
  }else{
    lb=paste0(round(mks/365,0),'-Years')
  }
  lb=paste0(mks,'-Weeks')
  legend('bottomright',paste0(lb,',AUC=',AUCs),col = cols[1:length(AUCs)],
         lty = c(1,1),border = '',title = 'AUC',bty="n")
  return(plots)
}
ggplotTimeROC=function(time,status,score,mks=c(1,3,5),pal=NULL){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC.temp(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  # mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lb=paste0(ROC.DSST$times,'-Weeks')
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    ########
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    ###############################
    if(is.null(pal)){
      pal=c(pal_npg(alpha =0.8)(9)[c(3,4,1,9)]) ## 4种
    }
    
    custome_theme=function (base_size = 12, base_family = "", border = FALSE, 
                            margin = TRUE, legend = c("top", "bottom", "left", 
                                                      "right", "none"), x.text.angle = 0) 
    {
      half_line <- base_size/2
      if (!is.numeric(legend)) 
        legend <- match.arg(legend)
      if (x.text.angle > 5) 
        xhjust <- 1
      else xhjust <- NULL
      if (border) {
        panel.border <- element_rect(fill = NA, colour = "black", 
                                     size = 0.7)
        axis.line <- element_blank()
      }
      else {
        panel.border <- element_blank()
        axis.line = element_line(colour = "black", size = 0.5)
      }
      if (margin) 
        plot.margin <- margin(half_line, half_line, half_line, 
                              half_line)
      else plot.margin <- unit(c(0.5, 0.3, 0.3, 0.3), "mm")
      .theme <- theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(panel.border = panel.border, panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), axis.line = axis.line, 
              axis.text = element_text(color = "black"), 
              legend.key = element_blank(), strip.background = element_rect(fill = "#F2F2F2", 
                                                                            colour = "black", size = 0.7), plot.margin = plot.margin, 
              legend.position = legend, complete = TRUE)
      if (x.text.angle != 0) 
        .theme <- .theme + theme(axis.text.x = element_text(angle = x.text.angle, 
                                                            hjust = xhjust))
      .theme
    }
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+custome_theme(base_size=12,base_family="Times",border=T,legend=c(0.7,0.2))
    p1=p1+xlab('False positive fraction')+ylab('True positive fraction')
    # p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    p1=p1+theme(legend.background = element_rect(fill = NA, colour = NA))
    p1=p1+scale_colour_manual(values = pal)
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}


## 无虚线
ggplotKMCox=function(dat,title=NULL,legend.title='Groups',labs=NULL,add_text=NULL,pal=NULL){
  library(ggplot2)
  library(ggsci)
  library(survival)
  library(ggpubr)
  library(survminer)
  ######
  colnames(dat)=c('time','status','groups')
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  if(is.null(pal)){
    pal=pal_lancet()(9)[c(2,4,3,1,5:6,9)]
  }
  surv=survminer::ggsurvplot(sf, data = dat
                             , palette = pal 
                             , pval = TRUE
                             , surv.median.line = 'hv'
                             , conf.int = T
                             # , linetype = "strata"
                             , title=title
                             , legend.title=legend.title
                             , xlab = "Time(years)"
                             , conf.int.style = 'step'
                             , pval.coord = c(0, 0.2)#Add p-value
                             , risk.table = TRUE
                             , ggtheme = theme_pubr(base_size=12,base_family='Times')
                             , font.family='Times'
                             , risk.table.y.text = FALSE
                             , legend.labs = labs)
  p1=surv$plot
  p2=surv$table
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(0.9,0.3),align = "v")
  return(g2)
}

################
gse65682_data=cbind(gse65682.cli,t(gse65682.exp[,rownames(gse65682.cli)]))
gse65682_data=gse65682_data[!is.na(gse65682_data$OS.time),]
dim(gse65682_data) ## 479*19271
table(gse65682_data$disease)
# sepsis 
# 479

gse65682.lasso.res=mg_lasso_cox_use(gse65682_data[,ppi.hub.genes]
                                    , time = gse65682_data$OS.time
                                    , event = gse65682_data$OS
                                    , nfolds = 5
                                    , lambda.min = T
                                    , figLabels = c('B', 'C'))
gse65682.lasso.res$Genes
### "LTF"     "CEACAM8" "PGLYRP1" "MAPK14"  "S100A12" "LCN2"
genes=gse65682.lasso.res$Genes
genes=gsub("-","_",genes)
genes
### "LTF"     "CEACAM8" "PGLYRP1" "MAPK14"  "S100A12" "LCN2"
fmla <- as.formula(paste0("Surv(OS.time, OS) ~",paste(ppi.hub.genes , collapse= "+")))
fmla 
######### 逐步多因素回归分析
# colnames(gse65682_data)=gsub("-","_",colnames(gse65682_data))
cox <- coxph(fmla, data =gse65682_data)
cox=step(cox)
pdf('analysis/6_Model/GSE65682.model.forest.pdf', height = 5, width = 5,onefile = F)
survminer::ggforest(cox,data=gse65682_data,noDigits = 3)
dev.off()
########
coef(cox)
lan=as.numeric(coef(cox))
lan
gene.coef=data.frame(Gene=names(coef(cox)),Coef=coef(cox))
gene.coef$Type=ifelse(gene.coef$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Protective','Risk'))
gene.coef$Coef=round(gene.coef$Coef,digits = 3)
table(gene.coef$Type)
head(gene.coef)
library(dplyr)
lasso.coef.barplot=ggbarplot(data = gene.coef,x='Gene',y='Coef',fill = 'Type',
                                  color = 'white',
                                  orientation = "horiz",   #横向显示
                                  palette = "aaas",    #配色方案
                                  legend = 'top',    #图例位置
                                  label = TRUE,
                                  lab.vjust = 0.5,
                                  lab.hjust=-0.1,
                                  sort.val = "asc",    #上升排序，区别于desc
                                  sort.by.groups = TRUE) +    #按组排序
  coord_flip()+xlab('')+ylab('LASSO cox coefficient')+
  scale_y_continuous(expand=c(0.05, 0.05)) + scale_x_discrete(expand=c(0,0))+theme_custom_3
lasso.coef.barplot

gse65682.lasso.model=mg_merge_plot(gse65682.lasso.res$plot,lasso.coef.barplot,nrow = 1,ncol = 2,widths = c(2,1))
gse65682.lasso.model

############### GSE65682
rs= as.matrix((gse65682_data[,names(coef(cox))])) %*% lan
rs=mosaic::zscore(rs)
score.GSE65682=data.frame(gse65682_data[,c('OS.time','OS')],rs=rs)
score.GSE65682$group=ifelse(score.GSE65682$rs>0,'High','Low')
dim(score.GSE65682) ## 479*4

########## 单因素cox分析
mg_Forestplot=function(df_m,outFile,width=6,height=3){
  colnames(df_m)=c('HR','HR.95L','HR.95H','pvalue')
  gene=rownames(df_m)
  # hr=sprintf("%.3f",df_m$"HR")
  # hrLow=sprintf("%.3f",df_m$"HR.95L")
  # hrHigh=sprintf("%.3f",df_m$"HR.95H")
  hr=format(df_m$"HR", digits = 3)
  hrLow=format(df_m$"HR.95L", digits = 3)
  hrHigh=format(df_m$"HR.95H",digits =3)
  Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
  # pVal=ifelse(df_m$pvalue<1e-8, "<1e-8", sprintf("%.3f", df_m$pvalue))
  pVal=format(df_m$pvalue, digits = 3)
  #########
  pdf(file=outFile, width = width, height =height,onefile = FALSE)
  n=nrow(df_m)
  nRow=n+1
  ylim=c(1,nRow)
  layout.show(layout(matrix(c(1,2),nc=2),width=c(2,1.2)))
  #森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,1,2,1),mpg=c(2,0.5,0))
  # par(mar=c(3,2,1.5,1.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mpg=c(2,0.5,0))
  # par(mar=c(3,1,1.5,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="black",lwd=2.5,lty=5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
  points(as.numeric(hr), n:1, pch = 20, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
unicox<-function(vars=c("T","N"),time=NULL,event=NULL,data=LUAD_clinical){
  library(survival)
  require(survminer)
  y<-Surv(as.numeric(time),
          as.numeric(event))
  
  ## for循环
  pvalue<-c()
  HR<-c()
  lower<-c()
  upper<-c()
  varname<-c()
  
  for(i in vars ){
    cox.fit_uni<-coxph(y~data[[i]])
    uni_res <-summary(cox.fit_uni)
    pvalue[[i]]<-uni_res$waldtest[[3]]#pvalue
    HR[[i]]<-uni_res$conf.int[[1]]# HR提取
    lower[[i]] <-uni_res$conf.int[[3]]#lower
    upper[[i]] <-uni_res$conf.int[[4]]# upper
  }
  ## 保存结果到数据框中
  univar_res<-data.frame(
    # varname<-vars,
    HR<-as.numeric(HR),
    # CI=paste(as.numeric(lower),"-",as.numeric(upper),sep = ""),
    HR.95L<-as.numeric(lower),
    HR.95H<-as.numeric(upper),
    pvalue<-as.numeric(pvalue)
  )
  colnames(univar_res)<-c("HR","HR.95L","HR.95H","pvalue")
  rownames(univar_res)=vars
  univar_res## 最后的函数运行结果
}
multicox<-function(vars=c("T","N","M","age"),time=NULL,event=NULL,data=LUAD_clinical,forest=T){
  library(survival)
  require(survminer)
  y<-Surv(as.numeric(time),
          as.numeric(event))
  ## 构建公式
  FM<-as.formula(paste0("y~",paste(vars,collapse = "+")))
  cox.fit_multi <- coxph(FM,data = data)
  munivar_res<-summary(cox.fit_multi)#cox的结果
  pvalue<-munivar_res$coefficients[,"Pr(>|z|)"]#pvalue
  HR<-munivar_res$coefficients[,"exp(coef)"]# HR
  lower<-munivar_res$conf.int[,3]
  upper<-munivar_res$conf.int[,4]
  ## 保存结果到数据框中
  munivar_res<-data.frame(
    # varname<-vars,
    HR<-as.numeric(HR),
    # CI=paste(as.numeric(lower),"-",as.numeric(upper),sep = ""),
    HR.95L<-as.numeric(lower),
    HR.95H<-as.numeric(upper),
    pvalue<-as.numeric(pvalue)
  )
  colnames(munivar_res)<-c("HR","HR.95L","HR.95H","pvalue")
  rownames(munivar_res)=vars
  
  ## 保存pdf到工作目录
  if (forest==T){
    ggforest(cox.fit_multi,data = data)
    ggsave(filename = "mult_cox.pdf")}
  munivar_res## 最后的函数运行结果
  
}

score.GSE65682_dt=cbind(score.GSE65682,gse65682.cli[rownames(score.GSE65682),])

univar_res<-unicox(vars=c("age","gender","rs"),time = score.GSE65682_dt$OS.time,event = score.GSE65682_dt$OS,data=score.GSE65682_dt)
univar_res
rownames(univar_res)[3]='riskscore'
univar_res[which(univar_res$pvalue<0.05),]
## 多因素
mutivar_res<-multicox(vars=c("age","rs"),time = score.GSE65682_dt$OS.time,event = score.GSE65682_dt$OS,data=score.GSE65682_dt,forest = F)
mutivar_res
rownames(mutivar_res)[2]='riskscore'
table(mutivar_res$pvalue<0.05)

mg_Forestplot(df_m = univar_res,outFile = 'analysis/6_Model/gse65682.univar.forestplot.pdf',height = 3.5,width = 5)
mg_Forestplot(df_m = mutivar_res,outFile = 'analysis/6_Model/gse65682.mutivar.forestplot.pdf',height = 3.5,width = 5)

####################################
dt=score.GSE65682_dt[!is.na(score.GSE65682_dt$OS),]
dt=reshape::rename(dt,c('OS'='status'))
write.table(dt[,c(1,2,3,7)],file = 'analysis/6_Model/nomogram_dat.txt',sep = "\t",row.names = F,quote = F)

library(rms)
pdf('analysis/6_Model/model.calibration.pdf',height = 5,width = 5)
mg_plot_calibrate(dt$OS.time,dt$status,data = data.frame(dt[,c('rs','age')]),cut.time = c(7,14,21),timeLabel = c("1-Weeks","2-Weeks",'3-Weeks'))
dev.off()

pdf('analysis/6_Model/model.DCA.pdf',height = 5,width = 5)
mg_plotDCA(dt$status,c('age','gender','rs','age+rs'),c('age','gender','riskscore','Nomogram'),dt)
dev.off()

#######################
head(gse65682.cli)
colnames(gse65682.cli)[c(18,4,15)]
p.all=list()
for(i in c(18,4,15)){
  print(i)
  dt=data.frame(type=gse65682.cli[rownames(score.GSE65682),i],riskscore=score.GSE65682$rs,stringsAsFactors = F)
  dt=dt[which(!is.na(dt$type)),]
  xlab=colnames(gse65682.cli)[i]
  if(xlab=='outcome'){
    xlab='survival'
  }else if (xlab=='Age1'){
    xlab='age'
  }
  p=mg_violin_use(data = dt
                  ,melt = T
                  ,group.col = mycolors
                  ,xlab = xlab
                  ,ylab = "riskscore"
                  ,jitter=F
                  ,test_method = 'kruskal.test'
                  ,cmp_test_method = 'wilcox.test'
                  ,legend.pos = NULL
                  ,show_compare = T)
  
  p
  p.all=c(p.all,list(p))
}
length(p.all)
gse65682.cli.riskscore.boxplot=mg_merge_plot(p.all,nrow = 1,ncol = 3)

#############################################
############### GSE57065
#############################################
gse57065_data=cbind(gse57065.cli,t(gse57065.exp[,rownames(gse57065.cli)]))
table(gse57065_data$disease)
gse57065_data=gse57065_data[which(gse57065_data$disease=='sepsis'),]

rs= as.matrix((gse57065_data[,names(coef(cox))])) %*% lan
rs=mosaic::zscore(rs)
score.GSE57065=data.frame(rs=rs)
score.GSE57065$group=ifelse(score.GSE57065$rs>0,'High','Low')
head(score.GSE57065)
dim(score.GSE57065)

head(gse57065.cli)
colnames(gse57065.cli)[c(8,5,6)]

p.all=list()
for(i in c(8,5,6)){
  print(i)
  dt=data.frame(type=gse57065.cli[rownames(score.GSE57065),i],riskscore=score.GSE57065$rs,stringsAsFactors = F)
  dt=dt[which(!is.na(dt$type)),]
  xlab=colnames(gse57065.cli)[i]
  if(xlab=='sapsii'){
    xlab='SAPS II'
  }else if(xlab=='Age1'){
    xlab='age'
  }
  p=mg_violin_use(data = dt
                  ,melt = T
                  ,group.col = mycolors
                  ,xlab = xlab
                  ,ylab = "riskscore"
                  ,jitter=F
                  ,test_method = 'kruskal.test'
                  ,cmp_test_method = 'wilcox.test'
                  ,legend.pos = NULL
                  ,show_compare = T)
  
  p
  p.all=c(p.all,list(p))
}
length(p.all)
gse57065.cli.riskscore.boxplot=mg_merge_plot(p.all,nrow = 1,ncol = 3)
gse57065.cli.riskscore.boxplot

#############################################
############### GSE95233
#############################################
gse95233_data=cbind(gse95233.cli,t(gse95233.exp[,rownames(gse95233.cli)]))
table(gse95233_data$disease)
# healthy  sepsis 
# 22     102
gse95233_data=gse95233_data[which(gse95233_data$disease=='sepsis'),]

rs= as.matrix((gse95233_data[,names(coef(cox))])) %*% lan
rs=mosaic::zscore(rs)
score.GSE95233=data.frame(gse95233_data[,c('OS.time','OS')],rs=rs)
score.GSE95233$group=ifelse(score.GSE95233$rs>0,'High','Low')
dim(score.GSE95233) ## 102*4

head(gse95233.cli)
colnames(gse95233.cli)[c(9,4,6)]
p.all=list()
for(i in c(9,4,6)){
  print(i)
  dt=data.frame(type=gse95233.cli[rownames(score.GSE95233),i],riskscore=score.GSE95233$rs,stringsAsFactors = F)
  dt=dt[which(!is.na(dt$type)),]
  xlab=colnames(gse95233.cli)[i]
  if(xlab=='Age1'){
    xlab='age'
  }
  p=mg_violin_use(data = dt
                ,melt = T
                ,group.col = mycolors
                ,xlab = xlab
                ,ylab = "riskscore"
                ,jitter=F
                ,test_method = 'kruskal.test'
                ,cmp_test_method = 'wilcox.test'
                ,legend.pos = NULL
                ,show_compare = T)
  
  p
  p.all=c(p.all,list(p))
}
length(p.all)
gse95233.cli.riskscore.boxplot=mg_merge_plot(p.all,nrow = 1,ncol = 3)

#############################################
############### GSE185263
#############################################
gse185263_data=cbind(gse185263.cli,t(gse185263.exp[,rownames(gse185263.cli)]))
table(gse185263_data$disease)
gse185263_data=gse185263_data[which(gse185263_data$disease=='sepsis'),]

rs= as.matrix((gse185263_data[,names(coef(cox))])) %*% lan
rs=mosaic::zscore(rs)
score.GSE185263=data.frame(gse185263_data[,c('OS.time','OS')],rs=rs)
score.GSE185263$group=ifelse(score.GSE185263$rs>0,'High','Low')
dim(score.GSE185263)

head(gse185263.cli)
colnames(gse185263.cli)[c(13,8,7)]

p.all=list()
for(i in c(13,8,7)){
  print(i)
  dt=data.frame(type=gse185263.cli[rownames(score.GSE185263),i],riskscore=score.GSE185263$rs,stringsAsFactors = F)
  dt=dt[which(!is.na(dt$type)),]
  xlab=colnames(gse185263.cli)[i]
  if(xlab=='sapsii'){
    xlab='SAPS II'
  }else if(xlab=='Age1'){
    xlab='age'
  }
  p=mg_violin_use(data = dt
                  ,melt = T
                  ,group.col = mycolors
                  ,xlab = xlab
                  ,ylab = "riskscore"
                  ,jitter=F
                  ,test_method = 'kruskal.test'
                  ,cmp_test_method = 'wilcox.test'
                  ,legend.pos = NULL
                  ,show_compare = T)
  
  p
  p.all=c(p.all,list(p))
}
length(p.all)
gse185263.cli.riskscore.boxplot=mg_merge_plot(p.all,nrow = 1,ncol = 3)

figure=mg_merge_plot(gse65682.cli.riskscore.boxplot,gse57065.cli.riskscore.boxplot,
                     gse95233.cli.riskscore.boxplot,gse185263.cli.riskscore.boxplot,
                     nrow = 4,ncol = 1,labels = LETTERS[1:4],font.label = list(size = 14, color = "black", face ="bold", family = 'Times'))
figure

savePDF('analysis/6_Model/bulk.cli.riskscore.pdf',figure,height = 14,width = 12)

cor_point(x=as.numeric(gse185263.cli[rownames(score.GSE185263),"Sofa score"]),
          y=score.GSE185263$rs)

###########################
ggplotTimeROC(time = score.GSE65682$OS.time,status = score.GSE65682$OS,score = as.numeric(as.factor(gse65682.cli[rownames(score.GSE65682),'gender'])),mks = c(7,14,21))
ggplotTimeROC(time = score.GSE65682$OS.time,status = score.GSE65682$OS,score = gse65682.cli[rownames(score.GSE65682),'age'],mks = c(7,14,21))

gse65682.TimeROC=ggplotTimeROC(time = score.GSE65682$OS.time,status = score.GSE65682$OS,score = score.GSE65682$rs,mks = c(7,14,21))
gse65682.TimeROC

RiskType.colors=alpha(mycolors[c(3,5)],alpha = 1)
gse65682.KM=ggplotKMCox(data.frame(time = score.GSE65682$OS.time
                                   , event = score.GSE65682$OS
                                   , groups=score.GSE65682$group)
                        , legend.title = 'RiskType'
                        , title = 'GSE65682'
                        , pal = RiskType.colors
                        , labs = c('High', 'Low')
                        , add_text = '')

gse95982.KM=ggplotKMCox(data.frame(time = score.GSE95233$OS.time
                                   , event = score.GSE95233$OS
                                   , groups=score.GSE95233$group)
                        , legend.title = 'RiskType'
                        , title = 'GSE95233'
                        , pal = RiskType.colors
                        , labs = c('High', 'Low')
                        , add_text = '')
gse95982.KM

figure=mg_merge_plot(mg_getplot_bank('A'),gse65682.TimeROC,gse65682.KM,gse95982.KM,
                     nrow = 2,ncol = 2,labels = LETTERS[1:4],font.label = list(size = 14, color = "black", face ="bold",family = 'Times'))

savePDF('analysis/6_Model/bulk.model.AUC.pdf',figure,height = 10,width = 10)

########################## 与 PCD 和 TME 的关系
########## riskscore与免疫细胞的相关性分析
############ 相关性分析
cor.res=psych::corr.test(x=data.frame(RiskScore=score.GSE65682$rs),
                         y = data.frame(gse65682.pcd.ssgsea[rownames(score.GSE65682),],check.names = F))
df_cor=cor.res$r
df_pval=cor.res$p.adj
########
library(tidyverse)
g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p
head(g)

################## 
df=data.frame(g,check.names = F,stringsAsFactors = F)
head(df)
write.table(df,file = 'analysis/6_Model/gse65682.model.PCD.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)
write.table(df,file = 'results/Files/gse65682.model.PCD.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)

model.pcd.cor.plot=df %>% 
  filter(p.adj<0.05) %>%
  ggplot(aes(x=cor, y=reorder(to, cor), color=cor))+
  # geom_point(aes(size=abs(cor)))+
  geom_point(aes(size=-log10(p.adj)))+
  geom_segment(aes(x=0, xend=cor, y=to, yend=to), color="black")+
  scale_colour_gradient2(high="#F97B72", mid = "white",low="#0072B5",midpoint = 0)+
  # coord_flip()+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=14))
model.pcd.cor.plot
############### Cibersort
load("raw_datas/Pathways_ssGSEA/gse65682.exp.cibersort.RData")
gse65682.exp.cibersort=get.IOBR.immu.format(gse65682.exp.cibersort)
gse65682.exp.cibersort=gse65682.exp.cibersort[,1:22]
dim(gse65682.exp.cibersort)

cor.res=psych::corr.test(x=data.frame(RiskScore=score.GSE65682$rs),
                         y = data.frame(gse65682.exp.cibersort[rownames(score.GSE65682),1:22],check.names = F))
df_cor=cor.res$r
df_pval=cor.res$p.adj
########
library(tidyverse)
g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p
head(g)
################## 
df=data.frame(g,check.names = F,stringsAsFactors = F)
head(df)
write.table(df,file = 'analysis/6_Model/gse65682.model.CIBERSORT.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)
write.table(df,file = 'results/Files/gse65682.model.CIBERSORT.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)
model.cibersort.cor.plot=df %>% 
  filter(p.adj<0.05) %>%
  ggplot(aes(x=cor, y=reorder(to, cor), color=cor))+
  geom_point(aes(size=-log10(p.adj)))+
  geom_segment(aes(x=0, xend=cor, y=to, yend=to), color="black")+
  scale_colour_gradient2(high="#F97B72", mid = "white",low="#0072B5",midpoint = 0)+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=14))
model.cibersort.cor.plot

############### Hallmark
load('raw_datas/Pathways_ssGSEA/gse65682.h.all.ssgsea.RData')
gse65682.h.all.ssgsea=t(gse65682.h.all.ssgsea)
gse65682.h.all.ssgsea=data.frame(gse65682.h.all.ssgsea,check.names = F,stringsAsFactors = F)
colnames(gse65682.h.all.ssgsea)=gsub("HALLMARK_","",colnames(gse65682.h.all.ssgsea))
dim(gse65682.h.all.ssgsea)

cor.res=psych::corr.test(x=data.frame(RiskScore=score.GSE65682$rs),
                         y = data.frame(gse65682.h.all.ssgsea[rownames(score.GSE65682),],check.names = F))
df_cor=cor.res$r
df_pval=cor.res$p.adj
########
library(tidyverse)
g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p
head(g)

################## 
df=data.frame(g,check.names = F,stringsAsFactors = F)
head(df)
write.table(df,file = 'analysis/6_Model/gse65682.model.HALLMARK.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)
write.table(df,file = 'results/Files/gse65682.model.HALLMARK.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)

model.h.all.cor.plot=df %>% filter(p.adj<0.05) %>% 
  ggplot(aes(x=cor, y=reorder(to, cor), color=cor))+
  geom_point(aes(size=-log10(p.adj)))+
  geom_segment(aes(x=0, xend=cor, y=to, yend=to), color="black")+
  scale_colour_gradient2(high="#F97B72", mid = "white",low="#0072B5",midpoint = 0)+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=10), axis.text.y=element_text(face="bold", size=10), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=14))
model.h.all.cor.plot

font.label =list(size = 14, color = "black", face ="bold",family = 'Times')

plot1=mg_merge_plot(model.pcd.cor.plot,model.cibersort.cor.plot,labels = LETTERS[1:2],font.label = font.label)

model.TME.corplot=mg_merge_plot(plot1,model.h.all.cor.plot,
                               nrow = 2,ncol = 1,heights = c(1,2),labels = c('','C'), font.label = font.label)
model.TME.corplot
savePDF('analysis/6_Model/model.TME.corplot.pdf',model.TME.corplot,height = 12,width = 12)

save.image(file = '20230907_Sepsis_NETosis_final.RData')
