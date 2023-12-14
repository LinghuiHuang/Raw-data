###############################
library(xgboost)
ppi.hub.genes
model.genes=c('CEACAM8','MAPK14','S100A9')

########## Construct xgb.DMatrix object
input_x=t(gse57065.exp[model.genes,])
# input_y=ifelse(gse57065.cli$disease=='healthy',0,1)
input_y=gse57065.cli$disease
input_y=factor(input_y,levels = c('healthy','sepsis'))

####### 基于caret
##### 模型调参
library(caret)
library(pROC)
##########################################
# # 模型评估标准
metric='Accuracy'
####### Step 1: Number of Iterations and the Learning Rate
nrounds <- 1000
# note to start nrounds from 200, as smaller learning rates result in errors so
# big with lower starting points that they'll mess the scales
tune_grid <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 10),
  eta = c(0.05,0.1,0.2,0.3),
  max_depth = seq(2,10,2),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 10, # with n folds
  savePredictions=TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_tune <- caret::train(
  x = input_x,
  y = input_y,
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  metric=metric,
  verbose = TRUE)

xgb_tune$bestTune
plot(xgb_tune)

####### Step 2: Maximum Depth and Minimum Child Weight
tune_grid2 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 10),
  eta = xgb_tune$bestTune$eta,
  max_depth = if(xgb_tune$bestTune$max_depth == 2)c(xgb_tune$bestTune$max_depth:4)else c(xgb_tune$bestTune$max_depth - 1:xgb_tune$bestTune$max_depth + 1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1, 2, 3),
  subsample = 1
)

xgb_tune2 <- caret::train(x = input_x,
                          y = input_y,
                          trControl = tune_control,
                          tuneGrid = tune_grid2,
                          method = "xgbTree",
                          metric=metric,
                          verbose = TRUE)
xgb_tune2$bestTune
plot(xgb_tune2)
######### Step 3: Column and Row Sampling
tune_grid3 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 10),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = seq(from = 0.6, to = 0.9, by = 0.1),
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample=seq(from = 0.6, to = 0.9, by = 0.1)
)

xgb_tune3 <- caret::train(x = input_x,
                          y = input_y,
                          trControl = tune_control,
                          tuneGrid = tune_grid3,
                          method = "xgbTree",
                          metric=metric,
                          verbose = TRUE)
xgb_tune3$bestTune
plot(xgb_tune3)

######### Step 4: Gamma
tune_grid4 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 10),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune4 <- caret::train(x = input_x,
                          y = input_y,
                          trControl = tune_control,
                          tuneGrid = tune_grid4,
                          method = "xgbTree",
                          metric=metric,
                          verbose = TRUE)

xgb_tune4$bestTune
plot(xgb_tune4)
######## Step 5: Reducing the Learning Rate
tune_grid5 <- expand.grid(
  nrounds = seq(from = 100, to = 10000, by = 100),
  eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = xgb_tune4$bestTune$gamma,
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune5 <- caret::train(x = input_x,
                          y = input_y,
                          trControl = tune_control,
                          tuneGrid = tune_grid5,
                          method = "xgbTree",
                          metric=metric,
                          verbose = TRUE)

plot(xgb_tune5)
xgb_tune5$bestTune

########## Fitting the Model
final_grid <- expand.grid(
  nrounds = xgb_tune5$bestTune$nrounds,
  eta = xgb_tune5$bestTune$eta,
  max_depth = xgb_tune5$bestTune$max_depth,
  gamma = xgb_tune5$bestTune$gamma,
  colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
  min_child_weight = xgb_tune5$bestTune$min_child_weight,
  subsample = xgb_tune5$bestTune$subsample)

train_control <- caret::trainControl(
  method = "none",
  savePredictions=TRUE,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary,
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_model <- caret::train(x = input_x,
                          y = input_y,
                          trControl = train_control,
                          tuneGrid = final_grid,
                          method = "xgbTree",
                          metric=metric,
                          verbose = TRUE)

xgb_model$bestTune

##########################################
ppi.hub.genes
########## Construct xgb.DMatrix object
input_x=t(gse57065.exp[model.genes,])
input_y=gse57065.cli$disease
input_y=factor(input_y,levels = c('healthy','sepsis'))

train=data.frame(label=as.character(input_y),input_x,stringsAsFactors = F)
# train$label[which(train$label=='healthy')]=0
# train$label[which(train$label=='sepsis')]=1
# train$label=factor(train$label,levels = c(0,1))
# train$label=as.numeric(train$label)
#### 设置拟合条件
fitControl <- trainControl(method = "cv",
                           number = 10,
                           verboseIter = FALSE,
                           savePredictions=TRUE,
                           classProbs = TRUE,                                                           # set to TRUE for AUC to be computed
                           summaryFunction = twoClassSummary,
                           seeds=123456,
                           search = "random")
# XGBoost 调参
caret_xgb <- train(label ~ .-1,
                   data = train,
                   method="xgbTree",# xgbLinear
                   metric=metric,
                   trainControl=fitControl)
####### 优化后的超参数
caret_xgb$bestTune
####### 模型预测
caret_best <- caret_xgb$finalModel
importance_matrix <- xgb.importance(model = caret_best)
print(importance_matrix)
#以ggplot风格绘图
xgb.ggplot.importance(importance_matrix = importance_matrix)
###### 查看单棵决策树
xgb.plot.tree(model = caret_best, trees = 1, plot_width = 500,plot_height = 500)


########################################################
tune_grid <- expand.grid(
  nrounds = 100,
  eta = c(0.001,0.005,0.01, 0.1, 0.3),
  max_depth = c(2, 3, 5, 10),
  gamma = c(0),
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

tune_control <- caret::trainControl(method = "cv",
                                    number = 10,
                                    verboseIter = FALSE,
                                    savePredictions=TRUE,
                                    classProbs = TRUE,                                                           # set to TRUE for AUC to be computed
                                    summaryFunction = twoClassSummary,
                                    search = "random")

xgb_tune <- train(label ~ .-1,
                  data = train,
                  trControl = tune_control,
                  tuneGrid = tune_grid,
                  method = "xgbTree",
                  metric='Accuracy',
                  verbose = TRUE,
                  nthreads = 4)

xgb_tune$bestTune
plot(xgb_tune)

#########################################################
# pred1 <- predict(caret_best,newdata = as.matrix(train[,-1]),type = "prob")
pred <- predict(caret_xgb,newdata = as.matrix(train[,-1]))
pred
auc(train$label,pred)

########## GSE57065
setdiff(ppi.hub.genes,rownames(gse57065.exp))
dt=data.frame(type=gse57065.cli$disease,t(gse57065.exp[ppi.hub.genes,]),stringsAsFactors = F)
# dt$type[which(dt$type=='healthy')]=0
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
##
pred <- predict(caret_xgb,as.matrix(dt[,-1]))
# pred <- as.numeric(as.character(pred))
auc(dt$type,pred)
##
truth=dt$type
xtab <- table(pred, truth)
caret::confusionMatrix(xtab,positive='1',mode = 'everything')

roc.res <- pROC::roc(response=truth,predictor=(pred), levels=c(0,1),direction='<',smooth = F, ci = T)
gse57065.ggcoc <- ggroc(roc.res) + theme_gray() + 
  annotate("text", x = 0.4, y = 0.4, 
           label = paste0('GSE57065 \n AUC = ', round(roc.res$auc, 3)))
gse57065.ggcoc

########## GSE65682
setdiff(ppi.hub.genes,rownames(gse65682.exp))
rownames(gse65682.exp)[grep("HIST2H2AA3",rownames(gse65682.exp))]
rownames(gse65682.exp)[which(rownames(gse65682.exp)=='HIST2H2AA3')]='H2AC18'
dt=data.frame(type=gse65682.cli$disease,t(gse65682.exp[ppi.hub.genes,]),stringsAsFactors = F)
table(dt$type)
# dt$type[which(dt$type=='healthy')]=0
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
###
pred <- predict(caret_xgb,as.matrix(dt[,-1]))
# pred <- as.numeric(as.character(pred))
auc(dt$type,pred) ## 0.7519737
##
truth=dt$type
xtab <- table(pred, truth)
caret::confusionMatrix(xtab,positive='1',mode = 'everything')

roc.res <- pROC::roc(response=truth,predictor=(pred), levels=c(0,1),direction='<',smooth = F, ci = T)
gse65682.ggcoc <- ggroc(roc.res) + theme_gray() + 
  annotate("text", x = 0.4, y = 0.4, 
           label = paste0('GSE65682 \n AUC = ', round(roc.res$auc, 3)))
gse65682.ggcoc

# ########## GSE54514
# setdiff(ppi.hub.genes,rownames(gse54514.exp))
# dt=data.frame(type=gse54514.cli$disease,t(gse54514.exp[ppi.hub.genes,]),stringsAsFactors = F)
# dt$type[dt$type %in% c('sepsis nonsurvivor','sepsis survivor')]='sepsis'
# table(dt$type)
# dt$type[which(dt$type=='healthy')]=0
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
# pred <- predict(xgb_model,as.matrix(dt[,-1]))
# auc(dt$type,pred)
# ########## GSE145227
# setdiff(ppi.hub.genes,rownames(gse145227.exp))
# dt=data.frame(type=gse145227.cli$Source,t(gse145227.exp[ppi.hub.genes,]),stringsAsFactors = F)
# table(dt$type)
# dt$type[which(dt$type=='control')]='healthy'
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
# pred <- predict(caret_best,as.matrix(dt[,-1]))
# auc(dt$type,pred) ## 0.6333

########## GSE95233
setdiff(ppi.hub.genes,rownames(gse95233.exp))
rownames(gse95233.exp)[which(rownames(gse95233.exp)=='HIST2H2AA3')]='H2AC18'
dt=data.frame(type=gse95233.cli$disease,t(gse95233.exp[model.genes,]),stringsAsFactors = F)
table(dt$type)
# dt$type[which(dt$type=='Control')]=0
# dt$type[which(dt$type=='Patient')]=1
# dt$type=as.numeric(dt$type)
#########
pred <- predict(caret_xgb,as.matrix(dt[,-1]))
# pred= as.numeric(as.character(pred))
auc(dt$type,pred) ## 0.9891

truth=dt$type
xtab <- table(pred, truth)
caret::confusionMatrix(xtab,positive='1',mode = 'everything')

roc.res <- pROC::roc(response=truth,predictor=(pred), levels=c(0,1),direction='<',smooth = F, ci = T)
gse95233.ggcoc <- ggroc(roc.res) + theme_gray() + 
  annotate("text", x = 0.4, y = 0.4, 
           label = paste0('GSE95233 \n AUC = ', round(roc.res$auc, 3)))
gse95233.ggcoc

########## GSE185263
setdiff(ppi.hub.genes,rownames(gse185263.exp))
dt=data.frame(type=gse185263.cli$disease,t(gse185263.exp[model.genes,]),stringsAsFactors = F)
table(dt$type)
# dt$type[which(dt$type=='healthy')]=0
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
#########
pred <- predict(caret_xgb,as.matrix(dt[,-1]))
# pred= as.numeric(as.character(pred))
auc(dt$type,pred) ## 0.9891

truth=dt$type
xtab <- table(pred, truth)
caret::confusionMatrix(xtab,positive='sepsis',mode = 'everything')

roc.res <- pROC::roc(response=truth,predictor=(pred), levels=c(0,1),direction='<',smooth = F, ci = T)
gse185263.ggcoc <- ggroc(roc.res) + theme_gray() + 
  annotate("text", x = 0.4, y = 0.4, 
           label = paste0('GSE185263 \n AUC = ', round(roc.res$auc, 3)))
gse185263.ggcoc

plot2=mg_merge_plot(gse57065.ggcoc,gse65682.ggcoc,gse95233.ggcoc,gse185263.ggcoc,nrow = 1,ncol = 4)
plot2
