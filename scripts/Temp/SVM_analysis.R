library(caret)
ppi.hub.genes
########## Construct xgb.DMatrix object
input_x=t(gse57065.exp[ppi.hub.genes,])
# input_y=ifelse(gse57065.cli$disease=='healthy',0,1)
input_y=gse57065.cli$disease
input_y=factor(input_y,levels = c('healthy','sepsis'))

train=data.frame(label=as.character(input_y),input_x,stringsAsFactors = F)
# train$label[which(train$label=='healthy')]=0
# train$label[which(train$label=='sepsis')]=1
# train$label=factor(train$label,levels = c(0,1))
# train$label=as.numeric(train$label)

head(train)

train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
svm1 <- train(label ~., data = train, method = "svmLinear", trControl = train_control,  preProcess = c("center","scale"))
# Fit the model 
svm2 <- train(label ~., data = train, method = "svmLinear", trControl = train_control,  preProcess = c("center","scale"), tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
svm2
# Plot model accuracy vs different values of Cost
plot(svm2)
# Print the best tuning parameter C that
# maximizes model accuracy
svm2$bestTune
res2<-as_tibble(svm2$results[which.min(svm2$results[,2]),])
res2


setdiff(ppi.hub.genes,rownames(gse65682.exp))
dt=data.frame(type=gse65682.cli$disease,t(gse65682.exp[ppi.hub.genes,]),stringsAsFactors = F)
table(dt$type)
# dt$type[which(dt$type=='healthy')]=0
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
###
pred <- predict(svm2,as.matrix(dt[,-1]))
# pred <- as.numeric(as.character(pred))
auc(dt$type,pred) ## 0.7519737

########## GSE95233
setdiff(ppi.hub.genes,rownames(gse95233.exp))
dt=data.frame(type=gse95233.cli$Source,t(gse95233.exp[ppi.hub.genes,]),stringsAsFactors = F)
table(dt$type)
dt$type[which(dt$type=='Control')]='healthy'
dt$type[which(dt$type=='Patient')]='sepsis'
# dt$type[which(dt$type=='Control')]=0
# dt$type[which(dt$type=='Patient')]=1
# dt$type=as.numeric(dt$type)
#########
pred <- predict(svm2,as.matrix(dt[,-1]))
# pred= as.numeric(as.character(pred))
auc(dt$type,pred) ## 0.9891

########## GSE185263
setdiff(ppi.hub.genes,rownames(gse185263.exp))
dt=data.frame(type=gse185263.cli$`disease state`,t(gse185263.exp[ppi.hub.genes,]),stringsAsFactors = F)
table(dt$type)
# dt$type[which(dt$type=='healthy')]=0
# dt$type[which(dt$type=='sepsis')]=1
# dt$type=as.numeric(dt$type)
#########
pred <- predict(svm2,as.matrix(dt[,-1]))
# pred= as.numeric(as.character(pred))
auc(dt$type,pred) ## 0.7816
