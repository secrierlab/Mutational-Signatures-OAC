####################################
#####################################
## The xgboost model

source("../../MutationTimer/shap-values-master/shap.R")
library(tidyverse)
library(xgboost)
library(caret)
library(reshape)

load("data/sigs.timed.annot.RData")
sigs.timed.annot[which(sigs.timed.annot$Category == "LymphNode"),]$Category <- "Metastasis"
# remove NAs:
#sigs.timed.annot <- sigs.timed.annot[which((!is.na(sigs.timed.annot$Timing))&
#                                       (!is.na(sigs.timed.annot$Clonality))),]
df.sigs.timed <- melt(sigs.timed.annot)
df.sigs.timed$Feature <- apply(df.sigs.timed[,c("variable","Clonality","Timing")],1,
                               function(x) paste(x,collapse="."))

sigs.reshaped <- data.frame(cast(Sample~Feature, fun.aggregate = "sum",
                      data=df.sigs.timed, value="value"))
sigs.reshaped <- merge(sigs.reshaped,
                       unique(df.sigs.timed[,c("Sample","Category")]),
                       by.x="Sample", by.y="Sample",
                       all.x=FALSE,all.y=FALSE)


set.seed(19875) 


## Create the xgboost model for primary vs mets

sigs.reshaped.mets <- sigs.reshaped[which(sigs.reshaped$Category %in% 
                                            c("Metastasis")),]
sigs.reshaped.prim <- sigs.reshaped[which(sigs.reshaped$Category %in% 
                                            c("PrimaryTumour")),]

randTrainingSet.mets <- sample(1:nrow(sigs.reshaped.mets), 
                               round(nrow(sigs.reshaped.mets)*70/100),
                               replace=F)
dat.train.mets <- sigs.reshaped.mets[randTrainingSet.mets,]
dat.test.mets <- sigs.reshaped.mets[setdiff(1:nrow(sigs.reshaped.mets),randTrainingSet.mets),]
randTrainingSet.prim <- sample(1:nrow(sigs.reshaped.prim), 
                               round(nrow(sigs.reshaped.prim)*70/100),
                               replace=F)
dat.train.prim <- sigs.reshaped.prim[randTrainingSet.prim,]
dat.test.prim <- sigs.reshaped.prim[setdiff(1:nrow(sigs.reshaped.prim),randTrainingSet.prim),]

dat.train <- rbind(dat.train.mets, dat.train.prim)
dat.test <- rbind(dat.test.mets, dat.test.prim)

dat.train$CategoryNumerical <- sapply(dat.train$Category,
                                      function(x) ifelse(x=="Metastasis",1,0))
dat.test$CategoryNumerical <- sapply(dat.test$Category,
                                     function(x) ifelse(x=="Metastasis",1,0))

df.input <- as.matrix(dat.train[,-c(1,ncol(dat.train)-1,ncol(dat.train))])
rownames(df.input) <- dat.train$Sample

model_stages = xgboost(data = df.input, 
                       nround = 20, 
                       objective="binary:logistic",
                       label= dat.train$CategoryNumerical)  

df.predict <- as.matrix(dat.test[,-c(1,ncol(dat.test)-1,ncol(dat.test))])
rownames(df.predict) <- dat.test$Sample
pred <- predict(model_stages, df.predict)
pred.binary <- as.numeric(pred > 0.5)

plot(pred,dat.test$CategoryNumerical)
err <- mean(as.numeric(pred > 0.5) != dat.test$CategoryNumerical)
print(paste("test-error=", err))
#"test-error= 0.0756972111553785"

table(pred.binary == dat.test$CategoryNumerical)
232/(232+19)
#92.4%3

importance_matrix <- xgb.importance(model = model_stages)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)

## Calculate shap values
shap_result_stages = shap.score.rank(xgb_model = model_stages, 
                                     X_train =df.predict,
                                     shap_approx = F
)

# `shap_approx` comes from `approxcontrib` from xgboost documentation. 
# Faster but less accurate if true. Read more: help(xgboost)

## Plot var importance based on SHAP
var_importance(shap_result_stages, top_n=28)

## Prepare data for top N variables
shap_long_stages = shap.prep(shap = shap_result_stages,
                             X_train = df.predict, 
                             top_n = 28
)

## Plot shap overall metrics
pdf("plots_supp/xgboost.shapvalues.MetsVsPrimary.moreFeatures.pdf")
plot.shap.summary(data_long = shap_long_stages)
dev.off()
write.csv(shap_long_stages, file="plots_supp/metModel.shap.csv")
