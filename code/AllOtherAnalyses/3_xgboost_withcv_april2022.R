####################################
#####################################
## The xgboost model

source("../../MutationTimer/shap-values-master/shap.R")
library(tidyverse)
library(xgboost)
library(caret)
library(reshape)
library(SHAPforxgboost)

load("data/sigs.timed.annot.RData")
sigs.timed.annot[which(sigs.timed.annot$Category == "LymphNode"),]$Category <- "Metastasis"
df.sigs.timed <- melt(sigs.timed.annot)
sigs.reshaped <- sigs.timed.annot
sigs.reshaped$Timing <- sapply(sigs.reshaped$Timing,
                               function(x) ifelse(is.na(x),1,ifelse(x=="Early",0,1)))
sigs.reshaped$Clonality <- sapply(sigs.reshaped$Clonality,
                                  function(x) ifelse(x=="Clonal",0,1))
sigs.reshaped.barrprim <- sigs.reshaped[which(sigs.reshaped$Category %in% 
                                                c("Barretts","PrimaryTumour")),]
sigs.reshaped.barr <- sigs.reshaped[which(sigs.reshaped$Category %in% 
                                            c("Barretts")),]
sigs.reshaped.prim <- sigs.reshaped[which(sigs.reshaped$Category %in% 
                                            c("PrimaryTumour")),]

set.seed(19875) 
randTrainingSet.barr <- sample(1:nrow(sigs.reshaped.barr), 
                               round(nrow(sigs.reshaped.barr)*70/100),
                               replace=F)
dat.train.barr <- sigs.reshaped.barr[randTrainingSet.barr,]
dat.test.barr <- sigs.reshaped.barr[setdiff(1:nrow(sigs.reshaped.barr),randTrainingSet.barr),]
randTrainingSet.prim <- sample(1:nrow(sigs.reshaped.prim), 
                               round(nrow(sigs.reshaped.prim)*70/100),
                               replace=F)
dat.train.prim <- sigs.reshaped.prim[randTrainingSet.prim,]
dat.test.prim <- sigs.reshaped.prim[setdiff(1:nrow(sigs.reshaped.prim),randTrainingSet.prim),]

dat.train <- rbind(dat.train.barr, dat.train.prim)
dat.test <- rbind(dat.test.barr, dat.test.prim)

dat.train$CategoryNumerical <- sapply(dat.train$Category,
                                      function(x) ifelse(x=="Barretts",0,1))
dat.test$CategoryNumerical <- sapply(dat.test$Category,
                                     function(x) ifelse(x=="Barretts",0,1))

df.input <- as.matrix(dat.train[,-c(1,ncol(dat.train)-2,ncol(dat.train)-1,ncol(dat.train))])
rownames(df.input) <- dat.train$Sample
df.predict <- as.matrix(dat.test[,-c(1,ncol(dat.test)-2,ncol(dat.test)-1,ncol(dat.test))])
rownames(df.predict) <- dat.test$Sample

dtrain <- xgb.DMatrix(df.input, label = dat.train$CategoryNumerical)
dtest <- xgb.DMatrix(df.predict, label = dat.test$CategoryNumerical)

model_stages = xgboost(data = df.input, 
                       nround = 100, 
                       objective="binary:logistic",
                       label= dat.train$CategoryNumerical)  

df.predict <- as.matrix(dat.test[,-c(1,ncol(dat.test)-2,ncol(dat.test)-1,ncol(dat.test))])
rownames(df.predict) <- dat.test$Sample
pred <- predict(model_stages, df.predict)
pred.binary <- as.numeric(pred > 0.5)

plot(pred,dat.test$CategoryNumerical)
err <- mean(as.numeric(pred > 0.5) != dat.test$CategoryNumerical)
print(paste("test-error=", err))
#"test-error= 0.147169811320755"

table(pred.binary == dat.test$CategoryNumerical)
243/(243+38)
#86.5% accuracy
655/(110+655)
#86%

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
var_importance(shap_result_stages, top_n=15)

## Prepare data for top N variables
shap_long_stages = shap.prep(shap = shap_result_stages,
                             X_train = df.predict, 
                             top_n = 15
)

## Plot shap overall metrics
pdf("plots.xgboost.april2022/xgboost.shapvalues.BarrettsVsPrimary.clontimingindep.pdf")
plot.shap.summary(data_long = shap_long_stages)
dev.off()

write.csv(shap_long_stages, file="plots_3/shap_long_stages.csv")

#########################

library("randomForest")
library("DALEX")

set.seed(1313)
dat.train$CategoryNumerical <- as.factor(dat.train$CategoryNumerical)
sigs_rf <- randomForest(CategoryNumerical ~ ., 
                        data = dat.train[,c(2:17,20)])
explain_rf <- DALEX::explain(model = sigs_rf,  
                             data = df.input,
                             y = dat.train$CategoryNumerical == 1, 
                             label = "Random Forest")
pred <- predict(explain_rf, dat.test[,c(2:17,20)])
pred.binary <- as.numeric(pred > 0.5)
table(pred.binary == dat.test$CategoryNumerical)
649/(649+116)
#0.848366 accuracy

library(pROC)
# predict(.., type = 'prob') returns a probability matrix
rf.roc <- roc(dat.test$CategoryNumerical, 
              pred)
plot(rf.roc)
auc(rf.roc)
#Area under the curve: 0.8458

xgb.roc <- roc(dat.test$CategoryNumerical, 
               predict(model_stages, df.predict))
plot(xgb.roc)
auc(xgb.roc)
# Area under the curve: 0.8138

plot(rf.roc,xgb.roc)

write.csv(cbind(rf.roc$sensitivities,rf.roc$specificities,rf.roc$thresholds,
                xgb.roc$sensitivities,xgb.roc$specificities,xgb.roc$thresholds), 
          file="plots_3/ROC.AUC.csv")

pdf("plots.xgboost.april2022/ROCcurves.pdf")
roc_rose <- plot(rf.roc, print.auc = TRUE, col = "darkblue")
roc_rose <- plot(xgb.roc, print.auc = TRUE, 
                 col = "brown", print.auc.y = .4, add = TRUE)
dev.off()

