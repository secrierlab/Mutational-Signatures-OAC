####################################
#####################################
## The xgboost model

source("../../MutationTimer/shap-values-master/shap.R")
library(tidyverse)
library(xgboost)
library(caret)
library(reshape)
library(SHAPforxgboost)

load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples[sigs.allSamples<0.05] <- 0
sigs.allSamples$Sample <- rownames(sigs.allSamples) 

## Add indels:
load("data/indelSignatures.fullCohort.RData")
sigs.indel[sigs.indel<0.05] <- 0
sigs.indel$Sample <- rownames(sigs.indel)

sigs <- merge(sigs.allSamples,sigs.indel,
              by.x = "Sample", by.y="Sample",
              all.x=FALSE, all.y = FALSE)

load("data/annotation.sampleIDs.RData")
annotation.sampleIDs[which(annotation.sampleIDs$Category == "LymphNode"),]$Category <- "Metastasis"
sigs.annot <- merge(sigs, annotation.sampleIDs[,c("Sample","Category")],
                    by.x="Sample", by.y="Sample",
                    all.x=FALSE, all.y=FALSE)


sigs.reshaped <- sigs.annot
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

df.input <- as.matrix(dat.train[,-c(1,ncol(dat.train)-1,ncol(dat.train))])
rownames(df.input) <- dat.train$Sample
df.predict <- as.matrix(dat.test[,-c(1,ncol(dat.test)-1,ncol(dat.test))])
rownames(df.predict) <- dat.test$Sample

dtrain <- xgb.DMatrix(df.input, label = dat.train$CategoryNumerical)
dtest <- xgb.DMatrix(df.predict, label = dat.test$CategoryNumerical)

model_stages = xgboost(data = df.input, 
                       nround = 100, 
                       objective="binary:logistic",
                       label= dat.train$CategoryNumerical)  

df.predict <- as.matrix(dat.test[,-c(1,ncol(dat.test)-1,ncol(dat.test))])
rownames(df.predict) <- dat.test$Sample
pred <- predict(model_stages, df.predict)
pred.binary <- as.numeric(pred > 0.5)

plot(pred,dat.test$CategoryNumerical)
err <- mean(as.numeric(pred > 0.5) != dat.test$CategoryNumerical)
print(paste("test-error=", err))
#"test-error= 0.156583629893238"

table(pred.binary == dat.test$CategoryNumerical)
237/(237+44)
#84%

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
var_importance(shap_result_stages, top_n=20)

## Prepare data for top N variables
shap_long_stages = shap.prep(shap = shap_result_stages,
                             X_train = df.predict, 
                             top_n = 15
)

## Plot shap overall metrics
pdf("plots_supp/xgboost.shapvalues.BarrettsVsPrimary.MutPlusIndel.pdf")
plot.shap.summary(data_long = shap_long_stages)
dev.off()

write.csv(shap_long_stages, file="plots_supp/shap_long_stages.mutplusindel.csv")
