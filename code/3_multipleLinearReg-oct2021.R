##############
#### Trying to model transitions between Barrett's, primaries and mets using HMMs.

library(depmixS4)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(cowplot)
library(reshape)

# Load timed mutational signature data:
load("data/sigs.timed.annot.RData")
sigs.timed.annot[which(sigs.timed.annot$Category == "LymphNode"),]$Category <- "Metastasis"
df.sigs.timed <- melt(sigs.timed.annot)
#df.sigs.timed[which(df.sigs.timed$Category == "LymphNode"),]$Category <- "Metastasis"

# mod = depmix(Category ~SBS1,#SBS1+SBS2+SBS3+SBS5+SBS17a+SBS17b+
#                #SBS18+SBS28+SBS30+SBS35+SBS40+SBS41+SBS44,
#              nstates = 3,
#              transition = ~Timing+Clonality,
#              family = gaussian(),
#              data=sigs.timed.annot[,-c(1,3:14,17)])


### Select train and test sets, dividing into 2/3 for train, 1/3 for test:
set.seed(19875) 
randTrainingSet <- sample(1:nrow(sigs.timed.annot), 662, replace=F)
dat.train <- sigs.timed.annot[randTrainingSet,]
dat.test <- sigs.timed.annot[setdiff(1:nrow(sigs.timed.annot),randTrainingSet),]
dat.train.x <- as.matrix(dat.train[,2:17])
dat.train.y <- dat.train$Category
dat.test.x <- as.matrix(dat.test[,2:17])
dat.test.y <- dat.test$Category

#############
# define model grid for best subset regression
# defines which predictors are on/off; all combinations presented
# from https://stackoverflow.com/questions/41061729/how-to-get-the-best-subset-for-a-multinomial-regression-in-r
model.grid <- function(n){
  n.list <- rep(list(0:1), n)
  expand.grid(n.list)
}

# function for best subset regression
# ranks predictor combos using 5 selection criteria

best.subset <- function(y, x.vars, data){
  # y       character string and name of dependent variable
  # xvars   character vector with names of predictors
  # data    training data with y and xvar observations
  
  require(dplyr)
  require(purrr)
  require(magrittr)
  require(forecast)
  
  length(x.vars) %>%
    model.grid %>%
    apply(1, function(x) which(x > 0, arr.ind = TRUE)) %>%
    map(function(x) x.vars[x]) %>%
    .[2:dim(model.grid(length(x.vars)))[1]] %>%
    map(function(x) fit=multinom(paste0(y, " ~ 0+", paste(x, collapse = "+")), data = data)) %>%
    map("AIC") %>%
    do.call(rbind, .) %>%
    cbind(model.grid(length(x.vars))[-1, ], .) 
}

library(fpp2)

# test the function
compareall <- best.subset("Category", colnames(dat.train)[c(2:15)], dat.train) 
colnames(compareall) <- c(colnames(dat.train)[c(2:15)],"AIC")
coptimum <- compareall[which(compareall$AIC == min(compareall$AIC)),]
# without clonality and timing:
#SBS17b SBS5 SBS3 SBS1 SBS30 SBS28 SBS44 SBS40 SBS18 SBS41 SBS2 SBS35 SBS17a SBS8      AIC
#1015      0    1    1    0     1     1     1     1     1     1    0     0      0    0 768.6652
save(compareall, file="compareall.RData")
save(coptimum, file="coptimum.RData")

#with clonality and timing:
#SBS17b SBS5 SBS3 SBS1 SBS30 SBS28 SBS44 SBS40 SBS18 SBS41 SBS2
# 19483      0    1    0    1     1     0     0     0     0     0    1
# SBS35 SBS17a SBS8 Clonality      AIC
# 19483     1      0    0         1 859.1959

# Reduced, optimal model with minimum AIC:
featuresToKeep <- names(coptimum)[1:(ncol(coptimum)-1)]#[which(coptimum==1)])
glm.fit=multinom(paste0("Category~ 0+", paste(featuresToKeep, collapse = "+")), 
                 data = dat.train)
summary(glm.fit)
z <- summary(glm.fit)$coefficients/summary(glm.fit)$standard.error
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p
#anova(glm.fit,glm.fit.null)

save(glm.fit, file="glm.fit.optimised.cltim.RData")
# AIC: 679.423 

## Predictions:
predicted=predict(glm.fit,dat.test,type="probs")
bpp=cbind(dat.test, predicted)
bpp2 <- melt(bpp[,c("Category",featuresToKeep,#glm.fit$coefnames,
                    "Barretts","PrimaryTumour","Metastasis")],
             id.vars = c("Category",featuresToKeep))#glm.fit$coefnames))
colnames(bpp2)[c(ncol(bpp2)-1,ncol(bpp2))] <- c("Prediction","Probability")
bpp.truth <- unique(bpp2[,c("Category",featuresToKeep)])#glm.fit$coefnames)])
bpp.truth$StagGroup <- sapply(bpp.truth$Category,
                         function(x) ifelse(x=="Barretts",0,
                                            ifelse(x=="PrimaryTumour",0.5,1)))
bpp.truth$Category <- factor(bpp.truth$Category,levels=c("Barretts","PrimaryTumour","Metastasis"))
bpp2$Category <- factor(bpp2$Category,levels=c("Barretts","PrimaryTumour","Metastasis"))

p.pred <- list()
p.truth <- list()
i <- 0 
for (r in featuresToKeep){#glm.fit$coefnames) {
  i <- i+1
  # p.pred.current <- ggplot(bpp2, aes_string(x = r, y = "Probability", colour = "Prediction")) +
  #   geom_line() + facet_grid(Prediction ~ ., scales="free")+
  #   ylim(c(0,1))+
  #   ylab("Probability")
  p.pred.current <- ggplot(bpp2, aes_string(x = r, y = "Probability", colour = "Prediction")) +
    geom_point(size=1) + 
    geom_smooth(method="loess",se=TRUE)+
    scale_color_manual(values=c("#6ba393", "#e8b44d","#94788c"))+
    facet_grid(Prediction ~ ., scales="free")+
    ylim(c(0,1))+
    ylab("Probability")+
    theme(legend.position = "none") # no legend
  p.truth.current <- ggplot(bpp.truth,aes_string(x=r,y="Category"))+
    geom_boxplot(width=5)+
    geom_point(size=2, alpha=0.5)+
    theme(plot.margin=unit(c(t=0,r=4.6,b=0,l=0),"cm"))+
    ylab("Category")
  p.pred[[i]] <- p.pred.current
  p.truth[[i]] <- p.truth.current
  
  pdf(paste0("plots.multinom.june2022/noclontim_",r,"predictedCategory.pdf"))
  print(plot_grid(p.pred.current,p.truth.current,nrow=2,align="hv", axis="l",
                  rel_heights = c(1,0.2)))
  dev.off()
}


pdf(paste0("plots.multinom.june2022/multinom.selected.together.pdf"),w=10,h=6)
print(plot_grid(p.pred[[11]],p.pred[[10]],p.truth[[11]],p.truth[[10]],
                nrow=2,align="hv", axis="l",
                rel_heights = c(1,0.2)))
dev.off()

shap_values=predict(glm.fit, dat.test, predcontrib = TRUE, approxcontrib = F)

bpp.truth$Category == bpp2$Category

write.csv(bpp2, file="plots_3/bpp2.csv")
