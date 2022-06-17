###################
#### This script evaluates the prognostic potential of individual mutational signatures.

library(ggpubr)

# Signatures:
load("../SigProfiler/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

# Clinical data:
clinical <- read.csv("../SigProfiler/fullset_occams_SA_SAZ_28052020.csv")

# Annotation:
load("processeddata/annotation.sampleIDs.RData")

## Merge signatures and IDs:

sigs.allSamples <- data.frame(sigs.allSamples)
sigs.allSamples$Sample <- rownames(sigs.allSamples)
sigs.allSamples$SBS17 <- sigs.allSamples$SBS17b+sigs.allSamples$SBS17a
sigs.allSamples$SBSddr <- sigs.allSamples$SBS3+sigs.allSamples$SBS8

sigs <- merge(sigs.allSamples, annotation.sampleIDs,
              by.x="Sample", by.y ="Sample", all.x=FALSE, all.y=FALSE)

# Keep only primary tumours:
sigs <- sigs[which(sigs$Category == "PrimaryTumour"),]

clinical$Stage <- sapply(clinical$PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
                         function(x) ifelse(x %in% c("T0","T1","T1a","T1b"),"T1",
                                                     ifelse(x %in% c("T4","T4a"),"T4",
                                                            ifelse(x == "Tx",NA,x))))

sigs.plusclin <- merge(sigs,
                       clinical[,c("Illumina.ID",
                                  "Weeks.Survival.c", 
                                  "DI.PatientDateOfDeath",
                                  "Stage","BarettsAdjacentToTumourMicroscopicGastricM",
                                  "EX.BarrettsOesophagusDiagnosed")],
                       by.x="TumourID",by.y="Illumina.ID",
                       all.x=FALSE, all.y=FALSE)
sigs.plusclin$censor <- sapply(sigs.plusclin$DI.PatientDateOfDeath,
                               function(x) ifelse(is.na(x),0,1))

### 409 samples

library(survminer)
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(sigs.plusclin, time = "Weeks.Survival.c", event = "censor",
                         variables = colnames(sigs.allSamples)[c(1:14,16:17)])

summary(res.cut)
# cutpoint statistic
# SBS17b 0.246777621 2.9337332
# SBS5   0.106412091 1.8542443
# SBS3   0.000000000 0.3404109
# SBS1   0.077088596 2.9175884
# SBS30  0.026458438 3.4936503
# SBS28  0.004215787 0.8492943
# SBS44  0.019553347 1.5086207
# SBS40  0.106490859 1.6651836
# SBS18  0.191566155 2.1094921
# SBS41  0.071867398 2.5132764
# SBS2   0.023909124 2.5106578
# SBS35  0.005581178 3.3900986
# SBS17a 0.128545665 2.5462258
# SBS8   0.057067940 2.6567942
# SBS17  0.382233522 2.5239336
# SBSddr 0.057222502 1.8024690

# 2. Plot cutpoint for S17a
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
for (s in colnames(sigs.allSamples)[c(1:14,16:17)]) {
  pdf(paste0("plots.prognostic/cutoff.",s,".pdf"),onefile=FALSE)
   print(plot(res.cut, s, palette = "npg"))
  dev.off()
}

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
# Weeks.Survival.c censor SBS17b SBS5 SBS3 SBS1 SBS30 SBS28 SBS44 SBS40 SBS18 SBS41
# 1              280      0   high  low  low  low   low  high   low  high   low  high
# 2              161      1    low  low  low  low   low   low   low  high   low  high
# 3              161      1    low  low  low  low   low   low   low  high  high  high
# 4              161      1   high  low  low  low   low   low   low  high   low  high
# 5               78      1    low high high high   low   low   low  high   low  high
# 6               78      1    low high high high   low   low   low  high   low  high
# SBS2 SBS35 SBS17a
# 1  low   low   high
# 2  low   low    low
# 3 high   low    low
# 4  low   low    low
# 5 high  high    low
# 6 high  high    low

# 4. Fit survival curves and visualize
library("survival")

for (s in colnames(sigs.allSamples)[c(1:14,16:17)]) {
  fit <- survfit(as.formula(sprintf("Surv(Weeks.Survival.c, censor) ~ %s", s)),data = res.cat)
  pdf(paste0("plots.prognostic/survival.",s,".pdf"),onefile=FALSE)
  print(ggsurvplot(fit, data = res.cat, risk.table = TRUE, 
                   conf.int = TRUE,pval = TRUE))
  dev.off()
}


##############
### Next, plot the survival curves adjusted by stage:

library(survminer)
sigs.plusclin$S17b_highlow <- sapply(sigs.plusclin$SBS17b, function(x) ifelse(x<0.25,
                                                                              "low","high"))
# high  low 
# 136  273

fit1 <- survfit( Surv(Weeks.Survival.c, censor) ~ S17b_highlow, data = sigs.plusclin )

fit2 <- coxph( Surv(Weeks.Survival.c, censor) ~ S17b_highlow + Stage, data = sigs.plusclin )
pdf("plots.prognostic/survival.S17b.stageAdjusted.pdf")
ggadjustedcurves(fit2, data = sigs.plusclin, 
                 variable = "S17b_highlow",
                    individual.curves=TRUE)
dev.off()
# no signif difference if cut-off 0.17 is used

fit3 <- coxph( Surv(Weeks.Survival.c, censor) ~ S17b_highlow + EX.BarrettsOesophagusDiagnosed, data = sigs.plusclin )
pdf("plots.prognostic/survival.S17b.BEadjusted.pdf")
ggadjustedcurves(fit3, data = sigs.plusclin, 
                 variable = "S17b_highlow",
                 individual.curves=TRUE)
dev.off()

pdf("plots.prognostic/survival.S17bplusBE.pdf")
ggsurvplot(
  fit1,
  data = sigs.plusclin,
  size = 1,                 # change line size
  #palette =
  # c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  group.by = "EX.BarrettsOesophagusDiagnosed", # you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

sigs.plusclin$S30_highlow <- sapply(sigs.plusclin$SBS30, function(x) ifelse(x<0.03,
                                                                              "low","high"))
fit3 <- coxph( Surv(Weeks.Survival.c, censor) ~ S30_highlow + Stage, data = sigs.plusclin )
pdf("plots.prognostic/survival.S30.stageAdjusted.pdf")
ggadjustedcurves(fit3, data = sigs.plusclin, 
                 variable = "S30_highlow",
                 individual.curves=TRUE)
dev.off()

### No stage adjustment, but different cut-off:
sigs.plusclin$S30_highlow <- sapply(sigs.plusclin$SBS30, function(x) ifelse(x<0.05,
                                                                            "low","high"))
fit3 <- survfit(Surv(Weeks.Survival.c, censor) ~ S30_highlow, data = sigs.plusclin )
pdf("plots.prognostic/survival.S30.cutoff005.pdf")
ggsurvplot(
  fit3,
  data = sigs.plusclin,
  size = 1,                 # change line size
  #palette =
  # c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

sigs.plusclin$S30_highlow <- sapply(sigs.plusclin$SBS30, function(x) ifelse(x==0,
                                                                            "low","high"))
fit3 <- survfit( Surv(Weeks.Survival.c, censor) ~ S30_highlow, data = sigs.plusclin )
pdf("plots.prognostic/survival.S30.cutoff0.pdf")
ggsurvplot(
  fit3,
  data = sigs.plusclin,
  size = 1,                 # change line size
  #palette =
  # c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

