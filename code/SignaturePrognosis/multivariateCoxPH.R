###################
#### This script performs multivariate CoxPH analysis.

library(survminer)

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

sigs.plusclin <- merge(sigs,
                       clinical,
                       by.x="TumourID",by.y="Illumina.ID",
                       all.x=FALSE, all.y=FALSE)
sigs.plusclin$censor <- sapply(sigs.plusclin$DI.PatientDateOfDeath,
                               function(x) ifelse(is.na(x),0,1))

require("survival")

model0 <- coxph( Surv(Weeks.Survival.c, censor) ~ SBS17a+
                  SBS17b+SBS2+SBS3+SBS8+
                  SBS18+SBS30+SBS44+SBS35+SBS41+SBS28+
                  SBS1+SBS5+SBS40,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.signaturesOnly.pdf",onefile = FALSE)
ggforest(model0)
dev.off()
 
model0.red <- coxph( Surv(Weeks.Survival.c, censor) ~ SBS8+
                   SBS30+SBS35+SBS41,
                 data = sigs.plusclin)
pdf("plots.prognostic/multivariateCox.signaturesOnlyReduced.pdf",onefile = FALSE)
ggforest(model0.red)
dev.off()


sigs.plusclin$Tstage <- sapply(sigs.plusclin$PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
                               function(x) ifelse(grepl("T1",x),"T1",
                                                         ifelse(grepl("T4",x),"T4",x)))
sigs.plusclin$TstageOverall <- sapply(sigs.plusclin$Tstage,
                              function(x) ifelse(x %in% c("T1","T2"),"Early",
                                                 ifelse(x %in% c("T3","T4"),"Late",NA)))
sigs.plusclin$Nstage <- sapply(sigs.plusclin$PS.NStage.PrimaryTumour.FinalPretreatmentStaging.TNM7,
                              function(x) ifelse(x =="N0","N0",
                                                 ifelse(x %in% c("N1","N2","N3"),"N1+",NA)))
sigs.plusclin$Mstage <- sapply(sigs.plusclin$PS.MStage.PrimaryTumour.FinalPretreatmentStaging,
                               function(x) ifelse(x =="M0","M0",
                                                  ifelse(x =="M1","M1",NA)))
sigs.plusclin$Treatment <- sapply(sigs.plusclin$TP.CurativeTreatmentModality,
                                  function(x) ifelse(x %in% c("surgery only",
                                                              "chemo-radiotherapy and surgery"),x,
                                                     "other"))
sigs.plusclin$EX.RefluxFrequencyOfSymptoms <- factor(sigs.plusclin$EX.RefluxFrequencyOfSymptoms)

model <- coxph( Surv(Weeks.Survival.c, censor) ~ SBS8+
                  SBS30+SBS35+SBS41+
                  Gender+
                  EX.PatientOnAcidSuppressant+
                  #EX.RefluxFrequencyOfSymptoms+ 
                  #EX.ProtonPumpInhibitor+
                  EX.IsSmoker+EX.HeavyDrinker+
                  EX.BarrettsOesophagusDiagnosed+
                  TstageOverall+
                  Nstage+
                  Mstage+
                  EX.CurrentBMI+
                  RD.SignetRingCellsPresent+
                  Treatment+
                  PS.SiewertClassification,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.full.pdf",onefile = FALSE)
ggforest(model)
dev.off()

model <- coxph( Surv(Weeks.Survival.c, censor) ~ 
                  SBS30+SBS41+
                  Gender+
                  EX.PatientOnAcidSuppressant+
                  #EX.RefluxFrequencyOfSymptoms+ 
                  #EX.ProtonPumpInhibitor+
                  EX.IsSmoker+EX.HeavyDrinker+
                  EX.BarrettsOesophagusDiagnosed+
                  TstageOverall+
                  Nstage+
                  Mstage+
                  EX.CurrentBMI+
                  Treatment+
                  PS.SiewertClassification,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.optimal.pdf",onefile = FALSE,w=10)
ggforest(model)
dev.off()

# This is not a very good model!
model.minimal <- coxph( Surv(Weeks.Survival.c, censor) ~ 
                  SBS30+
                  #EX.RefluxFrequencyOfSymptoms+ 
                  #EX.ProtonPumpInhibitor+
                  TstageOverall+
                  EX.CurrentBMI+
                  Treatment,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.minimal.pdf",onefile = FALSE,w=10)
ggforest(model.minimal)
dev.off()

model <- coxph( Surv(Weeks.Survival.c, censor) ~ 
                  SBS17b+
                  TstageOverall+
                  Treatment,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.S17badjusted.pdf",onefile = FALSE,w=10)
ggforest(model)
dev.off()

require("survival")
sigs.plusclin$SBS17bhighlow <- sapply(sigs.plusclin$SBS17b, function(x) ifelse(x<=0.25,"low","high"))
fit <- survfit(Surv(Weeks.Survival.c, censor) ~ SBS17bhighlow, data = sigs.plusclin[which(sigs.plusclin$TstageOverall=="Early"),])
pdf("plots.prognostic/KM.SBS17b.earlytumours.pdf",
    onefile=FALSE)
ggsurvplot(
  fit,
  data = sigs.plusclin[which(sigs.plusclin$TstageOverall=="Early"),],
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

fit <- survfit(Surv(Weeks.Survival.c, censor) ~ SBS17bhighlow, data = sigs.plusclin[which(sigs.plusclin$TstageOverall=="Late"),])
pdf("plots.prognostic/KM.SBS17b.latetumours.pdf",
    onefile=FALSE)
ggsurvplot(
  fit,
  data = sigs.plusclin[which(sigs.plusclin$TstageOverall=="Early"),],
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

### Load TMB:
load("../dndscv/vcf.primaries.plusextra.forCluster.RData")

## Calculate TMB:
tb <- table(vcf.primaries.plusextra$sampleID)
df.tb <- data.frame(tb)
df.tb$Sample <- names(tb)
df.tb$TMB <- log10(df.tb$Freq+1)

sigs.plusclin <- merge(sigs.plusclin, df.tb[,c("Sample","TMB")],
                       by.x="Sample", by.y="Sample",
                       all.x=FALSE, all.y=FALSE)

model <- coxph( Surv(Weeks.Survival.c, censor) ~ SBS8+
                  SBS30+SBS35+SBS41+
                  TMB+
                  Gender+
                  EX.PatientOnAcidSuppressant+
                  #EX.RefluxFrequencyOfSymptoms+ 
                  #EX.ProtonPumpInhibitor+
                  EX.IsSmoker+EX.HeavyDrinker+
                  EX.BarrettsOesophagusDiagnosed+
                  TstageOverall+
                  Nstage+
                  Mstage+
                  EX.CurrentBMI+
                  RD.SignetRingCellsPresent+
                  Treatment+
                  PS.SiewertClassification,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.full.plusTMB.pdf",onefile = FALSE)
ggforest(model)
dev.off()

model <- coxph( Surv(Weeks.Survival.c, censor) ~ 
                  SBS30+TMB+
                  Gender+
                  EX.PatientOnAcidSuppressant+
                  #EX.RefluxFrequencyOfSymptoms+ 
                  #EX.ProtonPumpInhibitor+
                  EX.IsSmoker+EX.HeavyDrinker+
                  EX.BarrettsOesophagusDiagnosed+
                  TstageOverall+
                  Nstage+
                  Mstage+
                  EX.CurrentBMI+
                  Treatment+
                  PS.SiewertClassification,
                data = sigs.plusclin )
pdf("plots.prognostic/multivariateCox.optimal.plusTMB.onlysbs30.pdf",onefile = FALSE,w=10)
ggforest(model)
dev.off()