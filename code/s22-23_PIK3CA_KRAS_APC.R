########
### PIK3CA, KRAS, APC

library(ggplot2)
library(ggmosaic)

# Load SNV and indel data:
load("data/snvsFull.nonsyn.all.RData")
load("data/indelsFull.all.RData")

# Select only consequential indels:
indels.nonsyn <- indels.all[which(indels.all$V7 %in%
                                    c("inframe_insertion","inframe_deletion",
                                      "frameshift_variant", "frameshift_variant,NMD_transcript_variant",
                                      "TF_binding_site_variant,TFBS_ablation",
                                      "frameshift_variant,splice_region_variant",
                                      "stop_gained,frameshift_variant",
                                      "stop_gained,frameshift_variant,NMD_transcript_variant",
                                      "frameshift_variant,splice_region_variant,NMD_transcript_variant",
                                      "inframe_insertion,NMD_transcript_variant",
                                      "splice_acceptor_variant,frameshift_variant",
                                      "frameshift_variant,initiator_codon_variant",
                                      "inframe_deletion,incomplete_terminal_codon_variant",
                                      "inframe_deletion,NMD_transcript_variant",
                                      "frameshift_variant,stop_lost,NMD_transcript_variant",
                                      "inframe_deletion,splice_region_variant",
                                      "frameshift_variant,stop_lost",
                                      "incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant",
                                      "inframe_deletion,splice_region_variant,NMD_transcript_variant",
                                      "frameshift_variant,stop_lost,splice_region_variant,NMD_transcript_variant",
                                      "stop_gained,frameshift_variant,splice_region_variant",
                                      "inframe_insertion,splice_region_variant",
                                      "frameshift_variant,stop_retained_variant",
                                      "stop_gained,inframe_insertion",
                                      "frameshift_variant,stop_lost,splice_region_variant",
                                      "inframe_insertion,incomplete_terminal_codon_variant,coding_sequence_variant",
                                      "frameshift_variant,stop_retained_variant,NMD_transcript_variant",
                                      "inframe_insertion,splice_region_variant,NMD_transcript_variant",
                                      "stop_gained,inframe_deletion",
                                      "initiator_codon_variant,inframe_deletion",
                                      "stop_gained,inframe_insertion,NMD_transcript_variant",
                                      "stop_gained,inframe_deletion,NMD_transcript_variant")),]

### Select only PIK3CA, APC, KRAS:
snvs.all.nonsyn$Type <- "SNV"
indels.nonsyn$Type <- "Indel"
muts <- rbind(snvs.all.nonsyn, indels.nonsyn)
muts.sel <- muts[which(muts$Gene %in% c("PIK3CA","KRAS","APC")),]

# Load annotation:
load("data/annotation.sampleIDs.RData")

muts.sel.annot <- merge (muts.sel, annotation.sampleIDs,
                         by.x="Sample", by.y="Sample",
                         all.x=TRUE, all.y-FALSE)
save(muts.sel.annot, file="processeddata/muts.sel.annot.APCKRASPIK3CA.RData")
tb <- table(muts.sel.annot[,c("Gene","Category")])

table(annotation.sampleIDs$Category)

sweep(tb, 2, table(annotation.sampleIDs$Category), FUN = '/')

pdf("plots.revision/mosaicPlot.APC_KRAS_PIK3CA.pdf",w=5,h=3)
ggplot(data = muts.sel.annot) +
  geom_mosaic(aes(x = product(Gene), fill=Category))+
  ylab("")+
  xlab("")
dev.off()

library(gridExtra)
pdf("plots.revision/table.APC_KRAS_PIK3CA.pdf")
grid.table(tb)
dev.off()

write.csv(tb, file="plots.revision/APC_KRAS_PIK3CA.table.csv")


### Next, check prognostic potential:

# Clinical data:
clinical <- read.csv("../../SigProfiler/fullset_occams_SA_SAZ_28052020.csv")
clinical$Stage <- sapply(clinical$PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
                         function(x) ifelse(x %in% c("T0","T1","T1a","T1b"),"T1",
                                            ifelse(x %in% c("T4","T4a"),"T4",
                                                   ifelse(x == "Tx",NA,x))))
clinical$APC <- sapply(clinical$Illumina.ID, function(x)
  ifelse(length(which(grepl(x,muts.sel.annot[which(muts.sel.annot$Gene == "APC"),]$Sample)))>0,
         "MUT","WT"))
clinical$KRAS <- sapply(clinical$Illumina.ID, function(x)
  ifelse(length(which(grepl(x,muts.sel.annot[which(muts.sel.annot$Gene == "KRAS"),]$Sample)))>0,
         "MUT","WT"))
clinical$PIK3CA <- sapply(clinical$Illumina.ID, function(x)
  ifelse(length(which(grepl(x,muts.sel.annot[which(muts.sel.annot$Gene == "PIK3CA"),]$Sample)))>0,
         "MUT","WT"))
clinical$censor <- sapply(clinical$DI.PatientDateOfDeath,
                               function(x) ifelse(is.na(x),0,1))

library(ggpubr)
library(survminer)
library(survival)
fitAPC <- survfit( Surv(Weeks.Survival.c, censor) ~ APC, data = clinical )
fitKRAS <- survfit( Surv(Weeks.Survival.c, censor) ~ KRAS, data = clinical )
fitPIK3CA <- survfit( Surv(Weeks.Survival.c, censor) ~ PIK3CA, data = clinical )

pdf("plots.revision/survival.APC.pdf",onefile = FALSE)
ggsurvplot(
  fitAPC,
  data = clinical,
  size = 1,                 # change line size
  palette =
   c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()


pdf("plots.revision/survival.KRAS.pdf",onefile = FALSE)
ggsurvplot(
  fitKRAS,
  data = clinical,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

pdf("plots.revision/survival.PIK3CA.pdf",onefile = FALSE)
ggsurvplot(
  fitPIK3CA,
  data = clinical,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

write.csv(clinical, file="plots_supp/clinical.PIK3CA_KRAS.csv")
