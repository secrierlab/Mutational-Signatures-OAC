##################
### BER signature survival + mutations

# Clinical data:
clinical <- read.csv("../../SigProfiler/fullset_occams_SA_SAZ_28052020.csv")
clinical$Stage <- sapply(clinical$PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
                         function(x) ifelse(x %in% c("T0","T1","T1a","T1b"),"T1",
                                            ifelse(x %in% c("T4","T4a"),"T4",
                                                   ifelse(x == "Tx",NA,x))))

## Load mutsig data:
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples[sigs.allSamples<0.05] <- 0

sigs.allSamples$Sample <- sapply(rownames(sigs.allSamples), function(x)
  strsplit(x,"_vs_")[[1]][1])

sigs.plusclin <- merge(sigs.allSamples, clinical,
              by.x="Sample", by.y ="Illumina.ID", all.x=FALSE, all.y=FALSE)


sigs.plusclin$BER <- sapply(sigs.plusclin$SBS30, function(x)
  ifelse(x<0.05,"<5%",">5%"))

sigs.plusclin$censor <- sapply(sigs.plusclin$DI.PatientDateOfDeath,
                          function(x) ifelse(is.na(x),0,1))

library(ggpubr)
library(survminer)
library(survival)
fit <- survfit( Surv(Weeks.Survival.c, censor) ~ BER, data = sigs.plusclin )

pdf("plots_7/survival.SBS30_BER_5percCutoff.pdf",onefile = FALSE)
ggsurvplot(
  fit,
  data = sigs.plusclin,
  size = 1,                 # change line size
  palette =
    c("#BE92A2","#575D90"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

write.csv(sigs.plusclin, file="plots_7/sigs.plusclin.ber.csv")

### Any specific BER/HR defects?

## Read DDR genes:
library(gdata)
ddr <- read.xls("~/Desktop/UCL/supervisedProjects/DanielJacobson/Rotation1-master/DDRpathways.xlsx")
ddr.ber <- ddr[which(ddr$Pathway.1=="BER"),]$Gene.ID
ddr.hrd <- ddr[which(ddr$Pathway.1=="HR (Homologous Recombination)"),]$Gene.ID

## Load snvs and indels:
load("vep.snvs.nonsyn")
load("processeddata/indelsFull.all.RData")
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
vep.snvs.nonsyn$Type <- "SNV"
indels.nonsyn$Type <- "Indel"
muts <- rbind(vep.snvs.nonsyn, indels.nonsyn)
muts.sel.ber <- muts[which(muts$Gene %in% ddr.ber),]
muts.sel.hrd <- muts[which(muts$Gene %in% ddr.hrd),]

## BER defects:
intersect(unique(muts.sel.ber$Sample),rownames(sigs.allSamples)[sigs.allSamples$SBS30>=0.05])
### 9 out of 59 samples with SNVs and indels in BER pathway

## HR defects:
load("processeddata/hrd.RData")
intersect(unique(muts.sel.hrd$Sample),hrd)
### 20 out of 71 samples with SNVs and indels in BER pathway

### What about CN changes? --dels?


