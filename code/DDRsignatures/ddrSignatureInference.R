##################################
### This is the DDR mutational signature analysis.

library(NMF)
library(reshape)
library(ggplot2)
library(ggpubr)

## Read DDR genes:
library(gdata)
ddr <- read.xls("~/Desktop/UCL/supervisedProjects/DanielJacobson/Rotation1-master/DDRpathways.xlsx")

## Load DDR clonality info:
load("gr.clonalityDDR.RData")

gr.clonalityDDR <- merge(gr.clonalityDDR, ddr[,c("Gene.ID",
                                                 "Pathway.1")],
                         by.x="Gene.name", by.y="Gene.ID",
                         all.x=FALSE, all.y=FALSE)

# Load dndscv functional consequence annotation:
load("processeddata/dnds_annot.full.2021.RData")

dnds_annot$chr <- paste0("chr",dnds_annot$chr)
gr.clonalityDDR.plusannot <- merge(gr.clonalityDDR[,c(1:13,19)],
                                   dnds_annot,
                                   by.x=c("Sample","seqnames","start","REF","ALT"),
                                   by.y=c("sampleID","chr","pos","ref","mut"),
                                   all.x=FALSE, all.y = FALSE)
gr.clonalityDDR.plusannot.nonsyn <- gr.clonalityDDR.plusannot[which(gr.clonalityDDR.plusannot$impact != "Synonymous"),]

gr.clonalityDDR.plusannot.nonsyn.pathways <- merge(gr.clonalityDDR.plusannot.nonsyn,
                                                   ddr[,c("Process","Pathway.1")],
                                                   all.x=FALSE, all.y = FALSE,
                                                   by.x="Pathway.1",by.y="Pathway.1")
save(gr.clonalityDDR.plusannot.nonsyn.pathways, file="gr.clonalityDDR.plusannot.nonsyn.pathways.RData")

#########################
#### NMF procedure for early/late patterns in DDR:

### Count how many early vs late mutations there are in each gene in ddr:

#ddr.clonality <- table(gr.clonalityDDR.plusannot.nonsyn[,c("Sample","Clonality")])
#ddr.timing <- table(gr.clonalityDDR.plusannot.nonsyn[,c("Sample","Timing")])
gr.clonalityDDR.plusannot.nonsyn.pathways$ClonalityTiming <- apply(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Clonality","Timing")],
                                                          1, function(x) paste0(x[1],":",x[2]))
gr.clonalityDDR.plusannot.nonsyn.pathways <- gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$Pathway.1!= ""),]
ddr.timing <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","ClonalityTiming")])
ddr.process <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","Process")])
ddr.forNMF <- cbind(ddr.timing,ddr.process)

estim.r <- nmf(data.matrix(t(ddr.forNMF)), 2:10,"brunet", 
               seed=214627,
               nrun=100)
save(estim.r, file="estim.r.nonsyn.process.RData")
pdf("plots.ddr.nonsyn/process.nmf.estimates.ranksurvey.ddr.pdf")
plot(estim.r)
dev.off()

### choose rank 5
res <- nmf(data.matrix(t(ddr.forNMF)), 4,"brunet", 
           seed=214627,
           nrun=100)
save(res, file="res.noclon.nonsyn.4.RData")

w <- basis(res)
# dim(w)
# w
# colnames(w) <- paste0("Sig",1:ncol(w))
# df.w <- melt(w)
# colnames(df.w) <- c("Feature","Signature","Value")
# df.w$Feature <- factor(df.w$Feature,
#                        levels=c("Clonal:Early","Clonal:Late","Subclonal:NA",sort(setdiff(df.w$Feature,c("Clonal:Early","Clonal:Late","Subclonal:NA")))))
# df.w[which(df.w$Feature %in% c("Clonal:Early","Clonal:Late","Subclonal:NA")),]$Value<- df.w[which(df.w$Feature %in% c("Clonal:Early","Clonal:Late","Subclonal:NA")),]$Value/10

df.w <- melt(scoef(w))
colnames(df.w) <- c("Feature","Signature","Proportion")
df.w$Feature <- factor(df.w$Feature,
                       levels=c("Clonal:Early","Clonal:Late","Subclonal:NA",sort(setdiff(df.w$Feature,c("Clonal:Early","Clonal:Late","Subclonal:NA")))))

pdf("plots.ddr.nonsyn/NMF.ddrplot.4sigs.pathway_normalised.2.pdf")
ggplot(df.w,#[which(!(df.w$Feature %in% c("Direct Repair (not in humans)",""))),], 
       aes(Proportion,Feature))+
  geom_bar(stat='identity')+
  facet_wrap(.~Signature,nrow=1)
dev.off()

library(pheatmap)
library(RColorBrewer)
library(viridis)

p <- res@fit@H
rownames(p) <- paste0("DDRsig",1:4)
p.norm <- apply(p, 2, function(x) x/sum(x))
pdf("plots.ddr.nonsyn/heatmap.nmf.4sigs.pathway_final.pdf",h=3)
pheatmap(p.norm,color= rev(cividis(100)),
         show_colnames = FALSE, cutree_cols = 4)
dev.off()

phm <- pheatmap(p.norm,color= rev(cividis(100)),
                show_colnames = FALSE, cutree_cols = 4)
df.ddrgroups <- data.frame(cutree(phm$tree_col, k=4))
df.ddrgroups$Sample <- rownames(df.ddrgroups)
colnames(df.ddrgroups) <- c("DDRgroup","Sample")
df.ddrgroups$G <- 0

pdf("plots.ddr.nonsyn/heatmap.nmf.4sigs.pathway_final_withgroups.pdf",h=3)
pheatmap(p.norm,color= rev(cividis(100)),
         show_colnames = FALSE, cutree_cols = 4,
         annotation = df.ddrgroups[,c(1,3)])
dev.off()

### Finally, check survival info:
# Clinical data:
clinical <- read.csv("../SigProfiler/fullset_occams_SA_SAZ_28052020.csv")
clinical$Stage <- sapply(clinical$PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
                         function(x) ifelse(x %in% c("T0","T1","T1a","T1b"),"T1",
                                            ifelse(x %in% c("T4","T4a"),"T4",
                                                   ifelse(x == "Tx",NA,x))))
df.ddrgroups$TumourID <- sapply(df.ddrgroups$Sample,
                                 function(x) strsplit(x,"_vs_")[[1]][1])
df.ddr.plusclin <- merge(df.ddrgroups, clinical[,c("Illumina.ID",
                                                      "Weeks.Survival.c", 
                                                      "DI.PatientDateOfDeath",
                                                      "Stage")],
                         by.x="TumourID",by.y="Illumina.ID",
                         all.x=FALSE, all.y=TRUE)
df.ddr.plusclin$censor <- sapply(df.ddr.plusclin$DI.PatientDateOfDeath,
                               function(x) ifelse(is.na(x),0,1))
df.ddr.plusclin[which(is.na(df.ddr.plusclin$DDRgroup)),]$DDRgroup <- "None"
require("survival")
library(survminer)
fit <- survfit(Surv(Weeks.Survival.c, censor) ~ DDRgroup, data = df.ddr.plusclin)

pdf("plots.ddr.nonsyn/survival.byDDRgroup.allSamples.2.pdf",onefile = FALSE,w=5,h=5)
ggsurvplot(
  fit,
  data = df.ddr.plusclin,
  size = 1,                 # change line size
  palette =
    c("#b2e2e2", "#66c2a4","#2ca25f","#006d2c","grey"),# custom color palettes
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()


df.ddr.plusclin$DDR2vsrest <- sapply(df.ddr.plusclin$DDRgroup,
                                     function(x) ifelse(x==2,"2","other"))
fit2 <- survfit(Surv(Weeks.Survival.c, censor) ~ DDR2vsrest, data = df.ddr.plusclin)

pdf("plots.ddr.nonsyn/survival.byDDRgroup_2vsrest.pdf",onefile = FALSE)
ggsurvplot(
  fit2,
  data = df.ddr.plusclin,
  size = 1,                 # change line size
  #palette =
  # c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

# Annotation:
load("processeddata/annotation.sampleIDs.RData")
df.sig.norm <- data.frame(t(p.norm))
df.sig.norm$Sample <- rownames(df.sig.norm)

df.norm.plusannot <- merge(df.sig.norm, annotation.sampleIDs,
              by.x="Sample", by.y ="Sample", all.x=FALSE, all.y=FALSE)
df.norm.plusannot[which(df.norm.plusannot$Category == "LymphNode"),]$Category <- "Metastasis"

df.norm.plusannot$Category <- factor(df.norm.plusannot$Category,
                                     levels=c("Barretts","PrimaryTumour","Metastasis"))
my_comparisons <- list(c("Barretts","PrimaryTumour"),
                  c("PrimaryTumour","Metastasis"),
                  c("Barretts","Metastasis"))
pdf("plots.ddr.nonsyn/DDRsig1.compareStages.pdf",w=3,h=4)
ggboxplot(df.norm.plusannot, x = "Category", y = "DDRsig1",
          color = "Category", palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("plots.ddr.nonsyn/DDRsig2.compareStages.pdf",w=3,h=4)
ggboxplot(df.norm.plusannot, x = "Category", y = "DDRsig2",
          color = "Category", palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("plots.ddr.nonsyn/DDRsig3.compareStages.pdf",w=3,h=4)
ggboxplot(df.norm.plusannot, x = "Category", y = "DDRsig3",
          color = "Category", palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("plots.ddr.nonsyn/DDRsig4.compareStages.pdf",w=3,h=4)
ggboxplot(df.norm.plusannot, x = "Category", y = "DDRsig4",
          color = "Category", palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

df.norm.plusannot$TumourID <- sapply(df.norm.plusannot$Sample,
                                     function(x) strsplit(x,"_vs_")[[1]][1])

### Merge clinical data and signatures too:
df.d <- merge(df.ddr.plusclin, df.norm.plusannot,
              by.x="TumourID", by.y="TumourID", all.x=FALSE, all.y=FALSE)

pdf("plots.ddr.nonsyn/DDRsig1.compareClinicalStages.pdf",w=3,h=4)
ggboxplot(df.d, x = "Stage", y = "DDRsig1",
          color = "Stage", #palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("plots.ddr.nonsyn/DDRsig2.compareClinicalStages.pdf",w=3,h=4)
ggboxplot(df.d, x = "Stage", y = "DDRsig2",
          color = "Stage", #palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("plots.ddr.nonsyn/DDRsig3.compareClinicalStages.pdf",w=3,h=4)
ggboxplot(df.d, x = "Stage", y = "DDRsig3",
          color = "Stage", #palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf("plots.ddr.nonsyn/DDRsig4.compareClinicalStages.pdf",w=3,h=4)
ggboxplot(df.d, x = "Stage", y = "DDRsig4",
          color = "Stage", #palette =c("#FC4E07","#00AFBB", "#2D8658"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()

### Check which signature is dominant:
df.d$Dominant <- apply(df.d[,c("DDRsig1","DDRsig2","DDRsig3","DDRsig4")],1,
                       function(x) c("DDRsig1","DDRsig2","DDRsig3","DDRsig4")[which(x==max(x))])

fit <- survfit(Surv(Weeks.Survival.c, censor) ~ Dominant, data = df.d)

pdf("plots.ddr.nonsyn/survival.byDDRdominantsignature.pdf",onefile = FALSE,w=5,h=5)
ggsurvplot(
  fit,
  data = df.d,
  size = 1,                 # change line size
  palette =
    c("#b2e2e2", "#66c2a4","#2ca25f","#006d2c","grey"),# custom color palettes
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

save(df.d, file="DDR.df.d.RData")
save(df.norm.plusannot, file="DDR.df.norm.plusannot.RData")

#### Finally, are there any expression patterns associated with this?

# Check expression of relevant DDR and downstream pathways:

# Next, load expression data and map against genomic changes:
load("../tracksig_barrettsmetsprimaries/oac.expr.final.RData")
load("../tracksig_barrettsmetsprimaries/mat.expr.RData")

ddrlist <- NULL
i<-0
for (c in unique(ddr$Pathway.1)) {
  i<-i+1
  ddrlist[[i]] <- unique(ddr[which(ddr$Pathway.1 == c),]$Gene.ID)
}
names(ddrlist) <- unique(ddr$Pathway.1)

library(GSVA)
m <- t(mat.expr[,-1])
rownames(m) <- colnames(mat.expr[,-1])
ddrscores <- gsva(as.matrix(m),ddrlist, 
                       method="gsva", kcdf="Gaussian",
                       mx.diff=TRUE, abs.ranking=FALSE)
colnames(ddrscores) <- mat.expr$DNAid

df.ddrscores <- melt(ddrscores)
colnames(df.ddrscores) <- c("Pathway","Sample","Score")
df.ddrgroups$DNAid <- sapply(df.ddrgroups$Sample,
                             function(x) strsplit(x,"_vs_")[[1]][1])
df.ddrscores <- merge(df.ddrscores,
                      df.ddrgroups,
                      by.x="Sample",by.y = "DNAid",
                      all.x=FALSE,all.y=FALSE)

pdf("plots.ddr.nonsyn/pathwayActivityCompared.betweenGroups.pdf",w=20,h=15)
ggboxplot(df.ddrscores, x = "DDRgroup", y = "Score",
          color = "DDRgroup", palette =c("#b2e2e2", "#66c2a4","#2ca25f","#006d2c"),
          add = "jitter")+
  stat_compare_means(label.y = 1)+
  facet_wrap(~Pathway)
dev.off()


pdf("plots.ddr.nonsyn/pathwayActivityCompared.forEachGroups.pdf",w=20,h=10)
ggboxplot(df.ddrscores, x = "Pathway", y = "Score",
          color = "Pathway", #palette =c("#b2e2e2", "#66c2a4","#2ca25f","#006d2c"),
          add = "jitter")+
  stat_compare_means(label.y = 1)+
  facet_wrap(~DDRgroup)
dev.off()

df.ddrscores <- df.ddrscores[which(df.ddrscores$Pathway !=""),]

library(dplyr)
df.select <- data.frame(df.ddrscores%>%
                                     group_by(Pathway,DDRgroup)%>% 
                                     summarise(Mean=mean(Score), Median=median(Score)))
#

# Set all values that are not particularly enriched or depleted to 0:

df.select.modi <- df.select
df.select.modi[which((df.select.modi$Median>((-1)*0.2))&
                       (df.select.modi$Median<0.2)),]$Median <- 0
df.select.modi[which(df.select.modi$Median<((-1)*0.5)),]$Median <- -0.5
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
pdf("plots.ddr.nonsyn/ballonplot.pathwaysComparedBetweenGroups.pdf")
ggballoonplot(df.select.modi,
              x = "DDRgroup", y = "Pathway",
              fill = "Median",color="white",
              show.zeros=TRUE, size.range = c(1, 8)) +
  #scale_fill_gradientn(colors = my_cols)+
  gradient_fill(c("blue", "white", "red"))+
  guides(size = "none")
dev.off()

## Next try with dominant sig:

df.norm.plusannot$Dominant <- apply(df.norm.plusannot[,c("DDRsig1","DDRsig2","DDRsig3","DDRsig4")],1,
                       function(x) c("DDRsig1","DDRsig2","DDRsig3","DDRsig4")[which(x==max(x))])

df.ddrscores <- merge(df.ddrscores,
                      df.norm.plusannot,
                      by.x="Sample",by.y = "TumourID",
                      all.x=FALSE,all.y=FALSE)

pdf("plots.ddr.nonsyn/pathwayActivityCompared.betweenDominantSig.pdf",w=20,h=15)
ggboxplot(df.ddrscores, x = "Dominant", y = "Score",
          color = "Dominant", palette =c("#b2e2e2", "#66c2a4","#2ca25f","#006d2c"),
          add = "jitter")+
  stat_compare_means(label.y = 1)+
  facet_wrap(~Pathway)
dev.off()


pdf("plots.ddr.nonsyn/pathwayActivityCompared.forEachDominantSig.pdf",w=20,h=10)
ggboxplot(df.ddrscores, x = "Pathway", y = "Score",
          color = "Pathway", #palette =c("#b2e2e2", "#66c2a4","#2ca25f","#006d2c"),
          add = "jitter")+
  stat_compare_means(label.y = 1)+
  facet_wrap(~Dominant)
dev.off()
df.ddrscores <- df.ddrscores[which(df.ddrscores$Pathway !=""),]

df.select <- data.frame(df.ddrscores%>%
                          group_by(Pathway,Dominant)%>% 
                          summarise(Mean=mean(Score), Median=median(Score)))
#

# Set all values that are not particularly enriched or depleted to 0:

df.select.modi <- df.select
df.select.modi[which((df.select.modi$Median>((-1)*0.2))&
                       (df.select.modi$Median<0.2)),]$Median <- 0
df.select.modi[which(df.select.modi$Median<((-1)*0.5)),]$Median <- -0.5
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
pdf("plots.ddr.nonsyn/ballonplot.pathwaysComparedBetweenDominantSigs.pdf")
ggballoonplot(df.select.modi,
              x = "Dominant", y = "Pathway",
              fill = "Median",color="white",
              show.zeros=TRUE, size.range = c(1, 8)) +
  #scale_fill_gradientn(colors = my_cols)+
  gradient_fill(c("blue", "white", "red"))+
  guides(size = "none")
dev.off()