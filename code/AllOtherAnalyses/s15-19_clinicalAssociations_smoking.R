####################
### This script links clinical data with signature prevalence.

library(reshape)
library(ggpubr)
library(stringr)

load("data/annotation.sampleIDs.RData")
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

annotation.sampleIDs$OCCAMS_ID <- sapply(annotation.sampleIDs$OCCAMS_ID,
                                         function(x) ifelse(grepl("OC/",x),str_replace(x,"OC/","OCCAMS/"),x))

sigs.allSamples$Sample <- rownames(sigs.allSamples)
sigs.allSamples$TumourID <- sapply(rownames(sigs.allSamples),
                                 function(x) strsplit(x,"_vs_")[[1]][1])
sigs <- merge(sigs.allSamples, annotation.sampleIDs[,c("Sample","Category","OCCAMS_ID")],
              all.x = FALSE, all.y=FALSE,
              by.x="Sample",by.y="Sample")

### Read in the clinical data:
barr <- read.csv("../Barrett's Clinical Data for Maria.csv")
oac <- read.csv("../Primary Tumour Clinical data for Maria.csv")

### Analysis for OAC:
sigs.oac <- merge(sigs, oac[,c("OCCAMS_ID","EX.IsSmoker")],
                  by.x="OCCAMS_ID", by.y="OCCAMS_ID",
                  all.x=FALSE, all.y=FALSE)
sigs.oac.primary <- sigs.oac[which(sigs.oac$Category == "PrimaryTumour"),]
sigs.barr <- merge(sigs, barr[,c("ID","EverSmoked","PatientGrade")],
                  by.x="TumourID", by.y="ID",
                  all.x=FALSE, all.y=FALSE)
sigs.oac.barr <- sigs.barr[which(sigs.barr$Category == "Barretts"),]


df.oac <- melt(sigs.oac.primary,id.vars = c("OCCAMS_ID","Sample","TumourID",
                                            "Category","EX.IsSmoker"))
df.oac$EverSmoked <- sapply(df.oac$EX.IsSmoker, function(x)
  ifelse(is.na(x),x,ifelse(x=="never","N","Y")))
df.barr <- melt(sigs.oac.barr,id.vars = c("OCCAMS_ID","Sample","TumourID",
                                            "Category","EverSmoked","PatientGrade"))
df.barr[which(df.barr$EverSmoked == "Unknown"),]$EverSmoked <- NA

df.oac$EX.IsSmoker <- factor(df.oac$EX.IsSmoker,
                             levels=c("never","former","current"))
my_comparisons <- list(c("N","Y"))
df.oac$variable <- factor(df.oac$variable,
                           levels=sort(unique(as.character(df.oac$variable))))

df.oac$EverSmoked <- factor(df.oac$EverSmoked,
                             levels=c("Y","N"))
df.barr$EverSmoked <- factor(df.barr$EverSmoked,
                            levels=c("Y","N"))

pdf("plots_supp/primaries.EverSmoked.pdf",w=10,h=8)
ggboxplot(df.oac[which(!is.na(df.oac$EverSmoked)),], x = "EverSmoked", y = "value",
          color = "EverSmoked", 
          add = "jitter")+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(comparisons = my_comparisons)+
  facet_wrap(~variable,nrow=2, scale="free")+
  xlab("Ever smoked")+
  ylab("Mutational signature contribution (%)")
dev.off()
  
pdf("plots_supp/barretts.EverSmoked.pdf",w=10,h=8)
ggboxplot(df.barr[which(!is.na(df.barr$EverSmoked)),], x = "EverSmoked", y = "value",
          color = "EverSmoked", 
          add = "jitter")+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(comparisons = my_comparisons)+
  facet_wrap(~variable,nrow=2, scale="free")+
  xlab("Ever smoked")+
  ylab("Mutational signature contribution (%)")
dev.off()

write.csv(df.barr[which(!is.na(df.barr$EverSmoked)),],
          file="plots_supp/smoking.sigs.barr.csv")
write.csv(df.oac[which(!is.na(df.oac$EverSmoked)),],
          file="plots_supp/smoking.sigs.prim.csv")

### Numbers:
table(unique(df.barr[,c("TumourID","EverSmoked")])$EverSmoked)
# Y  N 
# 79 31

table(unique(df.oac[,c("TumourID","EverSmoked")])$EverSmoked)
# Y   N 
# 365 143

# TMB:
load("../tmb.annot.RData")
df.oac.tmb <- merge(unique(df.oac[,c("TumourID","Category","EverSmoked")]), tmb.annot[,c("TumourID","TMB")],
                by.x="TumourID",by.y="TumourID", all.x=FALSE, all.y=FALSE)
df.barr.tmb <- merge(unique(df.barr[,c("TumourID","Category","EverSmoked")]), tmb.annot[,c("TumourID","TMB")],
                    by.x="TumourID",by.y="TumourID", all.x=FALSE, all.y=FALSE)
df.tmb <- rbind(df.oac.tmb, df.barr.tmb)

pdf("plots_supp/TMB.EverSmoked.pdf",w=4,h=6)
ggboxplot(df.tmb[which(!is.na(df.tmb$EverSmoked)),], x = "EverSmoked", y = "TMB",
          color = "EverSmoked", 
          add = "jitter")+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(comparisons = my_comparisons)+
  facet_wrap(~Category,nrow=1)+
  xlab("Ever smoked")+
  ylab("Mutational signature contribution (%)")
dev.off()

write.csv(df.tmb[which(!is.na(df.tmb$EverSmoked)),], file="plots_supp/tmb.eversmoked.csv")

### Difference in TMB by Barrett type: 
df.barr.tmb.type <- merge(unique(df.barr[,c("TumourID","Category","EverSmoked","PatientGrade")]), 
                          tmb.annot[,c("TumourID","TMB")],
                     by.x="TumourID",by.y="TumourID", all.x=FALSE, all.y=FALSE)

pdf("plots_supp/TMB.EverSmoked.BarrettCategory.pdf",w=4,h=6)
ggboxplot(df.barr.tmb.type[which(!is.na(df.barr.tmb.type$EverSmoked)),], x = "EverSmoked", y = "TMB",
          color = "EverSmoked", 
          add = "jitter")+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(comparisons = my_comparisons)+
  facet_wrap(~PatientGrade,nrow=1)+
  xlab("Ever smoked")+
  ylab("Mutational signature contribution (%)")
dev.off()

write.csv(df.barr.tmb.type[which(!is.na(df.barr.tmb.type$EverSmoked)),],
          file="plots_supp/tmb.barrtype.csv")

tb1 <- table(df.barr.tmb.type[,c("EverSmoked","PatientGrade")])
pvals <- NULL
or <- NULL
for (p in 1:ncol(tb1)) {
  f <- fisher.test(cbind(tb1[,p], rowSums(tb1[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb1)
names(or) <- colnames(tb1)
pvals[pvals<0.05]
# NDBO_NP 
# 0.008634826 
or[pvals<0.05]
# NDBO_NP 
# 0.2404003 

fisher.test(matrix(c(8,5,10,3),nrow=2))
# not signif when comparing NP with PP

df.barr.tmb.type$PatientGrade <- factor(df.barr.tmb.type$PatientGrade,
                                        levels=c("NDBO_NP","NDBO_PP",
                                                 "LGD","HGD","IMC","trios"))

library(ggmosaic)
pdf("plots_supp/mosaic.smoking.BarrettCateg.pdf")
ggplot(data = df.barr.tmb.type[which(!is.na(df.barr.tmb.type$EverSmoked)),]) +
  geom_mosaic(aes(x = product(EverSmoked, PatientGrade), fill=EverSmoked)) + 
  labs(title='f(EverSmoked | PatientGrade) f(PatientGrade)')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))
dev.off()

### DDR-associated mutations?
load("../df.pathwayAlterationCounts.RData")
df$TumourID <- sapply(df$Sample, function(x) strsplit(x,"_vs_")[[1]][1])
df.tmb.pluspathways <- merge(df.tmb,
                             df, by.x="TumourID", by.y="TumourID",
                             all.x=FALSE, all.y = FALSE)
tb <- table(df.tmb.pluspathways[,c("EverSmoked","PathType")])


pvals <- NULL
or <- NULL
for (p in 1:ncol(tb)) {
  f <- fisher.test(cbind(tb[,p], rowSums(tb[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb)
names(or) <- colnames(tb)
pvals[pvals<0.05]
# damage in S phase:AMP                          Direct Repair:SNV 
# 0.027202299                                0.005576643 
# G2-M checkpoint:LOH Ubiquitins and Ubiquitin-like proteins:AMP 
# 0.035366409                                0.037783322
or[pvals<0.05]
# damage in S phase:AMP                          Direct Repair:SNV 
# 0.4335236                                  0.1913635 
# G2-M checkpoint:LOH Ubiquitins and Ubiquitin-like proteins:AMP 
# 0.7817303                                  4.2219827
pvaladj <- p.adjust(pvals, method="BH")
pvaladj[pvaladj<0.1]
#nothing

tb <- table(df.tmb.pluspathways[which(df.tmb.pluspathways$Type=="SNV"),c("EverSmoked","PathType")])
pvals <- NULL
or <- NULL
for (p in 1:ncol(tb)) {
  f <- fisher.test(cbind(tb[,p], rowSums(tb[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb)
names(or) <- colnames(tb)
pvals[pvals<0.05]
# Direct Repair:SNV 
# 0.001862324 
or[pvals<0.05]
# Direct Repair:SNV 
# 0.1507208
pvaladj <- p.adjust(pvals, method="BH")
pvaladj[pvaladj<0.1]
#Direct Repair:SNV 
#0.04655811 

tb <- table(df.tmb.pluspathways[which(df.tmb.pluspathways$Type=="Indel"),c("EverSmoked","PathType")])
pvals <- NULL
or <- NULL
for (p in 1:ncol(tb)) {
  f <- fisher.test(cbind(tb[,p], rowSums(tb[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb)
names(or) <- colnames(tb)
pvals[pvals<0.05]
# TLS:Indel 
#0.01433138
or[pvals<0.05]
# TLS:Indel 
#0.2040062 
pvaladj <- p.adjust(pvals, method="BH")
pvaladj[pvaladj<0.1]
#nothing

tb <- table(df.tmb.pluspathways[which(df.tmb.pluspathways$Type=="AMP"),c("EverSmoked","PathType")])
pvals <- NULL
or <- NULL
for (p in 1:ncol(tb)) {
  f <- fisher.test(cbind(tb[,p], rowSums(tb[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb)
names(or) <- colnames(tb)
pvals[pvals<0.05]
# damage in S phase:AMP Ubiquitins and Ubiquitin-like proteins:AMP 
# 0.02616225                                 0.03686646
or[pvals<0.05]
# damage in S phase:AMP Ubiquitins and Ubiquitin-like proteins:AMP 
# 0.4275548                                  4.2555786 
pvaladj <- p.adjust(pvals, method="BH")
pvaladj[pvaladj<0.1]
# nothing

tb <- table(df.tmb.pluspathways[which(df.tmb.pluspathways$Type=="DEL"),c("EverSmoked","PathType")])
pvals <- NULL
or <- NULL
for (p in 1:ncol(tb)) {
  f <- fisher.test(cbind(tb[,p], rowSums(tb[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb)
names(or) <- colnames(tb)
pvals[pvals<0.05]
# nothing

tb <- table(df.tmb.pluspathways[which(df.tmb.pluspathways$Type=="LOH"),c("EverSmoked","PathType")])
pvals <- NULL
or <- NULL
for (p in 1:ncol(tb)) {
  f <- fisher.test(cbind(tb[,p], rowSums(tb[,-p])))
  pvals <- c(pvals, f$p.value)
  or <- c(or, f$estimate)
}
names(pvals) <- colnames(tb)
names(or) <- colnames(tb)
pvals[pvals<0.05]
# nothing

pdf("plots_supp/mosaic.DDRpathways.SNV.pdf",w=10)
ggplot(data = df.tmb.pluspathways[which((!is.na(df.tmb.pluspathways$EverSmoked))&
                                          (df.tmb.pluspathways$Type=="SNV")),]) +
  geom_mosaic(aes(x = product(EverSmoked, PathType), fill=EverSmoked)) + 
  labs(title='f(EverSmoked | PathType) f(PathType)')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plots_supp/mosaic.DDRpathways.AMP.pdf",w=10)
ggplot(data = df.tmb.pluspathways[which((!is.na(df.tmb.pluspathways$EverSmoked))&
                                          (df.tmb.pluspathways$Type=="AMP")),]) +
  geom_mosaic(aes(x = product(EverSmoked, PathType), fill=EverSmoked)) + 
  labs(title='f(EverSmoked | PathType) f(PathType)')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plots_supp/mosaic.DDRpathways.DEL.pdf",w=10)
ggplot(data = df.tmb.pluspathways[which((!is.na(df.tmb.pluspathways$EverSmoked))&
                                          (df.tmb.pluspathways$Type=="DEL")),]) +
  geom_mosaic(aes(x = product(EverSmoked, PathType), fill=EverSmoked)) + 
  labs(title='f(EverSmoked | PathType) f(PathType)')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plots_supp/mosaic.DDRpathways.LOH.pdf",w=10)
ggplot(data = df.tmb.pluspathways[which((!is.na(df.tmb.pluspathways$EverSmoked))&
                                          (df.tmb.pluspathways$Type=="LOH")),]) +
  geom_mosaic(aes(x = product(EverSmoked, PathType), fill=EverSmoked)) + 
  labs(title='f(EverSmoked | PathType) f(PathType)')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("plots_supp/mosaic.DDRpathways.Indel.pdf",w=10)
ggplot(data = df.tmb.pluspathways[which((!is.na(df.tmb.pluspathways$EverSmoked))&
                                          (df.tmb.pluspathways$Type=="Indel")),]) +
  geom_mosaic(aes(x = product(EverSmoked, PathType), fill=EverSmoked)) + 
  labs(title='f(EverSmoked | PathType) f(PathType)')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

write.csv(df.tmb.pluspathways, file="muts.DDR.smoking.csv")

###################################
### Treatment:

library(wesanderson)

# Treatment data:
clin <- read.csv("../OCCAMS-export_MS-cohort_SAZ_20220711.csv",
                 header=TRUE)
clin.sigs <- merge(unique(df.oac[,c("OCCAMS_ID","Category","EverSmoked","variable","value")]), 
                     clin,by.x="OCCAMS_ID",by.y="ID", all.x=FALSE, all.y=FALSE)
clin.sigs <- clin.sigs[which(!is.na(clin.sigs$EverSmoked)),]

clin.sigs$RP.Tumour.Growth.c <- factor(clin.sigs$RP.Tumour.Growth.c,
                                       levels=c("shrink","no change","grow"))
mycomp <- list(c("Y","N"))
pdf("plots_supp/cont.RP.Tumour.Growth.c.pdf",w=18,h=12)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.Tumour.Growth.c)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(RP.Tumour.Growth.c~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$RP.TumourResponse <- factor(clin.sigs$RP.TumourResponse,
                                      levels=c("0%","<20%","=20%","<50%","=50%"))
mycomp1 <- list(c("0%","<20%"),
                c("<20%","=20%"),
                c("=20%","<50%"),
                c("<50%","=50%"))
pdf("plots_supp/cont.RP.TumourResponse.c.pdf",w=18,h=16)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.TumourResponse)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(RP.TumourResponse~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()


clin.sigs$RP.MandardScoreForResponse <- factor(clin.sigs$RP.MandardScoreForResponse)
mycomp1 <- list(c("TRG1","TRG2"),
                c("TRG2","TRG3"),
                c("TRG3","TRG4"),
                c("TRG4","TRG5"))
pdf("plots_supp/cont.RP.MandardScoreForResponse.pdf",w=20,h=16)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.MandardScoreForResponse)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(RP.MandardScoreForResponse~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.PalChemo.Response <- factor(clin.sigs$TR.PalChemo.Response,
                                         levels=c("PR","SD","SD, PD","PD"))
mycomp1 <- list(c("PR","SD"),
                c("SD","SD, PD"),
                c("SD, PD","PD"))
pdf("plots_supp/cont.TR.PalChemo.Response.pdf",w=18,h=16)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.PalChemo.Response)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp, label="p.signif")+
  facet_grid(TR.PalChemo.Response~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.Adj.Response <- factor(clin.sigs$TR.Adj.Response,
                                    levels=c("CR","PR","PD"))
mycomp1 <- list(c("CR","PR"),
                c("PR","PD"),
                c("CR","PD"))
pdf("plots_supp/cont.TR.Adj.Response.pdf",w=20,h=16)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.Adj.Response)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp, label="p.signif")+
  facet_grid(TR.Adj.Response~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.RT.Response <- factor(clin.sigs$TR.RT.Response,
                                   levels=c("CR","PR","SD","SD, PD","PD"))
mycomp1 <- list(c("CR","PR"),
                c("PR","SD"),
                c("SD","SD, PD"),
                c("SD, PD","PD"))
pdf("plots_supp/cont.TR.RT.Response.pdf",w=21,h=16)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp, label="p.signif")+
  facet_grid(TR.RT.Response~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.RT.Response.RvsD <- sapply(clin.sigs$TR.RT.Response,
                                        function(x) ifelse (is.na(x),x,ifelse(x %in% c("CR","PR"),"Responder","Non-responder")))
mycomp1 <- list(c("Responder","Non-responder"))
pdf("plots_supp/cont.TR.RT.RespondervsNonResponder.2.pdf",w=18,h=10)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response.RvsD)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.RT.Response.RvsD~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()


pdf("plots_supp/cont.TR.RT.RespondervsNonResponder.2.Primaries.pdf",w=18,h=10)
ggboxplot(clin.sigs[which((!is.na(clin.sigs$TR.RT.Response.RvsD))&
                      (clin.sigs$Category == "PrimaryTumour")),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.RT.Response.RvsD~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

write.csv(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response.RvsD)),],
          file="smoker.sigs.respondernonresponder.csv")

clin.sigs$TR.RT.Response.CRvsrest <- sapply(clin.sigs$TR.RT.Response,
                                            function(x) ifelse (x %in% c("CR"),"Complete responder","Non-responder"))
mycomp1 <- list(c("Complete responder","Non-responder"))
pdf("plots_supp/cont.TR.RT.Response.CRvsrest.pdf",w=18,h=10)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response.CRvsrest)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.RT.Response.CRvsrest~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.NeoAdj.Response <- factor(clin.sigs$TR.NeoAdj.Response,
                                       levels=c("CR","PR","SD","PD"))
mycomp1 <- list(c("CR","PR"),
                c("PR","SD"),
                c("SD","PD"))
pdf("plots_supp/cont.TR.NeoAdj.Response.pdf",w=18,h=12)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp, label="p.signif")+
  facet_grid(TR.NeoAdj.Response~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

 clin.sigs$TR.NeoAdj.Response.StableVsProg <- sapply(clin.sigs$TR.NeoAdj.Response,
                                                    function(x) ifelse(is.na(x),x,ifelse(x=="PD","PD","Response/Stable")))
mycomp1 <- list(c("PD","Response/Stable"))
pdf("plots_supp/cont.TR.NeoAdj.Response.StableVsProg.2.pdf",w=18,h=10)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.NeoAdj.Response.StableVsProg~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.NeoAdj.Response.StableVsProg2 <- sapply(clin.sigs$TR.NeoAdj.Response,
                                                    function(x) ifelse(is.na(x),x,ifelse(x%in% c("PD","SD"),"Non-responder","Responder")))
mycomp1 <- list(c("Non-responder","Responder"))
pdf("plots_supp/cont.TR.NeoAdj.Response.StableVsProg.3.pdf",w=18,h=10)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg2)),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.NeoAdj.Response.StableVsProg2~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.NeoAdj.Response.StableVsProg.3.Barretts.pdf",w=18,h=10)
ggboxplot(clin.sigs[which((!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg2))&
                      (clin.sigs$Category == "Barretts")),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.NeoAdj.Response.StableVsProg2~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.NeoAdj.Response.StableVsProg.3.Primaries.pdf",w=18,h=10)
ggboxplot(clin.sigs[which((!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg2))&
                            (clin.sigs$Category == "PrimaryTumour")),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.NeoAdj.Response.StableVsProg2~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs.naive <- clin.sigs[which(clin.sigs$OCCAMS_ID %in% naive$OCCAMS_ID),]
clin.sigs.treated <- clin.sigs[which(!(clin.sigs$OCCAMS_ID %in% naive$OCCAMS_ID)),]

        
pdf("plots_supp/cont.TR.NeoAdj.Response.StableVsProg.3.naive.pdf",w=18,h=10)
ggboxplot(clin.sigs.naive[which((!is.na(clin.sigs.naive$TR.NeoAdj.Response.StableVsProg2))),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.NeoAdj.Response.StableVsProg2~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.NeoAdj.Response.StableVsProg.3.treated.pdf",w=18,h=10)
ggboxplot(clin.sigs.treated[which((!is.na(clin.sigs.treated$TR.NeoAdj.Response.StableVsProg2))),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.NeoAdj.Response.StableVsProg2~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.RT.Response.RvsD.3.naive.pdf",w=18,h=10)
ggboxplot(clin.sigs.naive[which((!is.na(clin.sigs.naive$TR.RT.Response.RvsD))),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.RT.Response.RvsD~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.RT.Response.RvsD.3.treated.pdf",w=18,h=10)
ggboxplot(clin.sigs.treated[which((!is.na(clin.sigs.treated$TR.RT.Response.RvsD))),], 
          x = "EverSmoked", y = "value",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(TR.RT.Response.RvsD~variable,scale="free")+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

write.csv(clin.sigs.naive, file="clin.sigs.naive.csv")
write.csv(clin.sigs.treated, file="clin.sigs.treated.csv")


selected <- clin.sigs[which(clin.sigs$variable == "SBS17a"),]
df.melt <- melt(unique(selected[,c("OCCAMS_ID","value",
                            "TR.RT.Response.RvsD",
                            "TR.NeoAdj.Response.StableVsProg2",
                            "EverSmoked")]),
                id.vars=c("OCCAMS_ID","EverSmoked","value"))
colnames(df.melt)[3] <- "Exposure"

pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17a.pdf",h=5)
ggboxplot(df.melt[((!is.na(df.melt$value))&
                    (df.melt$value=="Non-responder")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(value~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

### select only chemo-treated patients:

df.melt2 <- melt(unique(clin.sigs[,c("OCCAMS_ID","value",
                                   "TR.RT.Response.RvsD",
                                   "TR.NeoAdj.Response.StableVsProg2",
                                   "EverSmoked","variable","value")]),
                id.vars=c("OCCAMS_ID","EverSmoked","value","variable"))
colnames(df.melt2)[c(3,4)] <- c("Exposure","Signature")

naive <- read.csv("plots_1/naive.primmet.csv")
df.melt.naive <- df.melt[which(df.melt$OCCAMS_ID %in% naive$OCCAMS_ID),]
df.melt.treated <- df.melt[which(!(df.melt$OCCAMS_ID %in% naive$OCCAMS_ID)),]

df.melt2.naive <- df.melt2[which(df.melt2$OCCAMS_ID %in% naive$OCCAMS_ID),]
df.melt2.treated <- df.melt2[which(!(df.melt2$OCCAMS_ID %in% naive$OCCAMS_ID)),]


pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17a.naive.pdf",h=5)
ggboxplot(df.melt.naive[((!is.na(df.melt.naive$value))&
                     (df.melt.naive$value=="Non-responder")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(value~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

write.csv(df.melt.naive[((!is.na(df.melt.naive$value))&
                           (df.melt.naive$value=="Non-responder")),], file="plots_supp/neoadj_chemo.naive.csv")
write.csv(df.melt.treated[((!is.na(df.melt.treated$value))&
                           (df.melt.treated$value=="Non-responder")),], file="plots_supp/neoadj_chemo.treated.csv")



pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17a.treated.pdf",h=5)
ggboxplot(df.melt.treated[((!is.na(df.melt.treated$value))&
                     (df.melt.treated$value=="Non-responder")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(value~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.RT.allsigs.naive.pdf",h=5)
ggboxplot(df.melt2.naive[((!is.na(df.melt2.naive$value))&
                           (df.melt2.naive$variable=="TR.RT.Response.RvsD")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_wrap(Signature~Exposure+value,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()


selected <- clin.sigs[which(clin.sigs$variable == "SBS17b"),]
df.melt <- melt(unique(selected[,c("OCCAMS_ID","value",
                                   "TR.RT.Response.RvsD",
                                   "TR.NeoAdj.Response.StableVsProg2",
                                   "EverSmoked")]),
                id.vars=c("OCCAMS_ID","EverSmoked","value"))
colnames(df.melt)[3] <- "Exposure"

pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17b.pdf",h=5)
ggboxplot(df.melt[((!is.na(df.melt$value))&
                     (df.melt$value=="Non-responder")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(value~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

### select only chemo-treated patients:

naive <- read.csv("plots_1/naive.primmet.csv")
df.melt.naive <- df.melt[which(df.melt$OCCAMS_ID %in% naive$OCCAMS_ID),]
df.melt.treated <- df.melt[which(!(df.melt$OCCAMS_ID %in% naive$OCCAMS_ID)),]

pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17b.naive.pdf",h=5)
ggboxplot(df.melt.naive[((!is.na(df.melt.naive$value))&
                           (df.melt.naive$value=="Non-responder")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(value~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17b.treated.pdf",h=5)
ggboxplot(df.melt.treated[((!is.na(df.melt.treated$value))&
                             (df.melt.treated$value=="Non-responder")),],
          x = "EverSmoked", y = "Exposure",
          color = "EverSmoked",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)+
  facet_grid(value~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

mycomp5 <- list(c("Responder","Non-responder"))
pdf("plots_supp/cont.TR.NeoAdj.RT.Combined.SBS17a.respvsnonresp.pdf")
ggboxplot(df.melt[((!is.na(df.melt$value))),],
          x = "value", y = "Exposure",
          color = "value",palette = c("#E69F00", "#56B4E9"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp5)+
  facet_grid(EverSmoked~variable,scale="free")+
  xlab("")+
  ylab("SBS17a exposure (%)")
dev.off()

write.csv(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg2)),],
          file="smoker.sigs.progressivestable.csv")

