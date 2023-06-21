#############################
##### Analyse clonality output.

library(reshape)
library(ggpubr)
library(pheatmap)

# Load timed mutational signature data:
load("data/sigs.timed.annot.RData")
df.sigs.timed <- melt(sigs.timed.annot)
df.sigs.timed[which(df.sigs.timed$Category == "LymphNode"),]$Category <- "Metastasis"
colnames(df.sigs.timed)[6:7] <- c("Signature","Exposure")
df.sigs.timed$Category <- factor(df.sigs.timed$Category,
                                 levels=c("Barretts","PrimaryTumour",
                                          "Metastasis"))

### all in one pdf:
df.sigs.timed$Clonality <- factor(df.sigs.timed$Clonality,
                                 levels=c("Clonal","Subclonal"))
my_comparisons <- list(c("Clonal","Subclonal"))
pdf("plots.s17/clonalVsSubclonal.primbarmet.individualsigs.pdf",h=6,w=10)
for (s in unique(df.sigs.timed$Signature)) {
  print(ggviolin(df.sigs.timed[which(df.sigs.timed$Signature == s),], 
           x = "Clonality", y = "Exposure",
           fill = "Clonality", palette =c("#9F4A54","#D7BBA8"),
           add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
    #stat_compare_means(label.y = 1) +
    facet_grid(vars(Signature),vars(Category),scale="free"))+
    labs(s)
}
dev.off()


# Next, compare changes in tumour compared to respective adjacent Barrett's:
load("../pairs.RData")
df.sigs.timed.paired <- merge(df.sigs.timed,
                              unique(pairs[,c("Sample","OCCAMS_ID")]),
                              by.x="Sample",by.y="Sample",
                              all.x=FALSE,all.y=FALSE)

my_comparisonsbp <- list(c("Barretts","PrimaryTumour"))
pdf("plots.s17/paired.BarrettsPrimaries.pdf",w=12,h=8)
ggpaired(df.sigs.timed.paired[which(df.sigs.timed.paired$Category %in% 
                                      c("Barretts","PrimaryTumour")),], 
         x = "Category", y = "Exposure",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_grid(Clonality~Signature, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisonsbp, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

my_comparisonscl <- list(c("Clonal","Subclonal"))
pdf("plots.s17/allcateg.clonalvssubclonal.pdf",w=12,h=8)
ggviolin(df.sigs.timed, 
         x = "Clonality", y = "Exposure",
         palette ="npg",
         fill = "Clonality",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_grid(Category~Signature, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisonscl, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

####################
### Changes at subclones?

changes <- array(0,c(length(unique(sigs.timed.annot$Sample)),14))
rownames(changes) <- unique(sigs.timed.annot$Sample)
colnames(changes) <- colnames(sigs.timed.annot)[2:15]
for (sample in unique(sigs.timed.annot$Sample)) {
  for (s in colnames(sigs.timed.annot)[2:15]) {
    if ((length(which((sigs.timed.annot$Clonality == "Clonal")&
                     (sigs.timed.annot$Sample==sample)))>0)&
        (length(which((sigs.timed.annot$Clonality == "Subclonal")&
                      (sigs.timed.annot$Sample==sample)))>0)) {
    changes[sample,s] <- sum(sigs.timed.annot[which((sigs.timed.annot$Clonality == "Clonal")&
                                            (sigs.timed.annot$Sample==sample)),s])-
    sigs.timed.annot[which((sigs.timed.annot$Clonality == "Subclonal")&
                       (sigs.timed.annot$Sample==sample)),s]
    }
  }
}
df.changes <- melt(changes)
colnames(df.changes) <- c("Sample","Signature","Change")
df.changes$Change <- (-1)*df.changes$Change
save(changes,file="changes.RData")

df.changes.plusannot <- merge(df.changes,
                              annotation.sampleIDs[,c("Sample","Category")],
                              by.x="Sample",by.y="Sample",
                              all.x=FALSE, all.y=FALSE)
save(df.changes.plusannot, file="df.changes.plusannot.RData")

pdf("plots.s17/bottleneck.change.pdf",w=10)
ggboxplot(df.changes, x = "Signature", y = "Change",
          color = "Signature")+
  geom_hline(yintercept=0, linetype="dotted", 
             color = "black", size=0.7)#+
  #facet_wrap(~Category)
dev.off()

pdf("plots.s17/bottleneck.change.bycateg.pdf",w=10)
ggboxplot(df.changes.plusannot, x = "Signature", y = "Change",
          color = "Signature")+
  geom_hline(yintercept=0, linetype="dotted", 
             color = "black", size=0.7)+
facet_wrap(~Category)
dev.off()


library(pheatmap)
pdf("plots.s17/heatmap.bottlenecks.pdf",w=8,h=4)
pheatmap(t(changes),
         show_colnames = FALSE,
         #annotation_col =annotc,
         cutree_cols = 3)
dev.off()
# the heatmap is not clear




####### Bottlenecks:

bottles <- NULL
for (sample in unique(oac.sigClonality.annot$Sample)) {
  print(sample)
  current.oac <- oac.sigClonality.annot[which(oac.sigClonality.annot$Sample == sample),]
  if (nrow(current.oac)>1) {
    sortedccf <- rev(sort(current.oac$CCF))
    changes <- current.oac[which(current.oac$CCF==sortedccf[2]),2:14]-
      current.oac[which(current.oac$CCF==sortedccf[1]),2:14]
    bottles <- rbind(bottles,cbind(changes, sample,current.oac$Category))
  }
}
df.bottles <- melt(bottles, id.vars = c("sample","current.oac$Category"))
colnames(df.bottles) <- c("Sample","Category","Signature","Change")

library(pheatmap)
bottles <- unique(bottles)
pdf("plots/heatmap.bottlenecks.pdf",w=8,h=4)
annotc <- data.frame(XX=0,Category=bottles$`current.oac$Category`)
annotc$Category <- factor(annotc$Category)
rownames(annotc) <- bottles$sample
matx <- bottles[,c(1:13)]
rownames(matx) <- bottles$sample
pheatmap(t(matx),
         show_colnames = FALSE,
         annotation_col =annotc,
         cutree_cols = 4)
dev.off()

df.bottles$Category <- factor(df.bottles$Category,
                              levels=c("Barretts",
                                       "PrimaryTumour",
                                       "LymphNode","Metastasis"))
df.bottles$Signature <- factor(df.bottles$Signature,
                     levels=c("SBS17a",
                              "SBS17b",
                              "SBS2",
                              "SBS3",
                              "SBS30",
                              "SBS18",
                              "SBS44",
                              "SBS28",
                              "SBS35",
                              "SBS1",
                              "SBS5",
                              "SBS40",
                              "SBS41"))
save(df.bottles, file="df.bottles.RData")
save(bottles, file="bottles.RData")

# Changes during bottleneck:
pdf("plots/bottleneck.change.pdf",w=20)
ggboxplot(df.bottles, x = "Signature", y = "Change",
          color = "Signature")+
  geom_hline(yintercept=0, linetype="dotted", 
             color = "black", size=0.7)+
  facet_wrap(~Category)
dev.off()



#### Chemo naive vs treated:
clindat <- read.csv("../SigProfiler/fullset_occams_SA_SAZ_28052020.csv")

oac.sigClonality.clin <- merge(oac.sigClonality.annot,
                               clindat[,c("Illumina.ID",
                                          "EX.RefluxFrequencyOfSymptoms", 
                                          "EX.DurationOfRefluxSymptoms",
                                          "TR.NeoAdj.ChemotherapyNumberOfCyclesPlanned")],
                               by.x="TumourID",by.y="Illumina.ID",
                               all.x=FALSE, all.y=FALSE)
oac.sigClonality.clin$NeoAdjChemo <- sapply(oac.sigClonality.clin$TR.NeoAdj.ChemotherapyNumberOfCyclesPlanned,
                                            function(x) ifelse(is.na(x),"Naive","Chemo-treated"))

#163 naive vs 240 treated
sigclon.clin.melt <- melt(oac.sigClonality.clin,
                          id.vars = c("Sample","TumourID","CCF",
                                      "Cohort","Category",
                                      "Clones","EX.RefluxFrequencyOfSymptoms",
                                      "EX.DurationOfRefluxSymptoms",
                                      "TR.NeoAdj.ChemotherapyNumberOfCyclesPlanned",
                                      "NeoAdjChemo"))
colnames(sigclon.clin.melt)[11:12]<-c("Signature","Exposure")
sigclon.clin.melt$Category <- factor(sigclon.clin.melt$Category,
                              levels=c("Barretts",
                                       "PrimaryTumour",
                                       "LymphNode","Metastasis"))
sigclon.clin.melt$Signature <- factor(sigclon.clin.melt$Signature,
                               levels=c("SBS17a",
                                        "SBS17b",
                                        "SBS2",
                                        "SBS3",
                                        "SBS30",
                                        "SBS18",
                                        "SBS44",
                                        "SBS28",
                                        "SBS35",
                                        "SBS1",
                                        "SBS5",
                                        "SBS40",
                                        "SBS41"))
sigclon.clin.melt$Clonality <- sapply(sigclon.clin.melt$CCF,
                                 function(x) ifelse(x==1,"Clonal","Subclonal"))
sigclon.clin.melt$Clonality <- factor(sigclon.clin.melt$Clonality,
                                     levels=c("Clonal",
                                              "Subclonal"))
sigclon.clin.melt$NeoAdjChemo <- factor(sigclon.clin.melt$NeoAdjChemo,
                                      levels=c("Naive",
                                               "Chemo-treated"))

pdf("plots/clonalVsSubclonal.chemoNaiveVStreated.pdf",w=6,h=32)
ggviolin(sigclon.clin.melt, 
         x = "NeoAdjChemo", y = "Exposure",
         fill = "NeoAdjChemo", palette =c("#747C92", "#C36F09"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("Naive","Chemo-treated")),label = "p.signif")+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 1) +
  facet_grid(vars(Signature),vars(Clonality),scale="free")
dev.off()


pdf("plots/clonalVsSubclonal.chemoNaiveVStreated.otherway.pdf",w=6,h=32)
ggviolin(sigclon.clin.melt, 
         x = "Clonality", y = "Exposure",
         fill = "Clonality", palette =c("#9F4A54","#D7BBA8"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2,label = "p.signif")+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 1) +
  facet_grid(vars(Signature),vars(NeoAdjChemo),scale="free")
dev.off()

pdf("plots/clonalVsSubclonal.chemoNaiveVStreated.otherway.selected.pdf",w=6,h=18)
ggviolin(sigclon.clin.melt[which(sigclon.clin.melt$Signature
                                 %in% c("SBS17a","SBS17b",
                                        "SBS2","SBS3","SBS30",
                                        "SBS35","SBS41")),], 
         x = "Clonality", y = "Exposure",
         fill = "Clonality", palette =c("#9F4A54","#D7BBA8"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2,label = "p.signif")+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 1) +
  facet_grid(vars(Signature),vars(NeoAdjChemo),scale="free")
dev.off()


#### Finally, associate with expression:


## Catalogue cases with increase at bottleneck vs decrease at bottleneck:
load("../df.changes.plusannot.RData")

df.changes.plusannot$ChangeType <- sapply(df.changes.plusannot$Change,
                                          function(x) ifelse(x<0,"Decrease",ifelse(x>0,"Increase","No_change")))
df.changes.plusannot$TumourID <- sapply(df.changes.plusannot$Sample,
                                        function(x) strsplit(as.character(x),"_vs_")[[1]][1])
df.changes.plusannot$ChangeType <- factor(df.changes.plusannot$ChangeType,
                                          levels =c("Decrease","No_change","Increase"))

table(df.changes.plusannot[,c("Signature","ChangeType")])

# Next, load expression data and map against genomic changes:
load("../../tracksig_barrettsmetsprimaries/oac.expr.final.RData")
load("../../tracksig_barrettsmetsprimaries/mat.expr.RData")

cbioportal_invasion <- read.table("cbioportal_invasionmetastasis.txt",
                                  header=FALSE, stringsAsFactors = FALSE)$V1

rownames(mat.expr) <- mat.expr$DNAid 
mat.expr <- mat.expr[,-1]

library(ConsensusTME)

mat.transf <- t(as.matrix(mat.expr))
rownames(mat.transf) <- colnames(mat.expr)
colnames(mat.transf) <- rownames(mat.expr)

celldeconv <- consensusTMEAnalysis(as.matrix(mat.transf), 
                                                 cancer = "ESCA", 
                                   statMethod = "gsva")
df.celldeconv <- data.frame(t(celldeconv))
df.celldeconv$Sample <- rownames(df.celldeconv)

df.changes.plusdeconv <- merge(df.changes.plusannot,
                               df.celldeconv, 
                               by.x="TumourID", by.y="Sample",
                               all.x=FALSE, all.y=FALSE)
df.changes.plusdeconv.s17a <- df.changes.plusdeconv[which(df.changes.plusdeconv$Signature == "SBS17a"),]
df.changes.plusdeconv.s17b <- df.changes.plusdeconv[which(df.changes.plusdeconv$Signature == "SBS17b"),]

df.melt.s17a <- melt(df.changes.plusdeconv.s17a,
                     id.vars = c("TumourID", "Sample", "Signature", 
                                 "Category", "Change","ChangeType"))
df.melt.s17b <- melt(df.changes.plusdeconv.s17b,
                     id.vars = c("TumourID", "Sample", "Signature", 
                                 "Category", "Change","ChangeType"))

### Plot immune infiltrates compared:
my_comparisons <- list(c("Decrease","No_change"),
                       c("No_change","Increase"),
                       c("Decrease","Increase"))
pdf("plots.s17/s17a.TMEcompared.ssgsea.pdf",w=10,h=8)
ggviolin(df.melt.s17a, 
         x = "ChangeType", y = "value",
         palette ="npg",
         fill = "ChangeType",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Infiltration score")
dev.off()

pdf("plots.s17/s17b.TMEcompared.ssgsea.pdf",w=10,h=8)
ggviolin(df.melt.s17b, 
         x = "ChangeType", y = "value",
         palette ="npg",
         fill = "ChangeType",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Infiltration score")
dev.off()

df.melt.s17together <- merge(df.melt.s17a,
                             df.melt.s17b,
                             by.x=c("TumourID","Sample","Category",
                                    "variable","value"),
                             by.y=c("TumourID","Sample","Category",
                                    "variable","value"),
                             all.x=FALSE, all.y=FALSE)
df.melt.s17together$ChangeType <- paste0(df.melt.s17together$ChangeType.x,".",
                                        df.melt.s17together$ChangeType.y)
df.melt.s17together[which(df.melt.s17together$ChangeType %in% 
                            c("Decrease.Increase","Increase.Decrease")),]$ChangeType <- "Opposite"

mycomp2 <- list(c("Decrease.Decrease","Increase.Increase"),
                c("Decrease.Decrease","No_change.No_change"),
                c("Increase.Increase","No_change.No_change"),
                c("No_change.No_change","Opposite"),
                c("Decrease.Decrease","Opposite"),
                c("Increase.Increase","Opposite"))

pdf("plots_supp/s17abtogether.TMEcompared.gsva.pdf",w=10,h=14)
ggviolin(df.melt.s17together, 
         x = "ChangeType", y = "value",
         palette ="npg",
         fill = "ChangeType",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("Infiltration score")
dev.off()

write.csv(df.melt.s17together,file="plots_supp/df.melt.s17together.csv")

df.melt <- melt(df.changes.plusdeconv,
                     id.vars = c("TumourID", "Sample", "Signature", 
                                 "Category", "Change","ChangeType"))



pdf("plots.s17/allSigs.TMEcompared.gsva.pdf",w=20,h=20)
ggviolin(df.melt, 
         x = "ChangeType", y = "value",
         palette ="npg",
         fill = "ChangeType",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_grid(Signature~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Infiltration score")
dev.off()

## No significant immune changes observed; most of them in SBS40

### Now try hallmarks of cancer from CancerSEA:
emt <- read.table("EMT.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName
hypoxia <- read.table("Hypoxia.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName
metastasis <- read.table("Metastasis.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName
inflammation <- read.table("Inflammation.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName
invasion <- read.table("Invasion.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName
angiogenesis <- read.table("Angiogenesis.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName
stemness <- read.table("Stemness.txt", header=TRUE, stringsAsFactors = FALSE)$GeneName

hallmarks <- "NULL"
hallmarks <- list(metastasis, invasion, emt, inflammation,cbioportal_invasion,
                  hypoxia, angiogenesis, stemness)
names(hallmarks) <- c("metastasis", "invasion", "emt", 
                      "inflammation","cbio_invasion",
                      "hypoxia", "angiogenesis", "stemness")

library(GSVA)
hallmarkscores <- gsva(as.matrix(mat.transf),hallmarks, 
     method="gsva", kcdf="Gaussian",
     mx.diff=TRUE, abs.ranking=FALSE)
df.hallmark <- data.frame(t(hallmarkscores))
df.hallmark$Sample <- rownames(df.hallmark)

df.changes.plushallmark <- merge(df.changes.plusannot,
                                 df.hallmark, 
                               by.x="TumourID", by.y="Sample",
                               all.x=FALSE, all.y=FALSE)

df.changes.plushallmark.s17a <- df.changes.plushallmark[which(df.changes.plushallmark$Signature == "SBS17a"),]
df.changes.plushallmark.s17b <- df.changes.plushallmark[which(df.changes.plushallmark$Signature == "SBS17b"),]

df.melt.s17a <- melt(df.changes.plushallmark.s17a,
                     id.vars = c("TumourID", "Sample", "Signature", 
                                 "Category", "Change","ChangeType"))
df.melt.s17b <- melt(df.changes.plushallmark.s17b,
                     id.vars = c("TumourID", "Sample", "Signature", 
                                 "Category", "Change","ChangeType"))

### Plot hallmarks compared:
pdf("plots.s17/s17a.HallmarksCompared.gsva.pdf",w=10,h=8)
ggviolin(df.melt.s17a, 
         x = "ChangeType", y = "value",
         palette ="npg",
         fill = "ChangeType",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Infiltration score")
dev.off()

pdf("plots.s17/s17b.HallmarksCompared.gsva.pdf",w=10,h=8)
ggviolin(df.melt.s17b, 
         x = "ChangeType", y = "value",
         palette ="npg",
         fill = "ChangeType",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Infiltration score")
dev.off()


#### Make a heat map of only primary tumours:
mat.sigchanges <- cast(Signature~Sample, value="Change",
                       fun.aggregate="mean", 
                       data=df.changes.plushallmark[which(df.changes.plushallmark$ChangeType!="No_change"),])
df.mat.sigchanges <- as.matrix(mat.sigchanges[,-1])
colnames(df.mat.sigchanges) <- colnames(mat.sigchanges[,-1])
rownames(df.mat.sigchanges) <- mat.sigchanges$Signature

### remove least variable signatures:
sortedvars <- sort(apply(df.mat.sigchanges, 1, function(x) var(x,na.rm = TRUE)))
keep <- setdiff(names(sortedvars[sortedvars>0.01]),"SBS3")


pdf("plots.s17/heatmap.SignatureChanges.onlyPrimary.pdf",h=4)
pheatmap(df.mat.sigchanges[keep,], show_colnames = FALSE)
dev.off()