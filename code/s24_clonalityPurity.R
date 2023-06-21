#####################
### Test relationship between clonality and purity.

library(ggpubr)

load("data/sigs.timed.annot.RData")
load("data/annotation.sampleIDs.RData")

### Barrett's samples with subclonality:
withsub <- unique(sigs.timed.annot[which((sigs.timed.annot$Category == "Barretts")&
                   (sigs.timed.annot$Clonality == "Subclonal")),]$TumourID)

withsub.prim <- unique(sigs.timed.annot[which((sigs.timed.annot$Category == "PrimaryTumour")&
                                           (sigs.timed.annot$Clonality == "Subclonal")),]$TumourID)
#144
nosub <- setdiff(annotation.sampleIDs[which(annotation.sampleIDs$Category=="Barretts"),]$TumourID,
        unique(sigs.timed.annot[which((sigs.timed.annot$Category == "Barretts")&
                                        (sigs.timed.annot$Clonality == "Subclonal")),]$TumourID))
#17 samples

load("../../tracksig_barrettsmetsprimaries/barretts.purities.RData")
load("../../tracksig_barrettsmetsprimaries/primaries.purities.RData")
load("../../tracksig_barrettsmetsprimaries/extra.purities.RData")
barretts.purities.cohort1 <- barretts.purities[which(barretts.purities$samplename %in% annotation.sampleIDs[which(annotation.sampleIDs$Category == "Barretts"),]$Sample),]
barretts.purities.cohort2 <- extra.purities[which(extra.purities$samplename %in% annotation.sampleIDs[which(annotation.sampleIDs$Category == "Barretts"),]$Sample),]
barretts.purities.cohort <- rbind(barretts.purities.cohort1, barretts.purities.cohort2)

barretts.purities.cohort$TumourID <- sapply(barretts.purities.cohort$samplename,
                                            function(x) strsplit(as.character(x),"_vs_")[[1]][1])
barretts.purities.cohort$Subclonality <- sapply(barretts.purities.cohort$TumourID,
                                                function(x) ifelse(x %in% withsub, "yes","no"))
barretts.purities.cohort$purity <- as.numeric(as.character(barretts.purities.cohort$purity))

my_comparisons <- list(c("yes","no"))
pdf("plots.clonalitypurity/purityBarretts.withVSwithoutSubclonality.pdf",w=10,h=8)
ggboxplot(barretts.purities.cohort, x = "Subclonality", y = "purity",
          color = "Subclonality", 
          add = "jitter")+
  scale_color_manual(values=c("#BFAE48", "#3A0842"))+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Subclonality")+
  ylab("Purity")
dev.off()

primaries.purities.cohort1 <- primaries.purities[which(primaries.purities$samplename %in% annotation.sampleIDs[which(annotation.sampleIDs$Category == "PrimaryTumour"),]$Sample),]
primaries.purities.cohort2 <- extra.purities[which(extra.purities$samplename %in% annotation.sampleIDs[which(annotation.sampleIDs$Category == "PrimaryTumour"),]$Sample),]
primaries.purities.cohort <- rbind(primaries.purities.cohort1, primaries.purities.cohort2)

primaries.purities.cohort$TumourID <- sapply(primaries.purities.cohort$samplename,
                                            function(x) strsplit(as.character(x),"_vs_")[[1]][1])
primaries.purities.cohort$Subclonality <- sapply(primaries.purities.cohort$TumourID,
                                                function(x) ifelse(x %in% withsub.prim, "yes","no"))
primaries.purities.cohort$purity <- as.numeric(as.character(primaries.purities.cohort$purity))


pdf("plots.clonalitypurity/purityPrimaries.withVSwithoutSubclonality.pdf",w=10,h=8)
ggboxplot(primaries.purities.cohort, x = "Subclonality", y = "purity",
          color = "Subclonality", 
          add = "jitter")+
  scale_color_manual(values=c("#BFAE48", "#3A0842"))+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Subclonality")+
  ylab("Purity")
dev.off()

primaries.purities.cohort$Category <- "PrimaryTumour"
barretts.purities.cohort$Category <- "BarrettOesophagus"
together <- rbind(primaries.purities.cohort, barretts.purities.cohort)

pdf("plots.clonalitypurity/purityBarrettsPrimaries.withVSwithoutSubclonality.pdf",w=10,h=8)
ggboxplot(together, x = "Subclonality", y = "purity",
          color = "Category", 
          add = "jitter")+
  scale_color_manual(values=c("#BFAE48", "#3A0842"))+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Subclonality")+
  ylab("Purity")
dev.off()

my_comparisons3 <- list(c("PrimaryTumour","BarrettOesophagus"))
pdf("plots_6/purityBarrettsPrimaries.withVSwithoutSubclonality.2panels.pdf")
ggboxplot(together, x = "Category", y = "purity",
          color = "Category", 
          add = "jitter")+
  scale_color_manual(values=c("#BFAE48", "#3A0842"))+
  stat_compare_means(comparisons = my_comparisons3)+
  facet_wrap(~Subclonality)+
  xlab("")+
  ylab("Sample purity")
dev.off()

write.csv(together, file="purity.subclonality.together.csv")

together$Categ.Subclon <- paste(together$Category, together$Subclonality, sep=":")

my_comparisons2 <- list(c("BarrettOesophagus:yes","BarrettOesophagus:no"),
                        c("BarrettOesophagus:no","PrimaryTumour:no"),
                        c("BarrettOesophagus:yes","PrimaryTumour:yes"),
                        c("PrimaryTumour:yes","PrimaryTumour:no"))

pdf("plots.clonalitypurity/purityTogether.withVSwithoutSubclonality.pdf",w=10,h=8)
ggboxplot(together, x = "Categ.Subclon", y = "purity",
          color = "Categ.Subclon", 
          add = "jitter")+
  #scale_color_manual(values=c("#BFAE48", "#3A0842"))+
  stat_compare_means(comparisons = my_comparisons2)+
  xlab("")+
  ylab("Purity")
dev.off()

### What if I remove all samples below 30%?
together.keep <- together[which(together$purity>=0.30),]

pdf("plots.clonalitypurity/purityTogetherAbove30.withVSwithoutSubclonality.pdf",w=10,h=8)
ggboxplot(together.keep, x = "Categ.Subclon", y = "purity",
          color = "Categ.Subclon", 
          add = "jitter")+
  #scale_color_manual(values=c("#BFAE48", "#3A0842"))+
  stat_compare_means(comparisons = my_comparisons2)+
  xlab("")+
  ylab("Purity")
dev.off()

save(together, file="purities.Barretts.Primaries.RData")
