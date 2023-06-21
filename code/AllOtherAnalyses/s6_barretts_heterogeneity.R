##########
## Check mut sig and CN variation by BO type.

library(ggpubr)
library(wesanderson)

load("data/annotation.sampleIDs.RData")
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

sigs.allSamples$Sample <- rownames(sigs.allSamples)
sigs.allSamples$TumourID <- sapply(rownames(sigs.allSamples),
                                   function(x) strsplit(x,"_vs_")[[1]][1])
sigs <- merge(sigs.allSamples, annotation.sampleIDs[,c("Sample","Category","OCCAMS_ID")],
              all.x = FALSE, all.y=FALSE,
              by.x="Sample",by.y="Sample")
sigs$OCCAMS_ID <- str_replace(sigs$OCCAMS_ID,"OC/","OCCAMS/")

### Read in the clinical data:
barr <- read.csv("data/Barrett's Clinical Data for Maria.csv")

### Analysis for B0:
sigs.barr <- merge(sigs, barr[,c("ID","PatientGrade")],
                   by.x="TumourID", by.y="ID",
                   all.x=FALSE, all.y=FALSE)
sigs.oac.barr <- sigs.barr[which(sigs.barr$Category == "Barretts"),]

df.barr <- melt(sigs.oac.barr,id.vars = c("OCCAMS_ID","Sample","TumourID",
                                          "Category","PatientGrade"))

### Check mut sig prevalence by BO type:

df.barr$PatientGrade <- factor(df.barr$PatientGrade,
                                        levels=c("NDBO_NP","NDBO_PP",
                                                 "LGD","HGD","IMC","trios"))

pdf("plots_supp/compareSigsByBOtype.pdf",w=20,h=8)
ggboxplot(df.barr, 
          x = "PatientGrade", y = "value",
          color = "PatientGrade",palette = wes_palette("IsleofDogs1", n =6),
          add = "jitter")+
  stat_compare_means(label.x=2)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

write.csv(df.barr, file="plots_supp/df.barr.het.csv")


## Next, load CN data and compare against primary tumour CN:
load("../cns.ploid.RData")
cns.ploid$ID <- sapply(cns.ploid$sample, function(x) strsplit(as.character(x),"_vs_")[[1]][1])

cns.mean <- cns.ploid %>%
  group_by(sample, ID) %>%
  dplyr::summarize(MeanCN = mean(RelativeCN, na.rm=TRUE),
                   MeanTotalCN= mean(total_cn, na.rm=TRUE),
                   MeanMinorCN= mean(minor_cn, na.rm=TRUE),
                   MeanPloidy = mean(ploidy, na.rm=TRUE))
cns.mean <- data.frame(cns.mean)

# Merge with sample info:
cns.mean.barr <- merge(cns.mean, barr[,c("ID","PatientGrade")],
                       by.x="ID", by.y="ID",
                       all.x=FALSE, all.y=FALSE)
cns.mean.barr <- merge(cns.mean, barr[,c("ID","PatientGrade")],
                       by.x="ID", by.y="ID",
                       all.x=FALSE, all.y=FALSE)
colnames(cns.mean.barr)[7] <- "Category"
cns.mean.barr <- cns.mean.barr[,c(2,1,3,4,5,6,7)]
cns.mean.oac <- merge(cns.mean, annotation.sampleIDs[which(annotation.sampleIDs$Category == "PrimaryTumour"),
                      c("Sample","Category")],
                       by.x="sample", by.y="Sample",
                       all.x=FALSE, all.y=FALSE)
cns.mean.both <- rbind(cns.mean.barr, cns.mean.oac)

cns.mean.both$Category <- factor(cns.mean.both$Category,
                               levels=c("NDBO_NP","NDBO_PP",
                                        "LGD","HGD","IMC","trios","PrimaryTumour"))


mycomp <- list(c("NDBO_NP","NDBO_PP"),
               c("NDBO_PP","LGD"),
               c("LGD","HGD"),
               c("HGD","IMC"),
               c("IMC","trios"),
               c("trios","PrimaryTumour"))
pdf("plots_supp/compare.MeanCN.ByBOtype.pdf")
ggboxplot(cns.mean.both, 
          x = "Category", y = "MeanCN",
          color = "Category",palette = c(wes_palette("IsleofDogs1", n =6),"grey"),
          add = "jitter")+
  stat_compare_means(label.x=2, comparisons = mycomp)+
  xlab("")+
  ylab("Average copy number")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

write.csv(cns.mean.both, file="plots_supp/cns.mean.both.csv")

pdf("plots_supp/compare.MeanTotalCN.ByBOtype.pdf")
ggboxplot(cns.mean.both, 
          x = "Category", y = "MeanTotalCN",
          color = "Category",palette = c(wes_palette("IsleofDogs1", n =6),"grey"),
          add = "jitter")+
  stat_compare_means(label.x=2, comparisons = mycomp)+
  xlab("")+
  ylab("Average total copy number")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf("plots_supp/compare.MeanMinorCN.ByBOtype.pdf")
ggboxplot(cns.mean.both, 
          x = "Category", y = "MeanMinorCN",
          color = "Category",palette = c(wes_palette("IsleofDogs1", n =6),"grey"),
          add = "jitter")+
  stat_compare_means(label.x=2, comparisons = mycomp)+
  xlab("")+
  ylab("Average minor copy number")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf("plots_supp/compare.MeanPloidy.ByBOtype.pdf")
ggboxplot(cns.mean.both, 
          x = "Category", y = "MeanPloidy",
          color = "Category",palette = c(wes_palette("IsleofDogs1", n =6),"grey"),
          add = "jitter")+
  stat_compare_means(label.x=2, comparisons = mycomp)+
  xlab("")+
  ylab("Average ploidy")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


