########################
### Sigs compared plus heatmap

library(reshape)
library(ggpubr)
library(pheatmap)

load("data/indelSignatures.fullCohort.RData")
load("data/annotation.sampleIDs.RData")

## Merge signatures and IDs:

sigs.allSamples <- data.frame(sigs.indel)
sigs.allSamples$Sample <- rownames(sigs.allSamples)

sigs <- merge(sigs.allSamples, annotation.sampleIDs,
              by.x="Sample", by.y ="Sample", all.x=FALSE, all.y=FALSE)
df.sigs <- melt(sigs)
df.sigs$Category <- factor(df.sigs$Category, levels=c("Barretts","PrimaryTumour",
                                                  "LymphNode","Metastasis"))

### Now infer the trios:

tb <- table(annotation.sampleIDs$OCCAMS_ID)

pairs <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==2])),]
trios <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==3])),]
singles <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==1])),]
multi1 <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==16])),]
multi2 <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==17])),]


my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","LymphNode"))
pdf("plots.revision/indels.paired.pdf",w=18,h=10)
ggpaired(pairs, x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 0.025, paired = TRUE) +
  xlab("")+
  ylab("% contribution")
dev.off()

my_comparisonsbp <- list(c("Barretts","PrimaryTumour"))
pdf("plots.revision/indels.paired.BarrettsPrimaries.pdf",w=12,h=10)
ggpaired(pairs[which((pairs$Category != "LymphNode")),], 
         x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=3)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisonsbp, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

my_comparisonsbp <- list(c("Barretts","PrimaryTumour"))
pdf("plots.revision/indels.paired.BarrettsPrimaries.nostars.pdf",w=12,h=10)
ggpaired(pairs[which((pairs$Category != "LymphNode")),], 
         x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=3)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisonsbp)+
  xlab("")+
  ylab("% contribution")
dev.off()




### Need to add metastases and lymph nodes:
singles.allmets <- df.sigs[which(!(df.sigs$Sample %in% pairs$Sample)),]
singles.allmets[which(singles.allmets$Category == "LymphNode"),]$Category <- "Metastasis"
my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","Metastasis"))
pdf("plots.revision/indels.singles.pdf",w=12,h=9)
ggviolin(singles.allmets, x = "Category", 
         y = "value",
         palette ="npg",
         fill = "Category",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=3)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.revision/indels.singles.nostars.pdf",w=12,h=9)
ggviolin(singles.allmets, x = "Category", 
         y = "value",
         palette ="npg",
         fill = "Category",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=3)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()


pdf("plots.revision/indels.singles.density.pdf",w=18,h=20)
ggdensity(singles.allmets, x = "value",
         add = "mean",rug=TRUE,
         palette ="npg",
         color = "Category",fill="Category",alpha = 0.7)+
  facet_wrap(~variable,scales = "free",nrow=7)+
  xlab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

write.csv(singles.allmets, file="plots_supp/indels.singles.allmets.density.csv")

