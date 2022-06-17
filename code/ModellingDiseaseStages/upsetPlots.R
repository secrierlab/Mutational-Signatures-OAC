##############
#### This script generates upset plots in order to investigate signature combinations at different stages.

library(depmixS4)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(cowplot)

# Load mutational signature data:
load("../SigProfiler/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

# Load annotation of samples:
load("processeddata/annotation.sampleIDs.RData")

## Merge signatures and IDs:

sigs.allSamples <- data.frame(sigs.allSamples)
sigs.allSamples$Sample <- rownames(sigs.allSamples)

sigs.norm.annot <- merge(sigs.allSamples, annotation.sampleIDs,
              by.x="Sample", by.y ="Sample", all.x=FALSE, all.y=FALSE)
df.sigs <- melt(sigs.norm.annot)

### Next, check which signatures are dominant in each category:

library(dplyr)
df.sigs%>%
  group_by(Category,variable)%>% 
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value)) %>% 
  filter(Max>0.40)


#########
## Create an upset plot to reflect the number of samples with signatures greater than 10% prevalence:

library(ComplexHeatmap)
library(reshape)
library(UpSetR)
library("ggplotify")

sigs.norm.annot[,2:15] <- ifelse(sigs.norm.annot[,2:15]<0.05,0,1)
sigs.norm.annot$CategoryBin <- sapply(sigs.norm.annot$Category,
                                      function(x) ifelse(x=="Barretts",0,
                                                         ifelse(x=="PrimaryTumour",1,2)))
pdf("plots.upset/upset.primaries.pdf",h=5,onefile=FALSE)
g1 <- upset(sigs.norm.annot[which(sigs.norm.annot$Category == "PrimaryTumour"),], 
      sets = colnames(sigs.norm.annot)[2:15],order.by = "freq")
print(g1)
dev.off()
pdf("plots.upset/upset.barretts.pdf",h=5)
g2 <- upset(sigs.norm.annot[which(sigs.norm.annot$Category == "Barretts"),], 
      sets = colnames(sigs.norm.annot)[2:15],order.by = "freq")
print(g2)
dev.off()
pdf("plots.upset/upset.mets.pdf",h=5)
g3 <- upset(sigs.norm.annot[which(sigs.norm.annot$Category %in% c("LymphNode","Metastasis")),], 
      sets = colnames(sigs.norm.annot)[2:15],order.by = "freq")
print(g3)
dev.off()

pdf("plots.upset/upset.altogether.pdf", w=12)
plot_grid(as.grob(g2), as.grob(g1), as.grob(g3))
dev.off()


## This one does not seem to work:
pdf("plots.upset/upset.oneplotByGroup.pdf",h=5,onefile=FALSE)
upset(sigs.norm.annot, 
            sets = colnames(sigs.norm.annot)[2:15],order.by = c("Category"),
            group.by = "Category", col)
dev.off()

library(ComplexHeatmap)
m1 = make_comb_mat(sigs.norm.annot[which(sigs.norm.annot$Category == "Barretts"),2:15])
m2 = make_comb_mat(sigs.norm.annot[which(sigs.norm.annot$Category == "PrimaryTumour"),2:15])
m3 = make_comb_mat(sigs.norm.annot[which(sigs.norm.annot$Category %in% c("LymphNode","Metastasis")),2:15])

pdf("plots.upset/complexUpset.oneplotByGroup.pdf",w=12,h=5)
u1 <- UpSet(m1[comb_size(m1) >= 5], pt_size = unit(3, "mm"), lwd = 2,
            set_order = c("SBS17a","SBS17b","SBS2","SBS3","SBS8",
                          "SBS18","SBS30","SBS44","SBS35",
                          "SBS41","SBS28",
                          "SBS1","SBS5","SBS40"),
            comb_order = rev(order(comb_size(m1[comb_size(m1) >= 5]))))
u2 <- UpSet(m2[comb_size(m2) >= 5], pt_size = unit(3, "mm"), lwd = 2,
            set_order = c("SBS17a","SBS17b","SBS2","SBS3","SBS8",
                          "SBS18","SBS30","SBS44","SBS35",
                          "SBS41","SBS28",
                          "SBS1","SBS5","SBS40"),
            comb_order = rev(order(comb_size(m2[comb_size(m2) >= 5]))))
u3 <- UpSet(m3[comb_size(m3) >= 5], pt_size = unit(3, "mm"), lwd = 2,
            set_order = c("SBS17a","SBS17b","SBS2","SBS3","SBS8",
                          "SBS18","SBS30","SBS44","SBS35",
                          "SBS41","SBS28",
                          "SBS1","SBS5","SBS40"),
            comb_order = rev(order(comb_size(m3[comb_size(m3) >= 5]))))
plot_grid(as.grob(u1), as.grob(u2),as.grob(u3),nrow = 1,
          rel_widths=c(1,2,1))
dev.off()
