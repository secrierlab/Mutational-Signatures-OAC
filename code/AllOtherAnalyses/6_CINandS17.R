##############################
### This script checks the relation between CIN and S17.

library(reshape)
library(ggpubr)
library(gdata)

# Load mut sig data:
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

# Load CN data:
load("../../tracksig_barrettsmetsprimaries/barretts.cn.RData")
load("../../tracksig_barrettsmetsprimaries/mets.cn.RData")
load("../../tracksig_barrettsmetsprimaries/primaries.cn.RData")
load("../../tracksig_barrettsmetsprimaries/barretts.ploidies.RData")
load("../../tracksig_barrettsmetsprimaries/mets.ploidies.RData")
load("../../tracksig_barrettsmetsprimaries/primaries.ploidies.RData")
load("../../tracksig_barrettsmetsprimaries/extra.cn.RData")
load("../../tracksig_barrettsmetsprimaries/extra.ploidies.RData")

### First, check if S17 correlates with ploidy:
all.ploidies <- rbind(barretts.ploidies, primaries.ploidies, mets.ploidies,
                      extra.ploidies)

sigs.allSamples$Sample <- rownames(sigs.allSamples)
# Load annotation:
load("data/annotation.sampleIDs.RData")

sigs.allSamples <- merge(sigs.allSamples,
                         annotation.sampleIDs,
                         by.x="Sample",by.y="Sample",
                         all.x=FALSE, all.y=FALSE)
sigs.allSamples[which(sigs.allSamples$Category == "LymphNode"),]$Category <- "Metastasis"
sigs.plusploidy <- merge(sigs.allSamples,
                         all.ploidies,
                         by.x="Sample", by.y="samplename",
                         all.x=FALSE, all.y=FALSE)
sigs.plusploidy$ploidy <- as.numeric(as.character(sigs.plusploidy$ploidy))

### Iterate through every signature:
for (s in 2:15) {
  c <- cor.test(sigs.plusploidy[,s],sigs.plusploidy$ploidy)
  print(paste0(colnames(sigs.plusploidy)[s],": p=",c$p.value,"; R=",c$estimate))
}

### Barrett's only:
for (s in 2:15) {
  c <- cor.test(sigs.plusploidy[which(sigs.plusploidy$Category == "Barretts"),s],
                sigs.plusploidy[which(sigs.plusploidy$Category == "Barretts"),]$ploidy)
  print(paste0(colnames(sigs.plusploidy)[s],": p=",c$p.value,"; R=",c$estimate))
}

### Primaries only:
for (s in 2:15) {
  c <- cor.test(sigs.plusploidy[which(sigs.plusploidy$Category == "PrimaryTumour"),s],
                sigs.plusploidy[which(sigs.plusploidy$Category == "PrimaryTumour"),]$ploidy)
  print(paste0(colnames(sigs.plusploidy)[s],": p=",c$p.value,"; R=",c$estimate))
}

### Split into diploid and the rest:
sigs.plusploidy$PloidyCategory <- sapply(sigs.plusploidy$ploidy,
                                         function(x) ifelse(round(x)==2,"diploid",
                                                            ifelse(round(x)<2,"monoploid","polyploid")))
sigs.plusploidy.keep <- sigs.plusploidy[which(sigs.plusploidy$PloidyCategory != "monoploid"),]

## Compare sigs between ploidy categories:
sigs.melt <- melt(sigs.plusploidy.keep,
                  id.vars = c("Sample","TumourID","Category","OCCAMS_ID",
                              "ploidy","PloidyCategory"))

mycomp <- list(c("diploid","polyploid"))
pdf("plots_6/sigPrevalence.ComparingPloidies.pdf",w=10,h=8)
ggviolin(sigs.melt, 
         x = "PloidyCategory", y = "value",
         palette ="npg",
         fill = "PloidyCategory",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

pdf("plots_6/sigPrevalence.ComparingPloidies.bySampleType.pdf",w=10,h=8)
ggviolin(sigs.melt, 
         x = "PloidyCategory", y = "value",
         palette ="npg",
         fill = "PloidyCategory",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_grid(Category~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

### Now try to compare ploidy between S17 yes/no:
sigs.plusploidy$S17A_categ <- sapply(sigs.plusploidy$SBS17a,
                                     function(x) ifelse(x<=0.05,"no","yes"))
sigs.plusploidy$S17B_categ <- sapply(sigs.plusploidy$SBS17b,
                                     function(x) ifelse(x<=0.05,"no","yes"))
sigs.plusploidy$S17AB_categ <- apply(sigs.plusploidy[,c("SBS17a","SBS17b")],1,
                                     function(x) ifelse(x[1]+x[2]<=0.05,"no","yes"))

mycomp2 <- list(c("yes","no"))
pdf("plots_6/ploidy.bySBS17a.pdf")
ggviolin(sigs.plusploidy, 
         x = "S17A_categ", y = "ploidy",
         palette ="npg",
         fill = "S17A_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("Ploidy")
dev.off()

pdf("plots_6/ploidy.bySBS17b.pdf")
ggviolin(sigs.plusploidy, 
         x = "S17B_categ", y = "ploidy",
         palette ="npg",
         fill = "S17B_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("Ploidy")
dev.off()

write.csv(sigs.plusploidy, file="plots_6/sigs.plusploidy.csv")

pdf("plots_6/ploidy.bySBS17ab.pdf")
ggviolin(sigs.plusploidy, 
         x = "S17AB_categ", y = "ploidy",
         palette ="npg",
         fill = "S17AB_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("Ploidy")
dev.off()

###############
## Next, investigate CIN as defined by CN events:

all.cn <- rbind(barretts.cn,primaries.cn,mets.cn,extra.cn)

# Merge with ploidy:
all.cn <- merge(all.cn, all.ploidies,
                by.x="sample",by.y="samplename",
                all.x=FALSE, all.y=FALSE)
all.cn$ploidy <- as.numeric(as.character(all.cn$ploidy))
all.cn$RelativeCN <- all.cn$total_cn/all.cn$ploidy
all.cn$CNA <- sapply(all.cn$RelativeCN, 
                     function(x) ifelse(x>=2,"Gain",
                                        ifelse(x<=0.75,"Loss","Neutral")))
all.cn$LOH <- sapply(all.cn$minor_cn, function(x) ifelse(x==0,"yes","no"))                     
all.cn$SegmentLength <- apply(all.cn[,c("start","end")],1,
                              function(x) x[2]-x[1]+1)

# To calculate fraction of chromosome, I need to merge with chromosome lengths:
chrlen <- read.xls("../hg38_chrLengths.xlsx", header=FALSE)
colnames(chrlen)[c(1,2)] <- c("chr","length")

all.cn.merged <- merge(all.cn, chrlen[,c("chr","length")],
                       by.x="chromosome",by.y="chr",
                       all.x=FALSE, all.y=FALSE)
all.cn.merged$FractionChr <- all.cn.merged$SegmentLength/all.cn.merged$length
all.cn.merged$CINsegment <- apply(all.cn.merged[,c("CNA","LOH","FractionChr")],1,
                                  function(x) ifelse(((x[1]!="Neutral")|(x[2]=="yes"))&((as.numeric(x[3])>0.05)|(as.numeric(x[3])==0.05)),"yes","no"))

## Count the number of CIN segments that span >5% of a chromosome:
df.cin <- data.frame(table(all.cn.merged[,c("sample","CINsegment")]))
df.cin.yes <- df.cin[which(df.cin$CINsegment == "yes"),]
df.cin.yes$CINscore <- scale(df.cin.yes$Freq)
  
# Merge with signature info:
sigs.pluscin <- merge(sigs.allSamples,
                      df.cin.yes[,c("sample","Freq","CINscore")],
                         by.x="Sample", by.y="sample",
                         all.x=FALSE, all.y=FALSE)
cor.test(sigs.pluscin$SBS17b,sigs.pluscin$CINscore)
cor.test(sigs.pluscin$SBS17a,sigs.pluscin$CINscore)
cor.test(sigs.pluscin$SBS17b,sigs.pluscin$Freq)
cor.test(sigs.pluscin$SBS17a,sigs.pluscin$Freq)
# no obvious correlation here

sigs.pluscin$S17A_categ <- sapply(sigs.pluscin$SBS17a,
                                     function(x) ifelse(x<=0.05,"no","yes"))
sigs.pluscin$S17B_categ <- sapply(sigs.pluscin$SBS17b,
                                     function(x) ifelse(x<=0.05,"no","yes"))
sigs.pluscin$S17AB_categ <- apply(sigs.pluscin[,c("SBS17a","SBS17b")],1,
                                     function(x) ifelse(x[1]+x[2]<=0.05,"no","yes"))

mycomp2 <- list(c("yes","no"))
pdf("plots_6/CIN.bySBS17a.cutoff5perc.pdf")
ggviolin(sigs.pluscin, 
         x = "S17A_categ", y = "CINscore",
         palette ="npg",
         fill = "S17A_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("CIN score")
dev.off()

pdf("plots_6/CIN.bySBS17b.cutoff5perc.pdf")
ggviolin(sigs.pluscin, 
         x = "S17B_categ", y = "CINscore",
         palette ="npg",
         fill = "S17B_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("CIN score")
dev.off()

pdf("plots_6/CIN.bySBS17ab.cutoff5perc.pdf")
ggviolin(sigs.pluscin, 
         x = "S17AB_categ", y = "CINscore",
         palette ="npg",
         fill = "S17AB_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("CIN score")
dev.off()

## Also directly the frequency:
pdf("plots_6/CINfreq.bySBS17a.cutoff5perc.pdf")
ggviolin(sigs.pluscin, 
         x = "S17A_categ", y = "Freq",
         palette ="npg",
         fill = "S17A_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("CIN frequency")
dev.off()

pdf("plots_6/CINfreq.bySBS17b.cutoff5perc.pdf")
ggviolin(sigs.pluscin, 
         x = "S17B_categ", y = "Freq",
         palette ="npg",
         fill = "S17B_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("CIN frequency")
dev.off()

pdf("plots_6/CINfreq.bySBS17ab.cutoff5perc.pdf")
ggviolin(sigs.pluscin, 
         x = "S17AB_categ", y = "Freq",
         palette ="npg",
         fill = "S17AB_categ",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~Category, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = mycomp2)+
  xlab("")+
  ylab("CIN frequency")
dev.off()

write.csv(sigs.pluscin, file="plots_6/sigs.pluscin.sbs17b.csv")
