#####################
#### Here check relation between MMR high samples and immunity.

library(ggpubr)
library(reshape)
library(wesanderson)

load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

### Determine MMR high samples:

sigs.allSamples$Sample <- rownames(sigs.allSamples)
df.mmr <- sigs.allSamples[,c("Sample","SBS44")]
df.mmr$MMR <- 0
df.mmr[which(df.mmr$SBS44 > 0.05),]$MMR <- 1
table(df.mmr$MMR)
# 0   1 
# 961  36

### Merge with annotation file:
load("data/annotation.sampleIDs.RData")
df.mmr <- merge(df.mmr, annotation.sampleIDs,
                by.x="Sample",by.y="Sample",
                all.x=FALSE, all.y=FALSE)
table(df.mmr[,c("Category","MMR")])

mmr <- df.mmr[which(df.mmr$MMR==1),]$Sample
save(mmr, file="processeddata/MMR.RData")

### Calculate immunity info:
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
                                   statMethod = "ssgsea")
df.celldeconv <- data.frame(t(celldeconv))
df.celldeconv$Sample <- rownames(df.celldeconv)


## Now merge with immunity:
df.mmr$ID <- sapply(df.mmr$Sample, function(x) strsplit(x,"_vs_")[[1]][1])
df.mmr.tme <- merge(df.mmr, df.celldeconv,
                by.x="ID", by.y="Sample",
                all.x=FALSE, all.y=FALSE)
table(df.mmr.tme[,c("MMR")])
# 0   1 
# 191  12 

write.csv(df.mmr.tme, file="plots_supp/df.mmr.tme.csv")

df.melt <- melt(df.mmr.tme, id.vars = c("ID","Sample","SBS44","MMR","TumourID","Category","OCCAMS_ID"))

### Plot immune infiltrates compared:
my_comparisons <- list(c("0","1"))
pdf("plots_supp/MMR.TMEcompared.ssgsea.pdf",w=10,h=8)
ggviolin(df.melt, 
         x = "MMR", y = "value",
         palette =wes_palette("Zissou1",4, type = "continuous"),
         fill = "MMR",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable, scale="free")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("")+
  ylab("Infiltration score")
dev.off()

### Also link with TMB:
load("../tmb.annot.RData")
df.tmb <- merge(df.mmr, tmb.annot[,c("Sample","TMB")],
                by.x="Sample", by.y="Sample",
                all.x=FALSE, all.y=FALSE) 
pdf("plots_supp/MMR.TMBcompared.pdf",w=3,h=5)
ggviolin(df.tmb, 
         x = "MMR", y = "TMB",
         palette =wes_palette("Zissou1",4, type = "continuous"),
         fill = "MMR",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("")+
  ylab("Tumour mutational burden (log10)")
dev.off()

write.csv(df.tmb, file="plots_supp/df.tmb.csv")
