###########
### Check for relation between cancer drivers and signatures.

library(ggpubr)

load("vep.snvs.RData")
load("vep.indels.RData")

vep.snvs.missense <- vep.snvs[which(grepl("missense",vep.snvs$V7) |
                                    grepl("stop_gained",vep.snvs$V7)|
                                    grepl("stop_lost",vep.snvs$V7)),]
save(vep.snvs.missense, file="vep.snvs.missense.RData")

vep.indels.frameshift <- vep.indels[which(grepl("frameshift_variant",vep.indels$V7)|
                                        grepl("inframe_insertion",vep.indels$V7)|
                                        grepl("inframe_deletion",vep.indels$V7)),]
save(vep.indels.frameshift, file="vep.indels.frameshift.RData")

vep.nonsyn <- rbind(vep.snvs.missense, vep.indels.frameshift)
df.nonsyn <- unique(vep.nonsyn[,c("Sample","Gene")])
rev(sort(table(df.nonsyn$Gene)))[1:20]

#### Load sigs and compare between groups:
load("processeddata/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

  selected <- unique(df.nonsyn[which(df.nonsyn$Gene == "APC"),]$Sample)
wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% selected),]$SBS17a,
sigs.allSamples[which(!(rownames(sigs.allSamples) %in% selected)),]$SBS17a)
wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% selected),]$SBS17b,
            sigs.allSamples[which(!(rownames(sigs.allSamples) %in% selected)),]$SBS17b)
wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% selected),]$SBS30,
            sigs.allSamples[which(!(rownames(sigs.allSamples) %in% selected)),]$SBS30)
wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% selected),]$SBS41,
            sigs.allSamples[which(!(rownames(sigs.allSamples) %in% selected)),]$SBS41)



sigs.allSamples$Sample <- rownames(sigs.allSamples)
sbs41.APC <- merge(sigs.allSamples[,c("Sample","SBS41")],
                   df.nonsyn[which(df.nonsyn$Gene == "APC"),],
                   all.x=TRUE, all.y=FALSE)
sbs41.APC[which(is.na(sbs41.APC$Gene)),]$Gene <- "APC_wt"

my_comparisons <- c("APC","APC_wt")
pdf("plots.extra/SBS41byAPC.pdf")
ggboxplot(sbs41.APC[which(sbs41.APC$SBS41>0.05),], x = "Gene", y = "SBS41",
          color = "Gene", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means()
dev.off()

sigs.allSamples$Sample <- rownames(sigs.allSamples)
sbs2.PIK3CA <- merge(sigs.allSamples[,c("Sample","SBS2")],
                   df.nonsyn[which(df.nonsyn$Gene == "PIK3CA"),],
                   all.x=TRUE, all.y=FALSE)
sbs2.PIK3CA[which(is.na(sbs2.PIK3CA$Gene)),]$Gene <- "PIK3CA_wt"

pdf("plots.extra/SBS2byPIK3CA.pdf")
ggboxplot(sbs2.PIK3CA, x = "Gene", y = "SBS2",
          color = "Gene", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means()
dev.off()

sigs.allSamples$Sample <- rownames(sigs.allSamples)
sbs2.KRAS <- merge(sigs.allSamples[,c("Sample","SBS2")],
                     df.nonsyn[which(df.nonsyn$Gene == "KRAS"),],
                     all.x=TRUE, all.y=FALSE)
sbs2.KRAS[which(is.na(sbs2.KRAS$Gene)),]$Gene <- "KRAS_wt"

pdf("plots.extra/SBS2byKRAS.pdf")
ggboxplot(sbs2.KRAS, x = "Gene", y = "SBS2",
          color = "Gene", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means()
dev.off()

pdf("plots.extra/SBS41byAPC.dotplot.pdf")
ggplot(sbs41.APC, aes(x=Gene, y=SBS41)) + 
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, dotsize=1.2)+
  stat_summary(fun.y=mean_sdl, geom="point", shape=18,
               size=3, color="red")+
  stat_compare_means(comparisons = my_comparisons)
  
dev.off()

#W = 132848, p-value = 0.023  -TP53 mutants show greater SBS30; 
#W = 12323, p-value = 0.003508 - KRAS mutants show less SBS17b
#W = 31933, p-value = 0.003587 -ARID1A mutants show less SBS17b
#W = 4884.5, p-value = 0.02388 - EGFR mutants show less SBS30;
#W = 19799, p-value = 0.0008796 - APC mutants show less SBS17a;
#W = 16006, p-value = 0.007374 - PIK3CA mutants show less SBS17a;
#W = 3413, p-value = 0.001916 - PTEN, p-value = 0.01819 MET mutants show less SBS17b
#W = 21476, p-value = 0.006891; ABCB1 mutants show greater sbs17B

genes <- c("TP53","CDKN2A","KRAS","MYC","ERBB2","GATA4",
           "CCND1","GATA6","SMAD4","CDK6","ARID1A","EGFR",
           "CCNE1","CCND3","MUC6","MDM2","KCNQ3","APC",
           "SMARCA4","PIK3CA","ABCB1","PTEN","MET")

load("processeddata/df.bottles.tracksig.may2021.RData")
df.s17 <- df.bottles[which(df.bottles$Signature == "SBS17b"),]
df.s17$Type <- sapply(df.s17$Change,
                      function(x) ifelse(x<0,"Decrease","Increase"))
df.dec <- df.s17[which(df.s17$Type == "Decrease"),]
df.inc <- df.s17[which(df.s17$Type == "Increase"),]

for (g in genes) {
  print("=================")
  print(g)
  selected <- unique(df.nonsyn[which(df.nonsyn$Gene == g),]$Sample)
  dec.yes <- length(unique(df.dec[which(df.dec$Sample %in% selected),]$Sample))
  dec.no <- length(unique(df.dec[which(!(df.dec$Sample %in% selected)),]$Sample))
  inc.yes <- length(unique(df.inc[which(df.inc$Sample %in% selected),]$Sample))
  inc.no <- length(unique(df.inc[which(!(df.inc$Sample %in% selected)),]$Sample))
  m <- matrix(c(inc.yes,inc.no,dec.yes,dec.no),nrow=2)
  print(m)
  print(fisher.test(m))
}

g <- "TP53"
m <- matrix(c(dec.yes,dec.no,inc.yes,inc.no),ncol=2)
print(m)
fisher.test(m)
# dec inc
# yes  270   66
# no  206   97
# data:  m
# p-value = 0.0003854
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.321374 2.814642
# sample estimates:
#   odds ratio 
# 1.924327 
