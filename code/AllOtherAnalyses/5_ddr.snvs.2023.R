##################################
### This is the DDR mutational signature analysis.

library(NMF)
library(reshape)
library(ggplot2)

## Read DDR genes:
library(gdata)
ddr <- read.xls("~/Desktop/UCL/supervisedProjects/DanielJacobson/Rotation1-master/DDRpathways.xlsx")

## Load DDR clonality info:
load("data/gr.clonalityDDR.plusannot.nonsyn.pathways.june2021.RData")

#########################
#### NMF procedure for early/late patterns in DDR:

### Count how many early vs late mutations there are in each gene in ddr:

gr.clonalityDDR.plusannot.nonsyn.pathways$ClonalityTiming <- apply(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Clonality","Timing")],
                                                          1, function(x) paste0(x[1],":",x[2]))
gr.clonalityDDR.plusannot.nonsyn.pathways <- gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$Process!= "Direct Repair (not in humans)"),]
gr.clonalityDDR.plusannot.nonsyn.pathways <- gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$Pathway.1!= ""),]
ddr.process <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","Pathway.1")])

### load annot:
load("data/annotation.sampleIDs.RData")
table(annotation.sampleIDs[which(annotation.sampleIDs$Sample %in% rownames(ddr.process)),]$Category)
# Barretts PrimaryTumour 
# 1           355
# The BO is:
bo <- annotation.sampleIDs[which((annotation.sampleIDs$Sample %in% rownames(ddr.process))&
                             (annotation.sampleIDs$Category == "Barretts")),]$Sample

# Remove Barrett's sample:
gr.clonalityDDR.plusannot.nonsyn.pathways <- gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$Sample != bo),]
gr.clonalityDDR.plusannot.nonsyn.pathways <- unique(gr.clonalityDDR.plusannot.nonsyn.pathways)
ddr.process <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","Pathway.1")])

ddr.process.ce <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$ClonalityTiming == "Clonal:Early"),c("Sample","Pathway.1")])
ddr.process.cl <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$ClonalityTiming == "Clonal:Late"),c("Sample","Pathway.1")])
ddr.process.sb <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[which(gr.clonalityDDR.plusannot.nonsyn.pathways$ClonalityTiming == "Subclonal:NA"),c("Sample","Pathway.1")])

### Print 3 heatmaps:
library(pheatmap)
library(cowplot)
library(viridis)

p1 <- pheatmap(ddr.process.ce, cluster_cols = FALSE, show_rownames = FALSE, color = rev(magma(100)))
p2 <- pheatmap(ddr.process.cl, cluster_cols = FALSE, show_rownames = FALSE, color = rev(magma(100)))
p3 <- pheatmap(ddr.process.sb, cluster_cols = FALSE, show_rownames = FALSE, color = rev(magma(100)))

pdf("plots_5/ddr.snvs.heatmaps.clonalearly.pdf")
print(p1)
dev.off()

pdf("plots_5/ddr.snvs.heatmaps.clonallate.pdf")
print(p2)
dev.off()

pdf("plots_5/ddr.snvs.heatmaps.subclonal.pdf")
print(p3)
dev.off()

### Next, plot %samples altered per pathway:
ddr.pathwayTimings<- table(unique(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","Pathway.1","ClonalityTiming")])[,c("Pathway.1","ClonalityTiming")])
ddr.pathwayTimings.frac <- ddr.pathwayTimings/777

df.melt.ddr <- melt(ddr.pathwayTimings.frac)
df.melt.ddr$Pathway.1 <- factor(df.melt.ddr$Pathway.1,
                                levels = rev(unique(df.melt.ddr$Pathway.1)))

# Copy number variation:

load("../../tracksig_barrettsmetsprimaries/primaries.cn.RData")
load("../../tracksig_barrettsmetsprimaries/barretts.cn.RData")
load("../../tracksig_barrettsmetsprimaries/mets.cn.RData")
load("../../tracksig_barrettsmetsprimaries/extra.cn.RData")
cns <- rbind(primaries.cn,barretts.cn,mets.cn,extra.cn)

# Also load ploidies:
load("../../tracksig_barrettsmetsprimaries/barretts.ploidies.RData")
load("../../tracksig_barrettsmetsprimaries/primaries.ploidies.RData")
load("../../tracksig_barrettsmetsprimaries/mets.ploidies.RData")
load("../../tracksig_barrettsmetsprimaries/extra.ploidies.RData")
ploidies <- rbind(primaries.ploidies,barretts.ploidies,mets.ploidies,extra.ploidies)

cns.ploid <- merge(cns, ploidies,
                   by.x="sample",by.y="samplename",
                   all.x=FALSE, all.y=FALSE)
cns.ploid$LOH <- sapply(cns.ploid$minor_cn, function(x) ifelse(x==0,1,0))
cns.ploid$ploidy <- as.numeric(as.character(cns.ploid$ploidy))
cns.ploid$RelativeCN <- cns.ploid$total_cn/cns.ploid$ploidy
cns.ploid$CNchange <- sapply(cns.ploid$RelativeCN, 
                             function(x) ifelse(x>=2,"AMP",ifelse(x<=0.5,"DEL","NEUTRAL")))
cns.ploid[which((cns.ploid$CNchange == "NEUTRAL")&
                  (cns.ploid$LOH == 1)),]$CNchange <- "LOH"  
save(cns.ploid, file="cns.ploid.RData")

### Intersect CN with genes in DDR pathways:
hg38 <- read.csv("../hg38_geneLocations.txt")

library(GenomicRanges)
gr.hg38 <- GRanges(
  seqnames=hg38$Chromosome.scaffold.name,
  ranges=IRanges(start=hg38$Gene.start..bp., end=hg38$Gene.end..bp.),
  gene=hg38$Gene.name)

gr.cns <- GRanges(
  seqnames=cns.ploid$chromosome,
  ranges=IRanges(start=cns.ploid$start, end=cns.ploid$end),
  CN=cns.ploid$CNchange,
  sample=cns.ploid$sample)

ovlp <- data.frame(findOverlaps(gr.cns, gr.hg38, type="any"))
gr.overlap <- cbind(as.data.frame(gr.cns[as.matrix(ovlp)[,1]]),
                    as.data.frame(gr.hg38[as.matrix(ovlp)[,2]])) 

### Next, select only DDR genes:
gr.ddr <- merge(gr.overlap, ddr, by.x="gene", by.y="Gene.ID",
                all.x=FALSE, all.y=FALSE)
gr.ddr.cna <- gr.ddr[which(gr.ddr$CN != "NEUTRAL"),]
gr.ddr.cna <- gr.ddr.cna[which(gr.ddr.cna$Pathway.1 !=""),]
gr.ddr.cna <- gr.ddr.cna[which(gr.ddr.cna$sample %in% annotation.sampleIDs[which(annotation.sampleIDs$Category=="PrimaryTumour"),]$Sample),]

### Next, plot %samples altered per pathway:
ddr.cna<- table(unique(gr.ddr.cna[,c("Pathway.1","CN","sample")])[,c("Pathway.1","CN")])
ddr.cna <- ddr.cna/777

df.melt.cn <- melt(ddr.cna)
df.melt.cn$Pathway.1 <- factor(df.melt.cn$Pathway.1,
                                levels = rev(unique(df.melt.cn$Pathway.1)))
pdf("plots_5/cnFrequenciesByPathway.pdf",w=5)
ggplot(df.melt.cn,aes(x = CN, y = Pathway.1)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(CN, Pathway.1, label=round(value,2)), colour = "white", check_overlap = TRUE)  +
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  xlab("")
dev.off()

### Also indels:
load("data/indelsFull.all.RData")
indels.nonsyn <- indels.all[which(indels.all$V7 %in%
                                    c("inframe_insertion","inframe_deletion",
                                      "frameshift_variant", "frameshift_variant,NMD_transcript_variant",
                                      "TF_binding_site_variant,TFBS_ablation",
                                      "frameshift_variant,splice_region_variant",
                                      "stop_gained,frameshift_variant",
                                      "stop_gained,frameshift_variant,NMD_transcript_variant",
                                      "frameshift_variant,splice_region_variant,NMD_transcript_variant",
                                      "inframe_insertion,NMD_transcript_variant",
                                      "splice_acceptor_variant,frameshift_variant",
                                      "frameshift_variant,initiator_codon_variant",
                                      "inframe_deletion,incomplete_terminal_codon_variant",
                                      "inframe_deletion,NMD_transcript_variant",
                                      "frameshift_variant,stop_lost,NMD_transcript_variant",
                                      "inframe_deletion,splice_region_variant",
                                      "frameshift_variant,stop_lost",
                                      "incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant",
                                      "inframe_deletion,splice_region_variant,NMD_transcript_variant",
                                      "frameshift_variant,stop_lost,splice_region_variant,NMD_transcript_variant",
                                      "stop_gained,frameshift_variant,splice_region_variant",
                                      "inframe_insertion,splice_region_variant",
                                      "frameshift_variant,stop_retained_variant",
                                      "stop_gained,inframe_insertion",
                                      "frameshift_variant,stop_lost,splice_region_variant",
                                      "inframe_insertion,incomplete_terminal_codon_variant,coding_sequence_variant",
                                      "frameshift_variant,stop_retained_variant,NMD_transcript_variant",
                                      "inframe_insertion,splice_region_variant,NMD_transcript_variant",
                                      "stop_gained,inframe_deletion",
                                      "initiator_codon_variant,inframe_deletion",
                                      "stop_gained,inframe_insertion,NMD_transcript_variant",
                                      "stop_gained,inframe_deletion,NMD_transcript_variant")),]


indel.ddr <- merge(indels.nonsyn, ddr, by.x="Gene", by.y="Gene.ID",
                all.x=FALSE, all.y=FALSE)
indel.ddr <- indel.ddr[which(indel.ddr$Pathway.1 !=""),]
indel.ddr <- indel.ddr[which(indel.ddr$Sample %in% annotation.sampleIDs[which(annotation.sampleIDs$Category=="PrimaryTumour"),]$Sample),]

### Next, plot %samples altered per pathway:
ddr.indel<- table(unique(indel.ddr[,c("Pathway.1","Sample")])$Pathway.1)
ddr.indel <- ddr.indel/777

df.melt.indel <- melt(ddr.indel)
df.melt.indel$Var.1 <- factor(df.melt.indel$Var.1,
                               levels = rev(unique(df.melt.indel$Var.1)))
df.melt.indel$Type <- "Indel"
pdf("plots_5/indelFrequenciesByPathway.pdf",w=4)
ggplot(df.melt.indel,aes(x = Type, y = Var.1)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(Type, Var.1, label=round(value,2)), colour = "white", check_overlap = TRUE)  +
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  xlab("")
dev.off()



# Put everything together:
colnames(df.melt.ddr)[2] <- "Change"
colnames(df.melt.cn)[2] <- "Change"
colnames(df.melt.indel) <- c("Pathway.1","value","Change")
df.melt.indel <- df.melt.indel[,c(1,3,2)]

df.melt <- rbind(df.melt.ddr, df.melt.indel, df.melt.cn)
df.melt$Pathway.1 <- factor(df.melt$Pathway.1,
                               levels = rev(unique(df.melt$Pathway.1)))

pdf("plots_5/SNVandCNfrequenciesByPathway.All.pdf",w=8)
ggplot(df.melt,aes(x = Change, y = Pathway.1)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(Change, Pathway.1, label=round(value,2)), colour = "white", check_overlap = TRUE)  +
  scale_fill_gradient2(low="yellow", high="#1D3557", mid="#AD343E",
                       midpoint=0.25,na.value="black")+
  xlab("")
dev.off()

pdf("plots_5/SNVandCNfrequenciesByPathway.SNV.pdf",w=5)
ggplot(df.melt[which(df.melt$Change %in% unique(df.melt.ddr$Change)),],
       aes(x = Change, y = Pathway.1)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(Change, Pathway.1, label=round(value,2)), colour = "white", check_overlap = TRUE)  +
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  xlab("")
dev.off()

pdf("plots_5/SNVandCNfrequenciesByPathway.CN.pdf",w=5)
ggplot(df.melt[which(df.melt$Change %in% unique(df.melt.cn$Change)),],
       aes(x = Change, y = Pathway.1)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(Change, Pathway.1, label=round(value,2)), colour = "white", check_overlap = TRUE)  +
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  xlab("")
dev.off()

library(reshape2)
m <- acast(Pathway.1~Change, value="value", data=df.melt)
m[is.na(m)] <- 0

# cluster m:
p <- pheatmap(m)
h <- hclust(dist(m))
cls <- cutree(h,3)

# Reorder df.melt:
# df.melt$Pathway.1 <- factor(df.melt$Pathway.1,
#                             levels = c(rev(h$labels[h$order[cls == 2]]),
#                                        rev(h$labels[h$order[cls == 3]]),
#                                        rev(h$labels[h$order[cls == 1]])))

m1 <- m[order(m[,c("AMP")],decreasing=T),]
m2 <- m1[order(m1[,c("DEL")],decreasing=T),]

df.melt$Pathway.1 <- factor(df.melt$Pathway.1,
                            levels = rev(rownames(m2)))

pdf("plots_5/SNVandCNfrequenciesByPathway.All.clustered.pdf",w=8)
ggplot(df.melt,aes(x = Change, y = Pathway.1)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(Change, Pathway.1, label=round(value,2)), colour = "white", check_overlap = TRUE)  +
  scale_fill_gradient2(low="yellow", high="#1D3557", mid="#AD343E",
                       midpoint=0.25,na.value="black")+
  xlab("")
dev.off()
   
write.csv(df.melt, file="plots_5/ddr.alts.csv")

###############################################
### Next, check correlations with signatures:
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples$Sample <- rownames(sigs.allSamples)
load("data/indelSignatures.fullCohort.RData")
sigs.indel$Sample <- rownames(sigs.indel)

# For BER:
df1 <- unique(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","Pathway.1")])
df2 <- unique(gr.ddr.cna[,c("sample","Pathway.1")])
colnames(df2)[1] <- "Sample"
df3 <- unique(indel.ddr[,c("Sample","Pathway.1")])
df <- unique(rbind(df1,df2,df3))

df.ber <- df[which(df$Pathway.1 == "BER"),]
df.ber$BER <- 1

df.hr<- df[which(df$Pathway.1 == "HR (Homologous Recombination)"),]
df.hr$HR <- 1

df.mmr<- df[which(df$Pathway.1 == "MMR"),]
df.mmr$MMR <- 1

## Compare SBS30 sigs b/tgroups:

wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% df.ber$Sample),]$SBS30,
            sigs.allSamples[which(!(rownames(sigs.allSamples) %in% df.ber$Sample)),]$SBS30)

sigs.allSamples <- merge(sigs.allSamples, df.ber[,c("Sample","BER")],
                         by.x="Sample", by.y="Sample",
                         all.x=TRUE, all.y=FALSE)
sigs.allSamples[is.na(sigs.allSamples$BER),]$BER <- 0

library(ggpubr)

mycomp <- list(c("1","0"))

pdf("plots_5/BER.sigsCompared.pdf")
ggboxplot(sigs.allSamples, x = "BER", y = "SBS30",
          color = "BER", palette =c("#CC5803", "#1D3557"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)
dev.off()

sigs.allSamples <- merge(sigs.allSamples, df.hr[,c("Sample","HR")],
                         by.x="Sample", by.y="Sample",
                         all.x=TRUE, all.y=FALSE)
sigs.allSamples[is.na(sigs.allSamples$HR),]$HR <- 0

pdf("plots_5/HR.sigsCompared.pdf")
ggboxplot(sigs.allSamples, x = "HR", y = "SBS3",
          color = "HR", palette =c("#CC5803", "#1D3557"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)
dev.off()

sigs.indel <- merge(sigs.indel, df.hr[,c("Sample","HR")],
                         by.x="Sample", by.y="Sample",
                         all.x=TRUE, all.y=FALSE)
sigs.indel[is.na(sigs.indel$HR),]$HR <- 0

pdf("plots_5/HR.indelCompared.ID6.pdf")
ggboxplot(sigs.indel, x = "HR", y = "ID6",
         color = "HR", palette =c("#CC5803", "#1D3557"),
         add = "jitter")+
  stat_compare_means(comparisons = mycomp)
dev.off()

pdf("plots_5/HR.indelCompared.ID8.pdf")
ggboxplot(sigs.indel, x = "HR", y = "ID8",
         color = "HR", palette =c("#CC5803", "#1D3557"),
         add = "jitter")+
  stat_compare_means(comparisons = mycomp)
dev.off()

sigs.allSamples <- merge(sigs.allSamples, df.mmr[,c("Sample","MMR")],
                         by.x="Sample", by.y="Sample",
                         all.x=TRUE, all.y=FALSE)
sigs.allSamples[is.na(sigs.allSamples$MMR),]$MMR <- 0

pdf("plots_5/MMR.sigsCompared.pdf")
ggboxplot(sigs.allSamples, x = "MMR", y = "SBS44",
         color = "MMR", palette =c("#CC5803", "#1D3557"),
         add = "jitter")+
  stat_compare_means(comparisons = mycomp)
dev.off()

sigs.indel <- merge(sigs.indel, df.mmr[,c("Sample","MMR")],
                    by.x="Sample", by.y="Sample",
                    all.x=TRUE, all.y=FALSE)
sigs.indel[is.na(sigs.indel$MMR),]$MMR <- 0

pdf("plots_5/MMR.indelCompared.ID7.pdf")
ggboxplot(sigs.indel, x = "MMR", y = "ID7",
          color = "MMR", palette =c("#CC5803", "#1D3557"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp)
dev.off()

write.csv(sigs.allSamples, file="plots_5/ddr.sigs.compared.csv")
write.csv(sigs.indel, file="plots_5/ddr.sigsindel.compared.csv")


df1 <- table(gr.clonalityDDR.plusannot.nonsyn.pathways[,c("Sample","Pathway.1")])
colnames(df1)<- paste0(colnames(df1),":SNV")
df2 <- table(gr.ddr.cna[which(gr.ddr.cna$CN=="AMP"),c("sample","Pathway.1")])
colnames(df2)<- paste0(colnames(df2),":AMP")
df3 <- table(gr.ddr.cna[which(gr.ddr.cna$CN=="DEL"),c("sample","Pathway.1")])
colnames(df3)<- paste0(colnames(df3),":DEL")
df4 <- table(gr.ddr.cna[which(gr.ddr.cna$CN=="LOH"),c("sample","Pathway.1")])
colnames(df4)<- paste0(colnames(df4),":LOH")
df5 <- table(indel.ddr[,c("Sample","Pathway.1")])
colnames(df5)<- paste0(colnames(df5),":Indel")

df1 <- as.data.frame(df1)
df2 <- data.frame(df2)
df3 <- data.frame(df3)
df4 <- data.frame(df4)
df5 <- data.frame(df5)

df.all <- merge(df1,df2, by.x="Sample", by.y="Sample",
                all.x=TRUE, all.y=TRUE)

df$PathType <- paste0(df$Pathway.1,":",df$Type)
save(df, file="df.pathwayAlterationCounts.RData")

write.csv(df, file="plots_5/pathwayAlterations.csv")
