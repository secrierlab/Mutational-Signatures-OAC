#############
### Generate oncoprint for cancer drivers.

# load drivers:
drivers <- read.delim("data/OACdrivers_Frankell.txt", header=FALSE)$V1

load("data/cns.ploid.RData")

### Intersect CN with genes in DDR pathways:
hg38 <- read.csv("data/hg38_geneLocations.txt")

# Select only drivers:
hg38.drivers <- hg38[which(hg38$Gene.name %in% drivers),]

library(GenomicRanges)
gr.hg38 <- GRanges(
  seqnames=hg38.drivers$Chromosome.scaffold.name,
  ranges=IRanges(start=hg38.drivers$Gene.start..bp., end=hg38.drivers$Gene.end..bp.),
  gene=hg38.drivers$Gene.name)

gr.cns <- GRanges(
  seqnames=cns.ploid$chromosome,
  ranges=IRanges(start=cns.ploid$start, end=cns.ploid$end),
  CN=cns.ploid$CNchange,
  sample=cns.ploid$sample)

ovlp <- data.frame(findOverlaps(gr.cns, gr.hg38, type="any"))
gr.overlap <- cbind(as.data.frame(gr.cns[as.matrix(ovlp)[,1]]),
                    as.data.frame(gr.hg38[as.matrix(ovlp)[,2]])) 
gr.overlap$AMP <- sapply(gr.overlap$CN, function(x) ifelse(x=="AMP",1,0))
gr.overlap$DEL <- sapply(gr.overlap$CN, function(x) ifelse(x=="DEL",1,0))
gr.overlap$LOH <- sapply(gr.overlap$CN, function(x) ifelse(x=="LOH",1,0))

# Now make this into a matrix:
library(reshape2)
mat.cn.amp <- acast(gene~sample, value.var = "AMP",fun.aggregate = mean,
                data=gr.overlap)
mat.cn.del <- acast(gene~sample, value.var = "DEL",fun.aggregate = mean,
                    data=gr.overlap)
mat.cn.loh <- acast(gene~sample, value.var = "LOH",fun.aggregate = mean,
                    data=gr.overlap)

## Next, select snvs and indels:
load("data/snvsFull.nonsyn.all.RData")
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

# Select only genes of interest:
snvs.selected <- unique(snvs.all.nonsyn[which(snvs.all.nonsyn$Gene %in% drivers),])
library(stringr)
snvs.selected$Sample <- sapply(snvs.selected$Sample, function(x)
  ifelse(grepl("\\.",x),str_replace_all(x,"\\.","_"),x))
indels.selected <- unique(indels.nonsyn[which(indels.nonsyn$Gene %in% drivers),])

mat.snvs <- array(0,c(length(drivers),length(unique(c(snvs.selected$Sample,indels.selected$Sample,colnames(mat.cn.amp))))))
rownames(mat.snvs) <- drivers
colnames(mat.snvs) <- unique(c(snvs.selected$Sample,indels.selected$Sample,colnames(mat.cn.amp)))
for (i in 1:nrow(snvs.selected)) {
  print(i)
  mat.snvs[snvs.selected[i,]$Gene,snvs.selected[i,]$Sample] <- 1
}

mat.indels <- array(0,c(length(drivers),length(unique(c(snvs.selected$Sample,indels.selected$Sample,colnames(mat.cn.amp))))))
rownames(mat.indels) <- drivers
colnames(mat.indels) <- unique(c(snvs.selected$Sample,indels.selected$Sample,colnames(mat.cn.amp)))
for (i in 1:nrow(indels.selected)) {
  print(i)
  mat.indels[indels.selected[i,]$Gene,indels.selected[i,]$Sample] <- 1
}


# next, load annotation and plot one matrix for each cancer stage:
load("data/annotation.sampleIDs.RData")

b <- annotation.sampleIDs[which(annotation.sampleIDs$Category == "Barretts"),]$Sample
o <- annotation.sampleIDs[which(annotation.sampleIDs$Category == "PrimaryTumour"),]$Sample
m <- annotation.sampleIDs[which(annotation.sampleIDs$Category %in% c("LymphNode","Metastasis")),]$Sample

# also load mut sig data:
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples[sigs.allSamples<0.05] <- 0

mat_list.barr <- list()
genes <- intersect(rownames(mat.snvs),rownames(mat.cn.amp))
mat_list.barr$snv <- mat.snvs[genes,b]
mat_list.barr$indel <- mat.indels[genes,b]
mat.cn.amp[genes,b][is.na(mat.cn.amp[genes,b])] <- 0
mat_list.barr$amp <- mat.cn.amp[genes,b]
mat.cn.del[genes,b][is.na(mat.cn.del[genes,b])] <- 0
mat_list.barr$del <- mat.cn.del[genes,b]
names(mat_list.barr) <- c("SNV","Indel","AMP","DEL")

altered_nums<- Reduce( "+" , lapply(mat_list.barr, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.barr.sub<-  lapply(mat_list.barr, function(x) x[slice_indx,])

library(ComplexHeatmap)

alter_fun_list1 = list(
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#445E93", col = NA))
  },
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#7EB2DD", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#F93943", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#FFB563", col = NA))
  }
)

tb.barr <- sigs.allSamples[b,]
pdf("plots_4/oncoprint.drivers.barretts.pdf",w=15)
oncoPrint(mat_list.barr.sub, alter_fun = alter_fun_list1,
          col = c(SNV = "#445E93", Indel = "#7EB2DD",
                            AMP="#F93943", DEL="#FFB563", LOH="#191D32"),
          top_annotation = HeatmapAnnotation(SBS17a = anno_barplot(tb.barr[,"SBS17a"],ylim=c(0,0.6)),
                                             SBS17b = anno_barplot(tb.barr[,"SBS17b"],ylim=c(0,0.6)),
                                             SBS2 = anno_barplot(tb.barr[,"SBS2"],ylim=c(0,0.6)),
                                             SBS3 = anno_barplot(tb.barr[,"SBS3"],ylim=c(0,0.6)),
                                             SBS8 = anno_barplot(tb.barr[,"SBS8"],ylim=c(0,0.6)),
                                             SBS41 = anno_barplot(tb.barr[,"SBS41"],ylim=c(0,0.6)),
                                             SBS44 = anno_barplot(tb.barr[,"SBS44"],ylim=c(0,0.6)),
                                             SBS40 = anno_barplot(tb.barr[,"SBS40"],ylim=c(0,0.6)))
          )
dev.off()

write.csv(mat_list.barr.sub, file="plots_4/mat_list.barr.sub.csv")


mat_list.primary <- list()
genes <- intersect(rownames(mat.snvs),rownames(mat.cn.amp))
o <- intersect(o,colnames(mat.cn.amp))
mat_list.primary$SNV <- mat.snvs[genes,o]
mat_list.primary$Indel <- mat.indels[genes,o]
mat.cn.amp[genes,o][is.na(mat.cn.amp[genes,o])] <- 0
mat_list.primary$AMP <- mat.cn.amp[genes,o]
mat.cn.del[genes,o][is.na(mat.cn.del[genes,o])] <- 0
mat_list.primary$DEL <- mat.cn.del[genes,o]

altered_nums<- Reduce( "+" , lapply(mat_list.primary, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.primary.sub<-  lapply(mat_list.primary, function(x) x[slice_indx,])

tb.primary <- sigs.allSamples[o,]
pdf("plots_4/oncoprint.drivers.primary.pdf",w=15)
oncoPrint(mat_list.primary.sub, alter_fun = alter_fun_list1,
          col = c(SNV = "#445E93", Indel = "#7EB2DD",
                  AMP="#F93943", DEL="#FFB563"),
          top_annotation = HeatmapAnnotation(SBS17a = anno_barplot(tb.primary[,"SBS17a"],ylim=c(0,0.6)),
                                             SBS17b = anno_barplot(tb.primary[,"SBS17b"],ylim=c(0,0.6)),
                                             SBS2 = anno_barplot(tb.primary[,"SBS2"],ylim=c(0,0.6)),
                                             SBS3 = anno_barplot(tb.primary[,"SBS3"],ylim=c(0,0.6)),
                                             SBS8 = anno_barplot(tb.primary[,"SBS8"],ylim=c(0,0.6)),
                                             SBS41 = anno_barplot(tb.primary[,"SBS41"],ylim=c(0,0.6)),
                                             SBS44 = anno_barplot(tb.primary[,"SBS44"],ylim=c(0,0.6)),
                                             SBS40 = anno_barplot(tb.primary[,"SBS40"],ylim=c(0,0.6)))
)
dev.off()

write.csv(mat_list.primary.sub, file="plots_4/mat_list.primary.sub.csv")


mat_list.mets <- list()
genes <- intersect(rownames(mat.snvs),rownames(mat.cn.amp))
m <- intersect(m,colnames(mat.cn.amp))
mat_list.mets$SNV <- mat.snvs[genes,m]
mat_list.mets$Indel <- mat.indels[genes,m]
mat.cn.amp[genes,m][is.na(mat.cn.amp[genes,m])] <- 0
mat_list.mets$AMP <- mat.cn.amp[genes,m]
mat.cn.del[genes,m][is.na(mat.cn.del[genes,m])] <- 0
mat_list.mets$DEL <- mat.cn.del[genes,m]

altered_nums<- Reduce( "+" , lapply(mat_list.mets, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.mets.sub<-  lapply(mat_list.mets, function(x) x[slice_indx,])


tb.met <- sigs.allSamples[m,]

pdf("plots_4/oncoprint.drivers.mets.pdf",w=15)
oncoPrint(mat_list.mets.sub, alter_fun = alter_fun_list1,
          col = c(SNV = "#445E93", Indel = "#7EB2DD",
                  AMP="#F93943", DEL="#FFB563"),
          top_annotation = HeatmapAnnotation(SBS17a = anno_barplot(tb.met[,"SBS17a"],ylim=c(0,0.6)),
                                             SBS17b = anno_barplot(tb.met[,"SBS17b"],ylim=c(0,0.6)),
                                             SBS2 = anno_barplot(tb.met[,"SBS2"],ylim=c(0,0.6)),
                                             SBS3 = anno_barplot(tb.met[,"SBS3"],ylim=c(0,0.6)),
                                             SBS8 = anno_barplot(tb.met[,"SBS8"],ylim=c(0,0.6)),
                                             SBS41 = anno_barplot(tb.met[,"SBS41"],ylim=c(0,0.6)),
                                             SBS44 = anno_barplot(tb.met[,"SBS44"],ylim=c(0,0.6)),
                                             SBS40 = anno_barplot(tb.met[,"SBS40"],ylim=c(0,0.6)))
          #col = list(sigs=c("SBS17a" = "#961D4E", "SBS5" = "#A60067")))
          #                           "SBS2" = "#6153CC","SBS3"="#6C91C2","SBS8"="#C3C9E9",
          #                           "SBS44"="#44FFD1",
          #                           "SBS41"="#8B8982","SBS18"="#373F47","SBS35"="#F1D6B8",
          #                           "SBS28"="#FBACBE",
          #                           "SBS1"="#360A14","SBS5"="#B9D2B1","SBS40"="#A18276")))
)
dev.off()

write.csv(mat_list.mets.sub, file="plots_4/mat_list.mets.sub.csv")


### Now summarise analyses:

## Barretts:
mb <- mat_list.barr$SNV+mat_list.barr$Indel+mat_list.barr$AMP+mat_list.barr$DEL
mb[mb!=0] <- 1

sigs.p <- NULL
sigs.n <- NULL
m.assoc <- array(0,c(nrow(mb),ncol(sigs.allSamples)))
rownames(m.assoc) <- rownames(mb)
colnames(m.assoc) <- colnames(sigs.allSamples)
p.assoc <- array(1,c(nrow(mb),ncol(sigs.allSamples)))
rownames(p.assoc) <- rownames(mb)
colnames(p.assoc) <- colnames(sigs.allSamples)
for (g in rownames(mb)) {
  for (s in 1:ncol(sigs.allSamples)) {
    if (length(which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))))>0) {
      c <- wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))),s],
                       sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==0))),s])
      sigs.p <- c(sigs.p, c$p.value)
      sigs.n <- c(sigs.n, paste0(g,":",colnames(sigs.allSamples)[s]))
      p.assoc[g,colnames(sigs.allSamples)[s]] <- c$p.value
      m.assoc[g,colnames(sigs.allSamples)[s]] <- median(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))),s])-
        median(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==0))),s])
    }
    else {
      sigs.p <- c(sigs.p, 1)
      sigs.n <- c(sigs.n, paste0(g,":",colnames(sigs.allSamples)[s]))
    }
  }
}

padj.p <- p.adjust(sigs.p, method="BH")
padj.p[padj.p<0.1]
sigs.p[sigs.p<0.05]
sigs.n[sigs.p<0.05]
sigs.n[padj.p<0.05]
p.assoc.adj <- apply(p.assoc,2,p.adjust,method="fdr", n = length(p.assoc))#p.adjust(p.assoc, method="BH")

m.assoc[is.na(m.assoc)] <- 0
m.assoc[is.infinite(m.assoc)] <- 0
m.assoc[p.assoc>=0.05] <- 0

m.assoc.red <- m.assoc[which(rowSums(abs(m.assoc))>0.1),]
m.assoc.red2 <- m.assoc.red[,which(colSums(abs(m.assoc.red))>0.1)]

paletteLength <- 20
myColor <- colorRampPalette(c("#EBD2BE", "white", "#7F055F"))(paletteLength)
myBreaks <- c(seq(min(m.assoc.red2), 0, length.out=ceiling(paletteLength/2) ), 
              seq(0.01, max(m.assoc.red2), length.out=floor(paletteLength/2)))

pdf("plots_4/assoc_withdrivers.Barretts4.pdf",w=5,h=3)
pheatmap(t(m.assoc.red2), color=myColor, breaks = myBreaks)
dev.off()

write.csv(m.assoc.red2, file="plots_4/m.assoc.red2.barretts.csv")
 
## Primaries:
mb <- mat_list.primary$SNV+mat_list.primary$Indel+mat_list.primary$AMP+mat_list.primary$DEL
mb[mb!=0] <- 1

sigs.p <- NULL
sigs.n <- NULL
m.assoc <- array(0,c(nrow(mb),ncol(sigs.allSamples)))
rownames(m.assoc) <- rownames(mb)
colnames(m.assoc) <- colnames(sigs.allSamples)
p.assoc <- array(1,c(nrow(mb),ncol(sigs.allSamples)))
rownames(p.assoc) <- rownames(mb)
colnames(p.assoc) <- colnames(sigs.allSamples)
for (g in rownames(mb)) {
  for (s in 1:ncol(sigs.allSamples)) {
    if (length(which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))))>0) {
      c <- wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))),s],
                       sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==0))),s])
      sigs.p <- c(sigs.p, c$p.value)
      sigs.n <- c(sigs.n, paste0(g,":",colnames(sigs.allSamples)[s]))
      p.assoc[g,colnames(sigs.allSamples)[s]] <- c$p.value
      m.assoc[g,colnames(sigs.allSamples)[s]] <- median(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))),s])-
        median(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==0))),s])
    }
    else {
      sigs.p <- c(sigs.p, 1)
      sigs.n <- c(sigs.n, paste0(g,":",colnames(sigs.allSamples)[s]))
    }
  }
}

padj.p <- p.adjust(sigs.p, method="BH")
padj.p[padj.p<0.1]
sigs.p[sigs.p<0.05]
sigs.n[sigs.p<0.05]
sigs.n[padj.p<0.05]
p.assoc.adj <- apply(p.assoc,2,p.adjust,method="BH", n = length(p.assoc))#p.adjust(p.assoc, method="BH")
#p.adjust(p.assoc, method="BH")

m.assoc[is.na(m.assoc)] <- 0
m.assoc[is.infinite(m.assoc)] <- 0
m.assoc[p.assoc.adj>=0.05] <- 0

m.assoc.red <- m.assoc[which(rowSums(abs(m.assoc))>0.1),]
m.assoc.red2 <- m.assoc.red[,which(colSums(abs(m.assoc.red))>0.1)]

paletteLength <- 20
myColor <- colorRampPalette(c("#EBD2BE", "white", "#7F055F"))(paletteLength)
myBreaks <- c(seq(min(m.assoc.red2), 0, length.out=ceiling(paletteLength/2) ), 
              seq(0.01, max(m.assoc.red2), length.out=floor(paletteLength/2)))

pdf("plots_4/assoc_withdrivers.Primaries4.pdf",w=4,h=2.5)
pheatmap(t(m.assoc.red2), color=myColor, breaks=myBreaks)
dev.off()

write.csv(m.assoc.red2, file="plots_4/m.assoc.red2.primaries.csv")


## Mets:
mb <- mat_list.mets$SNV+mat_list.mets$Indel+mat_list.mets$AMP+mat_list.mets$DEL
mb[mb!=0] <- 1

sigs.p <- NULL
sigs.n <- NULL
m.assoc <- array(0,c(nrow(mb),ncol(sigs.allSamples)))
rownames(m.assoc) <- rownames(mb)
colnames(m.assoc) <- colnames(sigs.allSamples)
p.assoc <- array(1,c(nrow(mb),ncol(sigs.allSamples)))
rownames(p.assoc) <- rownames(mb)
colnames(p.assoc) <- colnames(sigs.allSamples)
for (g in rownames(mb)) {
  for (s in 1:ncol(sigs.allSamples)) {
    if (length(which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))))>0) {
      c <- wilcox.test(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))),s],
                       sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==0))),s])
      sigs.p <- c(sigs.p, c$p.value)
      sigs.n <- c(sigs.n, paste0(g,":",colnames(sigs.allSamples)[s]))
      p.assoc[g,colnames(sigs.allSamples)[s]] <- c$p.value
      m.assoc[g,colnames(sigs.allSamples)[s]] <- median(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==1))),s])-
        median(sigs.allSamples[which(rownames(sigs.allSamples) %in% names(which(mb[g,]==0))),s])
    }
    else {
      sigs.p <- c(sigs.p, 1)
      sigs.n <- c(sigs.n, paste0(g,":",colnames(sigs.allSamples)[s]))
    }
  }
}

padj.p <- p.adjust(sigs.p, method="BH")
padj.p[padj.p<0.1]
sigs.p[sigs.p<0.05]
sigs.n[sigs.p<0.05]
sigs.n[padj.p<0.05]
p.assoc.adj <- apply(p.assoc,2,p.adjust,method="BH", n = length(p.assoc))#p.adjust(p.assoc, method="BH")

#p.assoc.adj <- p.adjust(p.assoc, method="BH")

m.assoc[is.na(m.assoc)] <- 0
m.assoc[is.infinite(m.assoc)] <- 0
m.assoc[p.assoc.adj>=0.05] <- 0

m.assoc.red <- m.assoc[which(rowSums(abs(m.assoc))>0.1),]
m.assoc.red2 <- m.assoc.red[,which(colSums(abs(m.assoc.red))>0.1)]

paletteLength <- 20
myColor <- colorRampPalette(c("#EBD2BE", "white", "#7F055F"))(paletteLength)
myBreaks <- c(seq(min(m.assoc.red2), 0, length.out=ceiling(paletteLength/2) ), 
              seq(0.01, max(m.assoc.red2), length.out=floor(paletteLength/2)))

pdf("plots_4/assoc_withdrivers.Mets4.pdf",w=5,h=2.5)
pheatmap(t(m.assoc.red2), color=myColor, breaks = myBreaks)
dev.off()

write.csv(m.assoc.red2, file="plots_4/m.assoc.red2.mets.csv")





