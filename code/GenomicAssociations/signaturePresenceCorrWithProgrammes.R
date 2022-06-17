#############################
##### Associate the presence/absence of signatures with expression:

library(reshape)

## Catalogue cases with increase at bottleneck vs decrease at bottleneck:
load("processeddata/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

sigs.allSamples$Sample <- sapply(rownames(sigs.allSamples), 
                                 function(x) strsplit(x,"_vs_")[[1]][1])
df.sigs <- melt(sigs.allSamples)
colnames(df.sigs)[2:3] <- c("Signature","Exposure")
df.sigs$ExposureBinary <- sapply(df.sigs$Exposure,
                                 function(x) ifelse(x>0.05,"yes","no"))
table(df.sigs[,c("Signature","ExposureBinary")])


# Next, load expression data and map against genomic changes:
load("../tracksig_barrettsmetsprimaries/oac.expr.final.RData")
load("../tracksig_barrettsmetsprimaries/mat.expr.RData")

cbioportal_invasion <- read.table("cbioportal_invasionmetastasis.txt",
                                  header=FALSE, stringsAsFactors = FALSE)$V1

rownames(mat.expr) <- mat.expr$DNAid 
mat.expr <- mat.expr[,-1]

## Next, read all cbioportal programmes:
cbio <- read.xls("cbioportal_AllProgrammes.xlsx")
cbiolist <- NULL
i<-0
for (c in colnames(cbio)) {
  i<-i+1
  cbiolist[[i]] <- unique(setdiff(cbio[,c],""))
}
names(cbiolist) <- colnames(cbio)


library(ConsensusTME)

mat.transf <- t(as.matrix(mat.expr))
rownames(mat.transf) <- colnames(mat.expr)
colnames(mat.transf) <- rownames(mat.expr)

celldeconv <- consensusTMEAnalysis(as.matrix(mat.transf), 
                                                 cancer = "ESCA", 
                                   statMethod = "gsva")
df.celldeconv <- data.frame(t(celldeconv))
df.celldeconv$Sample <- rownames(df.celldeconv)
save(df.celldeconv, file="df.expr.celldeconv.RData")

df.sigs.plusdeconv <- merge(df.sigs,
                            df.celldeconv, 
                            by.x="Sample", by.y="Sample",
                            all.x=FALSE, all.y=FALSE)

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
save(df.hallmark, file="df.expr.hallmark.RData")

df.sigs.plushallmark <- merge(df.sigs,
                                 df.hallmark, 
                                 by.x="Sample", by.y="Sample",
                                 all.x=FALSE, all.y=FALSE)

# cbio scores:
cbioscores <- gsva(as.matrix(mat.transf),cbiolist, 
                       method="gsva", kcdf="Gaussian",
                       mx.diff=TRUE, abs.ranking=FALSE)
df.cbio <- data.frame(t(cbioscores))
df.cbio$Sample <- rownames(df.cbio)
save(df.cbio, file="df.expr.cbio.RData")

df.sigs.pluscbio <- merge(df.sigs,
                              df.cbio, 
                              by.x="Sample", by.y="Sample",
                              all.x=FALSE, all.y=FALSE)

### Iterate for each signature
for (s in unique(df.sigs$Signature)) {
  df.sigs.plusdeconv.selected <- df.sigs.plusdeconv[which(df.sigs.plusdeconv$Signature == s),]
  
  df.melt.selected.immune <- melt(df.sigs.plusdeconv.selected,
                       id.vars = c("Sample", "Signature", 
                                   "Exposure", "ExposureBinary"))
  
  df.sigs.plushallmark.selected <- df.sigs.plushallmark[which(df.sigs.plushallmark$Signature == s),]
  
  df.melt.selected.hallmark <- melt(df.sigs.plushallmark.selected,
                      id.vars = c("Sample", "Signature", 
                                  "Exposure", "ExposureBinary"))
  
  ### Plot immune infiltrates compared:
  my_comparisons <- list(c("yes","no"))
  pdf(paste0("plots.signatureCorrs/",s,".TMEcompared.gsva.pdf"),w=10,h=8)
  print(ggviolin(df.melt.selected.immune, 
           x = "ExposureBinary", y = "value",
           palette ="npg",
           fill = "ExposureBinary",alpha = 1,#color = NA,
           add = "boxplot", add.params = list(fill = "white"))+
    facet_wrap(~variable, scale="free")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank())+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
    xlab("")+
    ylab("Infiltration score"))
  dev.off()
  
  ### Plot hallmarks compared:
  pdf(paste0("plots.signatureCorrs/",s,".HallmarksCompared.gsva.pdf"),w=10,h=8)
  print(ggviolin(df.melt.selected.hallmark, 
           x = "ExposureBinary", y = "value",
           palette ="npg",
           fill = "ExposureBinary",alpha = 1,#color = NA,
           add = "boxplot", add.params = list(fill = "white"))+
    facet_wrap(~variable, scale="free")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank())+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
    xlab("")+
    ylab("Infiltration score"))
  dev.off()

}

### Iterate for each signature for cbio:
for (s in unique(df.sigs$Signature)) {
  df.sigs.pluscbio.selected <- df.sigs.pluscbio[which(df.sigs.pluscbio$Signature == s),]
  
  df.melt.selected.cbio <- melt(df.sigs.pluscbio.selected,
                                    id.vars = c("Sample", "Signature", 
                                                "Exposure", "ExposureBinary"))
  
  ### Plot immune infiltrates compared:
  my_comparisons <- list(c("yes","no"))
  pdf(paste0("plots.signatureCorrs/",s,".CBIOcompared.gsva.pdf"),w=10,h=8)
  print(ggviolin(df.melt.selected.cbio, 
                 x = "ExposureBinary", y = "value",
                 palette ="npg",
                 fill = "ExposureBinary",alpha = 1,#color = NA,
                 add = "boxplot", add.params = list(fill = "white"))+
          facet_wrap(~variable, scale="free")+
          theme(axis.text.x=element_blank(),
                axis.ticks.x = element_blank())+
          stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
          xlab("")+
          ylab("Activity score"))
  dev.off()
  
  ## Plot also correlations:
  pdf(paste0("plots.signatureCorrs/",s,".CBIOscatter.gsva.pdf"),w=10,h=8)
  print(ggscatter(df.melt.selected.cbio, x = "Exposure", y = "value",
                  color = "variable", shape = 21, size = 3, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.x = 0, label.sep = "\n")
  )+xlab("Signature exposure")+ylab("Activity score")+
    facet_wrap(~variable,scale="free",nrow=5))
  dev.off()
  
  
}

