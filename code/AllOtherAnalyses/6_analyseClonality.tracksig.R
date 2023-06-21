#############################
##### Analyse clonality output from TrackSig.

library(reshape)
library(ggpubr)

# Load clonality reconstruction from TrackSig and annotation:
load("data/oac.sigClonality.final.may2021.RData")
load("data/annotation.sampleIDs.RData")

oac.sigClonality.annot <- merge(oac.sigClonality.final,
                                annotation.sampleIDs,
                                by.x="Sample", by.y="Sample",
                                all.x=FALSE, all.y=FALSE)

## Also assign number of signatures:
oac.sigClonality.annot$Clones <- 1
oac.sigClonality.annot[which(oac.sigClonality.annot$Sample %in%
                               names(which(table(oac.sigClonality.annot$Sample)==2))),]$Clones <- 2
oac.sigClonality.annot[which(oac.sigClonality.annot$Sample %in%
                               names(which(table(oac.sigClonality.annot$Sample)==3))),]$Clones <- 3
oac.sigClonality.annot[which(oac.sigClonality.annot$Sample %in%
                               names(which(table(oac.sigClonality.annot$Sample)==4))),]$Clones <- 4

# Change LN to mets:
oac.sigClonality.annot[which(oac.sigClonality.annot$Category == "LymphNode"),]$Category <- "Metastasis"

sigclon.melt <- melt(oac.sigClonality.annot,
                     id.vars = c("Sample","TumourID","OCCAMS_ID",
                                 "Category","CCF","Clones","ChangePoint"))
colnames(sigclon.melt)[8:9] <- c("Signature","Exposure")


## Next, compare clonal prevalences across cohorts:
sigclon.melt$Category <- factor(sigclon.melt$Category,
                              levels=c("Barretts","PrimaryTumour","Metastasis"))
my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","Metastasis"),
                       c("Barretts","Metastasis"))
colorschemeoac <- c("#6BA292",##E4D6A7", 
                    "#E9B44C", "#94778B")
sigclon.melt$Signature <- factor(sigclon.melt$Signature,
                                 levels=c("SBS17a",
                                          "SBS17b",
                                          "SBS2",
                                          "SBS3",
                                          "SBS8",
                                          "SBS30",
                                          "SBS18",
                                          "SBS44",
                                          "SBS28",
                                          "SBS35",
                                          "SBS1",
                                          "SBS5",
                                          "SBS40",
                                          "SBS41"))
pdf("plots_6/clonalPrevalence.primbarmet.pdf",w=20,h=8)
ggviolin(sigclon.melt[which(sigclon.melt$CCF==1),], 
          x = "Category", y = "Exposure",
          fill = "Category", palette = colorschemeoac,
          add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 1) +
  facet_wrap(~Signature,nrow=2,scale="free")
dev.off()


pdf("plots_6/subclonalPrevalence.primbarmet.pdf",w=20,h=8)
ggviolin(sigclon.melt[which(sigclon.melt$CCF!=1),], 
         x = "Category", y = "Exposure",
         fill = "Category", palette =colorschemeoac,
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 1) +
  facet_wrap(~Signature,scale="free",nrow=2)
dev.off()



sigclon.melt$Clonality <- sapply(sigclon.melt$CCF,
                                 function(x) ifelse(x==1,"Clonal","Subclonal"))

my_comparisons2 <- list(c("Clonal","Subclonal"))
pdf("plots_6/clonalVsSubclonal.all.pdf",w=12,h=8)
ggviolin(sigclon.melt, 
         x = "Clonality", y = "Exposure",
         fill = "Clonality", palette =c("#D7BBA8", "#9F4A54"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 1) +
  facet_wrap(~Signature,scale="free",nrow=2)
dev.off()

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


####### Bottlenecks:
bottles <- NULL
for (sample in unique(oac.sigClonality.annot$Sample)) {
  print(sample)
  current.oac <- oac.sigClonality.annot[which(oac.sigClonality.annot$Sample == sample),]
  if (nrow(current.oac)>1) {
    sortedccf <- rev(sort(current.oac$CCF))
    changes <- current.oac[which(current.oac$CCF==sortedccf[2]),2:15]-
      current.oac[which(current.oac$CCF==sortedccf[1]),2:15]
    bottles <- rbind(bottles,cbind(changes, sample,current.oac$Category))
  }
}
df.bottles <- melt(bottles, id.vars = c("sample","current.oac$Category"))
colnames(df.bottles) <- c("Sample","Category","Signature","Change")

library(pheatmap)
bottles <- unique(bottles)
pdf("plots_6/heatmap.bottlenecks.pdf",w=8,h=4)
annotc <- data.frame(XX=0,Category=bottles$`current.oac$Category`)
annotc$Category <- factor(annotc$Category)
rownames(annotc) <- bottles$sample
matx <- bottles[,c(1:13)]
rownames(matx) <- bottles$sample
pheatmap(t(matx),
         show_colnames = FALSE,
         annotation_col =annotc,
         cutree_cols = 4)
dev.off()

df.bottles$Category <- factor(df.bottles$Category,
                              levels=c("Barretts",
                                       "PrimaryTumour","Metastasis"))
df.bottles$Signature <- factor(df.bottles$Signature,
                     levels=c("SBS17a",
                              "SBS17b",
                              "SBS2",
                              "SBS3",
                              "SBS8",
                              "SBS30",
                              "SBS18",
                              "SBS44",
                              "SBS28",
                              "SBS35",
                              "SBS1",
                              "SBS5",
                              "SBS40",
                              "SBS41"))
load("data/df.bottles.tracksig.may2021.RData")
load("data/bottles.tracksig.may2021.RData")

# Changes during bottleneck:
pdf("plots_6/bottleneck.change.pdf",w=10,h=6)
ggboxplot(df.bottles, x = "Signature", y = "Change",
          color = "Category", palette = colorschemeoac)+
  geom_hline(yintercept=0, linetype="dotted", 
             color = "black", size=0.7)
  #facet_wrap(~Signature,nrow=3)
dev.off()

write.csv(unique(df.bottles), file="plots_6/df.bottles.csv")


##########################################
#### Finally, associate with expression:

df.bottles$Direction <- sapply(df.bottles$Change, function(x)
  ifelse(x<0,"Decrease",ifelse(x>0,"Increase","NoChange")))

load("../df.expr.hallmark.RData")
load("../df.expr.celldeconv.RData")
load("../df.expr.cbio.RData")
df.bottles$TumourID <- sapply(df.bottles$Sample, function(x) strsplit(x,"_vs_")[[1]][1])

df.sigs.plusdeconv <- merge(df.bottles[which(df.bottles$Category == "PrimaryTumour"),],
                            df.celldeconv,
                            by.x="TumourID",by.y="Sample",
                            all.x=FALSE, all.y=FALSE)
df.sigs.plushallmark <- merge(df.bottles[which(df.bottles$Category == "PrimaryTumour"),],
                            df.hallmark,
                            by.x="TumourID",by.y="Sample",
                            all.x=FALSE, all.y=FALSE)
df.sigs.pluscbio <- merge(df.bottles[which(df.bottles$Category == "PrimaryTumour"),],
                              df.cbio,
                              by.x="TumourID",by.y="Sample",
                              all.x=FALSE, all.y=FALSE)

### Run for SBS17b:
s="SBS17b"
  df.sigs.plusdeconv.selected <- df.sigs.plusdeconv[which(df.sigs.plusdeconv$Signature == s),]
  
  df.melt.selected.immune <- melt(df.sigs.plusdeconv.selected,
                                  id.vars = c("Sample","TumourID","Signature","Category", 
                                              "Change", "Direction"))
  
  df.sigs.plushallmark.selected <- df.sigs.plushallmark[which(df.sigs.plushallmark$Signature == s),]
  
  df.melt.selected.hallmark <- melt(df.sigs.plushallmark.selected,
                                    id.vars = c("Sample","TumourID","Signature","Category", 
                                                "Change", "Direction"))
  
  ### Plot immune infiltrates compared:
  my_comparisons <- list(c("Increase","Decrease"))
  pdf(paste0("plots_6/",s,".TMEcompared.gsva.pdf"),w=10,h=8)
  print(ggviolin(df.melt.selected.immune, 
                 x = "Direction", y = "value",
                 palette ="npg",
                 fill = "Direction",alpha = 1,#color = NA,
                 add = "boxplot", add.params = list(fill = "white"))+
          facet_wrap(~variable, scale="free")+
          theme(axis.text.x=element_blank(),
                axis.ticks.x = element_blank())+
          stat_compare_means(comparisons = my_comparisons)+
          xlab("")+
          ylab("Infiltration score"))
  dev.off()
  
  write.csv(df.melt.selected.immune, file="plots_6/sbs17b.tme.csv")

  ### Plot hallmarks compared:
  pdf(paste0("plots_6/",s,".HallmarksCompared.gsva.pdf"),w=10,h=8)
  print(ggviolin(df.melt.selected.hallmark, 
                 x = "Direction", y = "value",
                 palette ="npg",
                 fill = "Direction",alpha = 1,#color = NA,
                 add = "boxplot", add.params = list(fill = "white"))+
          facet_wrap(~variable, scale="free")+
          theme(axis.text.x=element_blank(),
                axis.ticks.x = element_blank())+
          stat_compare_means(comparisons = my_comparisons)+
          xlab("")+
          ylab("Activity score"))
  dev.off()

write.csv(df.melt.selected.hallmark, file="plots_6/sbs17b.hallmarks.csv")



  df.sigs.pluscbio.selected <- df.sigs.pluscbio[which(df.sigs.pluscbio$Signature == s),]
  
  df.melt.selected.cbio <- melt(df.sigs.pluscbio.selected,
                                    id.vars = c("Sample","TumourID","Signature","Category", 
                                                "Change", "Direction"))
  
  ### Plot immune infiltrates compared:
  my_comparisons <- list(c("Increase","Decrease"))
  pdf(paste0("plots_6/",s,".CBIOcompared.gsva.pdf"),w=10,h=8)
  print(ggviolin(df.melt.selected.cbio, 
                 x = "Direction", y = "value",
                 palette ="npg",
                 fill = "Direction",alpha = 1,#color = NA,
                 add = "boxplot", add.params = list(fill = "white"))+
          facet_wrap(~variable, scale="free")+
          theme(axis.text.x=element_blank(),
                axis.ticks.x = element_blank())+
          stat_compare_means(comparisons = my_comparisons)+
          xlab("")+
          ylab("Activity score"))
  dev.off()
  
  write.csv(df.melt.selected.cbio, file="plots_6/sbs17b.cbio.csv")
  
  
