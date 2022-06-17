########################
### Sigs compared plus heatmap

library(reshape)
library(ggpubr)
library(pheatmap)

load("processeddata/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
load("processeddata/sigs.allSamples.absolute.202101127_cosmicref13sig_S3S8.RData")
load("processeddata/annotation.sampleIDs.RData")

## Merge signatures and IDs:

sigs.allSamples <- data.frame(sigs.allSamples)
sigs.allSamples$Sample <- rownames(sigs.allSamples)

sigs <- merge(sigs.allSamples, annotation.sampleIDs,
              by.x="Sample", by.y ="Sample", all.x=FALSE, all.y=FALSE)
df.sigs <- melt(sigs)
df.sigs$Category <- factor(df.sigs$Category, levels=c("Barretts","PrimaryTumour",
                                                  "LymphNode","Metastasis"))
df.sigs$variable <- factor(df.sigs$variable,
                           levels=c("SBS17a","SBS17b","SBS2","SBS3","SBS8",
                                    "SBS18","SBS30","SBS44","SBS35",
                                    "SBS41","SBS28",
                                    "SBS1","SBS5","SBS40"))

sigs.allSamples.absolute <- data.frame(sigs.allSamples.absolute)
sigs.allSamples.absolute$Sample <- rownames(sigs.allSamples.absolute)

sigs.abs <- merge(sigs.allSamples.absolute, annotation.sampleIDs,
              by.x="Sample", by.y ="Sample", all.x=FALSE, all.y=FALSE)
df.sigs.abs <- melt(sigs.abs)
df.sigs.abs$Category <- factor(df.sigs.abs$Category, levels=c("Barretts","PrimaryTumour",
                                                      "LymphNode","Metastasis"))
df.sigs.abs$variable <- factor(df.sigs.abs$variable,
                           levels=c("SBS17a","SBS17b","SBS2","SBS3","SBS8",
                                    "SBS18","SBS30","SBS44","SBS35",
                                    "SBS41","SBS28",
                                    "SBS1","SBS5","SBS40"))
df.sigs.abs$log10value <- log10(df.sigs.abs$value+1)

### Now infer the trios:

tb <- table(annotation.sampleIDs$OCCAMS_ID)

pairs <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==2])),]
trios <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==3])),]
singles <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==1])),]
save(pairs, file="pairs.RData")

multi1 <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==16])),]
multi2 <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==17])),]

pairs.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% names(tb[tb==2])),]
trios.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% names(tb[tb==3])),]
singles.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% names(tb[tb==1])),]

multi1.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% names(tb[tb==16])),]
multi2.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% names(tb[tb==17])),]


my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","LymphNode"))
pdf("plots.sigsCompared/paired.pdf",w=18,h=10)
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
pdf("plots.sigsCompared/paired.BarrettsPrimaries.pdf",w=12,h=4)
ggpaired(pairs[which((pairs$Category != "LymphNode")&
                 (pairs$variable %in% 
                    c("SBS2","SBS3","SBS18","SBS35",
                      "SBS1","SBS5","SBS40"))),], 
         x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisonsbp, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()


pdf("plots.sigsCompared/paired.BarrettsPrimaries.absolute.pdf",w=12,h=4)
ggpaired(pairs.abs[which((pairs.abs$Category != "LymphNode")&
                       (pairs.abs$variable %in% 
                          c("SBS2","SBS3","SBS18","SBS35",
                            "SBS1","SBS5","SBS40"))),], 
         x = "Category", y = "log10value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = get_palette(c("#6ba393", "#e8b44d"),2), id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisonsbp, label = "p.signif")+
  xlab("")+
  ylab("Total mutations (log10)")
dev.off()

pdf("plots.sigsCompared/paired.absolute.pdf",w=18,h=4)
ggpaired(pairs.abs, x = "Category", y = "log10value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  xlab("")+
  ylab("Total mutations (log10)")
dev.off()

pdf("plots.sigsCompared/trios.pdf",w=18,h=4)
ggpaired(trios, x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  xlab("")+
  ylab("% contribution")
dev.off()

### Need to add metastases and lymph nodes:
singles.allmets <- df.sigs[which(!(df.sigs$Sample %in% pairs$Sample)),]
singles.allmets[which(singles.allmets$Category == "LymphNode"),]$Category <- "Metastasis"
my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","Metastasis"))
pdf("plots.sigsCompared/singles.pdf",w=12,h=6)
ggviolin(singles.allmets, x = "Category", 
         y = "value",
         palette ="npg",
         fill = "Category",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

singles.allmets.abs <- df.sigs.abs[which(!(df.sigs.abs$Sample %in% pairs.abs$Sample)),]
singles.allmets.abs[which(singles.allmets.abs$Category == "LymphNode"),]$Category <- "Metastasis"
my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","Metastasis"))
pdf("plots.sigsCompared/singles.absolute.pdf",w=12,h=6)
ggviolin(singles.allmets.abs, x = "Category", 
         y = "log10value",
         palette =get_palette(c("#6ba393", "#e8b44d","#94788c"),3),
         fill = "Category",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Total mutations (log10)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.sigsCompared/singles.density.pdf",w=18,h=20)
ggdensity(singles.allmets, x = "value",
         add = "mean",rug=TRUE,
         palette ="npg",
         color = "Category",fill="Category",alpha = 0.7)+
  facet_wrap(~variable,scales = "free",nrow=7)+
  xlab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()


############################################
######### The next plots are not so useful

# Finally, make some heat maps for the multi-samples:
multi1[1:4,]

pdf("plots.sigsCompared/multi1.pdf",w=10,h=3)
multi1 %>%
group_by(Category) %>%
  mutate(rescale = scales::rescale(value)) %>%
  ggplot(., aes(x = factor(variable), y = Sample)) +
  geom_tile(aes(alpha = rescale,fill=Category), color = "white") +
  scale_alpha(range = c(0.1, 1)) +
    geom_text(aes(label=round(value,2)))+
  theme(
        #axis.text.y = element_text(color = scales::hue_pal()(length(levels(multi1$Category)))),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "", y = "Samples")
dev.off()

pdf("plots.sigsCompared/multi2.pdf",w=10,h=3)
multi2 %>%
  group_by(Category) %>%
  mutate(rescale = scales::rescale(value)) %>%
  ggplot(., aes(x = factor(variable), y = Sample)) +
  geom_tile(aes(alpha = rescale,fill=Category), color = "white") +
  scale_alpha(range = c(0.1, 1)) +
  geom_text(aes(label=round(value,2)))+
  theme(
        #axis.text.y = element_text(color = scales::hue_pal()(length(levels(multi1$Category)))),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "", y = "Samples")
dev.off()

## Finally, heatmap for all three categories:
sigs.keep <- sigs[which(!(sigs$Sample %in% pairs$Sample)),]

pdf("plots.sigsCompared/heatmap.barretts.pdf",w=6,h=3)
pheatmap(t(sigs.keep[which(sigs.keep$Category == "Barretts"),2:14]),
         show_colnames = FALSE)
dev.off()

pdf("plots.sigsCompared/heatmap.primaries.pdf",w=6,h=3)
pheatmap(t(sigs.keep[which(sigs.keep$Category == "PrimaryTumour"),2:14]),
         show_colnames = FALSE)
dev.off()

pdf("plots.sigsCompared/heatmap.mets.pdf",w=6,h=3)
pheatmap(t(sigs.keep[which(sigs.keep$Category %in% c("Metastasis","LymphNode")),2:14]),
         show_colnames = FALSE)
dev.off()

pdf("plots.sigsCompared/singles.density.pdf",w=18,h=20)
ggdensity(singles.allmets, x = "value",
          add = "mean",rug=TRUE,
          palette ="npg",
          color = "Category",fill="Category",alpha = 0.7)+
  facet_wrap(~variable,scales = "free",nrow=7)+
  xlab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()


#########################
### Compare chemo-naive and treated samples:

# Clinical data:
clinical <- read.csv("Sigs study cohort_clinical demographics.csv")
clinical$Chemotherapy.Treatment <- "Other"
clinical[which(clinical$TR.Chemotherapy.Treatment.Protocol %in% c("Cisplatin 5FU",
                                                                  "Cisplatin/5 FU",
                                                                  "Cisplatin/5-FU",
                                                                  "Cispltin + 5-FU",
                                                                  "5FU/Cisplatin",
                                                                  "C/5'FU","Cis/5-FU",
                                                                  "Cis/5FU","CIS/5FU",
                                                                  "CISP/5FU","Cisplatin + 5-FU",
                                                                  "Cisplatin and 5FU","Cisplatin/5'FU",
                                                                  "Cisplatin/5FU")),]$Chemotherapy.Treatment <- "Cisplatin/5FU"

clinical[which(clinical$TR.Chemotherapy.Treatment.Protocol %in% c("Cisplatin & Capecitabine",
                                                                  "Cisplatin/Capecitabine",
                                                                  "Cisplatin/Capecitabine (CX)",
                                                                  "Concurrent Cisplatin and Cabecitabine")),]$Chemotherapy.Treatment <- "Cisplatin/Capecitabine"

clinical[which(grepl("ECX",clinical$TR.Chemotherapy.Treatment.Protocol)|grepl("ECARBOX",clinical$TR.Chemotherapy.Treatment.Protocol)|
                 (clinical$TR.Chemotherapy.Treatment.Protocol %in% c("E Carbo X ","E-Carbo-X","ECarboX","(E)CX"))),]$Chemotherapy.Treatment <- "ECX"

clinical[which(grepl("EOX",clinical$TR.Chemotherapy.Treatment.Protocol)),]$Chemotherapy.Treatment <- "EOX"

df.sigs <- merge(df.sigs,clinical[which(clinical$Is.Chemo.Naive!=""),
                                  c("Sample","Is.Chemo.Naive","TR.Chemotherapy.Number.Of.Cycles.Given",
                                    "TR.Chemotherapy.Treatment.Protocol","Chemotherapy.Treatment",
                                    "TR.Response")],
                               by.x="Sample",by.y="Sample",
                               all.x=FALSE, all.y=FALSE)
df.sigs.abs <- merge(df.sigs.abs,clinical[which(clinical$Is.Chemo.Naive!=""),
                                          c("Sample","Is.Chemo.Naive","TR.Chemotherapy.Number.Of.Cycles.Given",
                                            "TR.Chemotherapy.Treatment.Protocol","Chemotherapy.Treatment",
                                            "TR.Response")],
                     by.x="Sample",by.y="Sample",
                     all.x=FALSE, all.y=FALSE)

my_comparisons.chemo <-list(c("true","false"))
pdf("plots.sigsCompared/chemoNaiveVsTreated_relative.pdf",w=12,h=6)
ggviolin(df.sigs, x = "Is.Chemo.Naive", 
         y = "value",
         #palette ="npg",
         fill = "Is.Chemo.Naive",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons.chemo, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.sigsCompared/chemoNaiveVsTreated_relative.selectedSignificant.pdf",w=4,h=4)
ggviolin(df.sigs[which(df.sigs$variable %in% c("SBS35","SBS44")),], 
         x = "Is.Chemo.Naive", 
         y = "value",
         #palette ="npg",
         fill = "Is.Chemo.Naive",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons.chemo, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.sigsCompared/chemoNaiveVsTreated_absolute.pdf",w=12,h=6)
ggviolin(df.sigs.abs, x = "Is.Chemo.Naive", 
         y = "log10value",
         #palette ="npg",
         fill = "Is.Chemo.Naive",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons.chemo, label = "p.signif")+
  xlab("")+
  ylab("Total mutations (log10)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

# Treatment protocol comparison:

my_comparisons.protocol <- list(c("Cisplatin/5FU","ECX"),
                                c("ECX","EOX"),
                                c("EOX","Other"))
pdf("plots.sigsCompared/chemoTreated_treatmentProtocolCompared.pdf",w=12,h=6)
ggviolin(df.sigs[which(df.sigs$Is.Chemo.Naive == "false"),], 
         x = "Chemotherapy.Treatment", 
         y = "value",
         #palette ="npg",
         fill = "Chemotherapy.Treatment",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons.protocol, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

# Cis/5fu vs rest:
df.sigs$Platinum <- "Other"
df.sigs[which(df.sigs$Chemotherapy.Treatment=="Cisplatin/5FU"),]$Platinum <-"Cis/5FU" 
my_comparisons.protocol <- list(c("Cis/5FU","Other"))
pdf("plots.sigsCompared/chemoTreated_treatmentProtocolComparedPlatinumOnly.pdf",w=12,h=6)
ggviolin(df.sigs[which(df.sigs$Is.Chemo.Naive == "false"),], 
         x = "Platinum", 
         y = "value",
         #palette ="npg",
         fill = "Platinum",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons.protocol, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

# Finally, compare treatment response:
my_comparisons.response <- list(
                                c("PR","SD")) # 19 PR, 54 SD
df.sigs$TR.Response <- factor(df.sigs$TR.Response,
                              levels = c("CR","PR","PD","SD","unknown"))
pdf("plots.sigsCompared/chemoTreated_treatmentResponse.pdf",w=12,h=6)
ggviolin(df.sigs[which((df.sigs$Is.Chemo.Naive == "false") &
                         (df.sigs$TR.Response %in% c("PR","SD"))),], 
         x = "TR.Response", 
         y = "value",
         #palette ="npg",
         fill = "TR.Response",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons.response, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

### Fisher's exact test for BER:

df.d <- df.sigs[which(df.sigs$Is.Chemo.Naive == "false"),]
df.d$BER <- sapply(df.d[which(df.d$variable == "SBS30"),]$value,
                   function(x) ifelse(x<0.05,"Null",">0.05"))
df.ber <- unique(df.d[,c("BER","TR.Response","Sample")])

tb <- table(df.ber[,c("BER","TR.Response")])
fisher.test(tb[,c(2,4)])
# not significant; equal numbers of BER>0 and BER=0 between the groups

# Compare treatment cycles:
pdf("plots.sigsCompared/chemoTreated_TR.Chemotherapy.Number.Of.Cycles.Given.pdf",w=12,h=6)
ggviolin(df.sigs[which((df.sigs$Is.Chemo.Naive == "false")),], 
         x = "TR.Chemotherapy.Number.Of.Cycles.Given", 
         y = "value",
         #palette ="npg",
         fill = "TR.Chemotherapy.Number.Of.Cycles.Given",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  #stat_compare_means(comparisons = my_comparisons.response, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()


pdf("plots.sigsCompared/chemoTreated_scatter.Number.Of.Cycles.Given.pdf",w=12,h=6)
ggscatter(df.sigs[which((df.sigs$Is.Chemo.Naive == "false")),], 
         x = "TR.Chemotherapy.Number.Of.Cycles.Given", 
         y = "value",
         #palette ="npg",
         add = "reg.line",
         fill = "TR.Chemotherapy.Number.Of.Cycles.Given",conf.int = TRUE, # Add confidence interval
         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
         cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  #stat_compare_means(comparisons = my_comparisons.response, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

#############
### Checking levels of signatures by category:
df.barr <- df.sigs[which(df.sigs$Category == "Barretts"),]
quantile(df.barr[which(df.barr$variable=="SBS28"),]$value)
quantile(df.barr[which(df.barr$variable=="SBS35"),]$value)
quantile(df.barr[which(df.barr$variable=="SBS17a"),]$value)
quantile(df.barr[which(df.barr$variable=="SBS3"),]$value)

## Set all sigs <5% to 0:
df.sigs[which(df.sigs$value<0.05),]$value <- 0

sink("sigsThresholds.txt")
for (s in unique(df.sigs$variable)) {
  print("======")
  print(s)
  df.sel <- df.sigs[which(df.sigs$variable==s),]
  print("Barretts")
  print(quantile(df.sel[which(df.sel$Category == "Barretts"),]$value))
  print("PT")
  print(quantile(df.sel[which(df.sel$Category == "PrimaryTumour"),]$value))
  print("Met")
  print(quantile(df.sel[which(df.sel$Category %in% c("LymphNode","Metastasis")),]$value))
}
sink()

length(which(df.sigs[which((df.sigs$variable=="SBS35")&
                             (df.sigs$Category=="PrimaryTumour")),]$value!=0))



######################
### Separation of chemo-naive and treated:
prim.naive <- clinical[which(clinical$Is.Chemo.Naive == "true"),]$Occams.Id
prim.treated <- clinical[which(clinical$Is.Chemo.Naive == "false"),]$Occams.Id

# read in the met data:
library(gdata)
mets <- read.xls("Copy of sujath_20220328_metastatic samples_chemoStatus_sample_ids_NGedits110422.xlsx")
table(mets$chemo_status)
mets.naive <- mets[which(mets$chemo_status == "Na\xefve"),]$occams_id
mets.treated <- mets[which(mets$chemo_status == "Treated"),]$occams_id

# replace OC with OCCAMS in df.sigs:
library(stringr)
df.sigs[which(grepl("OC/",df.sigs$OCCAMS_ID)),]$OCCAMS_ID <- str_replace(df.sigs[which(grepl("OC/",df.sigs$OCCAMS_ID)),]$OCCAMS_ID,
                                                                         "OC/","OCCAMS/")
annotation.sampleIDs[which(grepl("OC/",annotation.sampleIDs$OCCAMS_ID)),]$OCCAMS_ID <- str_replace(annotation.sampleIDs[which(grepl("OC/",annotation.sampleIDs$OCCAMS_ID)),]$OCCAMS_ID,
                                                                         "OC/","OCCAMS/")
df.sigs[which(df.sigs$Category == "LymphNode"),]$Category <- "Metastasis"

df.naive <- df.sigs[which(df.sigs$OCCAMS_ID %in% c(prim.naive,mets.naive)),]
df.treated <- df.sigs[which(df.sigs$OCCAMS_ID %in% c(prim.treated,mets.treated)),]
tb <- table(df.naive$OCCAMS_ID)

pairs <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==2])),]
trios <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==3])),]
singles <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==1])),]

my_comparisons <- list(
                       c("PrimaryTumour","Metastasis"))
pdf("plots.sigsCompared/all.Naive.PrimaryVsMets.2.pdf",w=12,h=6)
ggviolin(df.naive[which(df.naive$Category %in% c("PrimaryTumour","Metastasis")),],
         x = "Category", 
         y = "value",
         palette =get_palette(c( "#e8b44d","#94788c"),2),
         fill = "Category",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.sigsCompared/all.Treated.PrimaryVsMets.2.pdf",w=12,h=6)
ggviolin(df.treated[which(df.treated$Category %in% c("PrimaryTumour","Metastasis")),], 
         x = "Category", 
         y = "value",
         palette =get_palette(c( "#e8b44d","#94788c"),2),
         fill = "Category",alpha = 1,#color = NA,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.sigsCompared/paired.naive.PrimaryvsMet.2.pdf",w=12,h=6)
ggpaired(df.naive[which(df.naive$Category %in% c("PrimaryTumour","Metastasis")),], 
         x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

pdf("plots.sigsCompared/paired.treated.PrimaryvsMet.2.pdf",w=12,h=6)
ggpaired(df.treated[which(df.treated$Category!="Barretts"),], 
         x = "Category", y = "value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()


df.sigs.abs[which(grepl("OC/",df.sigs.abs$OCCAMS_ID)),]$OCCAMS_ID <- str_replace(df.sigs.abs[which(grepl("OC/",df.sigs.abs$OCCAMS_ID)),]$OCCAMS_ID,
                                                                         "OC/","OCCAMS/")
df.sigs.abs[which(df.sigs.abs$Category == "LymphNode"),]$Category <- "Metastasis"

df.naive.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% c(prim.naive,mets.naive)),]
df.treated.abs <- df.sigs.abs[which(df.sigs.abs$OCCAMS_ID %in% c(prim.treated,mets.treated)),]

pdf("plots.sigsCompared/paired.naive.PrimaryvsMet.abs.2.pdf",w=12,h=6)
ggpaired(df.naive.abs[which(df.naive.abs$Category %in% c("PrimaryTumour","Metastasis")),], 
         x = "Category", y = "log10value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

pdf("plots.sigsCompared/paired.treated.PrimaryvsMet.abs.2.pdf",w=12,h=6)
ggpaired(df.treated.abs[which(df.treated.abs$Category!="Barretts"),], 
         x = "Category", y = "log10value",
         color = "Category", line.color = "gray", line.size = 0.4,
         palette = "npg", id="OCCAMS_ID")+
  facet_wrap(.~variable, scale="free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("% contribution")
dev.off()

multi1 <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==16])),]
multi2 <- df.sigs[which(df.sigs$OCCAMS_ID %in% names(tb[tb==17])),]

pdf("plots.sigsCompared/multi1.pdf",w=10,h=3)
multi1 %>%
  group_by(Category) %>%
  mutate(rescale = scales::rescale(value)) %>%
  ggplot(., aes(x = factor(variable), y = Sample)) +
  geom_tile(aes(alpha = rescale,fill=Category), color = "white") +
  scale_alpha(range = c(0.1, 1)) +
  geom_text(aes(label=round(value,2)))+
  theme(
    #axis.text.y = element_text(color = scales::hue_pal()(length(levels(multi1$Category)))),
    axis.ticks = element_blank(),
    axis.text.y = element_blank()) +
  labs(x = "", y = "Samples")
dev.off()
