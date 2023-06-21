##########
## Generic DDR patterns.

library(pheatmap)
library(reshape)

load("data/df.pathwayAlterationCounts.RData")


### Check co-occurrence and mutual exclusivity, removing LOH:
df <- df[which(!(df$Type %in% c("LOH"))),]

#df <- df[which(df$Type %in% c("DEL")),]
tb <- table(df[,c(2,1)])
tb[tb>0] <- 1

pheatmap(tb,show_colnames = FALSE)

### Check only samples in BER group:
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples[sigs.allSamples<0.05] <- 0
ber <- rownames(sigs.allSamples[sigs.allSamples$SBS30>0.05,])


df.ber <- df[which(df$Sample %in% ber),]
tbber <- table(df.ber[,c(2,1)])
tbber[tbber>0] <- 1

pheatmap(tbber,show_colnames = FALSE)


# Next, load expression data and map against genomic changes:
load("../../tracksig_barrettsmetsprimaries/oac.expr.final.RData")
load("../../tracksig_barrettsmetsprimaries/mat.expr.RData")

library(gdata)
ddr <- read.xls("~/Desktop/UCL/supervisedProjects/DanielJacobson/Rotation1-master/DDRpathways.xlsx")


ddrlist <- NULL
i<-0
for (c in unique(ddr$Pathway.1)) {
  i<-i+1
  ddrlist[[i]] <- unique(ddr[which(ddr$Pathway.1 == c),]$Gene.ID)
}
names(ddrlist) <- unique(ddr$Pathway.1)

library(GSVA)
library(reshape)
m <- t(mat.expr[,-1])
rownames(m) <- colnames(mat.expr[,-1])
ddrscores <- gsva(as.matrix(m),ddrlist, 
                  method="gsva", kcdf="Gaussian",
                  mx.diff=TRUE, abs.ranking=FALSE)
colnames(ddrscores) <- mat.expr$DNAid

df.ddrscores <- melt(ddrscores)
colnames(df.ddrscores) <- c("Pathway","Sample","Score")

df.tb <- data.frame(tb)
df.tb$SampleID <- sapply(df.tb$Sample, function(x) strsplit(as.character(x),"_vs_")[[1]][1])

df.ddrscores.m <- merge(df.ddrscores,
                      df.tb,
                      by.x=c("Sample","Pathway"),by.y = c("SampleID","Pathway.1"),
                      all.x=FALSE,all.y=FALSE)

## try beanplots in ggplot()
## http://www.jstatsoft.org/v28/c01/paper
library(ggplot2)
library(beanplot)
data("singer", package="lattice")
## figure 4 from paper using beanplot():
beanplot(Score~Pathway, data=df.ddrscores.m, side="both", border=NA, 
         col=list("black",c("grey","white")), ll=.04)

pdf("plots_5/general.beanplots.onlyDEL.pdf",w=20,h=10)
ggplot(data=df.ddrscores.m)+
  geom_density(data=subset(df.ddrscores.m,Freq == 1), aes(x=Score,y=-..density..,fill=factor(Freq)), trim=F)+
  geom_density(data=subset(df.ddrscores.m,Freq == 0), aes(x=Score,y= ..density..,fill=factor(Freq)), trim=F)+
  coord_flip() + xlab("body height (inch)") + ylab("Smoothed Density") + 
  facet_wrap(~Pathway,nrow=3,scale="free") 
dev.off()

colnames(tb) <- sapply(colnames(tb), function(x) strsplit(x,"_vs_")[[1]][1])
int <- intersect(colnames(tb),colnames(ddrscores))
p <- pheatmap(ddrscores[which(rownames(ddrscores)!=""),int])
colorder <- colnames(ddrscores[which(rownames(ddrscores)!=""),int])[p$tree_col$order]

library(viridis)

pdf("plots_5/mutsandexpr.together.pdf")
pheatmap(tb[,colorder],cluster_cols = FALSE)
pheatmap(ddrscores[which(rownames(ddrscores)!=""),int],color=inferno(100))
dev.off()

write.csv(ddrscores[which(rownames(ddrscores)!=""),int], file="ddrscores.csv")

pdf("plots_5/mutsandexpr.together2.pdf")
pheatmap(ddrscores[which(rownames(ddrscores)!=""),int],color=inferno(100))
dev.off()

df.exp <- t(ddrscores)
colnames(df.exp) <- paste0(colnames(df.exp),".exp")
df.exp <- data.frame(df.exp)
df.mut <- t(tb)
colnames(df.mut) <- paste0(colnames(df.mut),".mut")
df.mut <- as.data.frame(df.mut, drop=FALSE)

df.exp$Sample <- rownames(df.exp)
df.mut$Sample <- rownames(df.mut)

