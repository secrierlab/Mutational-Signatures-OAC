########
### Mutect vs Strelka

library(ggpubr)

# Read SNVs:
files <- list.files("../calls_mutect/")
chrs <- c(paste0("chr",c(1:22)),"chrX","chrY")
vcfs.mutect <- NULL
for (f in files) {
  currentfile <- read.delim(paste0("../calls_mutect/",f),
                            comment.char ="#",header=FALSE)
  currentfile.keep <- currentfile[which((currentfile$V7 == "PASS")&
                                    (currentfile$V4 %in% c("A","C","G","T"))&
                                      (currentfile$V5 %in% c("A","C","G","T"))&
                                      (currentfile$V1 %in% chrs)),c(1,2,4,5)]
  currentfile.keep$Sample <- strsplit(f,"_filtered")[[1]][1] 
  vcfs.mutect <- rbind(vcfs.mutect,currentfile.keep)
}

# Next, compare to SNVs from Strelka from the same samples:
load("../../SigProfiler/all.vcf.20201214.RData")
vcfs.strelka <- all.vcf[which(all.vcf$Sample %in% unique(vcfs.mutect$Sample)),]

vcfs.intersect <- merge(vcfs.mutect, vcfs.strelka,
                        by.x=c("Sample","V1","V2","V4","V5"),
                        by.y=c("Sample","#CHROM","POS","REF","ALT"),
                        all.x=FALSE, all.y=FALSE)

# > dim(vcfs.intersect)
# [1] 1119919       9
# > nrow(vcfs.mutect)
# [1] 1251248
# > nrow(vcfs.strelka)
# [1] 1239201
# > 1251248-1119919
# [1] 131329
# > 131329/1251248
# [1] 0.1049584
# > 1239201-1119919
# [1] 119282
# > 119282/1239201
# [1] 0.09625719

## Infer mut sigs in mutect calls:
library(deconstructSigs)

sigs.input <- mut.to.sigs.input(mut.ref = vcfs.mutect, 
                                sample.id = "Sample", 
                                chr = "V1", 
                                pos = "V2", 
                                ref = "V4", 
                                alt = "V5")


## Read new sig probs:
sigprob <- read.delim("../SigProfiler/COSMIC_v3.2_SBS_GRCh38.txt")
sp <- sigprob[,-1]
sp <- t(sp)
colnames(sp) <- sigprob$Type
load("../deconstructSigs/data/signatures.genome.cosmic.v3.may2019.rda")

## Reorder sp according to previous order:
sp <- sp[,colnames(signatures.genome.cosmic.v3.may2019)]
sol13 <- read.table("../SigProfiler/run.december/sigextractor_output/SBS96/All_Solutions/SBS96_13_Signatures/Signatures/SBS96_S13_Signatures.txt",
                    sep="\t", stringsAsFactors = FALSE,
                    header=TRUE)
# Matching contexts:
rownames(sigprob) <- sigprob$Type
sbssigs.ordered <- sigprob[sol13$MutationType,]

load("processeddata/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")

sbssigs.keep <-  data.frame(t(sbssigs.ordered[,colnames(sigs.allSamples)]))
colnames(sbssigs.keep) <- rownames(sbssigs.ordered)
#save(signatures.genome.cosmic.v3.may2019.keep, file="signatures.genome.cosmic.v3.may2019.keep.RData")

sbssigs.all <- data.frame(t(sbssigs.ordered[,which(!(colnames(sbssigs.ordered) %in% c("Type","SubType","Context")))]))
colnames(sbssigs.all) <- rownames(sbssigs.ordered)

sigs <- NULL
i <- 0
mat.contexts.all <- NULL
for (sample in rownames(sigs.input)) {
  i <- i+1
  print(i)
  sigs_1 = whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = sbssigs.keep,
                           sample.id = sample, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0,
                           tri.counts.method = 'default')
  mat_contexts <- array(0,c(length(sigs_1$weights),length(sigs_1$tumor)))
  rownames(mat_contexts) <- names(sigs_1$weights)
  colnames(mat_contexts) <- names(sigs_1$tumor)
  for (sigsi in names(sigs_1$weights)) {
    for (context in names(sigs_1$tumor)) {
      mat_contexts[sigsi,context] <- sigs_1$weights[sigsi]*sigs_1$tumor[context]
    }
  }
  mat_contexts <- data.frame(mat_contexts)
  mat_contexts$Sample <- sample
  mat.contexts.all <- rbind(mat.contexts.all,
                            mat_contexts)
  sigs <- rbind(sigs,sigs_1$weights)
  
}
sigs.allSamples <- data.frame(sigs)
save(sigs.allSamples, file="processeddata/sigs.mutect.RData")

library(reshape)
sigs.allSamples$Sample <- rownames(sigs.allSamples)
df.mutect <- melt(sigs.allSamples)

load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.strelka <- sigs.allSamples[unique(vcfs.strelka$Sample),]
sigs.strelka$Sample <- rownames(sigs.strelka)
df.strelka <- melt(sigs.strelka)

df.both <- merge(df.mutect, df.strelka,
                 by.x=c("Sample","variable"), by.y=c("Sample","variable"),
                 all.x=FALSE, all.y=FALSE)

pdf("plots_supp/Mutect_vs_Strelka_signatureContributions.pdf",w=12,h=8)
ggscatter(df.both, x = "value.x", y = "value.y",
          fill = "#875053", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#2B061E", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 0, label.sep = "\n")
)+
  facet_wrap(~variable, scale="free",nrow=3)+
  xlab("Mutect signature contribution")+
  ylab("Strelka signature contribution")
dev.off()

write.csv(df.both, file="strelkavsmutect.csv")
