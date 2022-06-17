########################
##### This script processes SigProfiler output and runs 
##### deconstructSigs to infer signature contributions 
##### based on optimal solution

library(gdata)
library(ggplot2)
library(lsa)

## Read new sig probs:
sigprob <- read.delim("COSMIC_v3.2_SBS_GRCh38.txt")
sp <- sigprob[,-1]
sp <- t(sp)
colnames(sp) <- sigprob$Type
load("data/signatures.genome.cosmic.v3.may2019.rda")

## Reorder sp according to previous order:
sp <- sp[,colnames(signatures.genome.cosmic.v3.may2019)]

### Read in the signature stats:
stats1 <- read.csv("SigProfiler_solutions/sigextractor_output/SBS96/All_solutions_stat.csv",
                   stringsAsFactors = FALSE)
stats2 <- read.csv("SigProfiler_solutions/sigextractor_output_2/SBS96/All_solutions_stat.csv",
                   stringsAsFactors = FALSE)
stats3 <- read.csv("SigProfiler_solutions/sigextractor_output_3/SBS96/All_solutions_stat.csv",
                   stringsAsFactors = FALSE)

stats <- rbind(stats1,stats2,stats3)

### Plot avg stability against frobenius reconstruction error:

stats.overlay <- stats[,c("Total.Signatures","Matrix.Frobenius.")]
colnames(stats.overlay) <- "avgStability"
stats$Total.Signatures <- as.numeric(as.character(stats$Total.Signatures))
stats$avgStability <- as.numeric(as.character(stats$avgStability))
stats$Matrix.Frobenius. <- as.numeric(as.character(stats$Matrix.Frobenius.))

pdf("plots/statsStability.pdf",w=10,h=5)
qplot(factor(Total.Signatures), avgStability, geom = c("line"), group = 1,
      data= stats)
dev.off()

pdf("plots/statsError.pdf",w=10,h=5)
qplot(factor(Total.Signatures), Matrix.Frobenius., geom = c("line"), group = 1,
      col="red",
      data= stats)
dev.off()


### Read in the signature probabilities for the 13 sig solution:
sol13 <- read.table("SigProfiler_solutions/sigextractor_output/SBS96/All_Solutions/SBS96_13_Signatures/Signatures/SBS96_S13_Signatures.txt",
                    sep="\t", stringsAsFactors = FALSE,
                    header=TRUE)
#sbssigs <- read.csv("sigProfiler_SBS_signatures_2019_05_22.csv")
#sbssigs$Context <- apply(sbssigs[,c("Type","SubType")],
 #                        1,
#                         function(x) paste0(substring(x[2],1,1),
#                                            "[",x[1],"]",substring(x[2],3,3)))

# Matching contexts:
rownames(sigprob) <- sigprob$Type
sbssigs.ordered <- sigprob[sol13$MutationType,]

### Do this for a range of solutions:

for (chosensig in 2:16) {
  
  print(paste0("===Chosen signature: ",chosensig))

### Read in the signature probabilities for the  sig solution:
  if (chosensig <14) {
    currentfile <- "SigProfiler_solutions/sigextractor_output/SBS96/All_Solutions/SBS96_"
  } else {
    currentfile <- "SigProfiler_solutions/sigextractor_output_2/SBS96/All_Solutions/SBS96_"
  }
sol13 <- read.table(paste0(currentfile,chosensig,"_Signatures/Signatures/SBS96_S",chosensig,"_Signatures.txt"),
                    sep="\t", stringsAsFactors = FALSE,
                    header=TRUE)

# Next, calculate cosine similarity with every signature:
similarities.13 <- array(-1,c(ncol(sol13)-1,ncol(sbssigs.ordered[,-c(1:2,ncol(sbssigs.ordered))])))
rownames(similarities.13) <- colnames(sol13[,-1])
colnames(similarities.13) <- colnames(sbssigs.ordered[,-c(1:2,ncol(sbssigs.ordered))])
for (s in colnames(sol13[,-1])) {
  for (cosmics in colnames(sbssigs.ordered[,-c(1:2,ncol(sbssigs.ordered))])) {
    similarities.13[s,cosmics] <- cosine(sol13[,s],sbssigs.ordered[,cosmics])
  }
}

similarities.13 <- data.frame(similarities.13)
similarities.13$MaxSimilarity <- apply(similarities.13, 1,
                                    function(x) colnames(similarities.13)[which(x==max(x))])
similarities.13$MaxCosine <- apply(similarities.13, 1,
                                function(x) x[max(x)])

print(similarities.13$MaxSimilarity)
print(similarities.13$MaxCosine)
print("=========")

}

### Decided on sol13 because it is the first that includes the platinum signature, doesn't have repeated signatures and has good stability (including per signature)
# after s=13, signature 41 is doubled, which is not what we want to achieve

# After choosing chosensig=13 and running the above code again, here are the resulting signatures:
aliases <- similarities.13$MaxSimilarity
names(aliases) <- rownames(similarities.13)

# add sbs17a instead of the doubled sbs17b:
aliases.solution13 <- setdiff(c(unique(aliases),"SBS17a","SBS8"),"SBS3")

### Finally, run deconstructSigs to infer signature prevalence in every sample:
library(deconstructSigs)
load("primaries.vcf.RData")
load("barretts.vcf.RData")
load("mets.vcf.RData")
load("extra.vcf.RData")

all.vcf <- rbind(primaries.vcf, barretts.vcf, mets.vcf, extra.vcf)
save(all.vcf, file="all.vcf.20201214.RData")

### Calculate TMB for every samples:
all.TMB <- table(all.vcf$Sample)
save(all.TMB, file="all.TMB.20201214.RData")

# No need to remove samples for the moment, there was just 1 below 50, with 44 mutations which is roughly ok

sigs.input <- mut.to.sigs.input(mut.ref = all.vcf, 
                                sample.id = "Sample", 
                                chr = "#CHROM", 
                                pos = "POS", 
                                ref = "REF", 
                                alt = "ALT")
save(sigs.input, file="deconstructSigs.input.20201214.RData")

sbssigs.keep <- data.frame(t(sbssigs.ordered[,aliases.solution13]))
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
save(sigs.allSamples, file="sigs.allSamples.normalised.202101127_cosmicref13sig_S8only.RData")

sigs.absolute <- sigs.allSamples
for (s in rownames(sigs.allSamples)) {
  sigs.absolute[s,] = round(all.TMB[s]*sigs.allSamples[s,])
}
sigs.allSamples.absolute <- sigs.absolute
save(sigs.allSamples.absolute, file="sigs.allSamples.absolute.202101127_cosmicref13sig_S8only.RData")


#### All sigs deconstructsigs:
save(sigs, file="sigs.deconstructsigs.RData")

# Remove artefactual sigs:
artefacts <- paste0("SBS",c(27,43,45:60))
sigs <- sigs[,setdiff(colnames(sigs),artefacts)]
pdf("plots/deconstructSigs.allsigs.pdf",w=20)
barplot(rev(sort(apply(sigs, 2, function(x) median(x)))),las=2)
dev.off()

pdf("plots/deconstructSigs.allsigs.mean.pdf",w=20)
barplot(rev(sort(apply(sigs, 2, function(x) mean(x)))),las=2)
dev.off()



