#########
## Load and prepare input:

library(stringr)

load("../tracksig_barrettsmetsprimaries/primaries.vcf.RData")
vcf.primaries <- primaries.vcf[,c("Sample","#CHROM","POS","REF","ALT")]
colnames(vcf.primaries) <- c("sampleID","chr","pos","ref","mut")
vcf.primaries$chr <- sapply(vcf.primaries$chr,
                            function(x) str_replace(x,"chr",""))
save(vcf.primaries, file="vcf.primaries.RData")

load("../tracksig_barrettsmetsprimaries/barretts.vcf.RData")
vcf.barretts <- barretts.vcf[,c("Sample","#CHROM","POS","REF","ALT")]
colnames(vcf.barretts) <- c("sampleID","chr","pos","ref","mut")
vcf.barretts$chr <- sapply(vcf.barretts$chr,
                            function(x) str_replace(x,"chr",""))
save(vcf.barretts, file="vcf.barretts.RData")

load("../tracksig_barrettsmetsprimaries/mets.vcf.RData")
vcf.mets <- mets.vcf[,c("Sample","#CHROM","POS","REF","ALT")]
colnames(vcf.mets) <- c("sampleID","chr","pos","ref","mut")
vcf.mets$chr <- sapply(vcf.mets$chr,
                           function(x) str_replace(x,"chr",""))
save(vcf.mets, file="vcf.mets.RData")


library("dndscv")
dndsout.primaries = dndscv(vcf.primaries)
save(dndsout.primaries, file="dndsout.primaries.RData")


### Load signature info:
load("../SigProfiler/sigs.allSamples.norm.RData")
sigs.allSamples.norm$S17 <- sigs.allSamples.norm$SBS17a+sigs.allSamples.norm$SBS17b
sigs.allSamples.norm$MaxSig <- apply(sigs.allSamples.norm[,-1],1,
                                      function(x) paste(colnames(sigs.allSamples.norm[,-1])[which(x==max(x))],collapse=";"))
samples.s17 <- sigs.allSamples.norm[which(sigs.allSamples.norm$MaxSig %in% c("S17","SBS17b;S17")),]$Sample
samples.other <- setdiff(sigs.allSamples.norm$Samples,samples.s17)

samples.s1 <- sigs.allSamples.norm[which(sigs.allSamples.norm$MaxSig == "SBS1"),]$Sample
samples.s18 <- sigs.allSamples.norm[which(sigs.allSamples.norm$MaxSig == "SBS18"),]$Sample
samples.s3 <- sigs.allSamples.norm[which(sigs.allSamples.norm$MaxSig == "SBS3"),]$Sample
samples.s5 <- sigs.allSamples.norm[which(sigs.allSamples.norm$MaxSig == "SBS5"),]$Sample
samples.s44 <- sigs.allSamples.norm[which(sigs.allSamples.norm$MaxSig == "SB44"),]$Sample


load("vcf.primaries.RData")
load("vcf.barretts.RData")
load("vcf.mets.RData")

vcf.primaries.S17 <- vcf.primaries[which(vcf.primaries$sampleID %in%
                                           samples.s17),]
vcf.primaries.other <- vcf.primaries[which(vcf.primaries$sampleID %in%
                                           samples.other),]
vcf.barretts.S17 <- vcf.barretts[which(vcf.barretts$sampleID %in%
                                           samples.s17),]
vcf.barretts.other <- vcf.barretts[which(vcf.barretts$sampleID %in%
                                             samples.other),]
vcf.mets.S17 <- vcf.mets[which(vcf.mets$sampleID %in%
                                         samples.s17),]
vcf.mets.other <- vcf.mets[which(vcf.mets$sampleID %in%
                                           samples.other),]

vcf.s17 <- rbind(vcf.primaries.S17, vcf.barretts.S17,
                 vcf.mets.S17)
vcf.other <- rbind(vcf.primaries.other, vcf.barretts.other,
                 vcf.mets.other)
save(vcf.s17, file="vcf.s17.RData")
save(vcf.other, file="vcf.othersig.RData")

load("vcf.othersig.RData")
vcf.s1 <- vcf.other[which(vcf.other$sampleID %in% samples.s1),]
vcf.s18 <- vcf.other[which(vcf.other$sampleID %in% samples.s18),]
vcf.s3 <- vcf.other[which(vcf.other$sampleID %in% samples.s3),]
vcf.s5 <- vcf.other[which(vcf.other$sampleID %in% samples.s5),]
vcf.s44 <- vcf.other[which(vcf.other$sampleID %in% samples.s44),]
save(vcf.s1, file="vcf.s1.RData")
save(vcf.s18, file="vcf.s18.RData")
save(vcf.s3, file="vcf.s3.RData")
save(vcf.s5, file="vcf.s5.RData")
save(vcf.s44, file="vcf.s44.RData")
