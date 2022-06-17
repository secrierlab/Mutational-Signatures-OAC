###################
### This script infers mutational signatures for clonal and subclonal mutations.

### Load MutationTimer results:
#load("../MutationTimer/primaries.mutMT.cls.456.RData")
#load("../MutationTimer/barretts.mutMT.cls.174.RData")
load("../MutationTimer/mets.mutMT.cls.122.RData")
#load("../MutationTimer/extra.mutMT.cls.245.RData")

# Split into early/late:
mets.clonal.early <- mets.mutMT.cls[which(mets.mutMT.cls$CLS == "clonal [early]"),]
mets.clonal.late <- mets.mutMT.cls[which(mets.mutMT.cls$CLS == "clonal [late]"),]
mets.clonal <- mets.mutMT.cls[which(mets.mutMT.cls$CLS %in% c("clonal [early]",
                                                              "clonal [late]",
                                                              "clonal [NA]")),]
mets.subclonal <- mets.mutMT.cls[which(mets.mutMT.cls$CLS == "subclonal"),]

# Load signature info:
load("processeddata/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
setsigs <- colnames(sigs.allSamples)

sbssigs <- read.csv("../SigProfiler/sigProfiler_SBS_signatures_2019_05_22.csv")
sbssigs$Context <- apply(sbssigs[,c("Type","SubType")],
                         1,
                         function(x) paste0(substring(x[2],1,1),
                                            "[",x[1],"]",substring(x[2],3,3)))

# Matching contexts:
rownames(sbssigs) <- sbssigs$Context
sbssigs.keep <- data.frame(t(sbssigs[,setsigs]))
colnames(sbssigs.keep) <- rownames(sbssigs)


### Mutational signature inference with deconstructSigs:

library(deconstructSigs)
library("BSgenome.Hsapiens.UCSC.hg19")

######################
# Early signatures:
snvs.early <- mets.clonal.early[,c("Sample","Mutation")]
snvs.early$Chr <- sapply(mets.clonal.early$Mutation,
                     function(x) strsplit(x,":")[[1]][1])
snvs.early$Start <- sapply(mets.clonal.early$Mutation,
                         function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.early$REF <- sapply(mets.clonal.early$Mutation,
                           function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.early$ALT <- sapply(mets.clonal.early$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.early$Start <- as.numeric(snvs.early$Start)
snvs.early$Chr <- paste0("chr",snvs.early$Chr)
mets.snvs.clonal.early <- snvs.early
save(mets.snvs.clonal.early, file="processeddata/mets.snvs.clonal.early.RData")

sigs.input.early <- mut.to.sigs.input(mut.ref = snvs.early, 
                                sample.id = "Sample", 
                                chr = "Chr", 
                                pos = "Start", 
                                ref = "REF", 
                                alt = "ALT",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)

sigs.early <- NULL
i <- 0
for (sample in rownames(sigs.input.early)) {
  i <- i+1
  print(i)
  sigs_1 = whichSignatures(tumor.ref = sigs.input.early, 
                           signatures.ref = sbssigs.keep,
                           sample.id = sample, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0,
                           tri.counts.method = 'default')
  sigs.early <- rbind(sigs.early,sigs_1$weights)
  
}
mets.sigs.clonal.early <- data.frame(sigs.early)
save(mets.sigs.clonal.early, file="processeddata/mets.sigs.clonal.early.RData")


###################
### Now for the late signatures:
snvs.late <- mets.clonal.late[,c("Sample","Mutation")]
snvs.late$Chr <- sapply(mets.clonal.late$Mutation,
                         function(x) strsplit(x,":")[[1]][1])
snvs.late$Start <- sapply(mets.clonal.late$Mutation,
                           function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.late$REF <- sapply(mets.clonal.late$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.late$ALT <- sapply(mets.clonal.late$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.late$Start <- as.numeric(snvs.late$Start)
snvs.late$Chr <- paste0("chr",snvs.late$Chr)
mets.snvs.clonal.late <- snvs.late
save(mets.snvs.clonal.late, file="processeddata/mets.snvs.clonal.late.RData")


sigs.input.late <- mut.to.sigs.input(mut.ref = snvs.late, 
                                      sample.id = "Sample", 
                                      chr = "Chr", 
                                      pos = "Start", 
                                      ref = "REF", 
                                      alt = "ALT",
                                      bsg = BSgenome.Hsapiens.UCSC.hg19)

sigs.late <- NULL
i <- 0
for (sample in rownames(sigs.input.late)) {
  i <- i+1
  print(i)
  sigs_1 = whichSignatures(tumor.ref = sigs.input.late, 
                           signatures.ref = sbssigs.keep,
                           sample.id = sample, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0,
                           tri.counts.method = 'default')
  sigs.late <- rbind(sigs.late,sigs_1$weights)
  
}
mets.sigs.clonal.late <- data.frame(sigs.late)
save(mets.sigs.clonal.late, file="processeddata/mets.sigs.clonal.late.RData")

###################
### Now for all clonal signatures:
snvs.all <- mets.clonal[,c("Sample","Mutation")]
snvs.all$Chr <- sapply(mets.clonal$Mutation,
                        function(x) strsplit(x,":")[[1]][1])
snvs.all$Start <- sapply(mets.clonal$Mutation,
                          function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.all$REF <- sapply(mets.clonal$Mutation,
                        function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.all$ALT <- sapply(mets.clonal$Mutation,
                        function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.all$Start <- as.numeric(snvs.all$Start)
snvs.all$Chr <- paste0("chr",snvs.all$Chr)
mets.snvs.clonal.all <- snvs.all
save(mets.snvs.clonal.all, file="processeddata/mets.snvs.clonal.all.RData")


sigs.input.all <- mut.to.sigs.input(mut.ref = snvs.all, 
                                     sample.id = "Sample", 
                                     chr = "Chr", 
                                     pos = "Start", 
                                     ref = "REF", 
                                     alt = "ALT",
                                     bsg = BSgenome.Hsapiens.UCSC.hg19)

sigs.all <- NULL
i <- 0
for (sample in rownames(sigs.input.all)) {
  i <- i+1
  print(i)
  sigs_1 = whichSignatures(tumor.ref = sigs.input.all, 
                           signatures.ref = sbssigs.keep,
                           sample.id = sample, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0,
                           tri.counts.method = 'default')
  sigs.all <- rbind(sigs.all,sigs_1$weights)
  
}
mets.sigs.clonal.all <- data.frame(sigs.all)
save(mets.sigs.clonal.all, file="processeddata/mets.sigs.clonal.all.RData")

###################
### Now for the subclonal signatures:
snvs.subclonal <- mets.subclonal[,c("Sample","Mutation")]
snvs.subclonal$Chr <- sapply(mets.subclonal$Mutation,
                       function(x) strsplit(x,":")[[1]][1])
snvs.subclonal$Start <- sapply(mets.subclonal$Mutation,
                         function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.subclonal$REF <- sapply(mets.subclonal$Mutation,
                       function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.subclonal$ALT <- sapply(mets.subclonal$Mutation,
                       function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.subclonal$Start <- as.numeric(snvs.subclonal$Start)
snvs.subclonal$Chr <- paste0("chr",snvs.subclonal$Chr)
mets.snvs.subclonal <- snvs.subclonal
save(mets.snvs.subclonal, file="processeddata/mets.snvs.subclonal.RData")


sigs.input.subclonal <- mut.to.sigs.input(mut.ref = snvs.subclonal, 
                                    sample.id = "Sample", 
                                    chr = "Chr", 
                                    pos = "Start", 
                                    ref = "REF", 
                                    alt = "ALT",
                                    bsg = BSgenome.Hsapiens.UCSC.hg19)

sigs.subclonal <- NULL
i <- 0
for (sample in rownames(sigs.input.subclonal)) {
  i <- i+1
  print(i)
  sigs_1 = whichSignatures(tumor.ref = sigs.input.subclonal, 
                           signatures.ref = sbssigs.keep,
                           sample.id = sample, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0,
                           tri.counts.method = 'default')
  sigs.subclonal <- rbind(sigs.subclonal,sigs_1$weights)
  
}
mets.sigs.subclonal <- data.frame(sigs.subclonal)
save(mets.sigs.subclonal, file="processeddata/mets.sigs.subclonal.RData")
