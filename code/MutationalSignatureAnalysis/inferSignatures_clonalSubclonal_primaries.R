###################
### This script infers mutational signatures for clonal and subclonal mutations.

### Load MutationTimer results:
load("../MutationTimer/primaries.mutMT.cls.456.RData")
#load("../MutationTimer/barretts.mutMT.cls.174.RData")
#load("../MutationTimer/primaries.mutMT.cls.122.RData")
#load("../MutationTimer/extra.mutMT.cls.245.RData")

# Split into early/late:
primaries.clonal.early <- primaries.mutMT.cls[which(primaries.mutMT.cls$CLS == "clonal [early]"),]
primaries.clonal.late <- primaries.mutMT.cls[which(primaries.mutMT.cls$CLS == "clonal [late]"),]
primaries.clonal <- primaries.mutMT.cls[which(primaries.mutMT.cls$CLS %in% c("clonal [early]",
                                                              "clonal [late]",
                                                              "clonal [NA]")),]
primaries.subclonal <- primaries.mutMT.cls[which(primaries.mutMT.cls$CLS == "subclonal"),]

# Load signature info:
load("processeddata/sbssigs.keep.RData")

### Mutational signature inference with deconstructSigs:

library(deconstructSigs)
library("BSgenome.Hsapiens.UCSC.hg19")

######################
# Early signatures:
snvs.early <- primaries.clonal.early[,c("Sample","Mutation")]
snvs.early$Chr <- sapply(primaries.clonal.early$Mutation,
                     function(x) strsplit(x,":")[[1]][1])
snvs.early$Start <- sapply(primaries.clonal.early$Mutation,
                         function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.early$REF <- sapply(primaries.clonal.early$Mutation,
                           function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.early$ALT <- sapply(primaries.clonal.early$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.early$Start <- as.numeric(snvs.early$Start)
snvs.early$Chr <- paste0("chr",snvs.early$Chr)
primaries.snvs.clonal.early <- snvs.early
save(primaries.snvs.clonal.early, file="processeddata/primaries.snvs.clonal.early.RData")

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
primaries.sigs.clonal.early <- data.frame(sigs.early)
save(primaries.sigs.clonal.early, file="processeddata/primaries.sigs.clonal.early.RData")


###################
### Now for the late signatures:
snvs.late <- primaries.clonal.late[,c("Sample","Mutation")]
snvs.late$Chr <- sapply(primaries.clonal.late$Mutation,
                         function(x) strsplit(x,":")[[1]][1])
snvs.late$Start <- sapply(primaries.clonal.late$Mutation,
                           function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.late$REF <- sapply(primaries.clonal.late$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.late$ALT <- sapply(primaries.clonal.late$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.late$Start <- as.numeric(snvs.late$Start)
snvs.late$Chr <- paste0("chr",snvs.late$Chr)
primaries.snvs.clonal.late <- snvs.late
save(primaries.snvs.clonal.late, file="processeddata/primaries.snvs.clonal.late.RData")


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
primaries.sigs.clonal.late <- data.frame(sigs.late)
save(primaries.sigs.clonal.late, file="processeddata/primaries.sigs.clonal.late.RData")

###################
### Now for all clonal signatures:
primaries.clonal.full <-primaries.clonal
samples <- unique(primaries.clonal.full$Sample)
primaries.clonal <- primaries.clonal.full[which(primaries.clonal.full$Sample %in% samples[401:435]),] 

snvs.all <- primaries.clonal[,c("Sample","Mutation")]
snvs.all$Chr <- sapply(primaries.clonal$Mutation,
                        function(x) strsplit(x,":")[[1]][1])
snvs.all$Start <- sapply(primaries.clonal$Mutation,
                          function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.all$REF <- sapply(primaries.clonal$Mutation,
                        function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.all$ALT <- sapply(primaries.clonal$Mutation,
                        function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.all$Start <- as.numeric(snvs.all$Start)
snvs.all$Chr <- paste0("chr",snvs.all$Chr)
primaries.snvs.clonal.all.5 <- snvs.all
save(primaries.snvs.clonal.all.5, file="processeddata/primaries.snvs.clonal.all.5.RData")


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
primaries.sigs.clonal.all.5 <- data.frame(sigs.all)
save(primaries.sigs.clonal.all.5, file="processeddata/primaries.sigs.clonal.all.5.RData")

###################
### Now for the subclonal signatures:
snvs.subclonal <- primaries.subclonal[,c("Sample","Mutation")]
snvs.subclonal$Chr <- sapply(primaries.subclonal$Mutation,
                       function(x) strsplit(x,":")[[1]][1])
snvs.subclonal$Start <- sapply(primaries.subclonal$Mutation,
                         function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.subclonal$REF <- sapply(primaries.subclonal$Mutation,
                       function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.subclonal$ALT <- sapply(primaries.subclonal$Mutation,
                       function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.subclonal$Start <- as.numeric(snvs.subclonal$Start)
snvs.subclonal$Chr <- paste0("chr",snvs.subclonal$Chr)
primaries.snvs.subclonal <- snvs.subclonal
save(primaries.snvs.subclonal, file="processeddata/primaries.snvs.subclonal.RData")


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
primaries.sigs.subclonal <- data.frame(sigs.subclonal)
save(primaries.sigs.subclonal, file="processeddata/primaries.sigs.subclonal.RData")
