###################
### This script infers mutational signatures for clonal and subclonal mutations.

### Load MutationTimer results:
#load("barretts.mutMT.cls.456.RData")
load("../MutationTimer/barretts.mutMT.cls.174.RData")
#load("../MutationTimer/barretts.mutMT.cls.122.RData")
#load("../MutationTimer/extra.mutMT.cls.245.RData")

# Split into early/late:
barretts.clonal.early <- barretts.mutMT.cls[which(barretts.mutMT.cls$CLS == "clonal [early]"),]
barretts.clonal.late <- barretts.mutMT.cls[which(barretts.mutMT.cls$CLS == "clonal [late]"),]
barretts.clonal <- barretts.mutMT.cls[which(barretts.mutMT.cls$CLS %in% c("clonal [early]",
                                                              "clonal [late]",
                                                              "clonal [NA]")),]
barretts.subclonal <- barretts.mutMT.cls[which(barretts.mutMT.cls$CLS == "subclonal"),]

# Load signature info:
load("processeddata/sbssigs.keep.RData")

### Mutational signature inference with deconstructSigs:

library(deconstructSigs)
library("BSgenome.Hsapiens.UCSC.hg19")

######################
# Early signatures:
snvs.early <- barretts.clonal.early[,c("Sample","Mutation")]
snvs.early$Chr <- sapply(barretts.clonal.early$Mutation,
                     function(x) strsplit(x,":")[[1]][1])
snvs.early$Start <- sapply(barretts.clonal.early$Mutation,
                         function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.early$REF <- sapply(barretts.clonal.early$Mutation,
                           function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.early$ALT <- sapply(barretts.clonal.early$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.early$Start <- as.numeric(snvs.early$Start)
snvs.early$Chr <- paste0("chr",snvs.early$Chr)
barretts.snvs.clonal.early <- snvs.early
save(barretts.snvs.clonal.early, file="processeddata/barretts.snvs.clonal.early.RData")

sigs.input.early <- mut.to.sigs.input(mut.ref = barretts.snvs.clonal.early, 
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
barretts.sigs.clonal.early <- data.frame(sigs.early)
save(barretts.sigs.clonal.early, file="processeddata/barretts.sigs.clonal.early.RData")


###################
### Now for the late signatures:
snvs.late <- barretts.clonal.late[,c("Sample","Mutation")]
snvs.late$Chr <- sapply(barretts.clonal.late$Mutation,
                         function(x) strsplit(x,":")[[1]][1])
snvs.late$Start <- sapply(barretts.clonal.late$Mutation,
                           function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.late$REF <- sapply(barretts.clonal.late$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.late$ALT <- sapply(barretts.clonal.late$Mutation,
                         function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.late$Start <- as.numeric(snvs.late$Start)
snvs.late$Chr <- paste0("chr",snvs.late$Chr)
barretts.snvs.clonal.late <- snvs.late
save(barretts.snvs.clonal.late, file="processeddata/barretts.snvs.clonal.late.RData")


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
barretts.sigs.clonal.late <- data.frame(sigs.late)
save(barretts.sigs.clonal.late, file="processeddata/barretts.sigs.clonal.late.RData")

###################
### Now for all clonal signatures:
snvs.all <- barretts.clonal[,c("Sample","Mutation")]
snvs.all$Chr <- sapply(barretts.clonal$Mutation,
                        function(x) strsplit(x,":")[[1]][1])
snvs.all$Start <- sapply(barretts.clonal$Mutation,
                          function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.all$REF <- sapply(barretts.clonal$Mutation,
                        function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.all$ALT <- sapply(barretts.clonal$Mutation,
                        function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.all$Start <- as.numeric(snvs.all$Start)
snvs.all$Chr <- paste0("chr",snvs.all$Chr)
barretts.snvs.clonal.all <- snvs.all
save(barretts.snvs.clonal.all, file="processeddata/barretts.snvs.clonal.all.RData")


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
barretts.sigs.clonal.all <- data.frame(sigs.all)
save(barretts.sigs.clonal.all, file="processeddata/barretts.sigs.clonal.all.RData")

###################
### Now for the subclonal signatures:
snvs.subclonal <- barretts.subclonal[,c("Sample","Mutation")]
snvs.subclonal$Chr <- sapply(barretts.subclonal$Mutation,
                       function(x) strsplit(x,":")[[1]][1])
snvs.subclonal$Start <- sapply(barretts.subclonal$Mutation,
                         function(x) strsplit(strsplit(x,":")[[1]][2],"_")[[1]][1])
snvs.subclonal$REF <- sapply(barretts.subclonal$Mutation,
                       function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][1])
snvs.subclonal$ALT <- sapply(barretts.subclonal$Mutation,
                       function(x) strsplit(strsplit(strsplit(x,":")[[1]][2],"_")[[1]][2],"/")[[1]][2])
snvs.subclonal$Start <- as.numeric(snvs.subclonal$Start)
snvs.subclonal$Chr <- paste0("chr",snvs.subclonal$Chr)
barretts.snvs.subclonal <- snvs.subclonal
save(barretts.snvs.subclonal, file="processeddata/barretts.snvs.subclonal.RData")


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
barretts.sigs.subclonal <- data.frame(sigs.subclonal)
save(barretts.sigs.subclonal, file="processeddata/barretts.sigs.subclonal.RData")
