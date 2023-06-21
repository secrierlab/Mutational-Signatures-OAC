#################
### This scripts merges and annotates clonal and subclonal signatures from different cohorts.

# Load mutational signature data:
load("processeddata/mets.snvs.clonal.early.RData")
load("processeddata/mets.snvs.clonal.late.RData")
load("processeddata/mets.snvs.subclonal.RData")
load("processeddata/barretts.snvs.clonal.early.RData")
load("processeddata/barretts.snvs.clonal.late.RData")
load("processeddata/barretts.snvs.subclonal.RData")
load("processeddata/extra.snvs.clonal.early.RData")
load("processeddata/extra.snvs.clonal.late.RData")
load("processeddata/extra.snvs.subclonal.RData")
load("processeddata/primaries.snvs.clonal.early.RData")
load("processeddata/primaries.snvs.clonal.late.RData")
load("processeddata/primaries.snvs.subclonal.RData")

snvs.clonal.early <- rbind(mets.snvs.clonal.early,
                           primaries.snvs.clonal.early,
                           barretts.snvs.clonal.early,
                           extra.snvs.clonal.early)
snvs.clonal.late <- rbind(mets.snvs.clonal.late,
                           primaries.snvs.clonal.late,
                           barretts.snvs.clonal.late,
                           extra.snvs.clonal.late)
snvs.subclonal <- rbind(mets.snvs.subclonal,
                          primaries.snvs.subclonal,
                          barretts.snvs.subclonal,
                          extra.snvs.subclonal)

snvs.clonal.early$Timing <- "Early"
snvs.clonal.late$Timing <- "Late"
snvs.subclonal$Timing <- NA
snvs.clonal.early$Clonality <- "Clonal"
snvs.clonal.late$Clonality <- "Clonal"
snvs.subclonal$Clonality <- "Subclonal"

# Load annotation of samples:
load("../tracksig_barrettsmetsprimaries/annotation.sampleIDs.RData")

snvs.clonal.early.timed.annot <- merge(snvs.clonal.early,
                          annotation.sampleIDs[,c("Sample","Category")],
                          by.x="Sample",by.y="Sample",
                          all.x=FALSE,all.y=FALSE)
snvs.clonal.late.timed.annot <- merge(snvs.clonal.late,
                                       annotation.sampleIDs[,c("Sample","Category")],
                                       by.x="Sample",by.y="Sample",
                                       all.x=FALSE,all.y=FALSE)
snvs.subclonal.timed.annot <- merge(snvs.subclonal,
                                       annotation.sampleIDs[,c("Sample","Category")],
                                       by.x="Sample",by.y="Sample",
                                       all.x=FALSE,all.y=FALSE)

snvs.all.timed.annot <- rbind(snvs.clonal.early.timed.annot,
                              snvs.clonal.late.timed.annot,
                              snvs.subclonal.timed.annot)
save(snvs.all.timed.annot, file="processeddata/snvs.all.timed.annot.RData")

