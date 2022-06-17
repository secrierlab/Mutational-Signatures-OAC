#################
### This scripts merges and annotates clonal and subclonal signatures from different cohorts.

# Load mutational signature data:
load("processeddata/mets.sigs.clonal.early.RData")
load("processeddata/mets.sigs.clonal.late.RData")
load("processeddata/mets.sigs.subclonal.RData")
load("processeddata/primaries.sigs.clonal.early.RData")
load("processeddata/primaries.sigs.clonal.late.RData")
load("processeddata/primaries.sigs.subclonal.RData")
load("processeddata/barretts.sigs.clonal.early.RData")
load("processeddata/barretts.sigs.clonal.late.RData")
load("processeddata/barretts.sigs.subclonal.RData")
load("processeddata/extra.sigs.clonal.early.RData")
load("processeddata/extra.sigs.clonal.late.RData")
load("processeddata/extra.sigs.subclonal.RData")

sigs.clonal.early <- rbind(mets.sigs.clonal.early,
                           primaries.sigs.clonal.early,
                           barretts.sigs.clonal.early,
                           extra.sigs.clonal.early)
sigs.clonal.late <- rbind(mets.sigs.clonal.late,
                           primaries.sigs.clonal.late,
                           barretts.sigs.clonal.late,
                           extra.sigs.clonal.late)
sigs.subclonal <- rbind(mets.sigs.subclonal,
                          primaries.sigs.subclonal,
                          barretts.sigs.subclonal,
                          extra.sigs.subclonal)

sigs.clonal.early$Timing <- "Early"
sigs.clonal.late$Timing <- "Late"
sigs.subclonal$Timing <- NA
sigs.clonal.early$Clonality <- "Clonal"
sigs.clonal.late$Clonality <- "Clonal"
sigs.subclonal$Clonality <- "Subclonal"

sigs.all <- rbind(sigs.clonal.early,
                  sigs.clonal.late,
                  sigs.subclonal)
sigs.all$Sample <- rownames(sigs.all)

# Load annotation of samples:
load("../tracksig_barrettsmetsprimaries/annotation.sampleIDs.RData")

# Need to modify because annotation is weird:
sigs.all$TumourID <- sapply(sigs.all$Sample,
                            function(x) strsplit(x,"_vs_")[[1]][1])
sigs.all <- sigs.all[,which(!(colnames(sigs.all) == "Sample"))]
sigs.timed.annot <- merge(sigs.all,
                          annotation.sampleIDs[,c("TumourID","Sample","Category")],
                          by.x="TumourID",by.y="TumourID",
                          all.x=FALSE,all.y=FALSE)
save(sigs.timed.annot, file="processeddata/sigs.timed.annot.RData")

