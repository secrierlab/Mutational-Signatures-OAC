###############
### This script checks the relation between signature exposures and treatment response.

library(ggpubr)
library(reshape)
library(wesanderson)
library(stringr)

# Treatment data:
clin <- read.csv("data/OCCAMS-export_MS-cohort_SAZ_20220711.csv",
                 header=TRUE)
summary(clin)

clin$TR.all.Response <- apply(clin[,c("TR.PalChemo.Response",
                                           "TR.Adj.Response",
                                           "TR.NeoAdj.Response",
                                           "TR.RT.Response")],1,
                                   function(x) ifelse(!is.na(x[3]),x[3],
                                                      ifelse(!is.na(x[4]),x[4],
                                                             ifelse(!is.na(x[2]),x[2],x[1]))))

# Load signatures and annotation:
load("data/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples[sigs.allSamples<0.05] <- 0
sigs.allSamples$Sample <- rownames(sigs.allSamples)
load("data/annotation.sampleIDs.RData")

annotation.sampleIDs$OCCAMS_ID <- sapply(annotation.sampleIDs$OCCAMS_ID,
                                         function(x) ifelse(grepl("OC/",x),str_replace(x,"OC/","OCCAMS/"),x))

clin <- merge(clin, annotation.sampleIDs,
              by.x="ID", by.y="OCCAMS_ID",
              all.x=FALSE, all.y=FALSE)
df.sigs <- melt(sigs.allSamples)

clin.sigs <- merge(clin, df.sigs,
                   by.x="Sample", by.y="Sample",
                   all.x=FALSE, all.y=FALSE)
clin.sigs$value <- as.numeric(clin.sigs$value)

pal <- wes_palette("Royal2", 100, type = "continuous")

clin.sigs.all <- clin.sigs

clin.sigs <- clin.sigs[which(clin.sigs$Category=="PrimaryTumour"),]
clin.sigs <- clin.sigs[which(clin.sigs$value!=0),]

clin.sigs$RP.Tumour.Growth.c <- factor(clin.sigs$RP.Tumour.Growth.c,
                                       levels=c("shrink","no change","grow"))
mycomp1 <- list(c("shrink","no change"),
                c("no change","grow"),
                c("shrink","grow"))
pdf("plots_7/cont.RP.Tumour.Growth.c.pdf",w=18,h=6)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.Tumour.Growth.c)),], 
          x = "RP.Tumour.Growth.c", y = "value",
          color = "RP.Tumour.Growth.c",palette = wes_palette("Royal2", n = 3),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.Tumour.Growth.c.merged.pdf",h=5)
ggboxplot(clin.sigs.merged[which((!is.na(clin.sigs.merged$RP.Tumour.Growth.c))&
                                   (clin.sigs.merged$variable=="SBS17b")),], 
          x = "RP.Tumour.Growth.c", y = "value",
          color = "Type",
          palette = wes_palette("Royal2", n = 3),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  #facet_wrap(~RP.Tumour.Growth.c,scale="free",nrow=1)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

write.csv(unique(clin.sigs.merged[which((!is.na(clin.sigs.merged$RP.Tumour.Growth.c))&
                         (clin.sigs.merged$variable=="SBS17b")),]),
          file="plots_7/shrinkgrow.csv")

naive <- read.csv("plots_1/naive.primmet.csv")
clin.sigs.naive <- clin.sigs[which(clin.sigs$ID %in% naive$OCCAMS_ID),]
clin.sigs.treated <- clin.sigs[which(!(clin.sigs$ID %in% naive$OCCAMS_ID)),]

pdf("plots_7/cont.RP.Tumour.Growth.c.treated.pdf",w=18,h=6)
ggboxplot(clin.sigs.treated[which(!is.na(clin.sigs.treated$RP.Tumour.Growth.c)),], 
          x = "RP.Tumour.Growth.c", y = "value",
          color = "RP.Tumour.Growth.c",palette = wes_palette("Royal2", n = 3),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.Tumour.Growth.c.naive.pdf",w=18,h=6)
ggboxplot(clin.sigs.naive[which(!is.na(clin.sigs.naive$RP.Tumour.Growth.c)),], 
          x = "RP.Tumour.Growth.c", y = "value",
          color = "RP.Tumour.Growth.c",palette = wes_palette("Royal2", n = 3),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

write.csv(clin.sigs[which(!is.na(clin.sigs$RP.Tumour.Growth.c)),],
          file="plots_7/growshrink.csv")

clin.sigs$RP.TumourResponse <- factor(clin.sigs$RP.TumourResponse,
                                       levels=c("0%","<20%","=20%","<50%","=50%"))
mycomp1 <- list(c("0%","<20%"),
                c("<20%","=20%"),
                c("=20%","<50%"),
                c("<50%","=50%"))
pdf("plots_7/cont.RP.TumourResponse.c.pdf",w=18,h=6)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.TumourResponse)),], 
          x = "RP.TumourResponse", y = "value",
          color = "RP.TumourResponse",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()


clin.sigs$RP.MandardScoreForResponse <- factor(clin.sigs$RP.MandardScoreForResponse)
mycomp1 <- list(c("TRG1","TRG2"),
                c("TRG2","TRG3"),
                c("TRG3","TRG4"),
                c("TRG4","TRG5"))
pdf("plots_7/cont.RP.MandardScoreForResponse.pdf",w=20,h=7)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.MandardScoreForResponse)),], 
          x = "RP.MandardScoreForResponse", y = "value",
          color = "RP.MandardScoreForResponse",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.MandardScoreForResponse.naive.pdf",w=20,h=7)
ggboxplot(clin.sigs.naive[which(!is.na(clin.sigs.naive$RP.MandardScoreForResponse)),], 
          x = "RP.MandardScoreForResponse", y = "value",
          color = "RP.MandardScoreForResponse",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.MandardScoreForResponse.treated..pdf",w=20,h=7)
ggboxplot(clin.sigs.treated[which(!is.na(clin.sigs.treated$RP.MandardScoreForResponse)),], 
          x = "RP.MandardScoreForResponse", y = "value",
          color = "RP.MandardScoreForResponse",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$RP.MandardScoreForResponse.Binary <- sapply(clin.sigs$RP.MandardScoreForResponse,
                                                      function(x) ifelse(is.na(x), x, ifelse(x %in% c("TRG1","TRG2"),"TRG1-2","TRG3+")))
clin.sigs$RP.MandardScoreForResponse.Binary2 <- sapply(clin.sigs$RP.MandardScoreForResponse,
                                                      function(x) ifelse(is.na(x), x, ifelse(x %in% c("TRG1","TRG2","TRG3"),"TRG1-3","TRG4+")))

clin.sigs$RP.MandardScoreForResponse.Binary <- factor(clin.sigs$RP.MandardScoreForResponse.Binary,
                                                      levels=c("TRG1-2","TRG3+"))
clin.sigs$RP.MandardScoreForResponse.Binary2 <- factor(clin.sigs$RP.MandardScoreForResponse.Binary2,
                                                      levels=c("TRG1-3","TRG4+"))
mycomp4 <- list(c("TRG1-2","TRG3+"))
pdf("plots_7/cont.RP.MandardScoreForResponse.Binary.pdf",w=20,h=7)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.MandardScoreForResponse.Binary)),], 
          x = "RP.MandardScoreForResponse.Binary", y = "value",
          color = "RP.MandardScoreForResponse.Binary",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp4)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs.naive <- clin.sigs[which(clin.sigs$Sample %in% naive),]
clin.sigs.naive$Type <- "Naive"
clin.sigs.treated <- clin.sigs[which(!(clin.sigs$Sample %in% naive)),]
clin.sigs.treated$Type <- "Treated"

clin.sigs.merged <- rbind(clin.sigs.naive,clin.sigs.treated)

pdf("plots_7/cont.RP.MandardScoreForResponse.Binary.naive.pdf",w=20,h=7)
ggboxplot(clin.sigs.naive[which(!is.na(clin.sigs.naive$RP.MandardScoreForResponse.Binary)),], 
          x = "RP.MandardScoreForResponse.Binary", y = "value",
          color = "RP.MandardScoreForResponse.Binary",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp4)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()


write.csv(unique(clin.sigs.naive[which((!is.na(clin.sigs.naive$RP.MandardScoreForResponse.Binary))&
                                  (clin.sigs.naive$variable=="SBS17a")),
                          c("Sample","variable","value","RP.MandardScoreForResponse.Binary")]),
          file="plots_7/naive.Mandard.csv")

mycomp5 <- list(c("TRG1-3","TRG4+"))
pdf("plots_7/cont.RP.MandardScoreForResponse.Binary2.naive.pdf",w=20,h=7)
ggboxplot(clin.sigs.naive[which(!is.na(clin.sigs.naive$RP.MandardScoreForResponse.Binary2)),], 
          x = "RP.MandardScoreForResponse.Binary2", y = "value",
          color = "RP.MandardScoreForResponse.Binary2",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp5)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.MandardScoreForResponse.Binary2.treated.pdf",w=20,h=7)
ggboxplot(clin.sigs.treated[which(!is.na(clin.sigs.treated$RP.MandardScoreForResponse.Binary2)),], 
          x = "RP.MandardScoreForResponse.Binary2", y = "value",
          color = "RP.MandardScoreForResponse.Binary2",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp5)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.MandardScoreForResponse.Binary.treated.pdf",w=20,h=7)
ggboxplot(clin.sigs.treated[which(!is.na(clin.sigs.treated$RP.MandardScoreForResponse.Binary)),], 
          x = "RP.MandardScoreForResponse.Binary", y = "value",
          color = "RP.MandardScoreForResponse.Binary",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp4)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.RP.MandardScoreForResponse.Binary.marged.pdf")
ggboxplot(clin.sigs.merged[which((!is.na(clin.sigs.merged$RP.MandardScoreForResponse.Binary))&
                                   (clin.sigs.merged$variable=="SBS17b")),], 
          x = "RP.MandardScoreForResponse.Binary", y = "value",
          color = "RP.MandardScoreForResponse.Binary",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp4)+
  facet_wrap(Type~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

write.csv(clin.sigs[which(!is.na(clin.sigs$RP.MandardScoreForResponse)),],
          "plots_supp/mandard.csv")

clin.sigs$TR.PalChemo.Response <- factor(clin.sigs$TR.PalChemo.Response,
                                         levels=c("PR","SD","SD, PD","PD"))
mycomp1 <- list(c("PR","SD"),
                c("SD","SD, PD"),
                c("SD, PD","PD"))
pdf("plots_7/cont.TR.PalChemo.Response.pdf",w=18,h=7)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.PalChemo.Response)),], 
          x = "TR.PalChemo.Response", y = "value",
          color = "TR.PalChemo.Response",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.Adj.Response <- factor(clin.sigs$TR.Adj.Response,
                                    levels=c("CR","PR","PD"))
mycomp1 <- list(c("CR","PR"),
                c("PR","PD"),
                c("CR","PD"))
pdf("plots_7/cont.TR.Adj.Response.pdf",w=20,h=6)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.Adj.Response)),], 
          x = "TR.Adj.Response", y = "value",
          color = "TR.Adj.Response",palette = wes_palette("Royal2", n = 3),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.RT.Response <- factor(clin.sigs$TR.RT.Response,
                                    levels=c("CR","PR","SD","SD, PD","PD"))
mycomp1 <- list(c("CR","PR"),
                c("PR","SD"),
                c("SD","SD, PD"),
                c("SD, PD","PD"))
pdf("plots_7/cont.TR.RT.Response.pdf",w=21,h=7)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response)),], 
          x = "TR.RT.Response", y = "value",
          color = "TR.RT.Response",palette = wes_palette("Royal2", n = 5),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.RT.Response.RvsD <- sapply(clin.sigs$TR.RT.Response,
                                        function(x) ifelse (is.na(x),x,ifelse(x %in% c("CR","PR"),"Responder","Non-responder")))
mycomp1 <- list(c("Responder","Non-responder"))
pdf("plots_7/cont.TR.RT.RespondervsNonResponder.pdf",w=18,h=6)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response.RvsD)),], 
          x = "TR.RT.Response.RvsD", y = "value",
          color = "TR.RT.Response.RvsD",palette = wes_palette("Royal2", n = 2),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs.naive <- clin.sigs[which(clin.sigs$Sample %in% naive),]
clin.sigs.treated <- clin.sigs[which(!(clin.sigs$Sample %in% naive)),]

pdf("plots_7/cont.TR.RT.RespondervsNonResponder.naive.pdf",w=18,h=6)
ggboxplot(clin.sigs.naive[which(!is.na(clin.sigs.naive$TR.RT.Response.RvsD)),], 
          x = "TR.RT.Response.RvsD", y = "value",
          color = "TR.RT.Response.RvsD",palette = wes_palette("Royal2", n = 2),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.TR.RT.RespondervsNonResponder.treated.pdf",w=18,h=6)
ggboxplot(clin.sigs.treated[which(!is.na(clin.sigs.treated$TR.RT.Response.RvsD)),], 
          x = "TR.RT.Response.RvsD", y = "value",
          color = "TR.RT.Response.RvsD",palette = wes_palette("Royal2", n = 2),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.all.Response.Binary <- sapply(clin.sigs$TR.all.Response,
                                        function(x) ifelse (is.na(x),x,ifelse(x %in% c("CR"),"Responder","Non-responder")))
mycomp1 <- list(c("Responder","Non-responder"))
pdf("plots_7/cont.TR.all.Response.RespondervsNonResponder.pdf",w=18,h=6)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.all.Response.Binary)),], 
          x = "TR.all.Response.Binary", y = "value",
          color = "TR.all.Response.Binary",palette = wes_palette("Royal2", n = 2),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.RT.Response.CRvsrest <- sapply(clin.sigs$TR.RT.Response,
                                        function(x) ifelse (is.na(x), x, ifelse(x %in% c("CR"),"Complete responder","Non-responder")))
mycomp1 <- list(c("Complete responder","Non-responder"))
pdf("plots_7/cont.TR.RT.Response.CRvsrest.pdf",w=18,h=6)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.RT.Response.CRvsrest)),], 
          x = "TR.RT.Response.CRvsrest", y = "value",
          color = "TR.RT.Response.CRvsrest",palette = wes_palette("Royal2", n = 2),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.NeoAdj.Response <- factor(clin.sigs$TR.NeoAdj.Response,
                                   levels=c("CR","PR","SD","PD"))
mycomp1 <- list(c("CR","PR"),
                c("PR","SD"),
                c("SD","PD"))
pdf("plots_7/cont.TR.NeoAdj.Response.pdf",w=18,h=8)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response)),], 
          x = "TR.NeoAdj.Response", y = "value",
          color = "TR.NeoAdj.Response",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.NeoAdj.Response.StableVsProg <- sapply(clin.sigs$TR.NeoAdj.Response,
                                       function(x) ifelse(is.na(x),x,ifelse(x=="PD","PD","Response/Stable")))
mycomp1 <- list(c("PD","Response/Stable"))
pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg.pdf",w=18,h=8)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg)),], 
          x = "TR.NeoAdj.Response.StableVsProg", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.NeoAdj.Response.StableVsProg2 <- sapply(clin.sigs$TR.NeoAdj.Response,
                                                    function(x) ifelse(is.na(x),x,ifelse(x %in% c("PD","SD"),"Non-responder","Responder")))
mycomp2 <- list(c("Non-responder","Responder"))
pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg2.pdf",w=18,h=8)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg2)),], 
          x = "TR.NeoAdj.Response.StableVsProg2", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg2",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp2)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs.naive <- clin.sigs[which(clin.sigs$Sample %in% naive),]
clin.sigs.treated <- clin.sigs[which(!(clin.sigs$Sample %in% naive)),]


pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg2.naive.pdf",w=18,h=8)
ggboxplot(clin.sigs.naive[which(!is.na(clin.sigs.naive$TR.NeoAdj.Response.StableVsProg2)),], 
          x = "TR.NeoAdj.Response.StableVsProg2", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg2",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp2)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg2.treated.pdf",h=5)
ggboxplot(clin.sigs.treated[which(!is.na(clin.sigs.treated$TR.NeoAdj.Response.StableVsProg2)),], 
          x = "TR.NeoAdj.Response.StableVsProg2", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg2",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp2)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$TR.NeoAdj.Response.StableVsProg3 <- sapply(clin.sigs$TR.NeoAdj.Response,
                                                     function(x) ifelse(is.na(x),x,ifelse(x %in% c("CR"),"CR","Other")))
mycomp3 <- list(c("CR","Other"))
pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg3.pdf",w=18,h=8)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg3)),], 
          x = "TR.NeoAdj.Response.StableVsProg3", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg3",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp3)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

mycomp4 <- list(c("CR","PD"))
pdf("plots_7/cont.TR.NeoAdj.Response.CRvsPD.pdf",w=18,h=8)
ggboxplot(clin.sigs[which(clin.sigs$TR.NeoAdj.Response %in% c("CR","PD")),], 
          x = "TR.NeoAdj.Response", y = "value",
          color = "TR.NeoAdj.Response",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp4)+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

write.csv(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg)),],
          file="plots_7/stableprogressive.csv")

clin.sigs$RP.PS.Tumour.Diff.c <- factor(clin.sigs$RP.PS.Tumour.Diff.c,
                                       levels=sort(unique(clin.sigs$RP.PS.Tumour.Diff.c)))
pdf("plots_7/cont.RP.PS.Tumour.Diff.c.pdf",w=22,h=8)
ggboxplot(clin.sigs[which(!is.na(clin.sigs$RP.PS.Tumour.Diff.c)),], 
          x = "RP.PS.Tumour.Diff.c", y = "value",
          color = "RP.PS.Tumour.Diff.c",palette = wes_palette("Royal2", n = 14, type="continuous"),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

clin.sigs$RP.PS.Tumour.Diff.c <- as.numeric(as.character(clin.sigs$RP.PS.Tumour.Diff.c))
pdf("plots_7/cont.RP.PS.Tumour.Diff.c.scatter.pdf",w=22,h=8)
ggscatter(clin.sigs[which(!is.na(clin.sigs$RP.PS.Tumour.Diff.c)),], 
          x = "RP.PS.Tumour.Diff.c", y = "value",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))+
  facet_wrap(~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

fisher.test(matrix(c(1+1,24+9,4+22+23,25+49+9),nrow=2))
#4+22+23
#49
#25+49+9
# 83

# Join with primary tumour clinical data:
oac <- read.csv("../clin_data_20220609.csv")
table(oac$PS.TStage.PrimaryTumour.FinalPretreatmentStaging)
oac$Tsage.split <- sapply(oac$PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
                          function(x) ifelse(is.na(x), NA, ifelse(x %in% c("T3","T4","T4a","T4b"),"T3+","Tlow")))

clin.m <- merge(clin, oac,
                by.x="TumourID", by.y="TumourID",
                all.x=FALSE, all.y=FALSE)

table(clin.m[,c("Tsage.split","RP.Tumour.Growth.c")])

fisher.test(matrix(c(41,93,52,41),nrow=2))
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c(41, 93, 52, 41), nrow = 2)
# p-value = 0.0001962
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.1931788 0.6248621
# sample estimates:
# odds ratio 
#  0.3493195 

table(clin.m[,c("Tsage.split","RP.MandardScoreForResponse.x")])



clins.sigs.stage <- merge(clin.sigs,
                          clin.m[,c("Sample","Tsage.split")],
                          by.x='Sample', by.y="Sample",
                          all.x=FALSE, all.y=FALSE)

table(unique(clins.sigs.stage[,c("Sample","Tsage.split","TR.NeoAdj.Response.StableVsProg")])[,c("Tsage.split","TR.NeoAdj.Response.StableVsProg")])
fisher.test(matrix(c(17,216,6,58),nrow=2))

mycomp3<-list(c("T3+","Tlow"))
pdf("plots_7/sbs17bystage.pdf")
ggboxplot(clins.sigs.stage[which(clins.sigs.stage$variable %in% c("SBS17a","SBS17b")),], 
          x = "Tsage.split", y = "value",
          color = "Tsage.split",palette = wes_palette("Royal2", n = 2),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp3)+
  facet_wrap(~variable,scale="free",nrow=1)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()


# Clinical data:
clinical <- read.csv("../Sigs study cohort_clinical demographics.csv")
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

naive <- clinical[which(clinical$Is.Chemo.Naive=="true"),]$Sample
treated <- clinical[which(clinical$Is.Chemo.Naive=="false"),]$Sample
notannot <- clinical[which(clinical$Is.Chemo.Naive==""),]$Sample

length(unique(clin.sigs[which(clin.sigs$Sample %in% naive),]$Sample))
length(unique(clin.sigs[which(clin.sigs$Sample %in% treated),]$Sample))

clin.sigs$Is.Chemonaive <- "no"
clin.sigs[which(clin.sigs$Sample %in% naive),]$Is.Chemonaive <- "yes"
clin.sigs[which(clin.sigs$Sample %in% notannot),]$Is.Chemonaive <- ""

mycomp1 <- list(c("PD","Response/Stable"))
pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg.bychemostatus.pdf")
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg)&
                            (clin.sigs$variable %in% c("SBS17a","SBS17b"))&
                            (clin.sigs$Is.Chemonaive!="")),], 
          x = "TR.NeoAdj.Response.StableVsProg", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(Is.Chemonaive~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg.bychemostatus.pdf")
ggboxplot(clin.sigs[which(!is.na(clin.sigs$TR.NeoAdj.Response.StableVsProg)&
                            (clin.sigs$variable %in% c("SBS17a","SBS17b"))),], 
          x = "TR.NeoAdj.Response.StableVsProg", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(Is.Chemonaive~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

pdf("plots_7/cont.TR.NeoAdj.Response.StableVsProg.byTstage.pdf")
ggboxplot(clins.sigs.stage[which(!is.na(clins.sigs.stage$TR.NeoAdj.Response.StableVsProg)&
                            (clins.sigs.stage$variable %in% c("SBS17a","SBS17b"))),], 
          x = "TR.NeoAdj.Response.StableVsProg", y = "value",
          color = "TR.NeoAdj.Response.StableVsProg",palette = wes_palette("Royal2", n = 4),
          add = "jitter")+
  stat_compare_means(comparisons = mycomp1, label="p.signif")+
  facet_wrap(Tsage.split~variable,scale="free",nrow=2)+
  xlab("")+
  ylab("Signature exposure (%)")
dev.off()

c <- clin.sigs[which(!is.na(clin.sigs$RP.Tumour.Growth.c)&
                       (clin.sigs$variable %in% c("SBS17a","SBS17b"))),]
pd <- c[which(c$TR.NeoAdj.Response.StableVsProg=="PD"),]
grow <- c[which(c$RP.Tumour.Growth.c=="grow"),]
