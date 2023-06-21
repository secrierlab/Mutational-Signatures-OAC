################
## Load results and generate plot:

library(pheatmap)

load("data/mat.fc.oac.all.medians.RData")
m.primary <- mat.fc.oac.all
load("data/mat.fc.barretts.prim.medians.RData")
m.barr.prim <- mat.fc.oac.all
load("data/mat.fc.barr.separate.medians.RData")
m.barr.separate <- mat.fc.barr.all

# Correct barretts separate calculation:
load("data/mat.fc.barr.highestfuturegrade.RData")
m.barr.separate[3,] <- mat.fc.barr

rownames(m.primary) <- paste0("Primary.",rownames(m.primary))
rownames(m.barr.prim) <- paste0("Barretts.P.",rownames(m.barr.prim))
rownames(m.barr.separate) <- paste0("Barretts.",rownames(m.barr.separate))

m <- rbind(m.barr.separate,m.barr.prim,m.primary)
m <- m[,which(!(colnames(m) %in% c("SBS28","SBS3")))]

### Reorder rows:
m <- m[c(2,3,4,5,1,6,10,9,8,11,7),]
unlist(m)

pdf("plots_supp/clinicalCorrelations.Heatmap.medians.redone.pdf",w=5.8,h=4)
breaksList1 = seq(-round(max(m,na.rm=TRUE)), round(max(m,na.rm=TRUE)), by = 1)
breaksList2 <- c(seq(-1,0,by=0.05),seq(0.025,1,by=0.025),seq(1.1,3,by=0.1))
pheatmap(m, cluster_rows = FALSE,breaks = breaksList2,
        color=colorRampPalette(c("navy", "white", "#FEFFA5","#FC9E4F","#A63D40"))(length(breaksList2)))
dev.off()

write.csv(m, file="plots_supp/clin.corrs.csv")
