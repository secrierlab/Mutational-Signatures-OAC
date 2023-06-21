###################
### This script analyses dndscv output.

library(gdata)
library(ggpubr)
library(ggforce)

# Load cancer gene census:
cancergenes <- read.csv("Census_allMon Oct 19 13_46_22 2020.csv")$Gene.Symbol

# Load S17 output:
load("output.2021/dndsout.s17.RData")

sel_cv = dndsout.s17$sel_cv
print(head(sel_cv), digits = 3)
all_genes.s17<- sel_cv[,c("gene_name","qallsubs_cv",
                           "wmis_cv","wnon_cv","wspl_cv")]

signif_genes.s17 = sel_cv[sel_cv$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                     "wmis_cv","wnon_cv","wspl_cv")]
rownames(signif_genes.s17) = NULL
print(signif_genes.s17)
write.csv(signif_genes.s17, file="output.2021/signif_genes.s17.csv", row.names = FALSE)

signif_genes.s17$gene_short <- sapply(signif_genes.s17$gene_name,
                                  function(x) strsplit(x,"\\.")[[1]][1])
signif_genes.s17.onlycancer <- signif_genes.s17[which(toupper(signif_genes.s17$gene_short) %in% 
                                                toupper(cancergenes)),]
write.csv(signif_genes.s17.onlycancer, file="output.2021/signif_genes.s17.onlyCancerGenes.csv", row.names = FALSE)


# Load not-S17 output:
load("output.2021/dndsout.s17other.RData")

sel_cv = dndsout.s17other$sel_cv
print(head(sel_cv), digits = 3)
all_genes.nots17<- sel_cv[,c("gene_name","qallsubs_cv",
                          "wmis_cv","wnon_cv","wspl_cv")]


signif_genes.nots17 = sel_cv[sel_cv$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                    "wmis_cv","wnon_cv","wspl_cv")]
rownames(signif_genes.nots17) = NULL
print(signif_genes.nots17)
write.csv(signif_genes.nots17, file="output.2021/signif_genes.nots17.csv", row.names = FALSE)

signif_genes.nots17$gene_short <- sapply(signif_genes.nots17$gene_name,
                                      function(x) strsplit(x,"\\.")[[1]][1])
signif_genes.nots17.onlycancer <- signif_genes.nots17[which(toupper(signif_genes.nots17$gene_short) %in% 
                                                        toupper(cancergenes)),]
write.csv(signif_genes.nots17.onlycancer, file="output.2021/signif_genes.nots17.onlyCancerGenes.csv", row.names = FALSE)


### Now plot one against the other:

unionsigs <- unique(c(signif_genes.s17$gene_name,
                                 signif_genes.nots17$gene_name))
df.s17nots17 <- merge(all_genes.s17[which(all_genes.s17$gene_name %in%
                                            unionsigs),],
                     all_genes.nots17[which(all_genes.nots17$gene_name %in%
                                              unionsigs),],
                     by.x="gene_name",by.y="gene_name")
df.s17nots17$log10Pval_s17 <- -log10(df.s17nots17$qallsubs_cv.x)
df.s17nots17$log10Pval_nots17 <- -log10(df.s17nots17$qallsubs_cv.y)
df.s17nots17[which(is.infinite(df.s17nots17$log10Pval_s17)),]$log10Pval_s17 <- 12
df.s17nots17[which(is.infinite(df.s17nots17$log10Pval_nots17)),]$log10Pval_nots17 <- 12
df.s17nots17$Significance <- apply(df.s17nots17[,c("log10Pval_s17",
                                                 "log10Pval_nots17")],1,
                                  function(x) ifelse((x[1]>1)&(x[2]>1),"Common",ifelse(
                                    x[1]>1,"S17-specific",
                                    "Other signature-specific")))

# only cancer genes:
df.s17nots17$gene_short <- sapply(df.s17nots17$gene_name,
                                  function(x) strsplit(x,"\\.")[[1]][1])
df.s17nots17.onlycancer <- df.s17nots17[which(toupper(df.s17nots17$gene_short) %in% 
                                                toupper(cancergenes)),]

pdf("plots.output2021/signif.S17_vs_notS17.pdf",w=6,h=6)
ggscatter(df.s17nots17[which((df.s17nots17$log10Pval_s17>3)|
                         (df.s17nots17$log10Pval_nots17>3)),], 
          y = "log10Pval_s17", x = "log10Pval_nots17",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          color="Significance",
          label = "gene_name", repel = TRUE)+
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=0.3)+
  geom_vline(xintercept=1, linetype="dashed", 
             color = "red", size=0.3)+
  ylab("S17 dominant (-log10 p-value)")+
  xlab("Other signatures dominant (-log10 p-value)")
#facet_zoom(ylim = c(0, 5),xlim=c(0,5))
dev.off()

pdf("plots.output2021/signif.S17_vs_notS17.onlyCancerGenes.pdf",w=6,h=6)
ggscatter(df.s17nots17.onlycancer, 
          y = "log10Pval_s17", x = "log10Pval_nots17",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          color="Significance",
          label = "gene_name", repel = TRUE)+
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=0.3)+
  geom_vline(xintercept=1, linetype="dashed", 
             color = "red", size=0.3)+
  ylab("S17 dominant (-log10 p-value)")+
  xlab("Other signatures dominant (-log10 p-value)")
#facet_zoom(ylim = c(0, 5),xlim=c(0,5))
dev.off()

###############
#### Pathway enrichment:

library(pathfindR)

signif_genes.s17.input <- all_genes.s17[,1:2]
output_df.s17 <- run_pathfindR(signif_genes.s17.input, 
                                     p_val_threshold = 0.1,
                               min_gset_size=3,
                               max_gset_size=300,
                                     pin_name_path="GeneMania")
knitr::kable(head(output_df.s17, 10))

pdf("plots.output2021/S17.enrichedPathways.GeneMania.pdf")
enrichment_chart(result_df = output_df.s17, 
                 top_terms = 20)
dev.off()

signif_genes.nots17.input <- all_genes.nots17[,1:2]
output_df.nots17 <- run_pathfindR(signif_genes.nots17.input, 
                               p_val_threshold = 0.1,
                               min_gset_size=3,
                               max_gset_size=300,
                               pin_name_path="GeneMania")
knitr::kable(head(output_df.nots17, 10))

pdf("plots.output2021/notS17.enrichedPathways.GeneMania.pdf")
enrichment_chart(result_df = output_df.nots17, 
                 top_terms = 20)
dev.off()


### These are all the pathways enriched in both groups. What about the specific pathways?

output_df.nots17 <- run_pathfindR(signif_genes.nots17.input[which(signif_genes.nots17.input$gene_name %in%
                                                                 df.s17nots17[which(df.s17nots17$Significance == "Other signature-specific"),]$gene_name),1:2], 
                                  p_val_threshold = 0.1,
                                  min_gset_size=3,
                                  max_gset_size=300,
                                  pin_name_path="GeneMania")
pdf("plots.output2021/notS17_specific.enrichedPathways.GeneMania.pdf")
enrichment_chart(result_df = output_df.nots17, 
                 top_terms = 20)
dev.off()

pdf("plots.output2021/notS17_specific.graph.GeneMania.pdf",w=10,h=10)
term_gene_graph(result_df = output_df.nots17, use_description = TRUE)
dev.off()

