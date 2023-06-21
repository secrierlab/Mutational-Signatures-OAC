####################
### Checking TMB 

library(ggpubr)

load("data/tmb.annot.RData")
write.csv(tmb.annot,file="plots_supp/tmb.csv")

# Plot of TMB:

my_comparisons <- list(c("Barretts","PrimaryTumour"),
                       c("PrimaryTumour","Metastasis"),
                       c("Metastasis","Barretts"))
pdf("plots.revision/TMBcompared.pdf",w=4,h=5)
ggboxplot(tmb.annot, x = "Category", y = "TMB",
          color = "Category", palette ="npg",
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("")+
  ylab("Tumour mutational burden (log10)")
dev.off()

