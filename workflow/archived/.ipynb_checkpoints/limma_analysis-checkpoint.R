library(readxl)
library(sets)
library(tidyverse)
library(limma)
library(gplots)

################################################################################
# Load the datasets
################################################################################
fn = "results/comp_data/protein_levels.npq.tsv"
protein_data = read.table(fn, header=TRUE, sep = "\t")

fn = "results/comp_data/clinical_data.tsv"
clinical_data = read.table(fn, header=TRUE, sep="\t")

################################################################################
# Run Limma
################################################################################

# set groups and design
group <- factor(as.vector(clinical_data$ibd_diagnosis))
design <- model.matrix(~ group)

# # fit linear model for data
fit <- lmFit(protein_data, design)

# apply empirical bayes
fit <- eBayes(fit)

# get results
results <- topTable(fit, coef="groupIBD-U", number=Inf)

################################################################################
# Process the visualize the results 
################################################################################

# generate volcano plot
volcanoplot(fit, coef="groupIBD-U", highlight=10, names=results$Gene)

# generate ma plot
plotMA(fit, coef="groupIBD-U")

# # Assuming 'exprs' is your expression matrix and 'results' is the topTable output
# top_genes <- rownames(results_ibdu)[1:50]  # Select top 50 genes
# heatmap.2(t_data[top_genes, ], scale="row", trace="none", dendrogram="both", col=bluered(100))

# Assuming 'fit' is your model fit object
vennDiagram(decideTests(fit))



