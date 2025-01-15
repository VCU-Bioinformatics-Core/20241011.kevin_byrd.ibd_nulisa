library(readxl)
library(sets)
library(tidyverse)
library(limma)
library(gplots)
library(ggplot2)
library(dplyr)

################################################################################
# Load the datasets
################################################################################
fn = "results/comp_data/protein_levels.npq.tsv"
protein_data <- read.table(fn, header=TRUE, check.names=FALSE, sep = "\t")

fn = "results/comp_data/clinical_data.tsv"
clinical_data <- read.table(fn, header=TRUE, sep="\t")

# remove alamar samples completely
clinical_data <- clinical_data[clinical_data$ibd_diagnosis != "Alamar_Sample_Control",] 
protein_data <- protein_data[, clinical_data$sample_id]




################################################################################
# Run limma with control as the anchor
################################################################################

# set groups and design
group <- factor(clinical_data$ibd_diagnosis, levels=c("Control", "CD", "UC", "IBD-U"))
design <- model.matrix(~ group)

# # fit linear model for data
fit <- lmFit(protein_data, design)

# apply empirical bayes
fit <- eBayes(fit)

# get results
results_cd <- topTable(fit, adjust="BH", coef="groupCD", number=Inf)
results_uc <- topTable(fit, adjust="BH", coef="groupUC", number=Inf)
results_ibdu <- topTable(fit, adjust="BH", coef="groupIBD-U", number=Inf)

# get reorganized results results
binary_results <- decideTests(fit, method="global", adjust.method="BH")
binary_results <- as.data.frame(binary_results)

################################################################################
# Visualize the results (with control as the anchor)
################################################################################

# Define the function
create_volcano_plot <- function(results) {
  # Create a logical column for significance
  results$significant <- results$P.Value <= 0.05
  
  # Identify the top 10 genes based on P.Value
  top_genes <- results %>%
    arrange(P.Value) %>%
    slice_head(n = 10)
  
  # Create a new column for color based on logFC
  results$color <- ifelse(results$significant & results$logFC < 0, "blue", 
                          ifelse(results$significant & results$logFC > 0, "red", "grey"))
  
  # Create the volcano plot
  ggplot(results, aes(x = logFC, y = -log10(P.Value), color = color)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_identity() +  # Use the colors defined in the data
    theme_minimal() +
    labs(title = "Volcano Plot", x = "Log2 Fold Change (logFC)", y = "-Log10 P-value") +
    theme(legend.position = "none") +  # Hide legend if not needed
    geom_text(data = top_genes, aes(label = rownames(top_genes)), 
              vjust = -0.5, hjust = 0.5, size = 3, color = "black")  # Annotate top 10 genes
}

create_volcano_plot(results_cd)
create_volcano_plot(results_uc)
create_volcano_plot(results_ibdu)
# Draw venn diagram of results
#vennDiagram(decideTests(fit, method="global", adjust.method="BH"))



################################################################################
# Run limma with CD data only
################################################################################

crohn_clinical_data <- clinical_data[clinical_data$ibd_diagnosis == "CD",]
crohn_protein_data <- protein_data[, crohn_clinical_data$sample_id]

# set groups and design
group <- factor(crohn_clinical_data$disease_activity_indicator, levels=c("In-active Disease", "Active Disease"))
design <- model.matrix(~ group)

# # fit linear model for data
fit <- lmFit(crohn_protein_data, design)

# apply empirical bayes
fit <- eBayes(fit)

# get results
results_da <- topTable(fit, adjust="BH", coef="groupActive Disease", number=Inf)

# get reorganized results results
binary_results <- decideTests(fit, method="global", adjust.method="BH")
binary_results <- as.data.frame(binary_results)

################################################################################
# Visualize the results (with CD data only)
################################################################################

create_volcano_plot(results_da)

# Draw venn diagram of results
vennDiagram(decideTests(fit, method="global", adjust.method="BH"))



################################################################################
# Run Limma with UC data only
################################################################################

uc_clinical_data <- clinical_data[clinical_data$ibd_diagnosis == "UC",]
uc_protein_data <- protein_data[, uc_clinical_data$sample_id]

# set groups and design
group <- factor(uc_clinical_data$disease_activity_indicator, levels=c("In-active Disease", "Active Disease"))
design <- model.matrix(~ group)

# # fit linear model for data
fit <- lmFit(uc_protein_data, design)

# apply empirical bayes
fit <- eBayes(fit)

# get results
results_da <- topTable(fit, adjust="BH", coef="groupActive Disease", number=Inf)

# get reorganized results results
binary_results <- decideTests(fit, method="global", adjust.method="BH")
binary_results <- as.data.frame(binary_results)

################################################################################
# Visualize the results (with UC data only)
################################################################################

create_volcano_plot(results_da)

# Draw venn diagram of results
vennDiagram(decideTests(fit, method="global", adjust.method="BH"))
