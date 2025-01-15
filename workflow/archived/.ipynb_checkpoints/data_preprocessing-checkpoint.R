library(readxl)
library(tidyr)
library(tidyverse)


################################################################################
# Process the protein levels
################################################################################

# load assay data
fn <- "results/raw/P-000458_ADA_NULISAseq_Inflammation Panel_1-NPQ Counts_2-Target Detectability_3-Sample Information_2024_08_26.xlsx"
data <- read_excel(fn, sheet = "NPQ Counts")

# pivot the data to make a matrix
wide_data <- data %>% pivot_wider(id_cols = SampleName, names_from = Target, values_from = NPQ)

# transpose the data to be in the correct format for limma
# rows = genes
# columns = samples
t_data = t(wide_data)

# clean up the first row
colnames(t_data) <- t_data[1,]
t_data <- as.data.frame(t_data[-1,])

# get the column names
sample_names <- colnames(t_data)

# remove the pattern [A-Z]_[0-9]+
new_sample_names <- gsub("[A-Z]_[0-9]+_", "", sample_names)

# assign the new column names back to the data frame
colnames(t_data) <- new_sample_names

# # trim whitespace and convert to numeric
# t_data <- apply(t_data, 2, function(x) as.numeric(trimws(x)))

write.table(t_data, "results/comp_data/protein_levels.npq.tsv", sep = "\t", col.names = TRUE, quote=FALSE)


################################################################################
# Process the clinical data
################################################################################

# load clinical data
fn <- "results/raw/ADA_IBD_Saliva_Biospecimen Manifest Form-for-NULISA_241022.xlsx"
clinical_data <- read_excel(fn, sheet = "Aliquot Information")

# rename columns for programming use
transform_string <- function(x) {
  x <- tolower(x)            # Convert to lowercase
  x <- gsub(" ", "_", x)     # Replace spaces with "_"
  x <- gsub("\\(", "_", x)   # Replace "(" with "_"
  x <- gsub("\\)", "", x)    # Remove ")"
  return(x)
}

# Apply the function to the vector
transformed_cols <- sapply(as.vector(colnames(clinical_data)), transform_string)
colnames(clinical_data) <- transformed_cols

# remove nan samples
clinical_data <- clinical_data[!is.na(clinical_data$project_name),]
clinical_sample_names = clinical_data$original_subject_id

# Create a new row as a dataframe
new_row <- data.frame(
  original_subject_id = c("SC_Rep01", "SC_Rep02", "SC_Rep03"),
  sample_id = c("SC_Rep01", "SC_Rep02", "SC_Rep03"),
  ibd_diagnosis = c("Alamar_Sample_Control", "Alamar_Sample_Control", "Alamar_Sample_Control"),
  disease_activity = c("N/A", "N/A", "N/A")
)

# Add the new row to the dataframe
clinical_data <- bind_rows(clinical_data, new_row)

# order the dataframe based on the order of new sample names
clinical_data <- clinical_data[match(new_sample_names, clinical_data$original_subject_id), ]


# def indicator function
get_indicator <- function(x, check_list, categories) {
  
  if (x %in% check_list) {
    return(categories[[1]])
  } 
  else {
    return(categories[[2]])
  }
}

# add ibd disease indicator
ibd_check_list = c("CD", "UC", "IBD-U")
ibd_indicator <- sapply(as.vector(clinical_data$ibd_diagnosis), get_indicator, check_list=ibd_check_list, categories=c("IBD Super Group", "Control Super Group"))
clinical_data$ibd_indicator <- ibd_indicator

# add disease activity indicator
da_check_list = c("Moderate", "Mild")
da_indicator <- sapply(as.vector(clinical_data$disease_activity), get_indicator, check_list=da_check_list, categories=c("Active Disease", "In-active Disease"))
clinical_data$disease_activity_indicator <- da_indicator

write.table(clinical_data, "results/comp_data/clinical_data.tsv", sep = "\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
