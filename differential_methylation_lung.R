
# Install necessary packages
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# Set working directory
setwd("/Users/matanatmammadli/Desktop/Project/")

# Load TCGAbiolinks package
library(TCGAbiolinks)

# Read the LUSC aliquot information
LUSC_cohort <- read.delim('/Users/matanatmammadli/Desktop/Project/LUSC_aliquot.tsv', sep = '\t', header = TRUE, fill = TRUE)

# Extract unique submitter IDs
LUSC_submitter_ids <- unique(LUSC_cohort$aliquot_submitter_id)

# Query the LUSC methylation data
query_LUSC_methylation <- GDCquery(
  project = "TCGA-LUSC", 
  data.category = "DNA Methylation",
  data.type = "Masked Intensities",
  platform = "Illumina Human Methylation 450",
  barcode = LUSC_submitter_ids
)

# Download the LUSC methylation data
GDCdownload(query_LUSC_methylation, files.per.chunk = 10)

# Prepare the data
betas_LUSC <- GDCprepare(query_LUSC_methylation)



## Preprocess the Data

library(minfi)
library(ChAMP)



# Load the data from the "Masked Intensities" folder
data_directory <- "/Users/matanatmammadli/Desktop/Project/GDCdata/TCGA-LUSC/DNA_Methylation/Masked_Intensities/"


# Load data using ChAMP
myLoad <- champ.load(directory = data_directory,
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="450K")

# Check the loaded data
head(myLoad$beta)
head(myLoad$pd)


# List all IDAT files ending with _noid_Grn.idat or _noid_Red.idat across all subdirectories
idat_files <- list.files(data_directory, pattern = "_noid_Grn.idat$|_noid_Red.idat$", recursive = TRUE)

idat_files <- unique(idat_files)

idat_files <- list.files(data_directory, pattern = "\\.idat$", recursive = TRUE)

# Extract sample identifiers (assuming naming convention *_Grn.idat and *_Red.idat)
sample_ids <- unique(gsub("_Grn.idat|_Red.idat", "", basename(idat_files)))

# Group files by sample ID
idat_groups <- lapply(sample_ids, function(id) {
  grep(id, idat_files, value = TRUE)
})


set.seed(123)  # For reproducibility, set seed for random sampling
pairs_to_keep <- sample(idat_groups, size = length(idat_groups) / 2)


# Define the path to the new directory
new_directory <- "/Users/matanatmammadli/Desktop/Project/GDCdata/TCGA-LUSC/DNA_Methylation/Masked_Intensities_pairs_to_keep/"

# Create the directory if it doesn't exist
if (!file.exists(new_directory)) {
  dir.create(new_directory)
}

# Assuming idat_files and pairs_to_keep are defined as in the previous example



for (filename in pairs_to_keep) {
  file.copy(file.path(data_directory, filename), file.path(new_directory, filename))
  # Use file.rename() instead of file.copy() if you want to move the files instead of copying them
}


# Create a MethylSet object to read IDAT files
methylation_data <- read.metharray.exp(base = new_directory, recursive = TRUE)


# Load required packages
library(minfi)

# Quality control (QC) on methylation data
qc_lung <- qcReport(methylation_data)

# Data Preprocessing and Normalization
methylation_data_preprocessed <- preprocessQuantile(methylation_data)


# Extract the assay data (beta values)
beta_values <- assay(methylation_data_preprocessed)
print(dim(beta_values))  # Print the dimensions of the matrix: [1] 485512    190
print(colnames(beta_values))  # Print the column names

# Extract the sample metadata
sample_metadata <- colData(methylation_data_preprocessed)
print(sample_metadata)

# Verify the structure of the preprocessed data object
str(methylation_data_preprocessed)



# Load additional packages for statistical analysis
library(limma)
library(ggplot2)

betas <- getBeta(methylation_data_preprocessed)

# Design matrix: Define experimental groups (e.g., case vs control)
group <- factor(c(rep("Case", ncol(methylation_data_preprocessed) / 2), rep("Control", ncol(methylation_data_preprocessed) / 2)))
design <- model.matrix(~ group)

# Perform differential methylation analysis using limma
fit <- lmFit(betas, design)
fit <- eBayes(fit)

# Extract differential methylation results
results <- topTable(fit, coef=1, number=Inf)



# Density plot of beta values
densityPlot(betas, main="Density Plot of Beta Values")


# Step 1: Define Groups
# Example: Assuming the first half are cases and the second half are controls
group <- factor(c(rep("Case", ncol(betas) / 2), rep("Control", ncol(betas) / 2)))

# Step 2: Create the Design Matrix
design <- model.matrix(~ group)

# Step 3: Perform PCA on Beta Values
pca_results <- prcomp(t(betas), scale. = TRUE)

# Step 4: Plot the PCA Results
# Convert group factor to color
group_colors <- as.numeric(group)
color_map <- c("green", "red")  # Customize the colors as needed

plot(pca_results$x[, 1], pca_results$x[, 2], col=color_map[group_colors], pch=19, 
     xlab="PC1", ylab="PC2", main="PCA of Beta Values")
legend("topright", legend=levels(group), col=color_map, pch=19)




# Plot QC results
plotQC(qc)

# Volcano plot
volcanoPlot(fit, coef=1)

# MA plot
plotMA(fit, coef=1)

# Heatmap of top differentially methylated probes
top_probes <- head(results$ID, 100)  # Adjust as per your results object structure
betas_top <- getBeta(methylation_data)[top_probes, ]


# Example using pheatmap package
library(pheatmap)
betas_heatmap <- betas  # Replace 'betas' with your beta values
heatmap(as.matrix(betas_heatmap), scale="row", show_rownames=FALSE)

# Example using ggplot2
library(ggplot2)
df <- as.data.frame(betas)  # Replace 'betas' with your beta values
df$group <- group  # Replace 'group' with your group labels
ggplot(df, aes(x=group, y=value)) +
  geom_boxplot() +
  labs(x="Group", y="Beta Value", title="Boxplot of Beta Values")

# Example scatter plot
plot(betas$group1, betas$group2, xlab="Group 1", ylab="Group 2", main="Scatter Plot")

# Example using ggplot2
ggplot(df, aes(x=value, fill=group)) +
  geom_density(alpha=0.5) +
  labs(x="Beta Value", y="Density", title="Density Plot of Beta Values")

# Example using ggplot2 and corrplot packages
library(corrplot)
cor_matrix <- cor(t(betas))  # Compute correlation matrix
corrplot(cor_matrix, method="color")

# Example using limma MDSplot
MDSplot(fit, labels=group, col=c("blue", "red"))


----------------------
  


# Create a data.frame with Basename column containing the full paths of IDAT files
sample_sheet <- data.frame(Basename = unique(idat_files), stringsAsFactors = FALSE)



# Check loaded data
head(rgSet$beta)


# Check if sample_sheet is correctly populated
print(head(targets))  # Print the first few rows to verify



targets <- data.frame(Basename = sample_sheet$Basename, stringsAsFactors = FALSE)

rgSet <- read.metharray.exp(base = data_directory, targets = targets)

unique_targets <- targets[!duplicated(targets$Basename), ]

unique_targets <- targets  # Create a copy of targets dataframe
unique_targets$Basename <- make.unique(unique_targets$Basename)
rgSet <- read.metharray.exp(base = data_directory, targets = unique_targets)


library(ENmix)
# Load data
beta_values <- ENmixLoadData(directory = data_directory)



# Load necessary libraries
library(minfi)

# Define the base directory
base_directory <- "/Users/matanatmammadli/Desktop/Project/GDCdata/TCGA-HNSC/DNA_Methylation/Masked_Intensities"

# Step 1: List all IDAT files across subdirectories
idat_files <- list.files(base_directory, pattern = "_noid_Grn.idat$|_noid_Red.idat$", recursive = TRUE, full.names = TRUE)

# Step 2: Create sample sheet
sample_sheet <- data.frame(Basename = idat_files, stringsAsFactors = FALSE)

# Step 3: Read IDAT files into RGChannelSet
rgSet <- read.metharray(basenames = sample_sheet$Basename, extended = TRUE, verbose = TRUE)
