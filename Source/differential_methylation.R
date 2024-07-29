
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

install.packages("TCGAbiolinks")

setwd("/Users/matanatmammadli/Desktop/Project/")

# Load TCGAbiolinks package
library(TCGAbiolinks)

HNSC_cohort <- read.delim('/Users/matanatmammadli/Desktop/Project/HNSC_aliquot.tsv', sep = '\t', header = TRUE, fill = TRUE)
HNSC_submitter_ids <- unique(HNSC_cohort$aliquot_submitter_id)

query_HNSC_methylation <- GDCquery(
  project = "TCGA-HNSC", 
  data.category = "DNA Methylation",
  data.type = "Masked Intensities",
  platform = "Illumina Human Methylation 450"
)


GDCdownload(query_HNSC_methylation, files.per.chunk=10)
betas <- GDCprepare(query_HNSC_methylation)


## Preprocess the Data

library(minfi)
library(ChAMP)



# Load the data from the "Masked Intensities" folder
data_directory <- "/Users/matanatmammadli/Desktop/Project/GDCdata/TCGA-HNSC/DNA_Methylation/Masked_Intensities/"


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
new_directory <- "/Users/matanatmammadli/Desktop/Project/GDCdata/TCGA-HNSC/DNA_Methylation/Masked_Intensities_pairs_to_keep/"

# Create the directory if it doesn't exist
if (!file.exists(new_directory)) {
  dir.create(new_directory)
}

# Assuming idat_files and pairs_to_keep are defined as in the previous example

# Iterate over pairs_to_keep and copy/move the files to the new directory
for (pair in pairs_to_keep) {
  file.copy(pair, new_directory)  # Use file.move() if you want to move instead of copy
}

for (filename in pairs_to_keep) {
  file.copy(file.path(data_directory, filename), file.path(new_directory, filename))
  # Use file.rename() instead of file.copy() if you want to move the files instead of copying them
}


# Create a MethylSet object to read IDAT files
methylation_data <- read.metharray.exp(base = new_directory, recursive = TRUE)

# Load required packages
library(minfi)

# Quality control (QC) on methylation data
qc <- qcReport(methylation_data)

# Data Preprocessing and Normalization
methylation_data_preprocessed <- preprocessQuantile(methylation_data)


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


# Create a MethylSet object to read IDAT files
methylation_data <- read.metharray.exp(base = data_directory, recursive = TRUE)



methylation_data <- read.metharray.exp(targets = complete_idat_files, BPPARAM = MulticoreParam())



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



library(ima)
idat_files <- list.files(data_directory, pattern = "_noid_Grn.idat$|_noid_Red.idat$", recursive = TRUE, full.names = TRUE)
rgSet <- IMA(idatFiles = idat_files)

# Optionally, specify extended = TRUE if you need the output in RGChannelSetExtended format
# rgSet <- read.metharray.sheet(base = data_directory, extended = TRUE, verbose = TRUE)

# Check the structure of rgSet
print(rgSet)


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
