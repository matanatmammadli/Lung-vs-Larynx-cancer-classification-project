# Install necessary packages
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# Set working directory
setwd("/Users/matanatmammadli/Desktop/Project/")

# Load TCGAbiolinks package
library(TCGAbiolinks)

HNSC_cohort <-
  read.delim('/Users/matanatmammadli/Desktop/Project/HNSC_aliquot.tsv', sep
             = '\t', header = TRUE, fill = TRUE)
HNSC_submitter_ids <- unique(HNSC_cohort$aliquot_submitter_id)

query_HNSC_methylation <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "DNA Methylation",
  data.type = "Masked Intensities",
  platform = "Illumina Human Methylation 450",
  barcode = HNSC_submitter_ids
)


GDCdownload(query_HNSC_methylation, files.per.chunk=10)
betas <- GDCprepare(query_HNSC_methylation)


## Preprocess the Data

library(minfi)
library(ChAMP)



# Load the data from the "Masked Intensities" folder
data_directory <-
  "/Users/matanatmammadli/Desktop/Project/GDCdata/TCGA-HNSC/DNA_Methylation/Masked_Intensities/"




# List all IDAT files ending with _noid_Grn.idat or _noid_Red.idat across all subdirectories
#idat_files <- list.files(data_directory, pattern ="_noid_Grn.idat$|_noid_Red.idat$", recursive = TRUE)

idat_files <- unique(idat_files)

idat_files <- list.files(data_directory, pattern = "\\.idat$", recursive =
                           TRUE)

idat_files <- unique(idat_files)


# Create a MethylSet object to read IDAT files
methylation_data_larynx <- read.metharray.exp(base = data_directory,
                                              recursive = TRUE)

# Load required packages
library(minfi)

# Quality control (QC) on methylation data
qc <- qcReport(methylation_data_larynx)

# Data Preprocessing and Normalization
methylation_data_preprocessed_larynx <-
  preprocessQuantile(methylation_data_larynx)


# Load additional packages for statistical analysis
library(limma)
library(ggplot2)

betas <- getBeta(methylation_data_preprocessed_larynx)

densityPlot(betas, main="Density Plot of Beta Values for larynx dataset")


# Step 1: Define Groups
# Example: Assuming the first half are cases and the second half are controls
group <- factor(c(rep("Case", ncol(betas) / 2), rep("Control", ncol(betas)
                                                    / 2)))

# Step 2: Create the Design Matrix
design <- model.matrix(~ group)

# Step 3: Perform PCA on Beta Values
pca_results <- prcomp(t(betas), scale. = TRUE)

# Step 4: Plot the PCA Results
# Convert group factor to color
group_colors <- as.numeric(group)
color_map <- c("green", "red")  # Customize the colors as needed

plot(pca_results$x[, 1], pca_results$x[, 2], col=color_map[group_colors],
     pch=19,
     xlab="PC1", ylab="PC2", main="PCA of Beta Values larynx")
legend("topright", legend=levels(group), col=color_map, pch=19)


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
top_probes <- head(results$ID, 100)  # Adjust as per your results object
structure
betas_top <- getBeta(methylation_data)[top_probes, ]


#  using pheatmap package
library(pheatmap)
betas_heatmap <- betas

# Selecting the first 100 rows and all columns
subset_betas <- betas_heatmap[1:100, ]

# Plot subset data
heatmap(as.matrix(subset_betas), scale = "row", show_rownames = FALSE)


#correlation plot
library(corrplot)
cor_matrix <- cor(t(subset_betas))
corrplot(cor_matrix, method="color", tl.cex = 0.35)




----------------------
  
  ##Machine Learning analysis on methylation data
  
# Assuming 'betas' is your preprocessed methylation data and 'group' is the target variable
# Convert the 'betas' data frame to a matrix for compatibility with ML algorithms
betas_matrix <- as.matrix(betas)

# Create a sample group variable (e.g., Case and Control)
# Replace this with your actual labels
set.seed(123)
group <- factor(c(rep("Case", ncol(betas) / 2), rep("Control", ncol(betas)
                                                    / 2)))

# Split into training and test sets
train_index <- sample(seq_len(ncol(betas_matrix)), size = 0.7 *
                        ncol(betas_matrix))

train_data <- betas_matrix[, train_index]
test_data <- betas_matrix[, -train_index]

train_labels <- group[train_index]
test_labels <- group[-train_index]


# Feature selection using limma

fit <- lmFit(betas, design)
fit <- eBayes(fit)

topTable <- topTable(fit, number = nrow(train_data))

# Select top features
selected_features <- rownames(topTable)[1:1000]

# Subset train_data_selected and test_data_selected by selected features
train_data_selected <- train_data[selected_features, ]
test_data_selected <- test_data[selected_features, ]

train_data_selected <- as.matrix(train_data_selected)  # Convert to numeric matrix if necessary


library(randomForest)

# Example for training Random Forest with selected features
rf_model <- randomForest(t(train_data_selected), train_labels)


## DO NOT RUN THIS PART
library(doParallel)
library(randomForest)
library(foreach)

# Set up parallel backend
num_cores <- detectCores() - 1  # Leave one core for other tasks
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Train Random Forest with parallel processing
rf_model <- foreach(ntree = rep(100, num_cores), .combine = combine,
                    .packages = 'randomForest') %dopar% {
                      randomForest(t(train_data_selected), train_labels, ntree = ntree)
                    }

# Stop the cluster
stopCluster(cl)

## TILL HERE






# Predict using the trained model on test data
predictions <- predict(rf_model, newdata = t(test_data_selected))

# Assess accuracy
accuracy <- mean(predictions == test_labels)
cat("Accuracy on test set:", accuracy, "\n")

# Confusion matrix on test set
confusion <- table(Actual = test_labels, Predicted = predictions)
print(confusion)

# Calculate class-specific metrics
library(caret)
class_metrics <- confusionMatrix(predictions, test_labels)
print(class_metrics)



# Assuming predictions is a factor with levels "Case" and "Control"
# Convert predictions to numeric probabilities
prob_case <- ifelse(predictions == "Case", 1, 0)  # Assuming "Case" is the
positive class
names(prob_case) <- names(predictions)  # Set the names to match the
original data

# Check the structure of prob_case to ensure it's numeric
str(prob_case)

# Compute ROC curve
library(pROC)
roc_curve <- roc(test_labels, prob_case)

# Print the ROC curve
print(roc_curve)

# Plot the ROC curve
plot(roc_curve, main = "ROC Curve", col = "blue")

# Compute and print AUC
auc <- auc(roc_curve)
print(paste("AUC:", auc))

# Example: Assuming you need to create a binary target variable
# Let's say based on the value of a specific feature (e.g., cg27227537)


# Get feature importance from the model
importance <- importance(rf_model)
print(importance)

# Plot feature importance
varImpPlot(rf_model, cex = 0.6)


## HYPEROPARAMETER TUNING

# Define parameter grid for tuning
param_grid <- expand.grid(
  ntree = c(100, 200, 300),  # Adjust the number of trees
  mtry = c(10, 20, 30)       # Adjust the number of variables tried at
  each split
)

# Train using caret's train function with cross-validation
set.seed(123)
rf_model_tuned <- train(
  x = t(train_data_selected), y = train_labels,
  method = "rf", trControl = trainControl(method = "cv", number = 5),
  tuneGrid = param_grid
)

# Print tuned parameters
print(rf_model_tuned)

# Assess performance of the tuned model
predictions_tuned <- predict(rf_model_tuned, newdata = t(test_data_selected))
accuracy_tuned <- mean(predictions_tuned == test_labels)
cat("Tuned model accuracy on test set:", accuracy_tuned, "\n")







----------------------





