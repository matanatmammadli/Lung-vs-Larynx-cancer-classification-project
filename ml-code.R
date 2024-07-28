#install.packages("caret")
#install.packages("randomForest")

library(caret)
library(randomForest)
library(pROC)


# load feature table (and set first row with sampleIDs as rownames)
featureTable <- read.delim("/Users/matanatmammadli/Desktop/Project/featureTab.csv", sep=",")
rownames(featureTable) <- featureTable$X
featureTable <- featureTable[ , -1]

# factorize target column
#featureTable$cancersite <- ifelse(featureTable$cancersite == "larynx", 0, 1)
featureTable$cancersite <- as.factor(featureTable$cancersite)


# split into train and test set
trainIndex <- createDataPartition(y = featureTable$cancersite, p = 0.7, list = FALSE)
data_train <- featureTable[trainIndex, ]
data_test <- featureTable[-trainIndex, ]
data_test$cancersite <- as.factor(data_test$cancersite)
data_train$cancersite <- as.factor(data_train$cancersite)

# cross validation
ctrl <- trainControl(method = "cv", number = 10)#, classProbs = TRUE)
#ctrl <- trainControl(method = "oob", classProbs = TRUE)


# Random Forest
rf_model <- train(cancersite ~ ., data = data_train, method = "rf", trControl = ctrl)
predictions_rf <- predict(rf_model, newdata = data_test)
confusionMatrix(predictions_rf, data_test$cancersite)
roc_rf <- roc(ifelse(predictions_rf == "lung", 1, 0), ifelse(data_test$cancersite == "lung", 1, 0))
plot(roc_rf, main = "ROC Curve - Random Forest Model")

importance <- varImp(rf_model, scale = FALSE)
print(importance)


# Gradient Boosting Machines
gbm_model <- train(cancersite ~ ., data = data_train, method = "gbm", trControl = ctrl)
predictions_gbm <- predict(gbm_model, newdata = data_test)
confusionMatrix(predictions_gbm, data_test$cancersite)
roc_gbm <- roc(ifelse(predictions_gbm == "lung", 1, 0), ifelse(data_test$cancersite == "lung", 1, 0))
plot(roc_gbm, main = "ROC Curve - Gradient Boosting Machine")


# Support Vector Machines
svm_model <- train(cancersite ~ ., data = data_train, method = "svmLinear", trControl = ctrl)
predictions_svm <- predict(svm_model, newdata = data_test)
confusionMatrix(predictions_svm, data_test$cancersite)
roc_svm <- roc(ifelse(predictions_svm == "lung", 1, 0), ifelse(data_test$cancersite == "lung", 1, 0))
plot(roc_svm, main = "ROC Curve- Support Vector Machine")
