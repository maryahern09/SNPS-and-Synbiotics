############################################
############ SNPS Pipeline #################
#### Step 2: Regularization- Glucose #######
###### created by Mary 7-28-2025 ###########
############################################

#STEPS
#-1 Get data ready for PLINK
#0 Run SNP QC in PLINK 
#1.Data curation in R
#2 Regularization for coefficient reduction and feature selection using glmnet in R
#3 PRS calculation using PLINK through bigsnpr
#4 linear models using nlme or lme4 in R

#this was run on the HPC

#install packages ####
library(data.table)
library(tidyverse)
library(glmnet)
#install.packages("caret")
library(caret)
install.packages("car")
library(car)

#### Step 2: Regularization Models ####
#### Part 1: Parameter tuning with train/test/validate splits via outcomes (1 for glucose and 1 for A1c) ####

#load data
data_glucose <- fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/GlycogenesSNPs_qc_PLINKandR.csv")

#Split subset data into train, validate, and test splits
# Create initial split for training and remaining data (validation + test)
train_indices <- createDataPartition(y = data_glucose$glucose, 
                                     p = 0.6, list = FALSE) # 60% for training
train_data <- data_glucose[train_indices, ]
remaining_data <- data_glucose[-train_indices, ]

# Split the remaining data into validation and test sets
# Assuming a 50/50 split of the remaining data (15% validation, 15% test)
validation_indices <- createDataPartition(y = remaining_data$glucose, 
                                          p = 0.5, list = FALSE)
validation_data <- remaining_data[validation_indices, ]
test_data <- remaining_data[-validation_indices, ]

#Perform EN
# Define your predictor matrix (x) and response vector (y) for training
x_train <- train_data %>%
  select(matches("^(rs|PC)"))%>%
  scale()%>%
  as.matrix() 
y_train <- train_data$glucose

# Train the Elastic Net model using cross-validation on the training data
# cv.glmnet performs cross-validation to find optimal lambda (regularization parameter)
# Set up data
x_validation <- validation_data %>%
  select(matches("^(rs|PC)"))%>%
  scale() %>%
  as.matrix()
y_validation <- validation_data$glucose

#begin model tuning
alpha_values <- seq(0, 1, by = 0.05) # alpha = 1 for Lasso, alpha = 0 for Ridge. Values between 0 and 1 for Elastic Net.
cv_results <- list()
validation_errors <- numeric(length(alpha_values))

set.seed(42)
for (i in seq_along(alpha_values)) {
  # Train only on the training set
  fit <- cv.glmnet(x_train, y_train, alpha = alpha_values[i], nfolds = 10)
  cv_results[[i]] <- fit
  
  # Predict on the held-out validation set
  preds <- predict(fit, newx = x_validation, s = "lambda.min")
  validation_errors[i] <- mean((y_validation - preds)^2)
}

# Choose alpha with lowest validation MSE
best_alpha_index <- which.min(validation_errors)
best_alpha <- alpha_values[best_alpha_index]

# Combine train and validation data
full_train_data <- bind_rows(train_data, validation_data)

x_full_train <- full_train_data %>%
  select(matches("^(rs|PC)")) %>%
  scale() %>%
  as.matrix()
y_full_train <- full_train_data$glucose

# Re-train final model on full training data using best alpha
# Let glmnet automatically tune lambda across a sequence
final_model <- cv.glmnet(
  x_full_train,
  y_full_train,
  alpha = best_alpha,  # determined with tuning above
  family = "gaussian", 
  nfolds = 10           
)

#predict full_train model
predictions_full_train <- predict(final_model, newx = x_full_train, s = "lambda.min")

# Finally, evaluate the chosen model on the test set
x_test <- test_data %>%
  select(matches("^(rs|PC)")) %>%
  scale() %>%
  as.matrix()
y_test <- test_data$glucose

#predict test model
predictions_test <- predict(final_model, newx = x_test, s = "lambda.min")

# Print model performance measures
# Model tuning parameters
final_lambda <- final_model$lambda.min
print(paste("Best lambda:", final_lambda))

cat("Best alpha:", best_alpha, "\n") 

# On Full train (validation + train) set
SSE <- sum((y_full_train - predictions_full_train)^2)
SST <- sum((y_full_train - mean(y_full_train))^2)
R2_val <- 1 - SSE/SST

rmse_full_train <- sqrt(mean((y_full_train - predictions_full_train)^2))
mae_full_train <- mean(abs(y_full_train - predictions_full_train))
mse_full_train <- mean((y_full_train - predictions_full_train)^2)
accuracy_full_train <- sum(predictions_full_train == y_full_train)/length(y_full_train)

print(paste("Full Train MSE:", mse_full_train))
print(paste("Full Train RMSE:", round(rmse_full_train, 2)))
print(paste("Full Train MAE:", round(mae_full_train, 2)))
print(paste("Full Train Accuracy:", accuracy_full_train))
print(paste("Full Train R²:", round(R2_val, 3)))

# On test set
SSE_test <- sum((y_test - predictions_test)^2)
SST_test <- sum((y_test - mean(y_test))^2)
R2_test <- 1 - SSE_test/SST_test

rmse_test <- sqrt(mean((y_test - predictions_test)^2))
mae_test <- mean(abs(y_test - predictions_test))
mse_test <- mean((y_test - predictions_test)^2)
accuracy_test <- sum(predictions_test == y_test)/length(y_test)

print(paste("Test RMSE:", round(rmse_test, 2)))
print(paste("Test MAE:", round(mae_test, 2)))
print(paste("Test MSE:", mse_test))
print(paste("Test Accuracy:", accuracy_test))
print(paste("Test R²:", round(R2_test, 3)))

#### Part 2: Applying final model to full dataset ####
# Create final dataset
x_final <- data_glucose %>%
  select(matches("^(rs|PC)")) %>%
  scale() %>%
  as.matrix()
y_final <- data_glucose$glucose

# make predictions
predictions_final <- predict(final_model, newx = x_final, s = "lambda.min")

# Model outcomes on final full dataset
SSE_final <- sum((y_final - predictions_final)^2)
SST_final <- sum((y_final - mean(y_final))^2)
R2_final <- 1 - SSE_final/SST_final

rmse_final <- sqrt(mean((y_final - predictions_final)^2))
mae_final <- mean(abs(y_final - predictions_final))
mse_final <- mean((y_final - predictions_final)^2)
accuracy_final <- sum(predictions_final == y_final)/length(y_final)

print(paste("Final RMSE:", round(rmse_final, 2)))
print(paste("Final MAE:", round(mae_final, 2)))
print(paste("Final MSE:", mse_final))
print(paste("Full Train Accuracy:", accuracy_full_train))
print(paste("Test R²:", round(R2_final, 3)))

# Extract coefficients from the final model
coef_matrix <- as.matrix(coef(final_model))

# Convert to a tidy data frame
coef_df <- data.frame(
  SNP = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix)
)

# Filter to non-zero coefficients
nonzero_coef_df <- subset(coef_df, Coefficient != 0)

# Write all coefficients to file (including 0s)
write.table(coef_df, file = "EN_glucose.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Write only non-zero coefficients to file
write.table(nonzero_coef_df, file = "EN_filtered_glucose.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
