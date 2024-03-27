# Replace 'your_file_path.csv' with the path to your CSV file
SGspatial <- read.csv('C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/SGY22RSpatial.csv')

library(readxl)
library("PerformanceAnalytics")
library(car)
require(bestNormalize)
library(lme4)
library(sommer)
library(dplyr)
library(statgenGxE)
library(ggplot2)
library(metan)
library(emmeans)
library(caret)
install.packages('sommer')


devtools::install_github("dustinfife/flexplot", ref="development")
library(flexplot)
library(sjPlot)
library(glmmTMB)
library(nlme)
library(ASRtriala)
library(psych)
library(tidyr)

library(plyr)
library(MASS)
library(mice)
install.packages('C:/Users/BABlanchard/Downloads/ASRtriala_1.0.2.zip', repos = NULL, type = "win.binary")
library(ASRtriala)
library(asremlPlus)
#=============QC of the data===========
#This block will read in the master data, take out New Iberia year 2 data, calculate stalk weight, 
#replace missing stool values with survival, use individual data to estimate family data,
#remove certain families with missing data, make a csv dataset that is clean,
#and calculate phenotypic correlations in a chart


#set columns as categorical variables (factors)
SGspatial$Cross=as.factor(SGspatial$Cross)
SGspatial$Crop=as.factor(SGspatial$Crop)
SGspatial$Female=as.factor(SGspatial$Female)
SGspatial$Male=as.factor(SGspatial$Male)
SGspatial$Plot=as.factor(SGspatial$Plot)

SGspatial$Rep=as.factor(SGspatial$Rep)
SGspatial$Column=as.factor(SGspatial$Column)
SGspatial$Row=as.factor(SGspatial$Row)


#Set columns as numeric variables
SGspatial$PlotWeight=as.numeric(SGspatial$PlotWeight)

str(SGspatial)

# Load the 'car' package
table(SGspatial$Cross)
table(SGspatial$Rep)
table(SGspatial$Crop)
table(SGspatial$Row)
table(SGspatial$Column)

table(SGspatial$Row, SGspatial$Column)
table(SGspatial$Cross, SGspatial$Location, SGspatial$PlantYear)
table(SGspatial$Plot, SGspatial$Rep, SGspatial$Column, SGspatial$Row)

# Assuming 'NRspatial' is your dataframe
library(dplyr)

# Calculate the mean PlotWeight for each combination of Plot, Rep, Column, Row, and Crop


avgPlotWeight_aggregate <- fill.grid(
  data = SGspatial,
  row = "Row",
  col = "Column"
)

#set columns as categorical variables (factors)
avgPlotWeight_aggregate$Plot=as.factor(avgPlotWeight_aggregate$Plot)
avgPlotWeight_aggregate$Rep=as.factor(avgPlotWeight_aggregate$Rep)
avgPlotWeight_aggregate$Column=as.factor(avgPlotWeight_aggregate$Column)
avgPlotWeight_aggregate$Row=as.factor(avgPlotWeight_aggregate$Row)


#Set columns as numeric variables
avgPlotWeight_aggregate$PlotWeight=as.numeric(avgPlotWeight_aggregate$PlotWeight)

audit <- audit.single(data = avgPlotWeight_aggregate, gen = "Plot", rep = "Rep",
                      row = "Row", col = "Column",
                      resp = "PlotWeight", type.label = "none")
audit$trial.stats
audit$rep.stats
audit$trial.plot


anova_result <- aov(PlotWeight ~ Column, data = avgPlotWeight_aggregate)
# View the ANOVA results
summary(anova_result)

anova_result <- aov(PlotWeight ~ Row, data = avgPlotWeight_aggregate)
# View the ANOVA results
summary(anova_result)
anova_result <- aov(PlotWeight ~ Row*Column, data = avgPlotWeight_aggregate)
summary(anova_result)



# Create a data frame with unique values from TRS$Cross
unique_cross_values <- data.frame(Cross = unique(avgPlotWeight_aggregate$Cross))
# Create a numeric ID column called "gid"
unique_cross_values$gid <- seq_along(unique_cross_values$Cross)
# Merge the unique values and their IDs with the original data frame
PlotWeightgid <- merge(SGspatial, unique_cross_values, by = "Cross")
# The "gid" data frame will have a new column "gid" with numeric IDs for each value in TRS$Cross
print(PlotWeightgid)
PlotWeightgid$gid <- as.factor(PlotWeightgid$gid)
PlotWeightgid$Cross <- as.factor(PlotWeightgid$Cross)
PlotWeightgid$PlotWeight <- as.numeric(PlotWeightgid$PlotWeight)
PlotWeightgid$Row <- as.factor(PlotWeightgid$Row)
PlotWeightgid$Column <- as.factor(PlotWeightgid$Column)

PlotWeightgid <- PlotWeightgid[complete.cases(PlotWeightgid$Row, PlotWeightgid$Column), ]
PlotWeightgid$Inter <- interaction(PlotWeightgid$gid, 
                                   PlotWeightgid$Crop)
PlotWeightgid$SPAT <- interaction(PlotWeightgid$Row, 
                                  PlotWeightgid$Column)
PlotWeightgid$SPAT <- as.factor(PlotWeightgid$SPAT)
PlotWeightgid$Inter <- as.factor(PlotWeightgid$Inter)

# Running the model
solspat <-  mmer(fixed = PlotWeight ~ 1,
                 random = ~ gid + vsr(usr(Row)) + vsr(usr(Column)),
                 rcov = ~ vsr(dsr(Column), units),
                 data = PlotWeightgid)
summary(solspat)
BLUPs_Cross <- as.data.frame(solspat$U$'gid'$PlotWeight + mean(PlotWeightgid$PlotWeight))


sol <-  mmer(fixed = PlotWeight ~ 1,
             random = ~ gid + Rep,
             rcov = ~ vsr(units),
             data = PlotWeightgid)
summary(sol)
BLUPs_Cross2 <- as.data.frame(sol$U$'gid'$PlotWeight + mean(PlotWeightgid$PlotWeight))



cor(BLUPs_Cross[,1], BLUPs_Cross2[,1])
anova(sol, solspat)

# Rename the only column in BLUPs_Cross to "SPATBLUPs"
colnames(BLUPs_Cross) <- "SPATBLUPs"

# Rename the only column in BLUPs_Cross2 to "REPBLUPs"
colnames(BLUPs_Cross2) <- "REPBLUPs"

# Combine BLUPs_Cross with BLUPs_Cross2 by matching their rownames
combined_BLUPs <- merge(BLUPs_Cross, BLUPs_Cross2, by = "row.names", all = TRUE)

# Reshape the data into long format
combined_BLUPs_long <- combined_BLUPs %>%
  pivot_longer(cols = -Row.names, names_to = "Variable", values_to = "Value")

# Define a variable to represent varieties
combined_BLUPs_long$Variety <- as.factor(combined_BLUPs_long$Row.names)

# Plot parallel coordinate plot
ggplot(combined_BLUPs_long, aes(x = Variable, y = Value, group = Row.names, color = Variety)) +
  geom_line(size = 1.0) +  # Thicken the lines
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 8),  # Adjust y-axis text size
        legend.position = "none",  # Remove legend
        plot.title = element_text(face = "bold"),  # Bold plot title
        axis.title.x = element_text(face = "bold"),  # Bold x-axis title
        axis.title.y = element_text(face = "bold")) +  # Bold y-axis title
  labs(x = "Model", y = "BLUPs", color = NULL,  # Remove legend title
       title = "BLUPs per Family for PlotWeight") +  # Add plot title
  ggtitle("Family Ranks for Plot Weight") + 
  scale_y_continuous(breaks = seq(40, max(combined_BLUPs_long$Value), by = 20)) +
  scale_x_discrete(labels = c("SPATBLUPs" = "Model[5]", "REPBLUPs" = "Model[4]"))









solspat <- lme(fixed = PlotWeight ~ 1,  
               random = list(gid = pdDiag(~1), SPAT = pdDiag(~1)),
               data = PlotWeightgid,
               method = 'REML')
summary(solspat)
gid_BLUPs <- ranef(solspat)$gid + mean(PlotWeightgid$PlotWeight)


# Create the linear mixed model using lme function and specify heterogenous variance 

sol <- lme(fixed = PlotWeight ~ 1,  
           random = list(gid = pdDiag(~1)),
           data = PlotWeightgid,
           method = 'REML')
summary(sol)
gid_BLUPs <- ranef(sol)$gid + mean(PlotWeightgid$PlotWeight)

sol <- lmer(PlotWeight~1 + (1|gid), PlotWeightgid)
summary(sol)
coef(sol)$gid
solSPAT <- lmer(PlotWeight~1 + (1|gid) + (1|SPAT), PlotWeightgid)
summary(solSPAT)
coef(solSPAT)$gid




# Set the number of folds for cross-validation (e.g., 5-fold)
num_folds <- 5

# Create a data partition for cross-validation
set.seed(123)  # Set a random seed for reproducibility
cv_index <- createFolds(PlotWeightgid$gid, k = num_folds, list = TRUE)

# Initialize lists to store cross-validation results
cv_predictions <- vector("list", length = num_folds)
cv_actual <- vector("list", length = num_folds)

# Perform cross-validation
for (fold in 1:num_folds) {
  # Extract the indices for the training and validation sets
  train_indices <- unlist(cv_index[-fold])
  valid_indices <- cv_index[[fold]]
  
  # Split the data into training and validation sets
  train_data <- PlotWeightgid[train_indices, ]
  valid_data <- PlotWeightgid[valid_indices, ]
  
  # Fit the model on the training data
  fitMET <- mmer(fixed = PlotWeight ~ 1,
                 random = ~ gid + Rep,
                 rcov = ~ vsr(units),
                 data = PlotWeightgid)
  
  # Create a data frame with meaningful column names
  BLUPs_Cross <- as.data.frame(fitMET$U$'gid'$PlotWeight + mean(PlotWeightgid$PlotWeight))
  rownames(BLUPs_Cross) <- sub("^gid", "", rownames(BLUPs_Cross))
  #  rownames(BLUPs_Cross) <- sub(":Crop$", "", rownames(BLUPs_Cross))
  colnames(BLUPs_Cross) <- "prediction"  # Assign a column name to BLUPs_Cross
  
  # Initialize the prediction column in valid_data
  valid_data$prediction <- NA
  
  # Loop through each row in valid_data
  for (i in 1:nrow(valid_data)) {
    gid <- valid_data$gid[i]
    row_index <- which(rownames(BLUPs_Cross) == gid)
    
    if (length(row_index) > 0) {
      valid_data$prediction[i] <- BLUPs_Cross$prediction[row_index]
    }
  }
  
  # Store the predicted and actual values for this fold
  cv_predictions[[fold]] <- as.list(valid_data$prediction)
  cv_actual[[fold]] <- as.list(valid_data$PlotWeight)
}


# Calculate RMSE for each fold
rmse_results <- sapply(1:num_folds, function(fold) {
  sqrt(mean((unlist(cv_predictions[[fold]]) - unlist(cv_actual[[fold]]))^2))
})


# Print RMSE for each fold
print(rmse_results)

# Initialize a vector to store correlations for each fold
correlations <- numeric(length = num_folds)

# Loop through the folds and calculate correlations
for (fold in 1:num_folds) {
  x <- as.numeric(cv_predictions[[fold]])
  y <- as.numeric(cv_actual[[fold]])
  
  if (!anyNA(x) && !anyNA(y)) {
    correlation <- cor(x, y)
    correlations[fold] <- correlation
  }
  else {
    correlations[fold] <- NA  # Set to NA if there are missing values
  }
}

# View the correlations
correlations

# Assess prediction accuracy for each fold
accuracy_results <- sapply(1:num_folds, function(fold) {
  # Calculate RMSE
  rmse <- sqrt(mean((unlist(cv_predictions[[fold]]) - unlist(cv_actual[[fold]]))^2))
  
  # Calculate MAE
  mae <- mean(abs(unlist(cv_predictions[[fold]]) - unlist(cv_actual[[fold]])))
  
  # Calculate R-squared
  actual <- unlist(cv_actual[[fold]])
  predicted <- unlist(cv_predictions[[fold]])
  ss_residual <- sum((actual - predicted)^2)
  ss_total <- sum((actual - mean(actual))^2)
  r_squared <- 1 - (ss_residual / ss_total)
  
  # Combine the metrics into a list
  result <- list(
    RMSE = rmse,
    MAE = mae,
    R_squared = r_squared
  )
  
  return(result)
})

# View the accuracy results for each fold
accuracy_results <- as.data.frame(accuracy_results)



# Set the number of folds for cross-validation (e.g., 5-fold)
num_folds <- 5

# Create a data partition for cross-validation
set.seed(123)  # Set a random seed for reproducibility
cv_index <- createFolds(PlotWeightgid$gid, k = num_folds, list = TRUE)

# Initialize lists to store cross-validation results
cv_predictions <- vector("list", length = num_folds)
cv_actual <- vector("list", length = num_folds)

# Perform cross-validation
for (fold in 1:num_folds) {
  # Extract the indices for the training and validation sets
  train_indices <- unlist(cv_index[-fold])
  valid_indices <- cv_index[[fold]]
  
  # Split the data into training and validation sets
  train_data <- PlotWeightgid[train_indices, ]
  valid_data <- PlotWeightgid[valid_indices, ]
  
  # Fit the model on the training data
  fitMET2 <- mmer(fixed = PlotWeight ~ 1,
                  random = ~ gid + Row + Column,
                  rcov = ~ vsr(units),
                  data = PlotWeightgid)
  
  
  # Create a data frame with meaningful column names
  BLUPs_Cross <- as.data.frame(fitMET2$U$'gid'$PlotWeight + mean(PlotWeightgid$PlotWeight))
  rownames(BLUPs_Cross) <- sub("^gid", "", rownames(BLUPs_Cross))
  #  rownames(BLUPs_Cross) <- sub(":Crop$", "", rownames(BLUPs_Cross))
  colnames(BLUPs_Cross) <- "prediction"  # Assign a column name to BLUPs_Cross
  
  # Initialize the prediction column in valid_data
  valid_data$prediction <- NA
  
  # Loop through each row in valid_data
  for (i in 1:nrow(valid_data)) {
    gid <- valid_data$gid[i]
    row_index <- which(rownames(BLUPs_Cross) == gid)
    
    if (length(row_index) > 0) {
      valid_data$prediction[i] <- BLUPs_Cross$prediction[row_index]
    }
  }
  
  # Store the predicted and actual values for this fold
  cv_predictions[[fold]] <- as.list(valid_data$prediction)
  cv_actual[[fold]] <- as.list(valid_data$PlotWeight)
}


# Initialize a vector to store correlations for each fold
correlations2 <- numeric(length = num_folds)

# Loop through the folds and calculate correlations
for (fold in 1:num_folds) {
  x <- as.numeric(cv_predictions[[fold]])
  y <- as.numeric(cv_actual[[fold]])
  
  if (!anyNA(x) && !anyNA(y)) {
    correlation <- cor(x, y)
    correlations2[fold] <- correlation
  }
  else {
    correlations2[fold] <- NA  # Set to NA if there are missing values
  }
}

# View the correlations
correlations2


# Calculate RMSE for each fold
rmse_results2 <- sapply(1:num_folds, function(fold) {
  sqrt(mean((unlist(cv_predictions[[fold]]) - unlist(cv_actual[[fold]]))^2))
})


# Print RMSE for each fold
print(rmse_results2)


# Assess prediction accuracy for each fold
accuracy_results2 <- sapply(1:num_folds, function(fold) {
  # Calculate RMSE
  rmse <- sqrt(mean((unlist(cv_predictions[[fold]]) - unlist(cv_actual[[fold]]))^2))
  
  # Calculate MAE
  mae <- mean(abs(unlist(cv_predictions[[fold]]) - unlist(cv_actual[[fold]])))
  
  # Calculate R-squared
  actual <- unlist(cv_actual[[fold]])
  predicted <- unlist(cv_predictions[[fold]])
  ss_residual <- sum((actual - predicted)^2)
  ss_total <- sum((actual - mean(actual))^2)
  r_squared <- 1 - (ss_residual / ss_total)
  
  # Combine the metrics into a list
  result <- list(
    RMSE = rmse,
    MAE = mae,
    R_squared = r_squared
  )
  
  return(result)
})

# View the accuracy results for each fold
accuracy_results2 <- as.data.frame(accuracy_results2)
# Assuming accuracy_results2 is your dataframe
flipped_accuracy_results2 <- as.data.frame(t(accuracy_results2))
# Add a column called "results" with the specified value
flipped_accuracy_results2$model <- "with SPATIAL"


# View the accuracy results for each fold
accuracy_results <- as.data.frame(accuracy_results)
# Assuming accuracy_results2 is your dataframe
flipped_accuracy_results <- as.data.frame(t(accuracy_results))
# Add a column called "results" with the specified value
flipped_accuracy_results$model <- "w/o SPATIAL"

# Assuming flipped_accuracy_results and flipped_accuracy_results2 are dataframes
combined_results <- rbind(flipped_accuracy_results, flipped_accuracy_results2)

# Create a boxplot comparing RMSE for both levels of the model
boxplot((as.numeric(RMSE)) ~ model, data = combined_results, 
        ylab = "RMSE", main = "Comparison of RMSE by Model Level")

boxplot((as.numeric(R_squared)) ~ model, data = combined_results, 
        ylab = "R_squared", main = "Comparison of R_squared by Model Level")



correlations <- unlist(correlations)  # Convert to a numeric vector if not already
correlations2 <- unlist(correlations2)  # Convert to a numeric vector if not already

# Create a boxplot to compare correlations
boxplot(correlations, correlations2, 
        names = c("w/o Spatial", "with Spatial"), 
        ylab = "Correlation", main = "Comparison of Correlations")

anova(fitMET, fitMET2)