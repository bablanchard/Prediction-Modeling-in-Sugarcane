library(readxl)
library("PerformanceAnalytics")
library(car)
library(dplyr)
library(nlme)
library(ASRtriala)
library(sommer)
library(caret)



#=============QC of the data===========
#This block will read in the master data, take out New Iberia year 2 data, calculate stalk weight, 
#replace missing stool values with survival, use individual data to estimate family data,
#remove certain families with missing data, make a csv dataset that is clean,
#and calculate phenotypic correlations in a chart

#read master
df <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/Master.csv", header = TRUE)
#take out NI Y2
df_filtered <- df[!(df$PlantYear == 2021 & df$Location == "NI"), ]
#make Bustalks numeric with the number of stalks in the bundle
df_filtered$Bustalks <- as.numeric(gsub("^\\*([0-9]+)$", "\\1", df_filtered$Bustalks))
df_filtered$Bustalks[is.na(df_filtered$Bustalks) & is.na(df_filtered$Bustalks)] <- 10
#calcualte SW
df_filtered$SW <- df_filtered$Bu.Wt / df_filtered$Bustalks
#if Stool = missing then fill in with Survival
df_filtered$Stool[is.na(df_filtered$Stool) | df_filtered$Stool == ""] <- df_filtered$Survival[is.na(df_filtered$Stool) | df_filtered$Stool == ""]
#read individual data from New Roads
individ <- read_excel("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/SeedlingIndividual.xlsx", sheet = "NRdatasheet")
#remove observations in individ that contain NAs for multiple columns
individ <- individ[complete.cases(individ[, c("Stalks", "Height", "Dia1", "Dia2", "Brix")]), ]
#add location column and PlantYear according to crop
individ$Location <- "NR"
individ$PlantYear <- ifelse(individ$Crop == 2, 2020, ifelse(individ$Crop == 1, 2021, NA))
#average Dia1 and Dia2
individ$Diam <- rowMeans(individ[, c("Dia1", "Dia2")], na.rm = TRUE)
#calculate average family diameter from individuals and apply to each observation of family to input into master
famdiam <- aggregate(Diam ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famdiam, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Diam[is.na(individ$Diam)] <- individ$Diam_avg[is.na(individ$Diam)]
names(individ)[names(individ) == "Diam_avg"] <- "famdiam"
#calculate average family Brix from individuals and apply to each observation of family to input into master
famBrix <- aggregate(Brix ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famBrix, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Brix[is.na(individ$Brix)] <- individ$famBrix[is.na(individ$Brix)]
names(individ)[names(individ) == "Brix_avg"] <- "famBrix"
#calculate average family diameter from individuals and apply to each observation of family to input into master
famHt <- aggregate(Height ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famHt, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Height[is.na(individ$Height)] <- individ$famHt[is.na(individ$Height)]
names(individ)[names(individ) == "Height_avg"] <- "famHt"
#calculate average family stalks from individuals and apply to each observation of family to input into master
famStlk <- aggregate(Stalks ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famStlk, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Stalks[is.na(individ$Stalks)] <- individ$famStlk[is.na(individ$Stalks)]
names(individ)[names(individ) == "Stalks_avg"] <- "famStlk"
#remove observation for SG PlantYear2020 Plot 49
df_filtered <- df_filtered[!(df_filtered$PlantYear == 2020 & df_filtered$Location == "SG" & df_filtered$Plot == 49), ]
#remove observation for SG PlantYear2021 Plot 28
df_filtered <- df_filtered[!(df_filtered$PlantYear == 2021 & df_filtered$Location == "SG" & df_filtered$Plot == 28), ]
# Match combinations of PlantYear, Crop, Plot, Rep, and Location in famStlk_avg to df_filtered and add the famStlk column
df_filtered$famStlk <- individ$famStlk[match(paste(df_filtered$PlantYear, df_filtered$Crop, df_filtered$Plot, df_filtered$Rep, df_filtered$Location), 
                                             paste(individ$PlantYear, individ$Crop, individ$Plot, individ$Rep, individ$Location))]
#now multiple the famStalk column to the number of stools to get an estimated stalk count for those families where only individual data was taken (these will be inflated)
df_filtered$Stalk <- ifelse(is.na(df_filtered$famStlk), df_filtered$Stalk, df_filtered$famStlk * df_filtered$Stool)
# Write the dataframe to a CSV file
write.csv(df_filtered, file = "Master_1.2.csv", row.names = TRUE)

# phenotypic correlation between traits
sapply(df_filtered[, 11:22, 27], class)
df_filtered[, 11:22] <- lapply(df_filtered[, 11:22], as.numeric)
df_filtered[, 27] <- as.numeric(df_filtered[, 27])
round(cor(df_filtered[, c(11:22, 27)], use = "pairwise"), 2)
data_subset <- df_filtered[, c(11:18, 27)]
#install.packages("PerformanceAnalytics")
chart.Correlation(as.matrix(na.omit(data_subset)), histogram = TRUE, pch = 1)

# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS <- df_filtered[complete.cases(df_filtered$TRS) & df_filtered$TRS != 0,]
Stalk <- df_filtered[complete.cases(df_filtered$Stalk) & df_filtered$Stalk != 0,]
Height <- df_filtered[complete.cases(df_filtered$Height) & df_filtered$Height != 0,]
Dia <- df_filtered[complete.cases(df_filtered$Dia) & df_filtered$Dia != 0,]
PlotWeight <- df_filtered[complete.cases(df_filtered$PlotWeight) & df_filtered$PlotWeight != 0,]
Brix <- df_filtered[complete.cases(df_filtered$Brix) & df_filtered$Brix != 0,]
Fiber <- df_filtered[complete.cases(df_filtered$Fiber) & df_filtered$Fiber != 0,]
SW <- df_filtered[complete.cases(df_filtered$SW) & df_filtered$SW != 0,]
Sucrose <- df_filtered[complete.cases(df_filtered$Sucrose..) & df_filtered$Sucrose.. != 0,]

Dia <- Dia[Dia$PlantYear == 2020 & 
             Dia$Crop %in% c(1) & 
             Dia$Location %in% c("NI", "NR", "SG"), ]

TRS$TRS <- as.numeric(TRS$TRS)
TRS <- TRS[complete.cases(TRS$TRS) & TRS$TRS != 0,]


#now make a genotype id dataframe
# Create a data frame with unique values from TRS$Cross
unique_cross_values <- data.frame(Cross = unique(TRS$Cross))
# Create a numeric ID column called "gid"
unique_cross_values$gid <- seq_along(unique_cross_values$Cross)
# Merge the unique values and their IDs with the original data frame
TRSgid <- merge(TRS, unique_cross_values, by = "Cross")
# The "gid" data frame will have a new column "gid" with numeric IDs for each value in TRS$Cross
print(Diagid)
unique_values <- unique(TRSgid$Cross)
print(unique_values)

Acsv <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/FxEsubsetAmatrix.csv")
unique_values <- unique(Acsv$X)
print(unique_values)


# Assuming Acsv$X needs to be replaced with PlotWeightgid$gid based on matching values in PlotWeightgid$Cross

# Match the indices of the Cross values in Acsv$X with PlotWeightgid$Cross and replace X with corresponding gid values
matched_indices <- match(Acsv$X, TRSgid$Cross)

# Replace Acsv$X with corresponding values from PlotWeightgid$gid
Acsv$X <- TRSgid$gid[matched_indices]

unique_col_names <- unique(colnames(Acsv))
print(unique_col_names)
new_col_names <- gsub("\\.", "-", unique_col_names)
colnames(Acsv) <- new_col_names
# Get current column names
current_col_names <- colnames(Acsv)
print(current_col_names)
current_col_names <- colnames(Acsv)


# Find the index of the column to replace
index_to_replace <- which(current_col_names == "CP18-1901-POL")
# Check if the column is found before attempting replacement
if (length(index_to_replace) > 0) {
  # Replace the column name
  current_col_names[index_to_replace] <- "CP18-1901 POLY"
  
  # Update column names in the dataframe
  colnames(Acsv) <- current_col_names
  
  # Print or inspect the modified dataframe
  print(colnames(Acsv))
} else {
  print("Column 'CP18-1901-POL' not found.")
}
# Find the index of the column to replace
index_to_replace <- which(current_col_names == "CP18-1921-POL")
# Check if the column is found before attempting replacement
if (length(index_to_replace) > 0) {
  # Replace the column name
  current_col_names[index_to_replace] <- "CP18-1921 POLY"
  
  # Update column names in the dataframe
  colnames(Acsv) <- current_col_names
  
  # Print or inspect the modified dataframe
  print(colnames(Acsv))
} else {
  print("Column 'CP18-1921-POL' not found.")
}
# Find the index of the column to replace
index_to_replace <- which(current_col_names == "HB-17-3183")
# Check if the column is found before attempting replacement
if (length(index_to_replace) > 0) {
  # Replace the column name
  current_col_names[index_to_replace] <- "HB 17-3183"
  
  # Update column names in the dataframe
  colnames(Acsv) <- current_col_names
  
  # Print or inspect the modified dataframe
  print(colnames(Acsv))
} else {
  print("Column 'HB-17-3183' not found.")
}
# Find the index of the column to replace
index_to_replace <- which(current_col_names == "HB-17-3208")
# Check if the column is found before attempting replacement
if (length(index_to_replace) > 0) {
  # Replace the column name
  current_col_names[index_to_replace] <- "HB 17-3208"
  
  # Update column names in the dataframe
  colnames(Acsv) <- current_col_names
  
  # Print or inspect the modified dataframe
  print(colnames(Acsv))
} else {
  print("Column 'HB-17-3208' not found.")
}
# Find the index of the column to replace
index_to_replace <- which(current_col_names == "HB-19-3031")
# Check if the column is found before attempting replacement
if (length(index_to_replace) > 0) {
  # Replace the column name
  current_col_names[index_to_replace] <- "HB 19-3031"
  
  # Update column names in the dataframe
  colnames(Acsv) <- current_col_names
  
  # Print or inspect the modified dataframe
  print(colnames(Acsv))
} else {
  print("Column 'HB-19-3031' not found.")
}







# Loop through column names in Acsv
for (col in colnames(Acsv)) {
  # Match the column name with values in PlotWeightgid$Cross
  match_index <- match(col, TRSgid$Cross)
  
  # If a match is found, replace the column name with the corresponding PlotWeightgid$gid value
  if (!is.na(match_index)) {
    new_col_name <- TRSgid$gid[match_index]
    colnames(Acsv)[colnames(Acsv) == col] <- new_col_name
  }
}

# Replace row names with values in Acsv$X
rownames(Acsv) <- Acsv$X

# Remove the 'Acsv$X' column
Acsv <- Acsv[, -which(names(Acsv) == "X")]

# Convert Acsv to a matrix
Acsv_matrix <- as.matrix(Acsv)

# Apply the conversion function to all elements of Acsv_matrix
dim(Acsv_matrix)
length(unique(TRSgid$gid))

TRSgid$gid <- as.factor(TRSgid$gid)
TRSgid$Dia <- as.numeric(TRSgid$TRS)


PWfit <- mmer(fixed = TRS ~ 1,
              random = ~ gid + Rep + gid:Rep,
              rcov = ~ vsr(units),
              data = TRSgid)
summary(PWfit)



PWfit2 <- mmer(fixed = TRS ~ 1,
               random = ~ vsr(gid, Gu = Acsv_matrix) + Rep + gid:Rep,
               rcov = ~ vsr(units),
               data = TRSgid)
summary(PWfit2)


anova(PWfit, PWfit2)

# Set the number of folds for cross-validation (e.g., 5-fold)
num_folds <- 10

# Create a data partition for cross-validation
set.seed(123)  # Set a random seed for reproducibility
cv_index <- createFolds(TRSgid$gid, k = num_folds, list = TRUE)

# Initialize lists to store cross-validation results
cv_predictions <- vector("list", length = num_folds)
cv_actual <- vector("list", length = num_folds)

# Perform cross-validation
for (fold in 1:num_folds) {
  # Extract the indices for the training and validation sets
  train_indices <- unlist(cv_index[-fold])
  valid_indices <- cv_index[[fold]]
  
  # Split the data into training and validation sets
  train_data <- TRSgid[train_indices, ]
  valid_data <- TRSgid[valid_indices, ]
  
  # Fit the model on the training data
  fitMET <- mmer(fixed = TRS ~ 1,
                 random = ~ gid + Rep + gid:Rep,
                 rcov = ~ vsr(units),
                 data = train_data)
  
  # Create a data frame with meaningful column names
  BLUPs_Cross <- as.data.frame(fitMET$U$gid$TRS + mean(train_data$TRS))
  rownames(BLUPs_Cross) <- sub("^gid", "", rownames(BLUPs_Cross))
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
  cv_actual[[fold]] <- as.list(valid_data$TRS)
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
num_folds <- 10

# Create a data partition for cross-validation
set.seed(123)  # Set a random seed for reproducibility
cv_index <- createFolds(TRSgid$gid, k = num_folds, list = TRUE)

# Initialize lists to store cross-validation results
cv_predictions <- vector("list", length = num_folds)
cv_actual <- vector("list", length = num_folds)

# Perform cross-validation
for (fold in 1:num_folds) {
  # Extract the indices for the training and validation sets
  train_indices <- unlist(cv_index[-fold])
  valid_indices <- cv_index[[fold]]
  
  # Split the data into training and validation sets
  train_data <- TRSgid[train_indices, ]
  valid_data <- TRSgid[valid_indices, ]
  
  # Fit the model on the training data
  fitMET2 <- mmer(fixed = TRS ~ 1,
                  random = ~ vsr(gid, Gu = Acsv_matrix) + Rep + gid:Rep,
                  rcov = ~ vsr(units),
                  data = train_data) 
  
  
  # Create a data frame with meaningful column names
  BLUPs_Cross <- as.data.frame(fitMET2$U$"u:gid"$TRS + mean(train_data$TRS))
  rownames(BLUPs_Cross) <- sub("^gid", "", rownames(BLUPs_Cross))
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
  cv_actual[[fold]] <- as.list(valid_data$TRS)
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
flipped_accuracy_results2$model <- "with GRM"


# View the accuracy results for each fold
accuracy_results <- as.data.frame(accuracy_results)
# Assuming accuracy_results2 is your dataframe
flipped_accuracy_results <- as.data.frame(t(accuracy_results))
# Add a column called "results" with the specified value
flipped_accuracy_results$model <- "w/o GRM"

# Assuming flipped_accuracy_results and flipped_accuracy_results2 are dataframes
combined_results <- rbind(flipped_accuracy_results, flipped_accuracy_results2)

# Create a boxplot comparing RMSE for both levels of the model
boxplot((as.numeric(RMSE)) ~ model, data = combined_results, 
        ylab = "RMSE", main = "Comparison of RMSE by Model Level")
boxplot((as.numeric(R_squared)) ~ model, data = combined_results, 
        ylab = "R_squared", main = "Comparison of RMSE by Model Level")


correlations <- unlist(correlations)  # Convert to a numeric vector if not already
correlations2 <- unlist(correlations2)  # Convert to a numeric vector if not already

# Create a boxplot to compare correlations
boxplot(correlations, correlations2, 
        names = c("w/o Pedigree based GRM", "with Pedigree based GRM"), 
        ylab = "Correlation", main = "Comparison of Correlations")

# Define colors for the boxplot
boxplot_colors <- c("salmon", "seagreen2")
boxplot_outlines <- c("salmon4", "seagreen4")

par(mar = c(6, 6, 4, 4))  # Adjust the margin (bottom, left, top, right)

# Creating the boxplot with modified y-axis ticks
bp <- boxplot(correlations, correlations2, 
              names = c("w/o Pedigree based GRM", "with Pedigree based GRM"), 
              ylab = "Correlation of Predicted vs Actual Values", 
              col = boxplot_colors, border = boxplot_outlines, col.main = "black", font.main = 12, 
              cex.axis = 1.2, cex.lab = 1.6, las = 1)
title(main = "Prediction Accuracy of Mixed Models", cex.main = 1.8)  # Adjust cex.main for larger title text