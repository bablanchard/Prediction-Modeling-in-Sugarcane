library(readxl)
library(lme4)
library(sommer)
library(dplyr)
library(statgenGxE)
library(ggplot2)
library(metan)
library(broom)
library(tidyr)
df <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/LSU Breeding Program/2024ReleaseMeeting/active23metric_1.1.csv", header = TRUE)

#=========Outfield Subsetting=========
# Filter the dataframe
outfield <- df %>%
  filter((YEAR == 2021 & CROP == 0 & A_SERIES == 9999) |
           (YEAR == 2022 & CROP %in% c(0, 1) & A_SERIES == 9999) |
           (YEAR == 2023 & CROP %in% c(0, 1, 2) & A_SERIES == 9999))

# Create a new column $Env and set default as empty string
outfield$Env <- ""

# Assign value 'Sugar20PC' to $Env where conditions are met
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "ALLAINS" & outfield$CROP == 0] <- "AL210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "ALLAINS" & outfield$CROP == 0] <- "AL220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "ALLAINS" & outfield$CROP == 1] <- "AL221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "ALLAINS" & outfield$CROP == 0] <- "AL230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "ALLAINS" & outfield$CROP == 1] <- "AL231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "ALLAINS" & outfield$CROP == 2] <- "AL232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "ALMA" & outfield$CROP == 0] <- "ALM210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "ALMA" & outfield$CROP == 0] <- "ALM220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "ALMA" & outfield$CROP == 1] <- "ALM221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "ALMA" & outfield$CROP == 0] <- "ALM230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "ALMA" & outfield$CROP == 1] <- "ALM231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "ALMA" & outfield$CROP == 2] <- "ALM232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "BRUNSWICK" & outfield$CROP == 0] <- "BR210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "BRUNSWICK" & outfield$CROP == 0] <- "BR220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "BRUNSWICK" & outfield$CROP == 1] <- "BR221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "BRUNSWICK" & outfield$CROP == 0] <- "BR230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "BRUNSWICK" & outfield$CROP == 1] <- "BR231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "BRUNSWICK" & outfield$CROP == 2] <- "BR232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "GLENWOOD" & outfield$CROP == 0] <- "GL210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "GLENWOOD" & outfield$CROP == 0] <- "GL220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "GLENWOOD" & outfield$CROP == 1] <- "GL221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "GLENWOOD" & outfield$CROP == 0] <- "GL230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "GLENWOOD" & outfield$CROP == 1] <- "GL231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "GLENWOOD" & outfield$CROP == 2] <- "GL232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "HARPER" & outfield$CROP == 0] <- "HA210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "HARPER" & outfield$CROP == 0] <- "HA220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "HARPER" & outfield$CROP == 1] <- "HA221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "HARPER" & outfield$CROP == 0] <- "HA230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "HARPER" & outfield$CROP == 1] <- "HA231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "HARPER" & outfield$CROP == 2] <- "HA232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "LANAUX" & outfield$CROP == 0] <- "LANA210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "LANAUX" & outfield$CROP == 0] <- "LANA220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "LANAUX" & outfield$CROP == 1] <- "LANA221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "LANAUX" & outfield$CROP == 0] <- "LANA230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "LANAUX" & outfield$CROP == 1] <- "LANA231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "LANAUX" & outfield$CROP == 2] <- "LANA232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "LANDRY" & outfield$CROP == 0] <- "LAND210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "LANDRY" & outfield$CROP == 0] <- "LAND220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "LANDRY" & outfield$CROP == 1] <- "LAND221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "LANDRY" & outfield$CROP == 0] <- "LAND230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "LANDRY" & outfield$CROP == 1] <- "LAND231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "LANDRY" & outfield$CROP == 2] <- "LAND232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "MARY" & outfield$CROP == 0] <- "MA210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "MARY" & outfield$CROP == 0] <- "MA220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "MARY" & outfield$CROP == 1] <- "MA221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "MARY" & outfield$CROP == 0] <- "MA230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "MARY" & outfield$CROP == 1] <- "MA231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "MARY" & outfield$CROP == 2] <- "MA232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "RHEBERT" & outfield$CROP == 0] <- "RH210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "RHEBERT" & outfield$CROP == 0] <- "RH220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "RHEBERT" & outfield$CROP == 1] <- "RH221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "RHEBERT" & outfield$CROP == 0] <- "RH230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "RHEBERT" & outfield$CROP == 1] <- "RH231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "RHEBERT" & outfield$CROP == 2] <- "RH232"
outfield$Env[outfield$YEAR == 2021 & outfield$LOC == "STJOHN" & outfield$CROP == 0] <- "ST210"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "STJOHN" & outfield$CROP == 0] <- "ST220"
outfield$Env[outfield$YEAR == 2022 & outfield$LOC == "STJOHN" & outfield$CROP == 1] <- "ST221"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "STJOHN" & outfield$CROP == 0] <- "ST230"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "STJOHN" & outfield$CROP == 1] <- "ST231"
outfield$Env[outfield$YEAR == 2023 & outfield$LOC == "STJOHN" & outfield$CROP == 2] <- "ST232"
env_values <- unique(outfield$Env)

#============one example================
# Subset outfield by observations where $Env = ST232
ST232 <- outfield[outfield$Env == "ST232", ]
ST232$VARIETY <- as.factor(ST232$VARIETY)
ST232$REP <- as.factor(ST232$REP)


ST232fit <- mmer(fixed = T_SPACRE ~ 1 ,
               random = ~ VARIETY + REP,
               rcov = ~ vsr(units),
               data = ST232)

summary(ST232fit)


# 'VARIETY' is your fixed effect (variety), and 'replication' is your random effect
model <- lmer(T_SPACRE ~ 1 + (1|VARIETY) + (1|REP), data = ST232)

# Step 4: Model Diagnostics
summary(model)  # Check model summary
plot(model)     # Diagnostic plots

# Obtain variance components
var_components <- VarCorr(model)
varcomp <- as.data.frame(var_components)
# Obtain BLUPs (Best Linear Unbiased Predictors)
blups <- ranef(model)

G <- varcomp[varcomp$grp == "VARIETY", "vcov"]
P <- varcomp[varcomp$grp == "REP", "vcov"]
R <- varcomp[varcomp$grp == "Residual", "vcov"]
H <- G/(G+R/3)

mean_T_SPACRE <- mean(ST232$T_SPACRE)
sd_T_SPACRE <- sd(ST232$T_SPACRE)

CV <- (sd_T_SPACRE/mean_T_SPACRE * 100)


#==========loop and collect trial quality=============
sugaryield <- subset(outfield, !is.na(T_SPACRE) & T_SPACRE != 0)
sugaryield$VARIETY <- factor(sugaryield$VARIETY)
# Get unique values of outfield$Env
env_values <- unique(sugaryield$Env)
print(env_values)
# Initialize an empty list to store subsets
datasets_list <- list()

# Iterate over each unique value of outfield$Env
for (env_value in env_values) {
  # Subset dataframe for the current value of outfield$Env
  subset_df <- sugaryield[sugaryield$Env == env_value, ]
  
  # Assign the subset dataframe to a list element with the name of the value of outfield$Env
  datasets_list[[env_value]] <- subset_df
}


# Determine the length of env_values
n <- length(env_values)

# Create dataframe with Env values and empty columns for H and CV
outfieldtrials <- data.frame(
  Env = env_values,
  H = rep(NA, n),
  CV = rep(NA, n),
  stringsAsFactors = FALSE
)

# Loop through each dataset in datasets_list
for (env_value in names(datasets_list)) {
  # Check if there are valid cases for the response variable
  if (!any(is.na(datasets_list[[env_value]]$T_SPACRE))) {
    # Fit the linear mixed model
    model <- lmer(T_SPACRE ~ 1 + (1|VARIETY) + (1|REP), data = datasets_list[[env_value]])
    
    # Calculate CV using response variable
    mean_T_SPACRE <- mean(datasets_list[[env_value]]$T_SPACRE)
    sd_T_SPACRE <- sd(datasets_list[[env_value]]$T_SPACRE)
    CV <- (sd_T_SPACRE / mean_T_SPACRE) * 100
    
    # Obtain variance components
    var_components <- VarCorr(model)
    varcomp <- as.data.frame(var_components)
    
    # Calculate H
    G <- varcomp[varcomp$grp == "VARIETY", "vcov"]
    R <- varcomp[varcomp$grp == "REP", "vcov"]
    H <- G / (G + R / 3)
    
    # Find the corresponding row in outfieldtrials for the current env_value
    row_index <- match(env_value, outfieldtrials$Env)
    
    # Add the calculated values to the corresponding row in outfieldtrials
    outfieldtrials[row_index, c("CV", "H")] <- c(CV, H)
  }
}


# Print the resulting dataframe
print(outfieldtrials)


#===========loop and simulate heritability===========
# Initialize sugaryield_subsets list to store subsets
sugaryield_subsets <- list()

# Sort outfieldtrials by increasing values of H
outfieldtrials_sorted <- outfieldtrials[order(outfieldtrials$H), ]

# Iterate through each row of outfieldtrials_sorted
for (i in 1:nrow(outfieldtrials_sorted)) {
  # Get the Env value with the least H
  level_to_remove <- outfieldtrials_sorted$Env[i]
  
  # Print information about the removal
  cat("Removing Env level:", level_to_remove, "\n")
  
  # Create subset of sugaryield data removing current Env level
  sugaryield_subset <- sugaryield[!(sugaryield$Env == level_to_remove), , drop = FALSE]
  
  # Add the subset to the list
  sugaryield_subsets[[length(sugaryield_subsets) + 1]] <- sugaryield_subset
  
  # Update sugaryield dataframe for the next iteration
  sugaryield <- sugaryield_subset
  
  # Check if only two levels of Env remain
  if (length(unique(sugaryield$Env)) <= 2) {
    break
  }
}

# Print the final number of datasets created
print(length(sugaryield_subsets))


# Check the length of the list (should be 56)
length(sugaryield_subsets)

# Check the number of observations in each subset
subset_sizes <- sapply(sugaryield_subsets, nrow)
subset_sizes
if (nrow(sugaryield_subsets[[1]]) > 0) {
  print(unique(sugaryield_subsets[[1]]$Env))
} else {
  print("No observations remaining in the dataset.")
}

sugaryield <- subset(outfield, !is.na(T_SPACRE) & T_SPACRE != 0)
sugaryield$VARIETY <- factor(sugaryield$VARIETY)
# Fit the linear mixed model
modelall <- lmer(T_SPACRE ~ 1 + (1|VARIETY) + (1|Env) + (1|VARIETY:Env), data = sugaryield)



# Obtain variance components
var_components <- VarCorr(modelall)
varcomp <- as.data.frame(var_components)

# Calculate H
G <- varcomp[varcomp$grp == "VARIETY", "vcov"]
GE <- varcomp[varcomp$grp == "VARIETY:Env", "vcov"]
R <- varcomp[varcomp$grp == "Residual", "vcov"]
H <- G / (G + (GE/60) + (R/180))


blups <- ranef(modelall)$VARIETY
colnames(blups) <- paste0("modelall")
blups$VARIETY <- rownames(blups)
SYBLUPs <- blups



# Initialize dataframe to store correlations
correlation_df <- data.frame(Dataframe_Number = 1:57, Correlation = numeric(57))
# Create a new column for "H" with zeros
H_values <- rep(0, 57)

# Combine with the existing dataframe
correlation_df <- cbind(correlation_df, H = H_values)


# Loop through each dataframe in the sugaryield_subsets list
for (i in 1:length(sugaryield_subsets)) {
  
  # Fit the linear mixed model
  model <- lmer(T_SPACRE ~ 1 + (1|VARIETY) + (1|Env) + (1|VARIETY:Env), data = sugaryield_subsets[[i]])
  
  # Extract BLUPs for the dataframe
  blups <- ranef(model)$VARIETY
  colnames(blups) <- paste0("model_", i)
  blups$VARIETY <- rownames(blups)
  
  SYBLUPs <- merge(SYBLUPs, blups, by = "VARIETY", all.x = TRUE)
  
  # Calculate correlation with full model BLUPs
  cor_value <- cor(SYBLUPs[, paste0("model_", i)], SYBLUPs$modelall)
  
  # Obtain variance components
  var_components <- VarCorr(model)
  varcomp <- as.data.frame(var_components)
  
  # Calculate H
  G <- varcomp[varcomp$grp == "VARIETY", "vcov"]
  GE <- varcomp[varcomp$grp == "VARIETY:Env", "vcov"]
  R <- varcomp[varcomp$grp == "Residual", "vcov"]
  num_levels_Env <- length(unique(sugaryield_subsets[[i]]$Env))
  num_levels_reps <- num_levels_Env*3
  H <- G / (G + (GE/num_levels_Env) + (R/num_levels_reps))
  
  
  # Store correlation in dataframe
  correlation_df[i, "Correlation"] <- cor_value
  correlation_df[i, "H"] <- H
}

# Print correlation dataframe
print(correlation_df)

# Add a column called $Environments to correlation_df
correlation_df$Environments <- seq(57, 57 - nrow(correlation_df) + 1, by = -1)

# Convert Environments to factor
correlation_df$Environments <- as.numeric(correlation_df$Environments)
correlation_33 <- correlation_df %>%
  filter(Environments == 33)
correlation_33$Correlation <- round(correlation_33$Correlation, 3)

# Create the plot
ggplot(correlation_df, aes(x = -Environments, y = Correlation)) +
  geom_point(size = 3, color = "dodgerblue") +  # Increase point size and set color
  geom_line(size = 1.5, color = "dodgerblue") +  # Increase line size and set color
  geom_vline(xintercept = -33, size = 1.2, linetype = "dashed", color = "seagreen") +  # Add vertical line at Environment = 26
  geom_text(data = correlation_33, aes(x = -Environments, y = Correlation, label = Correlation), vjust = -0.5, hjust = -0.2, size = 8, color = "black") +  # Display Correlation value at Environment = 26
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  # Format y-axis ticks
  scale_x_continuous(breaks = c(-59, seq(-55, -10, by = 5), -4), labels = c(59, seq(55, 10, by = -5), 4)) +  # Set x-axis values at intervals of 5 and reverse the limits
  
  labs(x = "Number of Environment", y = "Correlation to Full Dataset BLUPs", title = "Simulation of Reduced Environments") +  # Set axis and plot titles
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 14, face = "bold"))

# Convert Environments to factor
correlation_df$H <- as.numeric(correlation_df$H)
H_17 <- correlation_df %>%
  filter(Environments == 17)
H_17$H <- round(H_17$H, 3)
# Create the plot
ggplot(correlation_df, aes(x = -Environments, y = H)) +
  geom_point(size = 3, color = "dodgerblue") +  # Increase point size and set color
  geom_line(size = 1.5, color = "dodgerblue") +  # Increase line size and set color
  geom_vline(xintercept = -17, size = 1.2, linetype = "dashed", color = "seagreen") +  # Add vertical line at Environment = 26
  geom_text(data = H_17, aes(x = -Environments, y = H, label = H), vjust = -0.5, hjust = -0.2, size = 8, color = "black") +  # Display Correlation value at Environment = 26
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  # Format y-axis ticks
  scale_x_continuous(breaks = c(-59, seq(-55, -10, by = 5), -4), labels = c(59, seq(55, 10, by = -5), 4)) +  # Set x-axis values at intervals of 5 and reverse the limits
  
  labs(x = "Number of Environment", y = "Broad Sense Heritability", title = "Simulation of Reduced Environments") +  # Set axis and plot titles
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 14, face = "bold"))




SYBLUPscorrelogram <- SYBLUPs[, c(2, 12, 22, 32, 42, 52)]
names(SYBLUPscorrelogram)[1] <- "ALLEnvs"
names(SYBLUPscorrelogram)[2] <- "50Envs"
names(SYBLUPscorrelogram)[3] <- "40Envs"
names(SYBLUPscorrelogram)[4] <- "30Envs"
names(SYBLUPscorrelogram)[5] <- "20Envs"
names(SYBLUPscorrelogram)[6] <- "10Envs"

library(psych)

pairs.panels(SYBLUPscorrelogram,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)          # If TRUE, adds confidence intervals

