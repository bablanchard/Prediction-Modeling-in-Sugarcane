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

#examine varieties
print(unique(outfield$VARIETY))
unique_combinations <- outfield %>% distinct(LOC, YEAR, CROP)
print(unique_combinations)
unique_env <- as.data.frame(unique_combinations)



#===========Making ERM=================
#loading packages
library(heatmaply)
require(foreach)
require(doParallel)
library(ggplot2)
require(EnvRtype)
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/Supplementary%20Source%20and%20Data/get_weather_v2')
library(SoilType)
library(caret)
library(stringr)
# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1) # type the number of cores you want to use
getDoParWorkers()

# matching positions
met <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/LSU Breeding Program/2024ReleaseMeeting/outfieldlocs_1.2.csv")

met$Name
dim(met)
head(met)
tail(met)
period <- 2020:2023
trial <- met$Abbreviation  
grid2 <- expand.grid(trial, period)
colnames(grid2) <- c("Abbreviation", "period")
grid2 <- merge(grid2, met)
#grid2$planting <- gsub("2020", "", grid2$planting, fixed = T)
#grid2$harvest <- gsub("2023", "", grid2$harvest, fixed = T)
dim(grid2)
head(grid2)

system.time(
  env.data <- foreach(i = 1:nrow(grid2), 
                      .packages = c("EnvRtype"), 
                      .combine = "rbind",
                      .export = c("get_weather"),
                      .multicombine = TRUE, 
                      .errorhandling = "remove",
                      .verbose = TRUE    
  ) %dopar% {
    
    
    sample <- grid2[i,]  
    
    output <- get_weather(env.id = sample$Abbreviation,
                          country = NULL,
                          lat = sample$Lat,
                          lon = sample$Long,
                          start.day = paste0(sample$planting),
                          end.day = paste0(sample$harvest),
                          parallel = F)
  }  
)  

dim(env.data)
head(env.data)
tail(env.data)
saveRDS(env.data, "OutFieldEnv.data")
unique(env.data$env)

# tuning for cardinal of temperature, and obtaining other traits
aux <- met[,c(2,5)] 
colnames(aux) <- c("env", "Alt")
aux2 <- merge(env.data, aux)
dim(aux2)
df.clim <- processWTH(env.data = aux2,Tbase1 = 4,Tbase2 = 40,Topt1 = 22.5,Topt2 = 35, Alt = aux2$Alt)
head(df.clim)
length(unique(df.clim$env))
dim(df.clim)
summary(df.clim)
sum(is.na(df.clim))

# defining stage interval (day of start), stage names to define de EC for T matrix 
(var.i = names(df.clim)[c(9:16,22:29)])
stages = c('Plant',"30", "61", "90", "121", "154", "184", "215", "245", "276", "306", "337")
interval <- c(0, 30, 61, 90, 121, 154, 184, 215, 245, 276, 306, 337) # do not include the last date

EC <- W_matrix(env.data = df.clim, 
               env.id = 'env', 
               var.id = var.i, 
               statistic = 'mean',
               scale = F, 
               center = F,
               by.interval = T, 
               sd.tol = 5,
               time.window = interval, 
               names.window = stages)


dim(EC)
W <- EC
sum(W == "NaN")
W[W == "NaN"] <- 0
W <- scale(W)
dim(W)
W[, 1:4]
# Export matrix to CSV
write.csv(W, file = "W_matrix.csv", row.names = FALSE)

#########################################################################
# Soil information
#########################################################################

# running
system.time(
  SC <- foreach(i = 1:nrow(met), 
                .packages = c("caret", "stringr"), 
                .combine = "rbind",
                .export = c("predict", "train", "trainControl", "str_replace_all"),
                .multicombine = TRUE, 
                .errorhandling = "stop",
                .verbose = TRUE    
  ) %do% {
    # subset the data  
    trial.MET <- droplevels.data.frame(met[i,])
    # retrieve the data
    output <- get_soil(env.id = trial.MET$Abbreviation, 
                       lat = trial.MET$Lat, 
                       long = trial.MET$Long, 
                       max.depth = 20, 
                       isric.data = soil.data)
  }   
)


head(SC)
tail(SC)
dim(SC)

# Soil Covariates visualization
SCov <- reshape2::dcast(SC, env ~ Trait, value.var = "Predicted.Value", mean)
head(SCov)
dim(SCov)
rownames(SCov) <- SCov[,1]
SCov <- SCov[,2:ncol(SCov)]

# eliminate soil traits without any information
SCov <- SCov[, apply(SCov, 2, function(x){!all(is.na(x))})]
SCov[is.na(SCov)] <- 0
SCov <- scale(SCov)
dim(SCov)

W <- W[match(rownames(SCov), rownames(W)),]
all(rownames(W) == rownames(SCov))

EC <- cbind(W, SCov)
dim(EC)

################# Quality Control and feature selection ###############
##### Correlation - eliminating those almost perfected correlated #####
library(caret)
# Find columns with missing values
cols_with_missing <- apply(EC, 2, function(x) any(is.na(x)))
ECc <- EC[, !cols_with_missing]
EC <- ECc

highlyCorrelated <- findCorrelation(cor(EC), cutoff = 0.98)
print(highlyCorrelated)
EC.clean <- EC[,-highlyCorrelated]
dim(EC.clean)
EC.clean[EC.clean == "NaN"] <- 0
EC.clean[EC.clean == "NA"] <- 0




install.packages("RColorBrewer")
library(RColorBrewer)

# Choose a color palette
color_palette <- brewer.pal(9, "YlOrRd")
# plot the EC - like environmental markers
heatmaply(EC.clean, 
          fontsize_row = 6,
          fontsize_col = 6,
          col = color_palette,
          file = c("EC_heatmapoutfield.html", "EC_heatmap.png"))


#===============Extract BLUPs per Env using CYfit2===========
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


sugaryield <- subset(outfield, !is.na(T_SPACRE) & T_SPACRE != 0)
sugaryield$VARIETY <- factor(sugaryield$VARIETY)

SYfit2 <- mmer(fixed = T_SPACRE ~ 1 ,
               random = ~ VARIETY + Env + VARIETY:Env,
               rcov = ~ vsr(units),
               data = sugaryield)


# Create a data frame with meaningful column names
SYBLUPs <- as.data.frame(SYfit2$U$Env$T_SPACRE + mean(sugaryield$T_SPACRE))
colnames(SYBLUPs) <- "SYfit2"  # Assign a column name to BLUPs_Cross
SYBLUPs$Env <- rownames(SYBLUPs)
SYBLUPs$Env <- gsub("^Env", "", SYBLUPs$Env)

#=======attach all ECs to the SYBLUPs dataframe==================
# Assuming 'Env' column is character type in SYBLUPs
SYBLUPs$Env <- as.character(SYBLUPs$Env)

# Find unique values in Env column of SYBLUPs
unique_env <- unique(SYBLUPs$Env)

# Create an empty list to store merged dataframes
merged_data_list <- list()

for (env_value in unique_env) {
  # Subset SYBLUPs for the current env value
  subset_SYBLUPs <- subset(SYBLUPs, Env == env_value)
  
  # Find matching row names in EC.clean for the current env value
  matching_row_names <- rownames(EC.clean)[rownames(EC.clean) %in% subset_SYBLUPs$Env]
  
  # Subset EC.clean based on matching row names
  subset_EC_clean <- EC.clean[matching_row_names, , drop = FALSE]  # Ensure it remains a matrix
  
  # Create a dataframe from the subset of EC.clean with row names preserved
  subset_EC_clean_df <- data.frame(subset_EC_clean, row.names = matching_row_names)
  
  # Merge subset_SYBLUPs and subset_EC_clean_df based on row names
  merged_data <- merge(subset_SYBLUPs, subset_EC_clean_df, by.x = "Env", by.y = "row.names", all.x = TRUE)
  
  # Store merged dataframe in the list
  merged_data_list[[env_value]] <- merged_data
}




# Combine all merged dataframes into one dataframe
final_merged_data <- do.call(rbind, merged_data_list)
write.csv(final_merged_data, file = "SYBLUPswECs.csv", row.names = FALSE)


# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1) # type the number of cores you want to use
getDoParWorkers()
control <- rfeControl(functions = rfFuncs, method = "repeatedcv", number = 5, repeats = 5, allowParallel = TRUE)
opt.reg <- rfe(x = final_merged_data[, 3:ncol(final_merged_data)], y = final_merged_data$SYfit2, sizes = c(2:154),
               metric = "Rsquared", maximize = TRUE, rfeControl = control)
opt.reg
RFSresults <- as.data.frame(opt.reg$results)
write.csv(RFSresults, file = "RFSresults.csv", row.names = FALSE)

predictors(opt.reg)

opt.reg$bestSubset

(traits.reg <- opt.reg$optVariables[1:which(round(opt.reg$results$Rsquared, 2) == max(round(opt.reg$results$Rsquared, 2)))[1]])

ggplot(opt.reg)

sink("the most important EC to explain SY via ML.txt")

opt.reg

sink()

sink("EC.predictors.txt")

data.frame(EC = predictors(opt.reg))

sink()

# Cleaning W again using RFE

EC.clean <- EC.clean[, predictors(opt.reg)]

dim(EC.clean)

#=================determine most important ECs for SY=================




# plot the EC - like environmental markers
heatmaply(EC.clean, 
          fontsize_row = 6,
          fontsize_col = 6,
          col = color_palette,
          file = c("EC_heatmapoutfield.html", "EC_heatmap.png"))

# Obtain the environmental relationship matrix - ERM
kinship <- env_kernel(env.data = as.matrix(EC.clean), is.scaled = T, gaussian = T, sd.tol = 5)[[2]]
dim(kinship)
kinship
saveRDS(kinship, "W_kinshipoutfieldenv.")
heatmaply(kinship, 
          fontsize_row = 6,
          fontsize_col = 6,
          col = color_palette,
          file = c("EC_kinship_heatmapoutfieldenv.html", "EC_kinship_heatmap.png"))

#=============test ERM to see if it provides increased fitness to the prediction model============
SYfitERM <- mmer(fixed = T_SPACRE ~ 1 ,
                 random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
                 rcov = ~ vsr(units),
                 data = sugaryield)

summary(SYfit2)
summary(SYfitERM)

anova(SYfitERM, SYfit2)

#==============clustering========================
# clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

# or we can use this function to determine the optimal number of clusters
set.seed(123)
dev.off()
fviz_nbclust(EC.clean, 
             FUNcluster = kmeans, 
             k.max = nrow(EC.clean)-1, 
             nboot = nrow(EC.clean)/2, 
             method = "wss")

# K-Means Clustering
k <- kmeans(scale(EC.clean), centers = 6, nstart = 25)
p <- fviz_cluster(k, data = EC.clean)
p

ggsave(filename = paste0("./", 'outfieldenvcluster.tiff'),
       plot = p,
       device = 'tiff',
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)

k$size
sum(k$size)

(clusters <- data.frame(Name = names(k$cluster), 
                        cluster = k$cluster))

dim(clusters)
write.table(clusters, "outfieldenvclusters.txt")







#================Simulation========================
# Create an empty list to store subsets of sugaryield dataframe
sugaryield_subsets <- list()

# Initialize vector to keep track of removed Env levels
removed_envs <- character()

# Initialize sugaryield dataframe for the first iteration
sugaryield_current <- sugaryield

# Loop until all levels of Env are removed
removed_count <- 0
while (removed_count < 54) {
  # Loop through each cluster
  for (clust in 1:6) {
    # Get levels of Env for current cluster
    env_levels <- unique(clusters[clusters$cluster == clust, "Name"])
    
    # Filter out already removed Env levels for the current cluster
    available_env <- setdiff(env_levels, removed_envs)
    
    # Check if all Env levels for the current cluster are removed
    if (length(available_env) == 0) {
      next
    }
    
    # Randomly select an Env level
    level_to_remove <- sample(available_env, 1)
    
    # Add the removed Env level to the list
    removed_envs <- c(removed_envs, level_to_remove)
    
    # Print information about the removal
    cat("Removing Env level:", level_to_remove, "from cluster", clust, "\n")
    
    # Create subset of sugaryield data removing current Env level
    sugaryield_subset <- sugaryield_current[!(sugaryield_current$Env == level_to_remove), , drop = FALSE]
    
    # Add the subset to the list
    sugaryield_subsets[[length(sugaryield_subsets) + 1]] <- sugaryield_subset
    
    # Update sugaryield dataframe for the next iteration
    sugaryield_current <- sugaryield_subset
    
    # Increment counter
    removed_count <- removed_count + 1
    
    # Check if all levels of Env have been removed
    if (removed_count >= 54) {
      break
    }
  }
}

# Check the length of the list (should be 56)
length(sugaryield_subsets)

# Check the number of observations in each subset
subset_sizes <- sapply(sugaryield_subsets, nrow)
subset_sizes



SYfitERM <- mmer(fixed = T_SPACRE ~ 1 ,
                 random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
                 rcov = ~ vsr(units),
                 data = sugaryield)

# Create a data frame with meaningful column names
SYBLUPs <- as.data.frame(SYfitERM$U$"VARIETY"$T_SPACRE)
rownames(SYBLUPs) <- sub("^VARIETY", "", rownames(SYBLUPs))
colnames(SYBLUPs) <- "SYfit2"  # Assign a column name to BLUPs_Cross
SYBLUPs$VARENV <- rownames(SYBLUPs)





# Initialize dataframe to store correlations
correlation_df <- data.frame(Dataframe_Number = 1:54, Correlation = numeric(54))

# Loop through each dataframe in the sugaryield_subsets list
for (i in 1:length(sugaryield_subsets)) {
  # Perform the mixed model prediction on the dataframe
  SYfit <- mmer(fixed = T_SPACRE ~ 1,
                random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
                rcov = ~ vsr(units),
                data = sugaryield_subsets[[i]])
  
  # Extract BLUPs for the dataframe
  SYBLUPs_reduced <- as.data.frame(SYfit$U$"VARIETY"$T_SPACRE)
  rownames(SYBLUPs_reduced) <- sub("^VARIETY", "", rownames(SYBLUPs_reduced))
  colnames(SYBLUPs_reduced) <- paste0("SYfit_", i)
  SYBLUPs_reduced$VARENV <- rownames(SYBLUPs_reduced)
  
  SYBLUPs <- merge(SYBLUPs, SYBLUPs_reduced, by = "VARENV", all.x = TRUE)
  
  # Calculate correlation with full model BLUPs
  cor_value <- cor(SYBLUPs[, paste0("SYfit_", i)], SYBLUPs$SYfit2)
  
  # Store correlation in dataframe
  correlation_df[i, "Correlation"] <- cor_value
}

# Print correlation dataframe
print(correlation_df)

# Add a column called $Environments to correlation_df
correlation_df$Environments <- seq(59, 59 - nrow(correlation_df) + 1, by = -1)

# Convert Environments to factor
correlation_df$Environments <- as.numeric(correlation_df$Environments)
correlation_28 <- correlation_df %>%
  filter(Environments == 28)
correlation_28$Correlation <- round(correlation_28$Correlation, 3)

# Create the plot
ggplot(correlation_df, aes(x = -Environments, y = Correlation)) +
  geom_point(size = 3, color = "dodgerblue") +  # Increase point size and set color
  geom_line(size = 1.5, color = "dodgerblue") +  # Increase line size and set color
  geom_vline(xintercept = -28, size = 1.2, linetype = "dashed", color = "seagreen") +  # Add vertical line at Environment = 26
  geom_text(data = correlation_28, aes(x = -Environments, y = Correlation, label = Correlation), vjust = -0.5, hjust = -0.2, size = 8, color = "black") +  # Display Correlation value at Environment = 26
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  # Format y-axis ticks
  scale_x_continuous(breaks = c(-59, seq(-55, -10, by = 5), -4), labels = c(59, seq(55, 10, by = -5), 4)) +  # Set x-axis values at intervals of 5 and reverse the limits
  
  labs(x = "Number of Environment", y = "Correlation to Full Dataset BLUPs", title = "Simulation of Reduced Environments") +  # Set axis and plot titles
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


# Merge the dataframes based on matching values in data$Env and met$Abbreviation
clusters <- clusters %>% 
  rename(Env = Name)
merged_data <- merge(met, clusters[, c("Env", "cluster")], by.x = "Abbreviation", by.y = "Env", all.x = TRUE)
write.csv(merged_data, "outfieldlocwithcluster.csv", row.names = FALSE)

merged_data <- read.csv("outfieldlocwithcluster.csv")


#=============map================
# Assuming 'merged_data' is your dataset
library(leaflet)

# Define a function to add noise to coordinates to prevent overlap
add_noise <- function(coord) {
  noise <- runif(length(coord), min = -0.006, max = 0.006)  # Adjust noise range as needed
  return(coord + noise)
}


# Apply noise to longitude and latitude
merged_data$Long <- add_noise(merged_data$Long)
merged_data$Lat <- add_noise(merged_data$Lat)

# Create a color palette with distinct colors for each cluster
cluster_colors <- rainbow(length(unique(merged_data$"cluster")))
library(xts)


map <- leaflet(data = merged_data) %>%
  addTiles()
# Create a leaflet map
map <- leaflet(data = merged_data) %>%
  addTiles() %>%
  addLegend(
    position = "bottomright",
    colors = cluster_colors,
    labels = unique(merged_data$"cluster"),
    title = "Cluster"
  )

# Assign colors to markers based on cluster values
for (i in unique(merged_data$"cluster")) {
  map <- addCircleMarkers(
    map,
    data = merged_data[merged_data$"cluster" == i, ],
    lng = ~Long,
    lat = ~Lat,
    radius = 7,
    fillColor = cluster_colors[i],
    fillOpacity = 0.5,
    stroke = FALSE,
    color = cluster_colors[i],
    group = paste0("cluster_", i)  # Group markers by cluster
  )
}

# Add layer control to toggle marker groups
map <- addLayersControl(
  map,
  baseGroups = "markers",
  overlayGroups = paste0("cluster_", unique(merged_data$"cluster")),
  options = layersControlOptions(collapsed = FALSE)
)

# Display the map
map




