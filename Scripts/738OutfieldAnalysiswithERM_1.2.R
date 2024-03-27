library(readxl)
library(lme4)
library(sommer)
library(dplyr)
library(statgenGxE)
library(ggplot2)
library(metan)
library(broom)
library(tidyr)
df <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/LSU Breeding Program/2024ReleaseMeeting/active23_1.2.csv", header = TRUE)

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
met <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/LSU Breeding Program/2024ReleaseMeeting/outfieldlocs.csv")

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
k <- kmeans(scale(EC.clean), centers = 4, nstart = 25)
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

######## the end #################
kinship <- readRDS("W_kinshipoutfieldenv.")


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


#For CY
caneyield <- subset(outfield, !is.na(TCA) & TCA != 0)
caneyield$VARIETY <- factor(caneyield$VARIETY)

library(sommer)
CYfit <- mmer(fixed = TCA ~ 1,
              random = ~ VARIETY + LOC + CROP + YEAR + VARIETY:LOC,
              rcov = ~ vsr(units),
              data = caneyield)
summary(CYfit)

CYfit3 <- mmer(fixed = TCA ~ 1,
               random = ~ VARIETY + Env + VARIETY:Env,
               rcov = ~ vsr(units),
               data = caneyield)
summary(CYfit3)


CYfit2 <- mmer(fixed = TCA ~ 1 ,
               random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
               rcov = ~ vsr(units),
               data = caneyield)
summary(CYfit2)

anova(CYfit, CYfit2)


# Extract unique values of VARIETY
unique_varieties <- unique(caneyield$VARIETY)

# Create a dataframe GVs with VARIETY column
GVs <- data.frame(VARIETY = unique_varieties)


# Create a data frame with meaningful column names
BLUPs_VARIETY <- as.data.frame(CYfit$U$"VARIETY"$TCA + mean(caneyield$TCA))
rownames(BLUPs_VARIETY) <- sub("^VARIETY", "", rownames(BLUPs_VARIETY))
colnames(BLUPs_VARIETY) <- "CYfit"  # Assign a column name to BLUPs_Cross

# Create a data frame with meaningful column names
BLUPs_CYfit3 <- as.data.frame(CYfit3$U$"VARIETY"$TCA + mean(caneyield$TCA))
rownames(BLUPs_CYfit3) <- sub("^VARIETY", "", rownames(BLUPs_CYfit3))
colnames(BLUPs_CYfit3) <- "CYfit3"  # Assign a column name to BLUPs_Cross

# Create a data frame with meaningful column names
BLUPs_CYfit2 <- as.data.frame(CYfit2$U$"VARIETY"$TCA + mean(caneyield$TCA))
rownames(BLUPs_CYfit2) <- sub("^VARIETY", "", rownames(BLUPs_CYfit2))
colnames(BLUPs_CYfit2) <- "CYfit2"  # Assign a column name to BLUPs_Cross

GVs$CYfit <- BLUPs_VARIETY[match(GVs$VARIETY, rownames(BLUPs_VARIETY)), "CYfit"]
GVs$CYfit3 <- BLUPs_CYfit3[match(GVs$VARIETY, rownames(BLUPs_CYfit3)), "CYfit3"]
GVs$CYfit2 <- BLUPs_CYfit2[match(GVs$VARIETY, rownames(BLUPs_CYfit2)), "CYfit2"]

# Compute the average TCA for each value of GVs$VARIETY in caneyield
average_TCA <- aggregate(TCA ~ VARIETY, data = caneyield, FUN = mean)

# Merge the average TCA values to GVs by matching VARIETY values
GVs <- merge(GVs, average_TCA, by = "VARIETY", suffixes = c("", ".mean"))

# Rename the merged column to GVs$RawMean
names(GVs)[ncol(GVs)] <- "RawMean"
# Arrange plots side by side
library(gridExtra)
# Create density plot overlaying GVs$CYfit and GVs$RawMean
plot_overlay <- ggplot(GVs, aes(x = CYfit)) +
  geom_density(fill = "blue", alpha = 0.5) +
  geom_density(aes(x = RawMean), fill = "red", alpha = 0.5) +
  geom_density(aes(x = CYfit3), fill = "purple", alpha = 0.5) +
  geom_density(aes(x = CYfit2), fill = "seagreen", alpha = 0.5) +
  labs(title = "Density Plot Overlay of GVs$CYfit and GVs$RawMean",
       x = "Values",
       y = "Density") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

# Display the overlay plot
print(plot_overlay)


data_subset <- GVs[, c(2:5)]

library(psych)

pairs.panels(data_subset,
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


#===========now reduce the number of environments based on the ERM==============
data <- read.table("outfieldenvclusters.txt", header = TRUE, sep = "\t")
data <- read.csv("outfieldenvclusters.txt")
data$cluster <- as.numeric(substring(data$Name.cluster, nchar(data$Name.cluster), nchar(data$Name.cluster)))
data$Env <- sub("\\s.*", "", data$Name.cluster)
data <- subset(data, select = -Name.cluster)

# Group data by cluster and sample two random Env values per cluster
library(dplyr)
sampled_env <- data %>%
  group_by(cluster) %>%
  sample_n(6) %>%
  select(cluster, Env)

# Merge sampled_env with caneyield to subset caneyield by sampled Env values
subset_caneyield <- caneyield[caneyield$Env %in% sampled_env$Env, ]

CYfit4 <- mmer(fixed = TCA ~ 1,
              random = ~ VARIETY + LOC + CROP + YEAR + VARIETY:LOC,
              rcov = ~ vsr(units),
              data = subset_caneyield)
summary(CYfit4)

CYfit5 <- mmer(fixed = TCA ~ 1,
               random = ~ VARIETY + Env + VARIETY:Env,
               rcov = ~ vsr(units),
               data = subset_caneyield)
summary(CYfit5)


CYfit6 <- mmer(fixed = TCA ~ 1 ,
               random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
               rcov = ~ vsr(units),
               data = subset_caneyield)
summary(CYfit6)

# Create a data frame with meaningful column names
BLUPs_CYfit4 <- as.data.frame(CYfit4$U$"VARIETY"$TCA + mean(subset_caneyield$TCA))
rownames(BLUPs_CYfit4) <- sub("^VARIETY", "", rownames(BLUPs_CYfit4))
colnames(BLUPs_CYfit4) <- "CYfit4"  # Assign a column name to BLUPs_Cross

# Create a data frame with meaningful column names
BLUPs_CYfit5 <- as.data.frame(CYfit5$U$"VARIETY"$TCA + mean(subset_caneyield$TCA))
rownames(BLUPs_CYfit5) <- sub("^VARIETY", "", rownames(BLUPs_CYfit5))
colnames(BLUPs_CYfit5) <- "CYfit5"  # Assign a column name to BLUPs_Cross

# Create a data frame with meaningful column names
BLUPs_CYfit6 <- as.data.frame(CYfit6$U$"VARIETY"$TCA + mean(subset_caneyield$TCA))
rownames(BLUPs_CYfit6) <- sub("^VARIETY", "", rownames(BLUPs_CYfit6))
colnames(BLUPs_CYfit6) <- "CYfit6"  # Assign a column name to BLUPs_Cross

GVs$CYfit4 <- BLUPs_CYfit4[match(GVs$VARIETY, rownames(BLUPs_CYfit4)), "CYfit4"]
GVs$CYfit5 <- BLUPs_CYfit5[match(GVs$VARIETY, rownames(BLUPs_CYfit5)), "CYfit5"]
GVs$CYfit6 <- BLUPs_CYfit6[match(GVs$VARIETY, rownames(BLUPs_CYfit6)), "CYfit6"]

data_subset <- GVs[, c(2:8)]

library(psych)

pairs.panels(data_subset,
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



# Create a boxplot of TCA for every VARIETY
ggplot(caneyield, aes(x = VARIETY, y = TCA)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Boxplot of TCA for every VARIETY",
       x = "VARIETY",
       y = "TCA")


#===============Extract BLUPs per Env using CYfit2===========
sugaryield <- subset(outfield, !is.na(T_SPACRE) & T_SPACRE != 0)
sugaryield$VARIETY <- factor(sugaryield$VARIETY)

SYfit2 <- mmer(fixed = T_SPACRE ~ 1 ,
               random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
               rcov = ~ vsr(units),
               data = sugaryield)

TRS <- subset(outfield, !is.na(TRS_TON) & TRS_TON != 0)
TRS$VARIETY <- factor(TRS$VARIETY)

TRSfit2 <- mmer(fixed = TRS_TON ~ 1 ,
               random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
               rcov = ~ vsr(units),
               data = TRS)

# Create a data frame with meaningful column names
CYBLUPs <- as.data.frame(CYfit2$U$"VARIETY:Env"$TCA + mean(caneyield$TCA))
rownames(CYBLUPs) <- sub("^VARIETY", "", rownames(CYBLUPs))
colnames(CYBLUPs) <- "CYfit2"  # Assign a column name to BLUPs_Cross
CYBLUPs$VARENV <- rownames(CYBLUPs)

# Create a data frame with meaningful column names
TRSBLUPs <- as.data.frame(TRSfit2$U$"VARIETY:Env"$TRS_TON + mean(TRS$TRS_TON))
rownames(TRSBLUPs) <- sub("^VARIETY", "", rownames(TRSBLUPs))
colnames(TRSBLUPs) <- "TRSfit2"  # Assign a column name to BLUPs_Cross
TRSBLUPs$VARENV <- rownames(TRSBLUPs)

# Create a data frame with meaningful column names
SYBLUPs <- as.data.frame(SYfit2$U$"VARIETY:Env"$T_SPACRE + mean(sugaryield$T_SPACRE))
rownames(SYBLUPs) <- sub("^VARIETY", "", rownames(SYBLUPs))
colnames(SYBLUPs) <- "SYfit2"  # Assign a column name to BLUPs_Cross
SYBLUPs$VARENV <- rownames(SYBLUPs)

# Merge CYBLUPs, TRSBLUPs, and SYBLUPs based on the 'VARENV' column
merged_df <- merge(merge(CYBLUPs, TRSBLUPs, by = "VARENV", all = TRUE), SYBLUPs, by = "VARENV", all = TRUE)
GVsEnv <- merged_df
# Split the 'VARENV' column by the ':' character
split_values <- strsplit(as.character(GVsEnv$VARENV), ":")
# Extract the first and second parts into separate columns
GVsEnv$VARIETY <- sapply(split_values, `[`, 1)
GVsEnv$Env <- sapply(split_values, `[`, 2)
# Remove the last character from 'Env' column
GVsEnv$CROP <- substr(GVsEnv$Env, nchar(GVsEnv$Env), nchar(GVsEnv$Env))
# Update 'Env' column to remove the last character
GVsEnv$Env <- substr(GVsEnv$Env, 1, nchar(GVsEnv$Env) - 1)
# Remove the last two characters from 'Env' column
trimmed_chars <- substr(GVsEnv$Env, 1, nchar(GVsEnv$Env) - 2)
# Attach "20" to the removed characters
year <- paste0("20", substr(GVsEnv$Env, nchar(GVsEnv$Env) - 1, nchar(GVsEnv$Env)))
# Create 'YEAR' column
GVsEnv$YEAR <- year
GVsEnv$Env <- substr(GVsEnv$Env, 1, nchar(GVsEnv$Env) - 2)
# Remove "Env" from the beginning of each value in 'Env' column
GVsEnv$LOC <- gsub("^Env", "", GVsEnv$Env)
# Remove the 'Env' column
GVsEnv <- subset(GVsEnv, select = -c(Env))
# Convert specific columns to factors
GVsEnv$VARIETY <- factor(GVsEnv$VARIETY)
GVsEnv$CROP <- factor(GVsEnv$CROP)
GVsEnv$YEAR <- factor(GVsEnv$YEAR)
GVsEnv$LOC <- factor(GVsEnv$LOC)

# Aggregate data by unique combination of VARIETY and CROP and calculate mean
GVsCrop <- aggregate(cbind(SYfit2, CYfit2, TRSfit2) ~ VARIETY + CROP, data = GVsEnv, FUN = mean)
# Rename columns
names(GVsCrop)[names(GVsCrop) == "SYfit2"] <- "T_SPACRE"
names(GVsCrop)[names(GVsCrop) == "CYfit2"] <- "TCA"
names(GVsCrop)[names(GVsCrop) == "TRSfit2"] <- "TRS_TON"

# Save the dataframe as a CSV file
write.csv(GVsCrop, "OutfieldGVsCrop.csv", row.names = FALSE)


#==============Assess the number of environments that provides the same estimates=========
SYfit2 <- mmer(fixed = T_SPACRE ~ 1 ,
               random = ~ VARIETY + vsr(Env, Gu = kinship) + VARIETY:Env,
               rcov = ~ vsr(units),
               data = sugaryield)

# Create a data frame with meaningful column names
SYBLUPs <- as.data.frame(SYfit2$U$"VARIETY"$T_SPACRE)
rownames(SYBLUPs) <- sub("^VARIETY", "", rownames(SYBLUPs))
colnames(SYBLUPs) <- "SYfit2"  # Assign a column name to BLUPs_Cross
SYBLUPs$VARENV <- rownames(SYBLUPs)


# Create an empty list to store subsets of sugaryield dataframe
sugaryield_subsets <- list()

# Initialize vector to keep track of removed Env levels
removed_envs <- character()

# Initialize sugaryield dataframe for the first iteration
sugaryield_current <- sugaryield

# Loop until all levels of Env are removed
removed_count <- 0
while (removed_count < 56) {
  # Loop through each cluster
  for (clust in 1:4) {
    # Get levels of Env for current cluster
    env_levels <- unique(data[data$cluster == clust, "Env"])
    
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
    if (removed_count >= 56) {
      break
    }
  }
}

# Check the length of the list (should be 56)
length(sugaryield_subsets)

# Check the number of observations in each subset
subset_sizes <- sapply(sugaryield_subsets, nrow)
subset_sizes



# Initialize dataframe to store correlations
correlation_df <- data.frame(Dataframe_Number = 1:56, Correlation = numeric(56))

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
correlation_26 <- correlation_df %>%
  filter(Environments == 26)
correlation_26$Correlation <- round(correlation_26$Correlation, 3)

# Create the plot
ggplot(correlation_df, aes(x = -Environments, y = Correlation)) +
  geom_point(size = 3, color = "dodgerblue") +  # Increase point size and set color
  geom_line(size = 1.5, color = "dodgerblue") +  # Increase line size and set color
  geom_vline(xintercept = -26, size = 1.2, linetype = "dashed", color = "seagreen") +  # Add vertical line at Environment = 26
  geom_text(data = correlation_26, aes(x = -Environments, y = Correlation, label = Correlation), vjust = -0.5, hjust = -0.2, size = 8, color = "black") +  # Display Correlation value at Environment = 26
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
merged_data <- merge(met, data[, c("Env", "cluster")], by.x = "Abbreviation", by.y = "Env", all.x = TRUE)
write.csv(merged_data, "outfieldlocwithcluster.csv", row.names = FALSE)

