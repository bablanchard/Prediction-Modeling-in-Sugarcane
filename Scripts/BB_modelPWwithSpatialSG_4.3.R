# Replace 'your_file_path.csv' with the path to your CSV file
SGspatial <- read.csv('C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/SGY22RSpatial.csv')

library(readxl)
library(car)
require(bestNormalize)
library(lme4)
library(sommer)
library(dplyr)
library(ggplot2)
library(caret)

library(nlme)
library(psych)
library(tidyr)
library(plyr)
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
#this model specifies unstructured variance structure with Row and Column as covariates and residual variances are computed for each level of column
solspat <-  mmer(fixed = PlotWeight ~ 1,
                 random = ~ gid + vsr(usr(Row)) + vsr(usr(Column)),
                 rcov = ~ vsr(dsr(Column), units),
                 data = PlotWeightgid)
summary(solspat)
BLUPs_Cross <- as.data.frame(solspat$U$'gid'$PlotWeight + mean(PlotWeightgid$PlotWeight))

#this model represents the model initially planned for the experiment (as a RCBD)
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


