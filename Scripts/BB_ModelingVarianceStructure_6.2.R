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



#============make CY from PlotWeight and convert TRS to kg per Mg and make SY=========
# Convert plot area from square feet to hectares
plot_area_hectares <- 252 * 0.000009290304

# Convert weight from pounds to metric tons (megagrams)
df_filtered$PlotWeight_megagrams <- df_filtered$PlotWeight * 0.00045359237

# Convert weight to megagrams per hectare
df_filtered$CY <- df_filtered$PlotWeight_megagrams / plot_area_hectares

# Convert lbs to kg
df_filtered$TRS <- df_filtered$TRS * 453.592

# Convert tons to Megagrams
df_filtered$TRS <- df_filtered$TRS * (1 / 907.185)

# Rename column header if needed
colnames(df_filtered)[colnames(df_filtered) == "TRS"] <- "TRS_kg_per_Mg"

# Calculate SY (kilograms per hectare)
df_filtered$SY <- df_filtered$TRS_kg_per_Mg * df_filtered$CY

# Convert pounds to kilograms (1 pound = 0.453592 kilograms)
df_filtered$SW_kg <- df_filtered$SW * 0.453592
df_filtered$SW <- df_filtered$SW_kg
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
CY <- df_filtered[complete.cases(df_filtered$CY) & df_filtered$CY != 0,]


#==============BLUPs=================
#Create empty data frame for BLUP output
DataOutput <- data.frame(matrix(vector(),50,1, dimnames=list(c(), c("Entry"))))

#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutput$Entry <- unique(df_filtered[,2]) #fill in Entry numbers
DataOutput$Row <- c(1:50)
DataOutput$Cross <- DataOutput$Entry
DataOutput <- subset(DataOutput, select = -Entry)


#=========sommer=============
CY$Env <- with(CY, interaction(PlantYear, Location, Crop, sep = "_"))

# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
fitMET <- mmer(fixed = CY ~ 1,
               random = ~Cross + Location + Crop + PlantYear + Cross:Location,
               rcov = ~ vsr(units),
               data = CY)

summary(fitMET)
# Fitting genotype by environment models - unstructured model (US)
fitMET.US <- mmer(fixed = CY ~ 1,
                  random = ~Cross + Location + Crop + PlantYear + vsr(usr(Location), Cross),
                  rcov = ~ vsr(units),
                  data = CY)
summary(fitMET.US)

# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
fitMET.US.H <- mmer(fixed = CY ~ 1,
                    random = ~Cross + Location + Crop + PlantYear + vsr(usr(Location), Cross),
                    rcov = ~ vsr(dsr(Location), units),
                    data = CY)

summary(fitMET.US.H)
anova(fitMET, fitMET.US.H)


# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
fitMET <- mmer(fixed = CY ~ 1,
               random = ~Cross + Env + Cross:Env,
               rcov = ~ vsr(units),
               data = CY)

summary(fitMET)

# Fitting genotype by environment models - unstructured model (US)
fitMET.US <- mmer(fixed = CY ~ 1,
                  random = ~Cross + Env + vsr(usr(Env), Cross),
                  rcov = ~ vsr(units),
                  data = CY)
summary(fitMET.US)

# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
fitMET.US.H <- mmer(fixed = CY ~ 1,
                    random = ~Cross + Env,
                    rcov = ~ vsr(dsr(Env), units),
                    data = CY)

summary(fitMET.US.H)


anova(fitMET, fitMET.US)
anova(fitMET.US, fitMET.US.H)

