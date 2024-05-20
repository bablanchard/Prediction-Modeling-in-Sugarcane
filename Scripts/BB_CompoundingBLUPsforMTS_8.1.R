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




#==============BLUPs=================
#Create empty data frame for BLUP output
DataOutput <- data.frame(matrix(vector(),50,1, dimnames=list(c(), c("Entry"))))

#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutput$Entry <- unique(df_filtered[,2]) #fill in Entry numbers
DataOutput$Row <- c(1:50)
DataOutput$Cross <- DataOutput$Entry
DataOutput <- subset(DataOutput, select = -Entry)

#============TRS=============
#Combine columns to create the 'Env' factor
TRS <- TRS %>%
  mutate(Env = interaction(PlantYear, Crop))
table(TRS$Env)

TRS$FxL <- interaction(TRS$Cross, 
                         TRS$Location)
TRS$Inter <- interaction(TRS$Cross, 
                       TRS$Env)


#now make a genotype id dataframe
# Create a data frame with unique values from TRS$Cross
unique_cross_values <- data.frame(Cross = unique(TRS$Cross))
# Create a numeric ID column called "gid"
unique_cross_values$gid <- seq_along(unique_cross_values$Cross)
# Merge the unique values and their IDs with the original data frame
TRSgid <- merge(TRS, unique_cross_values, by = "Cross")
# The "gid" data frame will have a new column "gid" with numeric IDs for each value in TRS$Cross
print(TRSgid)
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

# Extract unique combinations of Cross and gid from TRSgid
ID <- unique(TRSgid[, c("Cross", "gid")])

# Print the ID dataframe
print(ID)
#LMM with specified heterogenous variance per kinship
fitMET2 <- mmer(fixed = TRS ~ 1,
                random = ~ vsr(gid, Gu = Acsv_matrix) + Env + Location + FxL + Inter,
                rcov = ~ vsr(units),
                data = TRSgid)
summary(fitMET2)


BLUPs_Cross <- as.data.frame(fitMET2$U$"u:gid"$TRS + mean(TRSgid$TRS))
rownames(BLUPs_Cross) <- sub("^gid", "", rownames(BLUPs_Cross))
colnames(BLUPs_Cross) <- "prediction"

# Merge ID dataframe with BLUPs_Cross based on gid
ID_with_predictions <- merge(ID, BLUPs_Cross, by.x = "gid", by.y = "row.names", all.x = TRUE)

# Print the resulting dataframe
print(ID_with_predictions)

# Merge DataOutput with ID_with_predictions based on Cross
MergedData <- merge(DataOutput, ID_with_predictions, by = "Cross", all.x = TRUE)

# Add a new column TRSBLUPSALL with the predicted values
MergedData$TRSBLUPSALL <- MergedData$prediction

# Drop unnecessary columns if needed
DataOutput <- MergedData[, c("Cross", "TRSBLUPSALL", "gid")]



TRSBLUPSperloc <- as.data.frame(fitMET2$U$FxL$TRS)
# Assuming TRSBLUPSperloc is your dataframe
# Extract row names and create new columns
TRSBLUPSperloc$RowNames <- rownames(TRSBLUPSperloc)
TRSBLUPSperloc$Location <- sub(".*\\.", "", TRSBLUPSperloc$RowNames)
TRSBLUPSperloc$Crosssub <- sub("^...|\\..*", "", TRSBLUPSperloc$RowNames)
TRSBLUPSperloc$Cross <- sub("\\..*", "", TRSBLUPSperloc$Crosssub)


# Drop the original row names column if not needed
TRSBLUPSperloc <- TRSBLUPSperloc[, c("RowNames", "Location", "Cross", colnames(TRSBLUPSperloc)[-ncol(TRSBLUPSperloc)])]

colnames(TRSBLUPSperloc)[colnames(TRSBLUPSperloc) == "fitMET2$U$FxL$TRS"] <- "prediction"
# Match values and create new column
DataOutput$TRSBLUPSNI <- DataOutput$TRSBLUPSALL + TRSBLUPSperloc$prediction[TRSBLUPSperloc$Cross %in% DataOutput$Cross & TRSBLUPSperloc$Location == "NI"]
DataOutput$TRSBLUPSNR <- DataOutput$TRSBLUPSALL + TRSBLUPSperloc$prediction[TRSBLUPSperloc$Cross %in% DataOutput$Cross & TRSBLUPSperloc$Location == "NR"]
DataOutput$TRSBLUPSSG <- DataOutput$TRSBLUPSALL + TRSBLUPSperloc$prediction[TRSBLUPSperloc$Cross %in% DataOutput$Cross & TRSBLUPSperloc$Location == "SG"]




#============PlotWeight=============
#Combine columns to create the 'Env' factor
PlotWeight <- PlotWeight %>%
  mutate(Env = interaction(PlantYear, Crop))
table(PlotWeight$Env)

PlotWeight$FxL <- interaction(PlotWeight$Cross, 
                              PlotWeight$Location)
PlotWeight$Inter <- interaction(PlotWeight$Cross, 
                                PlotWeight$Env)


#now make a genotype id dataframe
# Create a data frame with unique values from TRS$Cross
unique_cross_values <- data.frame(Cross = unique(PlotWeight$Cross))
# Create a numeric ID column called "gid"
unique_cross_values$gid <- seq_along(unique_cross_values$Cross)
# Merge the unique values and their IDs with the original data frame
PlotWeightgid <- merge(PlotWeight, unique_cross_values, by = "Cross")
# The "gid" data frame will have a new column "gid" with numeric IDs for each value in TRS$Cross
unique_values <- unique(PlotWeightgid$Cross)
print(unique_values)

Acsv <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/FxEsubsetAmatrix.csv")
unique_values <- unique(Acsv$X)
print(unique_values)


# Assuming Acsv$X needs to be replaced with PlotWeightgid$gid based on matching values in PlotWeightgid$Cross

# Match the indices of the Cross values in Acsv$X with PlotWeightgid$Cross and replace X with corresponding gid values
matched_indices <- match(Acsv$X, PlotWeightgid$Cross)

# Replace Acsv$X with corresponding values from PlotWeightgid$gid
Acsv$X <- PlotWeightgid$gid[matched_indices]

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
  match_index <- match(col, PlotWeightgid$Cross)
  
  # If a match is found, replace the column name with the corresponding PlotWeightgid$gid value
  if (!is.na(match_index)) {
    new_col_name <- PlotWeightgid$gid[match_index]
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
length(unique(PlotWeightgid$gid))

PlotWeightgid$gid <- as.factor(PlotWeightgid$gid)
PlotWeightgid$PlotWeight <- as.numeric(PlotWeightgid$PlotWeight)

# Extract unique combinations of Cross and gid from TRSgid
ID <- unique(PlotWeightgid[, c("Cross", "gid")])

#LMM with specified heterogenous variance per kinship
fitMET2 <- mmer(fixed = PlotWeight ~ 1,
                random = ~ vsr(gid, Gu = Acsv_matrix) + Env + Location + FxL + Inter,
                rcov = ~ vsr(units),
                data = PlotWeightgid)
summary(fitMET2)


BLUPs_Cross <- as.data.frame(fitMET2$U$"u:gid"$PlotWeight + mean(PlotWeightgid$PlotWeight))
rownames(BLUPs_Cross) <- sub("^gid", "", rownames(BLUPs_Cross))
colnames(BLUPs_Cross) <- "prediction"

# Merge ID dataframe with BLUPs_Cross based on gid
ID_with_predictions <- merge(ID, BLUPs_Cross, by.x = "gid", by.y = "row.names", all.x = TRUE)

# Merge DataOutput with ID_with_predictions based on Cross
MergedData <- merge(DataOutput, ID_with_predictions, by = "Cross", all.x = TRUE)

# Add a new column TRSBLUPSALL with the predicted values
MergedData$PWBLUPSALL <- MergedData$prediction

# Drop unnecessary columns if needed
DataOutput <- MergedData[, c("Cross", "PWBLUPSALL", "gid.x", "TRSBLUPSALL", "TRSBLUPSNI",
                             "TRSBLUPSNR", "TRSBLUPSSG")]
colnames(DataOutput)[colnames(DataOutput) == "gid.x"] <- "gid"



PWBLUPSperloc <- as.data.frame(fitMET2$U$FxL$PlotWeight)
# Assuming TRSBLUPSperloc is your dataframe
# Extract row names and create new columns
PWBLUPSperloc$RowNames <- rownames(PWBLUPSperloc)
PWBLUPSperloc$Location <- sub(".*\\.", "", PWBLUPSperloc$RowNames)
PWBLUPSperloc$Crosssub <- sub("^...|\\..*", "", PWBLUPSperloc$RowNames)
PWBLUPSperloc$Cross <- sub("\\..*", "", PWBLUPSperloc$Crosssub)


# Drop the original row names column if not needed
PWBLUPSperloc <- PWBLUPSperloc[, c("RowNames", "Location", "Cross", colnames(PWBLUPSperloc)[-ncol(PWBLUPSperloc)])]

colnames(PWBLUPSperloc)[colnames(PWBLUPSperloc) == "fitMET2$U$FxL$PlotWeight"] <- "prediction"
# Match values and create new column
DataOutput$PWBLUPSNI <- DataOutput$PWBLUPSALL + PWBLUPSperloc$prediction[PWBLUPSperloc$Cross %in% DataOutput$Cross & PWBLUPSperloc$Location == "NI"]
DataOutput$PWBLUPSNR <- DataOutput$PWBLUPSALL + PWBLUPSperloc$prediction[PWBLUPSperloc$Cross %in% DataOutput$Cross & PWBLUPSperloc$Location == "NR"]
DataOutput$PWBLUPSSG <- DataOutput$PWBLUPSALL + PWBLUPSperloc$prediction[PWBLUPSperloc$Cross %in% DataOutput$Cross & PWBLUPSperloc$Location == "SG"]




#============model visuals for multi-trait selection==========
#modify the BLUPS dataframe
# Calculate mean and standard deviation of BLUPs
blup_mean <- mean(DataOutput$TRSBLUPSALL, na.rm = TRUE)
blup_sd <- sd(DataOutput$TRSBLUPSALL, na.rm = TRUE)
# Create a new column for highlighting
DataOutput$TRSALLDiff <- ifelse(DataOutput$TRSBLUPSALL > (blup_mean + blup_sd), "Greater", ifelse(DataOutput$TRSBLUPSALL < (blup_mean - blup_sd), "Less", "Same"))
# Order the dataframe by BLUPs
DataOutput <- DataOutput[order(DataOutput$TRSBLUPSALL),]



# Create a bar plot with highlighted crosses and their BLUPs
base_fig <- ggplot(DataOutput, aes(x=reorder(Cross, TRSBLUPSALL), y=TRSBLUPSALL, fill=TRSALLDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of TRS Across All Locations", x="Family", y="BLUP",
       fill="Significantly Different than 1 SD") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 10, face = "bold")) +
  coord_cartesian(ylim = c(200, 240))

base_fig



# Specify the columns of interest
columns_of_interest <- c("TRSBLUPSALL", "TRSBLUPSNI", "TRSBLUPSNR", "TRSBLUPSSG", 
                         "PWBLUPSALL", "PWBLUPSNI", "PWBLUPSNR", "PWBLUPSSG")


# Iterate over each column
for (column in columns_of_interest) {
  # Calculate mean and standard deviation of the column
  column_mean <- mean(DataOutput[[column]])
  column_sd <- sd(DataOutput[[column]])
  
  # Create a new column for highlighting
  DataOutput[[paste0(column, "Diff")]] <- ifelse(DataOutput[[column]] > (column_mean + column_sd),
                                                 "Greater",
                                                 ifelse(DataOutput[[column]] < (column_mean - column_sd),
                                                        "Less",
                                                        "Same"))
  
  # Order the dataframe by the column values
  DataOutput <- DataOutput[order(DataOutput[[column]]), ]
  
  # Calculate the ylim range for the plot
  # Calculate the ylim range for the plot
  ylim_min <- min(DataOutput[[column]]) - (min(DataOutput[[column]])*0.10)
  ylim_max <- max(DataOutput[[column]]) + (max(DataOutput[[column]])*0.10)
  
  # Create the bar plot with highlighted values and adjusted ylim
  plot <- ggplot(DataOutput, aes(x = reorder(Cross, !!rlang::sym(column)), y = !!rlang::sym(column), fill = !!rlang::sym(paste0(column, "Diff")))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = column_mean, color = "royalblue3", linetype = "dashed", size = 1) +
    labs(title = paste(column, "BLUPs of Crosses Across All Locations"), x = "Cross", y = "BLUP") +
    scale_fill_manual(values = c("Greater" = "seagreen3", "Less" = "tomato3", "Same" = "grey")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_cartesian(ylim = c(ylim_min, ylim_max))
  
  # Print the plot
  print(plot)
}


#Multi-"Trait" Selection for BLUPs in NI
rank_PWBLUPALL <- rank(DataOutput$PWBLUPSALL)
rank_PWBLUPNI <- rank(DataOutput$PWBLUPSNI)
rank_cutoff_PWBLUPALL <- quantile(rank_PWBLUPALL, 0.85)
rank_cutoff_PWBLUPNI <- quantile(rank_PWBLUPNI, 0.85)
top_15_percent <- DataOutput[(rank_PWBLUPALL >= rank_cutoff_PWBLUPALL) & (rank_PWBLUPNI >= rank_cutoff_PWBLUPNI), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- DataOutput$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_PWBLUPALL >= rank_cutoff_PWBLUPALL & rank_PWBLUPNI >= rank_cutoff_PWBLUPNI, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupall_mean <- mean(DataOutput$PWBLUPSALL)
blupall_sd <- sd(DataOutput$PWBLUPSALL)
blupNI_mean <- mean(DataOutput$PWBLUPSNI)
blupNI_sd <- sd(DataOutput$PWBLUPSNI)

#find high outliers
blupall_outliers <- DataOutput$PWBLUPSALL > (blupall_mean + blupall_sd)
blupNI_outliers <- DataOutput$PWBLUPSNI > (blupNI_mean + blupNI_sd)
blupall_outlier_cross <- DataOutput$Cross[blupall_outliers]
blupNI_outlier_cross <- DataOutput$Cross[blupNI_outliers]
print(blupall_outlier_cross)
print(blupNI_outlier_cross)

point_colors <- ifelse(DataOutput$Cross %in% top_15_percent$Cross, "#1B4F72",
                       ifelse(DataOutput$Cross %in% DataOutput$Cross[blupall_outliers | blupNI_outliers], "#E67E22", "#2E86C1"))

# Scatterplot with colors modified for outliers and a gridded background
ggplot(DataOutput, aes(x = PWBLUPSALL, y = PWBLUPSNI, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Plot weight BLUPs across all locations and within NI",
       x = "BLUPs overall",
       y = "BLUPs in NI") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)





#Multi-"Trait" Selection for BLUPs in NI
rank_PWBLUPALL <- rank(DataOutput$PWBLUPSALL)
rank_PWBLUPNR <- rank(DataOutput$PWBLUPSNR)
rank_cutoff_PWBLUPALL <- quantile(rank_PWBLUPALL, 0.85)
rank_cutoff_PWBLUPNR <- quantile(rank_PWBLUPNR, 0.85)
top_15_percent <- DataOutput[(rank_PWBLUPALL >= rank_cutoff_PWBLUPALL) & (rank_PWBLUPNR >= rank_cutoff_PWBLUPNR), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- DataOutput$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_PWBLUPALL >= rank_cutoff_PWBLUPALL & rank_PWBLUPNR >= rank_cutoff_PWBLUPNR, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupall_mean <- mean(DataOutput$PWBLUPSALL)
blupall_sd <- sd(DataOutput$PWBLUPSALL)
blupNR_mean <- mean(DataOutput$PWBLUPSNR)
blupNR_sd <- sd(DataOutput$PWBLUPSNR)

#find high outliers
blupall_outliers <- DataOutput$PWBLUPSALL > (blupall_mean + blupall_sd)
blupNR_outliers <- DataOutput$PWBLUPSNR > (blupNR_mean + blupNR_sd)
blupall_outlier_cross <- DataOutput$Cross[blupall_outliers]
blupNR_outlier_cross <- DataOutput$Cross[blupNR_outliers]
print(blupall_outlier_cross)
print(blupNR_outlier_cross)

point_colors <- ifelse(DataOutput$Cross %in% top_15_percent$Cross, "#1B4F72",
                       ifelse(DataOutput$Cross %in% DataOutput$Cross[blupall_outliers | blupNR_outliers], "#E67E22", "#2E86C1"))

# Scatterplot with colors modified for outliers and a gridded background
ggplot(DataOutput, aes(x = PWBLUPSALL, y = PWBLUPSNR, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Plot weight BLUPs across all locations and within NR",
       x = "BLUPs overall",
       y = "BLUPs in NR") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)


#Multi-"Trait" Selection for BLUPs in NI
rank_PWBLUPALL <- rank(DataOutput$PWBLUPSALL)
rank_PWBLUPSG <- rank(DataOutput$PWBLUPSSG)
rank_cutoff_PWBLUPALL <- quantile(rank_PWBLUPALL, 0.85)
rank_cutoff_PWBLUPSG <- quantile(rank_PWBLUPSG, 0.85)
top_15_percent <- DataOutput[(rank_PWBLUPALL >= rank_cutoff_PWBLUPALL) & (rank_PWBLUPSG >= rank_cutoff_PWBLUPSG), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- DataOutput$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_PWBLUPALL >= rank_cutoff_PWBLUPALL & rank_PWBLUPSG >= rank_cutoff_PWBLUPSG, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupall_mean <- mean(DataOutput$PWBLUPSALL)
blupall_sd <- sd(DataOutput$PWBLUPSALL)
blupSG_mean <- mean(DataOutput$PWBLUPSSG)
blupSG_sd <- sd(DataOutput$PWBLUPSSG)

#find high outliers
blupall_outliers <- DataOutput$PWBLUPSALL > (blupall_mean + blupall_sd)
blupSG_outliers <- DataOutput$PWBLUPSSG > (blupSG_mean + blupSG_sd)
blupall_outlier_cross <- DataOutput$Cross[blupall_outliers]
blupSG_outlier_cross <- DataOutput$Cross[blupSG_outliers]
print(blupall_outlier_cross)
print(blupSG_outlier_cross)

point_colors <- ifelse(DataOutput$Cross %in% top_15_percent$Cross, "#1B4F72",
                       ifelse(DataOutput$Cross %in% DataOutput$Cross[blupall_outliers | blupSG_outliers], "#E67E22", "#2E86C1"))

# Scatterplot with colors modified for outliers and a gridded background
ggplot(DataOutput, aes(x = PWBLUPSALL, y = PWBLUPSSG, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Plot weight BLUPs across all locations and within SG",
       x = "BLUPs overall",
       y = "BLUPs in SG") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)




# Scatterplot with colors modified for outliers and a gridded background
ggplot(DataOutput, aes(x = PWBLUPSALL, y = PWBLUPSNI, label = Cross)) +
  geom_point(size = 3, alpha = 0.8, color = "#2E86C1") +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  labs(title = "Plot weight BLUPs across all locations and within NI",
       x = "BLUPs overall",
       y = "BLUPs in NI") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "none",  # Remove the legend
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(color = "#2E86C1", size = 3, nudge_x = 0.3, nudge_y = 0.3)

#=================combine BLUPS ALL and NI===============
# Fit the linear regression model
lm_fit <- lm(PWBLUPSALL ~ PWBLUPSNI, data = DataOutput)
plot(lm_fit)
# Calculate the perpendicular distances and the X-coordinates of the feet of the perpendiculars
calculate_perpendicular_distance <- function(x, y, coef) {
  m <- coef[2] # Slope of the regression line
  b <- coef[1] # Intercept of the regression line
  
  x_foot <- (x + m * (y - b)) / (1 + m^2) # X-coordinate of the foot of the perpendicular
  y_foot <- m * x_foot + b # Y-coordinate of the foot of the perpendicular
  
  return(sqrt((x - x_foot)^2 + (y - y_foot)^2)) # Perpendicular distance
}

DataOutput$Perpendicular_Distance <- calculate_perpendicular_distance(DataOutput$PWBLUPSNI, DataOutput$PWBLUPSALL, coef(lm_fit))
DataOutput$X_foot <- (DataOutput$PWBLUPSNI + coef(lm_fit)[2] * (DataOutput$PWBLUPSALL - coef(lm_fit)[1])) / (1 + coef(lm_fit)[2]^2)
DataOutput$Y_foot <- (lm_fit$coefficients[1] + DataOutput$X_foot * lm_fit$coefficients[2])
# Rank the data points based on the X-coordinates of the feet of the perpendiculars
ranked_data <- DataOutput[order(DataOutput$X_foot), ]

# Print the ranked data frame
print(ranked_data)

# Optionally, plot the ranked data points on the scatter plot
plot(PWBLUPSALL ~ PWBLUPSNI, data = DataOutput, main = "BLUPs across all locations and within NI",
     xlab = "BLUPs within NI", ylab = "BLUPs Across All Locations", col = "#2E86C1", pch=16)
points(X_foot ~ Y_foot, data = DataOutput, col = "seagreen3", pch = 16)
abline(lm_fit, col = "black")
# Add line segments connecting each data point to its foot of the perpendicular
for (i in 1:nrow(DataOutput)) {
  segments(
    DataOutput$PWBLUPSNI[i],                # x-coordinate of the point
    DataOutput$PWBLUPSALL[i],               # y-coordinate of the point
    DataOutput$X_foot[i],                   # x-coordinate of the foot of the perpendicular
    lm_fit$coefficients[1] + DataOutput$X_foot[i] * lm_fit$coefficients[2],  # y-coordinate of the foot of the perpendicular
    col = "gray",                           # color of the line
    lty = 1                                # line type (dashed)
  )
}
legend("topleft", legend = c("Original BLUPs", "Ranked BLUPs", "Line of Best Fit", "Perpendicular Distance to Line of Best Fit"),
       col = c("#2E86C1", "seagreen3", "black", "gray"), pch = c(16, 16, NA, NA),
       lty = c(0, 0, 1, 1), # Use 1 for solid line
       bty = "n" # Don't draw a box around the legend
)


#Multi-"Trait" Selection for BLUPs in NI
rank_X_foot <- rank(DataOutput$X_foot)
rank_TRSBLUPSALL <- rank(DataOutput$TRSBLUPSALL)
rank_cutoff_X_foot <- quantile(rank_X_foot, 0.85)
rank_cutoff_TRSBLUPSALL <- quantile(rank_TRSBLUPSALL, 0.85)
top_15_percent <- DataOutput[(rank_X_foot >= rank_cutoff_X_foot) & (rank_TRSBLUPSALL >= rank_cutoff_TRSBLUPSALL), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- DataOutput$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_X_foot >= rank_cutoff_X_foot & rank_TRSBLUPSALL >= rank_cutoff_TRSBLUPSALL, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupPW_mean <- mean(DataOutput$X_foot)
blupPW_sd <- sd(DataOutput$X_foot)
blupTRS_mean <- mean(DataOutput$TRSBLUPSALL)
blupTRS_sd <- sd(DataOutput$TRSBLUPSALL)

#find high outliers
blupPW_outliers <- DataOutput$X_foot > (blupPW_mean + blupPW_sd)
blupTRS_outliers <- DataOutput$TRSBLUPSALL > (blupTRS_mean + blupTRS_sd)
blupPW_outlier_cross <- DataOutput$Cross[blupPW_outliers]
blupTRS_outlier_cross <- DataOutput$Cross[blupTRS_outliers]
print(blupPW_outlier_cross)
print(blupTRS_outlier_cross)

point_colors <- ifelse(DataOutput$Cross %in% top_15_percent$Cross, "seagreen3",
                       ifelse(DataOutput$Cross %in% DataOutput$Cross[blupTRS_outliers | blupPW_outliers], "#E67E22", "#2E86C1"))

# Scatterplot with colors modified for outliers and a gridded background
ggplot(DataOutput, aes(x = X_foot, y = TRSBLUPSALL, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Multi-trait selection for plot weight and TRS",
       x = "Plot weight ranked BLUPs",
       y = "TRS BLUPs over all locations") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)




point_colors <- case_when(
  DataOutput$Cross %in% top_15_percent$Cross ~ "seagreen3",
  DataOutput$Cross %in% DataOutput$Cross[blupTRS_outliers | blupPW_outliers] ~ "#E67E22",
  TRUE ~ "#2E86C1"
)

# Scatterplot with colors modified for outliers and a gridded background
ggplot(DataOutput, aes(x = X_foot, y = TRSBLUPSALL, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(
    values = c("seagreen3" = "seagreen3", "#2E86C1" = "#2E86C1", "#E67E22" = "#E67E22"),
    name = "Quantile:",
    labels = c("Bottom 15%", "High Outliers (one trait)", "Top 15% (both traits)")
  ) +
  labs(title = "Multi-trait selection for plot weight and TRS",
       x = "Plot weight ranked BLUPs",
       y = "TRS BLUPs over all locations") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)




merged_data <- merge(TRSBLUPS, ranked_data[c("Cross", "X_foot")], by = "Cross", suffixes = c("", "_ranked"))
colnames(merged_data)[colnames(merged_data) == "X_foot_ranked"] <- "ALLSG"
print(merged_data)
columns_to_remove <- c("Perpendicular_Distance", "X_on_Line", "X_foot")
merged_data <- merged_data[, !(names(merged_data) %in% columns_to_remove)]
TRSBLUPS <- merged_data