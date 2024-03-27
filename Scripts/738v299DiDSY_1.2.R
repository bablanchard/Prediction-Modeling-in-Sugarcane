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


# Filter data to include only observations of specified varieties
filtered_varieties <- df %>%
  filter(VARIETY %in% c(2017738, 2001299) &
           YEAR %in% 2021:2023 &
           STAGE == "OUTFIELD" &
           !(LOC %in% "ARDOYNE") &
           CROP %in% c(0, 1, 2))

# Group by specified variables and filter groups that contain both varieties
filtered_combinations <- filtered_varieties %>%
  group_by(YEAR, STAGE, LOC, CROP) %>%
  filter(all(c(2017738, 2001299) %in% VARIETY)) %>%
  distinct(YEAR, STAGE, LOC, CROP)

# Subset the original dataframe based on the filtered combinations
table738 <- filtered_varieties %>%
  semi_join(filtered_combinations, by = c("YEAR", "STAGE", "LOC", "CROP"))

# Assuming 'table730' is the name of your dataframe
table738 <- table738 %>%
  mutate(CROP = case_when(
    CROP == 0 ~ 1,
    CROP == 1 ~ 2,
    CROP == 2 ~ 3,
    TRUE ~ as.numeric(CROP)  # Keep other values unchanged
  ))
table738 <- table738[complete.cases(table738$T_SPACRE), ]

# Replace 'T_SPACRE' with the actual column name if different
ggplot(table738, aes(x = VARIETY, y = T_SPACRE)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, width = 0.3) +  # Add jittered points
  labs(x = "VARIETY", y = "T_SPACRE") +
  ggtitle("Boxplot of T_SPACRE by VARIETY")

table738$T_SPACRE <- as.numeric(table738$T_SPACRE)
table738$CROP <- as.factor(table738$CROP)
table738$VARIETY <- as.factor(table738$VARIETY)

# Fit the difference-in-differences model with an intercept
model <- lm(T_SPACRE ~ VARIETY * CROP, data = table738)

# Extract model summary
model_summary <- summary(model)
coefficients <- coef(model)

# Extract model coefficients and tidy the output
tidy_output <- tidy(model_summary)
# Print the model coefficients and other statistics
print(tidy_output)
output <- as.data.frame(tidy_output)
# Convert CROP to a factor if it's a numeric variable
table738$CROP <- as.factor(table738$CROP)

# Create the time series plot
ggplot(table738, aes(x = CROP, y = T_SPACRE, group = VARIETY, color = as.factor(VARIETY))) +
  geom_line() +
  labs(x = "Time Points (CROP)", y = "T_SPACRE", color = "VARIETY") +
  theme_minimal()


average_T_SPACREM_2001299 <- mean(table738[table738$VARIETY == "2001299", "T_SPACRE"], na.rm = TRUE)
average_T_SPACREM_2017738 <- mean(table738[table738$VARIETY == "2017738", "T_SPACRE"], na.rm = TRUE)



estimates <- as.data.frame(tidy_output)
# Extract the coefficient for the variety main effect
variety_coef <- estimates$estimate[estimates$term == "VARIETY2017738"]
# Extract the coefficient for the variety interaction effects
varietycrop2_coefs <- estimates$estimate[estimates$term == "VARIETY2017738:CROP2"]
varietycrop3_coefs <- estimates$estimate[estimates$term == "VARIETY2017738:CROP3"]

# Sum all the VARIETY coefficients
sum_variety_coef <- sum(variety_coef, varietycrop2_coefs, varietycrop3_coefs, varietycrop4_coefs)

# Divide by the average T_SPACREM for VARIETY2001299
percentage_change <- (sum_variety_coef / average_T_SPACREM_2001299) * 100

# Print the percentage change
print(percentage_change)


