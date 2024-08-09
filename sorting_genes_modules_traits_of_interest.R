setwd("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/results/downstream_fasting_NGT/")

getwd()

library(WGCNA)
library(readr)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
allowWGCNAThreads()

library(openxlsx)

gc() #to free unused memory in R

#-----reading data for fasting glucose-------------------------------------------------------------------------
data_fasting_ngt_gs_mm_fasting_glucose=read.csv('geneInfo_fasting_NGT.csv',row.names=1)

# Create separate data frames based on moduleColor column
blue_module_fasting_glucose <- subset(data_fasting_ngt_gs_mm_fasting_glucose, moduleColor == 'blue')
purple_module_fasting_glucose <- subset(data_fasting_ngt_gs_mm_fasting_glucose, moduleColor == 'purple')
yellow_module_fasting_glucose <- subset(data_fasting_ngt_gs_mm_fasting_glucose, moduleColor == 'yellow')

# Check the first few rows of each new data frame
head(blue_module_fasting_glucose)
head(purple_module_fasting_glucose)
head(yellow_module_fasting_glucose)

# # Filter blue_module based on the given conditions
# # filtered_blue_module_fasting_glucose <- subset(blue_module_fasting_glucose, 
# abs(GS.Fasting.Glucose) >= 0.5 & 
#   p.GS.Fasting.Glucose < 0.05 &
#   abs(MM.blue) >= 0.5 &
#   p.MM.blue < 0.05)

# Filter blue_module based on the given conditions
filtered_blue_module_fasting_glucose <- subset(blue_module_fasting_glucose, 
                               abs(GS.Fasting.Glucose) >= 0.5 & 
                               abs(MM.blue) >= 0.5)

filtered_blue_module_fasting_glucose  <- filtered_blue_module_fasting_glucose %>%
  select(1:9, MM.blue, p.MM.blue)

# Filter purple_module based on the given conditions
filtered_purple_module_fasting_glucose <- subset(purple_module_fasting_glucose, 
                                 abs(GS.Fasting.Glucose) >= 0.5 & 
                                  abs(MM.purple) >= 0.5)

filtered_purple_module_fasting_glucose  <- filtered_purple_module_fasting_glucose %>%
  select(1:9, MM.purple, p.MM.purple)

# Filter yellow_module based on the given conditions
filtered_yellow_module_fasting_glucose <- subset(yellow_module_fasting_glucose, 
                                 abs(GS.Fasting.Glucose) >= 0.5 & 
                                   abs(MM.yellow) >= 0.5)

filtered_yellow_module_fasting_glucose  <- filtered_yellow_module_fasting_glucose %>%
  select(1:9, MM.yellow, p.MM.yellow)

# Check the first few rows of each filtered data frame
head(filtered_blue_module_fasting_glucose)
head(filtered_purple_module_fasting_glucose)
head(filtered_yellow_module_fasting_glucose)




#-----reading data for fasting insulin-------------------------------------------------------------------------
data_fasting_ngt_gs_mm_fasting_insulin=read.csv('geneInfo_fasting_insulin_NGT.csv',row.names=1)#gene significance with fasting insulin

# Create separate data frames based on moduleColor column
blue_module_fasting_insulin <- subset(data_fasting_ngt_gs_mm_fasting_insulin, moduleColor == 'blue')
purple_module_fasting_insulin <- subset(data_fasting_ngt_gs_mm_fasting_insulin, moduleColor == 'purple')
yellow_module_fasting_insulin <- subset(data_fasting_ngt_gs_mm_fasting_insulin, moduleColor == 'yellow')

# Check the first few rows of each new data frame
head(blue_module_fasting_insulin)
head(purple_module_fasting_insulin)
head(yellow_module_fasting_insulin)

# # Filter blue_module based on the given conditions
# filtered_blue_module <- subset(blue_module, 
#                                GS.Fasting.Insulin >= 0.5 & 
#                                  p.GS.Fasting.Insulin < 0.05 & 
#                                  MM.blue >= 0.5 & 
#                                  p.MM.blue < 0.05)

# Filter blue_module based on the given conditions
filtered_blue_module_fasting_insulin <- subset(blue_module_fasting_insulin, 
                               abs(GS.Fasting.Insulin) >= 0.5 & 
                                 abs(MM.blue) >= 0.5)
filtered_blue_module_fasting_insulin  <- filtered_blue_module_fasting_insulin %>%
  select(1:9, MM.blue, p.MM.blue)

# Filter purple_module based on the given conditions
filtered_purple_module_fasting_insulin <- subset(purple_module_fasting_insulin, 
                                 abs(GS.Fasting.Insulin) >= 0.5 & 
                                   abs(MM.purple) >= 0.5)

filtered_purple_module_fasting_insulin  <- filtered_purple_module_fasting_insulin %>%
  select(1:9, MM.purple, p.MM.purple)

# Filter yellow_module based on the given conditions
filtered_yellow_module_fasting_insulin <- subset(yellow_module_fasting_insulin, 
                                abs (GS.Fasting.Insulin) >= 0.5 & 
                                   abs (MM.yellow) >= 0.5)

filtered_yellow_module_fasting_insulin  <- filtered_yellow_module_fasting_insulin %>%
  select(1:9, MM.yellow, p.MM.yellow)

# Check the first few rows of each filtered data frame
head(filtered_blue_module_fasting_insulin)
head(filtered_purple_module_fasting_insulin)
head(filtered_yellow_module_fasting_insulin)



#-----reading data for HOMA IR-------------------------------------------------------------------------
data_fasting_ngt_gs_mm_fasting_homa_ir=read.csv('geneInfo_HOMA_IR_NGT.csv',row.names=1)#gene significance with fasting insulin

# Create separate data frames based on moduleColor column
blue_module_fasting_homa_ir <- subset(data_fasting_ngt_gs_mm_fasting_homa_ir, moduleColor == 'blue')
purple_module_fasting_homa_ir <- subset(data_fasting_ngt_gs_mm_fasting_homa_ir, moduleColor == 'purple')
yellow_module_fasting_homa_ir <- subset(data_fasting_ngt_gs_mm_fasting_homa_ir, moduleColor == 'yellow')

# Check the first few rows of each new data frame
head(blue_module_fasting_homa_ir)
head(purple_module_fasting_homa_ir)
head(yellow_module_fasting_homa_ir)

# # Filter blue_module based on the given conditions
# filtered_blue_module <- subset(blue_module, 
#                                GS.Fasting.Insulin >= 0.5 & 
#                                  p.GS.Fasting.Insulin < 0.05 & 
#                                  MM.blue >= 0.5 & 
#                                  p.MM.blue < 0.05)

# Filter blue_module based on the given conditions
filtered_blue_module_fasting_homa_ir <- subset(blue_module_fasting_homa_ir, 
                                             abs(GS.HOMA.IR) >= 0.5 & 
                                                 abs(MM.blue) >= 0.5)

filtered_blue_module_fasting_homa_ir  <- filtered_blue_module_fasting_homa_ir %>%
  select(1:9, MM.blue, p.MM.blue)

# Filter purple_module based on the given conditions
filtered_purple_module_fasting_homa_ir <- subset(purple_module_fasting_homa_ir, 
                                                 abs(GS.HOMA.IR) >= 0.5 & 
                                                   abs(MM.purple) >= 0.5)

filtered_purple_module_fasting_homa_ir  <- filtered_purple_module_fasting_homa_ir %>%
  select(1:9, MM.purple, p.MM.purple)

# Filter yellow_module based on the given conditions
filtered_yellow_module_fasting_homa_ir <- subset(yellow_module_fasting_homa_ir, 
                                                 abs (GS.HOMA.IR) >= 0.5 & 
                                                   abs (MM.yellow) >= 0.5)

filtered_yellow_module_fasting_homa_ir  <- filtered_yellow_module_fasting_homa_ir %>%
  select(1:9, MM.yellow, p.MM.yellow)

# Check the first few rows of each filtered data frame
head(filtered_blue_module_fasting_homa_ir)
head(filtered_purple_module_fasting_homa_ir)
head(filtered_yellow_module_fasting_homa_ir)


# Create a new workbook
wb <- createWorkbook()

# Add sheets to the workbook
addWorksheet(wb, "Filtered Blue Fasting Glucose")
addWorksheet(wb, "Filtered Blue Fasting Insulin")
addWorksheet(wb, "Filtered Blue HOMA-IR")

# Write the data to the sheets
writeData(wb, sheet = "Filtered Blue Fasting Glucose", filtered_blue_module_fasting_glucose)
writeData(wb, sheet = "Filtered Blue Fasting Insulin", filtered_blue_module_fasting_insulin)
writeData(wb, sheet = "Filtered Blue HOMA-IR", filtered_blue_module_fasting_homa_ir)

# Save the workbook
saveWorkbook(wb, file = "Filtered_Blue_Modules.xlsx", overwrite = TRUE)



# Create a new workbook
wb <- createWorkbook()

# Add sheets to the workbook
addWorksheet(wb, "Filtered Purple Fasting Glucose")
addWorksheet(wb, "Filtered Purple Fasting Insulin")
addWorksheet(wb, "Filtered Purple HOMA-IR")

# Write the data to the sheets
writeData(wb, sheet = "Filtered Purple Fasting Glucose", filtered_purple_module_fasting_glucose)
writeData(wb, sheet = "Filtered Purple Fasting Insulin", filtered_purple_module_fasting_insulin)
writeData(wb, sheet = "Filtered Purple HOMA-IR", filtered_purple_module_fasting_homa_ir)

# Save the workbook
saveWorkbook(wb, file = "Filtered_Purple_Modules.xlsx", overwrite = TRUE)




# Create a new workbook
wb <- createWorkbook()

# Add sheets to the workbook
addWorksheet(wb, "Filtered Yellow Fasting Glucose")
addWorksheet(wb, "Filtered Yellow Fasting Insulin")
addWorksheet(wb, "Filtered Yellow HOMA-IR")

# Write the data to the sheets
writeData(wb, sheet = "Filtered Yellow Fasting Glucose", filtered_yellow_module_fasting_glucose)
writeData(wb, sheet = "Filtered Yellow Fasting Insulin", filtered_yellow_module_fasting_insulin)
writeData(wb, sheet = "Filtered Yellow HOMA-IR", filtered_yellow_module_fasting_homa_ir)

# Save the workbook
saveWorkbook(wb, file = "Filtered_Yellow_Modules.xlsx", overwrite = TRUE)



#change the variable that goes into the df1, df2, and df3 variable

# Assuming your data frames are named as mentioned
df1 <- filtered_yellow_module_fasting_glucose
df2 <- filtered_yellow_module_fasting_insulin
df3 <- filtered_yellow_module_fasting_homa_ir

# Select only the desired columns including 'geneSymbol'
df1_selected <- df1 %>% select(,c(1:9))
df2_selected <- df2 %>% select(,c(2,8,9))
df3_selected <- df3 %>% select(,c(2,8:11))

# Perform the merge based on the common 'geneSymbol' column using dplyr
merged_df <- df1_selected %>%
  inner_join(df2_selected, by = "geneSymbol") %>%
  inner_join(df3_selected, by = "geneSymbol")

write.csv(merged_df,file='fasting_NGT_yellow_module_gs_mm_fasting_glu_insu_homa_ir.csv')

# View the merged dataframe
print(merged_df)



# Save the current workspace
save.image(file = "sorting_genes_modules_traits_of_interest.RData")

# Load the saved workspace
load("sorting_genes_modules_traits_of_interest.RData")
