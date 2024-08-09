setwd('/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/results/downstream_fasting_NGT')
getwd()

# Load the dplyr package for easier data manipulation 
library(dplyr)

# Define the file path
# Change filepath and the file as needed
input_file_path <- "CytoscapeInput-edges-yellow.txt"
output_file_path <- "CytoscapeInput-edges-yellow_cutoff_0.2.txt"

# Read the file into a data frame
data <- read.table(input_file_path, header = TRUE, sep = "\t")

# Filter the rows where weight is greater than or equal to 0.3
filtered_data <- data %>%
  filter(weight >= 0.2) #only 14 with cutoff of 0.3, at cutoff 0.2 616 pairs

# Save the filtered data to a new file
write.table(filtered_data, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Print a message indicating the file has been saved
cat("Filtered data saved to", output_file_path)

