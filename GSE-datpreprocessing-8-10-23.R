#setwd('C:/Users/ADMIN/Desktop/GSE-datpreprocessing-8-10-23') #HOME Workstation directory
setwd('/home/admin/R_data/GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23')

# Loading all workspace env data objects from GSE-datpreprocessing-8-10-23.RData
load(file ="/home/admin/R_data/GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23/GSE-datpreprocessing-8-10-23.RData") 

# Load required libraries-----
library(Biobase)
library(GEOquery)
library(limma)
library(tibble)
library(dplyr)
library(ggplot2)
library(WGCNA)



#------------Downloading processed data from GEO using GEOquery package-----
# Set the GEO series accession number
geo_accession <- "GSE153837"

# Download the data from GEO
gse <- getGEO(geo_accession, GSEMatrix = TRUE)

# Download the data from GEO supplementary files
#supplementary_gse <- getGEOSuppFiles("GSE153837")

# Extract the expression data
gsm_data <- exprs(gse[[1]])

# Access the sample IDs from the metadata (replace 'metadata_column' with the actual column name)
sample_ids <- phenoData(gse[[1]])$title

# Convert the expression data to a data frame
exprs_df <- as.data.frame(gsm_data)

#replace GSM IDs with Sample Names
colnames(exprs_df) <- sample_ids




#--------raw data shared by collaborator with p-val which was used further for filtering with p-val <0.05-----
raw_data <- read.csv('raw_data_GSE153837.csv',row.names = 1)

# Create a condition to filter rows where all p-values are less than 0.05
condition <- rowSums(raw_data[, grepl("PVALUE", names(raw_data))] < 0.05) == ncol(raw_data[, grepl("PVALUE", names(raw_data))])

# Filter rows based on the condition
filtered_df <- raw_data[condition, ]
filtered_df_pvalue <-  rownames_to_column(filtered_df, var = "Probe_ID")



#------ Data frame with filtered protein coding genes with NM tag and ensemble gene ID------
#filteredespr in file 'preprocessing-split-data-28-3-23' @ getwd('/home/admin/R_data/GSE153837_normalized_22-10-2023')
filteredexpr <- read.csv('filtered.exprData.csv')
protein_coding_filtered <-  as.data.frame(filteredexpr[,2])

# Change the name of the 'filteredexpr[,2]' column to 'ILMN_ID_PCG'
colnames(protein_coding_filtered)[1] <- "ILMN_ID_PCG"



#------ merging cols of protein coding genes with p-val < 0.05------
merged_filtered_pval_pcg <- merge(filtered_df_pvalue, protein_coding_filtered, by.x = 'Probe_ID', by.y = "ILMN_ID_PCG")
write.csv(merged_filtered_pval_pcg, "C:/Users/ADMIN/Desktop/GSE-datpreprocessing-8-10-23/merged_filtered_pval_pcg.csv", row.names=FALSE)

# Illumina probe list for protein coding genes with pval < 0.05 column name is 'PROBE_ID'
merged_filtered_pval_pcg_ILMN_list <- merged_filtered_pval_pcg[1] 



#-----merging log2normalized data with protein coding genes and p value <0.05----
#GSM expression data with first column as ILMN_ID
exprs_df_ILMN_ID <-  rownames_to_column(exprs_df, var = "ILMN_ID")

#merging log2 normalized data with probe list of protein coding genes and pval < 0.05
merge_log2norm_data_pcg_pval <-  merge(merged_filtered_pval_pcg_ILMN_list, exprs_df_ILMN_ID, by.x= 'Probe_ID', by.y='ILMN_ID')

#list of illuminaIDs for log2 normalized values for protein coding genes at cutoff of pval <0.05
#merge_log2norm_data_pcg_pval_ILMN_list <-  merge_log2norm_data_pcg_pval[1]

# Assuming your dataframe is named 'df'
# Move the first column to rownames
rownames(merge_log2norm_data_pcg_pval) <- merge_log2norm_data_pcg_pval[, 1]

# Remove the first column from the dataframe
merge_log2norm_data_pcg_pval <- merge_log2norm_data_pcg_pval[, -1]

#saving csv file for 10647 log2normalized protein coding genes with pval<0.05 
write.csv(merge_log2norm_data_pcg_pval, 'C:/Users/ADMIN/Desktop/GSE-datpreprocessing-8-10-23/merge_log2norm_data_pcg_pval.csv',row.names=TRUE)

#finding missing data point: there was 1 probe for which there was no value in LOG2normalized data from GEO 
# Assuming you have merged_filtered_pval_pcg_ILMN_list and exprs_df_ILMN_ID dataframes
missing_data_point <- anti_join(merged_filtered_pval_pcg_ILMN_list, merge_log2norm_data_pcg_pval, by = c("Probe_ID" = "Probe_ID"))

# 'missing_data_point' will contain the rows from merged_filtered_pval_pcg_ILMN_list that have no matching 'ILMN_ID' in exprs_df_ILMN_ID



#------plotting 10647 log2normalized protein coding genes with p-val<0.05 data-----

#clean and tidy the dataset for exploratory graphics
col_sel=names(merge_log2norm_data_pcg_pval)

mdata <- merge_log2norm_data_pcg_pval %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )

# ==== Plot groups (Sample Groups vs avg signal intensity) to identify outliers
(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "Avg Signal Intensity") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour
)



#-----reading annotation data for gene symbol and remove the alternative splice variants of probes----
annodata <- as.data.frame(read.csv("Annotation.csv",header = T))
annodata <- annodata[,c(2,3)]

#merging log2norm pcg with pval<0.05 with annodata (GENE SYMBOL and ILMN_ID)
merge_alternative_splice_log2norm_data_pcg_pval <-  merge(merge_log2norm_data_pcg_pval, annodata, by.x= 'Probe_ID', by.y='ILMN_ID')

# Remove duplicates and store unique values in another dataframe
unique_merge_removed_alternative_splice_log2norm_data_pcg_pval <- merge_alternative_splice_log2norm_data_pcg_pval %>%
  distinct(SYMBOL, .keep_all = TRUE) #8596 probes unique

# Move the first column to rownames
rownames(unique_merge_removed_alternative_splice_log2norm_data_pcg_pval) <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval[, 1]

# Remove the first column and gene symbol column from the dataframe
unique_merge_removed_alternative_splice_log2norm_data_pcg_pval <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval[, -c(1,30)]

#saving csv file for 8596 with no alternatively spliced duplicates, log2normalized protein coding genes with pval<0.05 
write.csv(unique_merge_removed_alternative_splice_log2norm_data_pcg_pval, 'C:/Users/ADMIN/Desktop/GSE-datpreprocessing-8-10-23/unique_merge_removed_alternative_splice_log2norm_data_pcg_pval.csv',row.names=TRUE)

# Find and store non-unique (duplicate) values in another dataframe# Find duplicate values in the specified column
duplicates <- merge_removed_alternative_splice_log2norm_data_pcg_pval[duplicated(merge_removed_alternative_splice_log2norm_data_pcg_pval$SYMBOL) | duplicated(merge_removed_alternative_splice_log2norm_data_pcg_pval$SYMBOL, fromLast = TRUE), ]
#'duplicates' contains the rows with duplicate values from the 'column_name' column



#-----plotting: violin, clustering, PCA 8596 with no alternatively spliced duplicates, log2normalized protein coding genes with pval<0.05 data-----

#clean and tidy the dataset for exploratory graphics
col_sel=names(unique_merge_removed_alternative_splice_log2norm_data_pcg_pval)

mdata <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )


# ==== Plot groups (Sample Groups vs avg signal intensity) to identify outliers
# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23/exploratory_plots/violin_scatter_plot_preprocessed.svg")

(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90,size=12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12),          # Rotate treatment text
      axis.title = element_text(size = 12)  # Increase font size for axis titles
    ) + #change font size by changing the value of size= 112, 14, etc.
    labs(x = "Treatment Groups", y = "Avg Signal Intensity") #+
    #facet_grid(cols = vars(group), drop = FALSE, scales = "free_x")      # Facet by hour
)

# Close the graphics device
dev.off() 


# ==== samples plot using WGCNA
datExpr0 <- as.data.frame(t(unique_merge_removed_alternative_splice_log2norm_data_pcg_pval))
sampleTree=hclust(dist(datExpr0),method = "average")

## Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# change the dimensions if the window is too large or too small.
sizeGrWindow(12,10)

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);

# Set the width and height for the SVG output (in inches)
#width_in_inches <- 10
#height_in_inches <- 6
# SVG graphics device
#svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23/exploratory_plots/all-sample-clustering.svg")

par(cex=0.6)
par(mar=c(9,9,9,9))
plot(sampleTree, main = "Sample clustering"
     ,sub="Outlier Detection",
     xlab="Samples",
     cex.lab=1.5,
     cex.axis=1,
     cex.main=2)

# Close the graphics device
dev.off() 

#creating a new dataset for making a combined plot of expression data and trait data
data_for_traits= read.csv('GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23/unique_merge_removed_alternative_splice_log2norm_data_pcg_pval.csv')
rownames(data_for_traits) <- data_for_traits[,1]
data_for_traits <- data_for_traits[,-c(1)]
data_for_traits=as.data.frame(t(data_for_traits))

#determining sample clusters
clust=data_for_traits
#table(clust)

#keeping clust 1
datExpr=(clust)
#datExpr=datExpr0[keepSamples,]
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

#----trait data-----------------------------------------------------------
#playing with Trait data
allTraits=read.csv("/home/admin/R_data/GSE153837_normalized_22-10-2023/traitdata.csv")
#removing unwanted columns
# Form a data frame analogous to expression data that will hold the clinical traits.
sample_ID=rownames(datExpr)
traitRows=match(sample_ID,allTraits$samples)
datTraits=allTraits[traitRows,]
rownames(datTraits)=datTraits[,1]
datTraits=datTraits[,-c(1:3)]


collectGarbage()


#------sample dendrogram and the colors underneath----------------------------
#reclustering samples
sampleTree2=hclust(dist(datExpr),method = 'average')
traitColors=numbers2colors(datTraits,signed = F)


# Plot the sample dendrogram and the colors underneath.

#png(file="Sample_dendrogram_and_trait_heatmap.png",
#width=1000, height=750)
# Set the width and height for the SVG output (in inches)
width_in_inches <- 12
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23/exploratory_plots/sample_clustering_and_heatmap_all_sample.svg")


plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),marAll = c(1, 12, 3, 1),# Bottom, left, top, right margins
                    cex.colorLabels = 1.5,
                    cex.dendroLabels = 1.5,
                    cex.lab=1.5,
                    cex.axis=1.5)



dev.off()

                    

# ==== sample PCA
pca <- prcomp(datExpr0)
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)


# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/GSE153837-datapreprocessing-8-10-23/exploratory_plots/all-PCA.svg",
    width = width_in_inches, height = height_in_inches)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat),alpha=0.6,size=3) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# Close the graphics device
dev.off() 




#-----separating data files ----
#separating data files for 8596 (preprocessed data) probes across different phenotypes: NGT fast, PD fast, NGT 2 hr, PD 2 hr (refer to preprocessing-split-data-28-3-23)

#split_data-directory: 
#/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data

#folder names and description (no. denote S.no. from analyses ppt slide 'Meeting: 19-5-23'): 


#fasting_NGT_PD  : fasting samples NGT and PD (1)
#glucose_load_NGT_PD  : 2 hour post glucose load NGT  and PD (2)
#PD_fasting_glucose_load : PD ONLY fasting and 2 hour post glucose load (3)
#NGT_fasting_glucose_load  : NGT ONLY fasting and 2 hour post glucose load (4)
#fasting_NGT  : fasting samples NGT ONLY (5)
#glucose_load_NGT  : 2 hour post glucose load NGT ONLY (6)
#fasting_PD  : fasting samples PD ONLY (7)
#glucose_load_PD  : 2 hour post glucose load PD ONLY (8)



# using datasets (5) and (7) do (1) analysis
 
# using datasets (7) and (8) do (3) analysis

# using datasets (5) and (6) do (4) analysis


#fasting_NGT_PD  : fasting samples NGT and PD (1)
fasting_NGT_PD <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("Fast"))
write.csv(fasting_NGT_PD,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT_PD/fasting_NGT_PD.csv",row.names=TRUE)

#glucose_load_NGT_PD  : 2 hour post glucose load NGT and PD (2)
glucose_load_NGT_PD <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("hour"))
write.csv(glucose_load_NGT_PD,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/glucose_load_NGT_PD/glucose_load_NGT_PD.csv",row.names=TRUE)

#PD_fasting_glucose_load : PD ONLY fasting and 2 hour post glucose load (3)
PD_fasting_glucose_load <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("PD"))
write.csv(PD_fasting_glucose_load,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/PD_fasting_glucose_load/PD_fasting_glucose_load.csv",row.names=TRUE)

#NGT_fasting_glucose_load  : NGT ONLY fasting and 2 hour post glucose load (4)
NGT_fasting_glucose_load <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("NGT"))
write.csv(NGT_fasting_glucose_load,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/NGT_fasting_glucose_load/NGT_fasting_glucose_load.csv",row.names=TRUE)

#fasting_NGT  : fasting samples NGT ONLY (5)
fasting_NGT <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("NGT_Fast"))
write.csv(fasting_NGT,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/fasting_NGT.csv")

#glucose_load_NGT  : 2 hour post glucose load NGT ONLY (6)
glucose_load_NGT <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("NGT_2"))
write.csv(glucose_load_NGT,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/glucose_load_NGT/glucose_load_NGT.csv")

#fasting_PD  : fasting samples PD ONLY (7)
fasting_PD <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("PD_Fast"))
write.csv(fasting_PD,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_PD/fasting_PD.csv")

#glucose_load_PD  : 2 hour post glucose load PD ONLY (8)
glucose_load_PD <- unique_merge_removed_alternative_splice_log2norm_data_pcg_pval %>% select(contains("PD_2"))
write.csv(glucose_load_PD,file="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/glucose_load_PD/glucose_load_PD.csv")



# ----Save all workspace env data objects in GSE-datpreprocessing-8-10-23.RData ----
