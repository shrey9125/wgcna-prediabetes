setwd("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT")

getwd()

library(WGCNA)
library(readr)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
allowWGCNAThreads()

# ---- All results are saved in results folder

# ----Save and load all workspace env data objects .RData file----
#save.image(file ="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/fasting_NGT.RData")


# Loading all workspace env data objects from GSE-datpreprocessing-8-10-23.RData
load(file ="/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/fasting_NGT.RData")

#imp follow this link for more info:
#https://statisticsglobe.com/stringsasfactors-argument-of-data-frame-function-in-r
#https://simplystatistics.org/posts/2015-07-24-stringsasfactors-an-unauthorized-biography/#:~:text=The%20argument%20'stringsAsFactors'%20is%20an,or%20as%20just%20plain%20strings.
options(stringsAsFactors = F)

#-----reading data-------------------------------------------------------------------------
Data=read.csv('fasting_NGT.csv',row.names=1)
 
# rownames(Data) <- Data [,1] Use if needed check df first
# Data <- Data[,-c(1)]
datExpr0=as.data.frame(t(Data))


#----Violin plot of the data----
violin_plot <- Data
col_sel = names(violin_plot)     # Get all but first column name
mdata <- violin_plot %>%
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
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/violin_scatter_fasting_NGT.svg")

(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "Avg Signal Intensity")# +
    #facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour
)

# Close the file
dev.off()

#----good sample genes check--------------------------------------------------------------------------
gsg=goodSamplesGenes(datExpr0,verbose=10)#2460 NA values
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples
# from the data:
#   if (!gsg$allOK)
#   {
#     # Optionally, print the gene and sample names that were removed:
#     if (sum(!gsg$goodGenes)>0)
#       printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#     if (sum(!gsg$goodSamples)>0)
#       printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
#     # Remove the offending genes and samples from the data:
#     datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
#   }

#-------samples plot-----------------------------------------------------------------------
#samples plot
sampleTree=hclust(dist(datExpr0),method = "average")

## Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,10)

# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/sample_clustering_fasting_NGT.svg")

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);

par(cex=0.6)
par(mar=c(1,5,3,1))

plot(sampleTree, main = "Sample clustering"
     ,sub="Outlier Detection",
     xlab="Samples",
     cex.lab=1.5,
     cex.axis=1,
     cex.main=2)

# Close the file
dev.off()

#dealing with outlier samples
#plot a cutoff line
#abline(h=15,col='green')

#determining sample clusters
clust=datExpr0
#table(clust)

#keeping clust 1
datExpr=(clust)
#datExpr=datExpr0[keepSamples,]
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)



#------sample PCA----
pca <- prcomp(datExpr)
pca.dat <- pca$x


pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/PCA_fasting_NGT.svg")

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# 3. Close the file
dev.off()

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
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/sample_clustering_and_heatmap_fasting_NGT.svg")


plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),marAll = c(1, 9, 3, 1),# Bottom, left, top, right margins
                    main = "Sample dendrogram and trait heatmap")



dev.off()
#---------------set of powers chosen--------------------------------

#set of powers chosen
powers=c(c(1:10),seq(from=12,to=60,by=2))
sft=pickSoftThreshold(datExpr,powerVector = powers,verbose = 10)

# Plot the results:
# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/pick_sft_fasting_NGT.svg")

par(mfrow = c(2,1),mar = c(5, 6, 4, 4))




cex1 = 1.5

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = 'Soft threshold (power)',
     ylab='Scale Free Topology Model Fit, signed R^2',type='n',
     main=paste('Scale Independence'),
     cex.lab = 1.5,  # Increase font size for axis labels
     cex.main = 1.5)  # Increase font size for main title

text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col='blue')

# this line corresponds to using an R^2 cut-off of h
#abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),
     cex.lab = 1.5,  # Increase font size for axis labels
     cex.main = 1.5)  # Increase font size for main title
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

dev.off()

#----Gyaan: Picking Power (refer to WGCNA FAQ)----

# no. of samples      unsigned          signed
# Less than 20 	        9 	              18
# 20-30               	8 	              16
# 30-40 	              7 	              14
# more than 40 	        6 	              12

#picked power: 18 (signed network i.e., negative correlations set to zero and positives raised to the soft threshold power)

#----WGCNA Run---------------------------------------------
picked_power = 18
temp_cor <- cor


# Force it to use WGCNA cor function
cor <- WGCNA::cor   


#blockwisemodule autoMATIC NETWORK CONSTRUCTION

net = blockwiseModules(datExpr,  # == Adjacency Function ==
                       power = picked_power,                # <= power here
                       networkType = "signed",
                       
                       # == Tree and Block Options ==
                       deepSplit = 2,
                       pamRespectsDendro = F,
                       #detectCutHeight = 0.75,
                       minModuleSize = 100,#can try to change minModuleSize for better distribution
                       maxBlockSize = 15000,
                       
                       # == Module Adjustments ==
                       reassignThreshold = 0,
                       mergeCutHeight = 0.5,#change the distance at which closely relate modules should merge
                       
                       # == TOM == Archive the run results in TOM file (saves time)
                       saveTOMs = T,
                       saveTOMFileBase = "TOM_data",
                       
                       # == Output Options
                       numericLabels = T,
                       verbose = 3)
cor <- temp_cor     # Return cor function to original n

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/cluster_dendrogram_fasting_NGT.svg")


# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",main='Cluster Dendrogram for NGT Fasting State',
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.main=1.2,
                    cex.lab=1.2)

dev.off()

table(net$colors)

#0    1    2    3    4    5    6    7    8    9   10 
#9 2471 1082  974  962  888  727  524  430  368  161 
 
moduleLabels=net$colors
moduleColors=labels2colors(net$colors)
MEs=net$MEs
geneTree=net$dendrograms[[1]]

#--------------associations across modules and treatments-------------------------------
#associations across modules and treatments
association <- moduleEigengenes(datExpr,mergedColors)$eigengenes

association <- orderMEs(association)
association_order=names(association)%>% gsub ("ME","",.)

association$Treatment=row.names(association) #if treatment was not there, it created a new column for it

# tidy & plot data
mME = association %>%
  pivot_longer(-Treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = association_order)
  )

library(ggplot2)
mME %>% ggplot(., aes(x=Treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Sample Relationships", y = "Modules", fill="Association (a.u.)")

module.table <- as.data.frame (table(moduleColors))


#----Generate csv file with ILMN IDs and module color and save----
module_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

#------plot for gene expression in each module (from Bioinformatics workbook)-------------------------------

# tan_of_interest = c("magenta")
# tan.submod=module_df %>% subset (colors %in% tan_of_interest)
# row.names(module_df) = module_df$gene_id
# 
# tan.subexpr = Data[tan.submod$gene_id,]
# #tan.subexpr=tan.subexpr[,-c(14,15)]
# 
# tan.submod_df = data.frame(tan.subexpr) %>%
#   mutate(
#     gene_id = row.names(.)
#   ) %>%
#   pivot_longer(-gene_id) %>%
#   mutate(
#     module = module_df[gene_id,]$colors
#   )
# 
# tan.color <- c("magenta")
# tan.submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
#   geom_line(aes(color = module),
#             alpha = 0.2,) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90)
#   ) + scale_color_manual(values = c('magenta4'))+
#   facet_grid(rows = vars(module)) +
#   labs(x = "Treatment",
#        y = "Expression")
# 

#----------Module-Trait Relationships------------------------------

MEs0=moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,datTraits,use='p')
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

heatmap.data <- merge(MEs, datTraits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

#https://github.com/kevinblighe/CorLevelPlot



# Set the width and height for the SVG output (in inches)
width_in_inches <- 10
height_in_inches <- 6

# SVG graphics device  https://r-coder.com/save-plot-r/
svg("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/plots_fasting_NGT/module_trait_relationship_fasting_NGT.svg")


CorLevelPlot(heatmap.data,
             y = names(heatmap.data)[1:11],
             x = names(heatmap.data)[12:22],
             titleX = "Trait Data",
             titleY = "Modules",
             rotTitleY = 90,
             corFUN='pearson',
             cexLabX = 1.0,
             rotLabX = 90,
             main = "Module-Trait Relationships",
             signifSymbols = c("***", "**", "*", ""),
             signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
             col = c("blue1", "skyblue", "white", "pink", "red"))

dev.off()

#----Gene module membership and trait significance and save GeneInfo CSV file-----
# Define a variable containing the trait column of datTrait
fasting_glucose = as.data.frame(datTraits$Fast_glucose);
names(fasting_glucose) = "Fasting Glucose"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership=as.data.frame(cor(datExpr,MEs,use='p'))
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership)=paste("MM",modNames,sep = "")
names(MMPvalue)=paste("p.MM",modNames,sep="")

geneTraitSignificance=as.data.frame(cor(datExpr,fasting_glucose,use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(fasting_glucose), sep="");
names(GSPvalue) = paste("p.GS.", names(fasting_glucose), sep="")

# dataframe with Most significant genes associated with Trait of interest
#Gyaan: # # Using the gene significance you can identify genes that have a high significance for trait of interest 
# # Using the module membership measures you can identify genes with high module membership in interesting modules.

annot=read.csv(file = "/home/admin/R_data/GSE153837_normalized_22-10-2023/Anno.csv")#annotation Data

gene_trait_significance_pval <- as.data.frame(GSPvalue)#here trait of interest was fasting glucose
gene_trait_significance_pval$ILMN_ID <- rownames(gene_trait_significance_pval)
sorted_gene_trait_significance_pval <- gene_trait_significance_pval[order(gene_trait_significance_pval$`p.GS.Fasting Glucose`),]
x=match(sorted_gene_trait_significance_pval$ILMN_ID,annot$ILMN_ID)
y=data.frame(ILMN_ID=sorted_gene_trait_significance_pval$ILMN_ID,
             GSP_val=sorted_gene_trait_significance_pval$`p.GS.Fasting Glucose`,
             gene_symbol=annot$SYMBOL[x],
             definition=annot$DEFINITION[x])
write.csv(y,file="genes_significantly_associated_with_fasting_glu.csv")

# sanity check
names(datExpr)
names(datExpr)[moduleColors=="purple"]

#dataframe with information about annotation and modules
#annot=read.csv(file = "C:/Users/ADMIN/Desktop/dataset_20percent/Anno.csv")
dim(annot)
names(annot)
probes=names(datExpr)
probes2annot=match(probes,annot$ILMN_ID)
sum(is.na(probes2annot))

#making a dataframe with information about annotation and modules
geneInfo0=data.frame(ILMN_ID=probes,
                     geneSymbol=annot$SYMBOL[probes2annot],
                     entrezID=annot$ENTREZ_GENE_ID[probes2annot],moduleColor=moduleColors,
                     ontology_process=annot$ONTOLOGY_PROCESS[probes2annot]
                     ,ontology_function=annot$ONTOLOGY_FUNCTION[probes2annot]
                     ,ontology_component=annot$ONTOLOGY_COMPONENT[probes2annot],
                     geneTraitSignificance,
                     GSPvalue)

#ordering module by their significance for choice of trait
modOrder=order(-abs(cor(MEs,fasting_glucose,use='p')))

#adding module membership information in the chosen order in previous line
for (mod in 1:ncol(geneModuleMembership))
{
  oldnames=names(geneInfo0)
  geneInfo0=data.frame(geneInfo0,
                       geneModuleMembership[,modOrder[mod]],
                       MMPvalue[,modOrder[mod]])
  names(geneInfo0)=c(oldnames,
                     paste("MM.",modNames[modOrder[mod]],sep=""),
                     paste("p.MM.",modNames[modOrder[mod]],sep=""))
}

#rearranging the data and saving a csv file

geneOrder=order(geneInfo0$moduleColor,-abs(geneInfo0$GS.Fasting.Glucose))
geneInfo=geneInfo0[geneOrder,]

write.csv(geneInfo,file="geneInfo.csv")





#----Plot for Gene Module Membership and Gene Trait Significance for Trait of Interest----
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Fasting Glucose Levels",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




#exporting network to cytoscape (we need to select the module we want to export)-----

TOM=TOMsimilarityFromExpr(datExpr,power=18)

#reading the annotation file
#annot=read.csv(file = "C:/Users/ADMIN/Desktop/dataset_20percent/Anno.csv")

#selecting modules
modules=c("yellow") #module of choice
probes=names(datExpr)
inModule=is.finite(match(moduleColors,modules))
modProbes=probes[inModule]
modGenes=annot$SYMBOL[match(modProbes,annot$ILMN_ID)]

#select corresponding topoogical overlap
modTOM=TOM[inModule,inModule]
dimnames(modTOM)=list(modGenes,modGenes)#changing names to gene symbols

setwd("/home/admin/R_data/GSE153837_normalized_22-10-2023/split-data/fasting_NGT/results/downstream_fasting_NGT/")

#exporting to cytoscape
cyt=exportNetworkToCytoscape(modTOM,
                             edgeFile = paste("CytoscapeInput-edges-",paste(modules,collapse="-"),".txt",sep=""),
                             nodeFile = paste("CytoscapeInput-nodes-",paste(modules,collapse="-"),".txt",sep=""),
                             weighted = TRUE,
                             threshold = 0.02,
                             nodeNames = modGenes,
                             altNodeNames = modProbes,
                             nodeAttr = moduleColors[inModule])





# #----------------------------------E X T R A-------------------------------------
# #----------(NOT Needed)generate and export network from Bioinformatics workbook---------#
# tan_of_interest = c("tan")
# tan.submod=module_df %>% subset (colors %in% tan_of_interest)
# row.names(module_df) = module_df$gene_id
# 
# tan.subexpr = femData[tan.submod$gene_id,]
# tan.subexpr=tan.subexpr[,-c(14,15)]
# 
# # Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
# TOM = TOMsimilarityFromExpr(t(tan.subexpr),
#                             power = picked_power)
# 
# row.names(TOM) = row.names(tan.subexpr)
# colnames(TOM) = row.names(tan.subexpr)
# 
# edge_list = data.frame(TOM) %>%
#   mutate(
#     gene1 = row.names(.)
#   ) %>%
#   pivot_longer(-gene1) %>%
#   dplyr::rename(gene2 = name, correlation = value) %>%
#   unique() %>%
#   subset(!(gene1==gene2)) %>%
#   mutate(
#     module1 = module_df[gene1,]$colors,
#     module2 = module_df[gene2,]$colors
#   )
# 
# # Export Network file to be read into Cytoscape, VisANT, etc
# write_delim(edge_list,
#             file = "edgelist.tsv",
#             delim = "\t")
# 
# 
# #the number of genes are equivalent to no. of columns in DatExpr
# ##the number of samples are equivalent to no. of rows in DatExpr
# 
# MEs0=moduleEigengenes(datExpr,moduleColors)$eigengenes
# MEs=orderMEs(MEs0)
# moduleTraitCor=cor(MEs,datTraits,use='p')
# moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)
# 
# #
# #------------(NOT needed)---module trait corr WGCNA Tutorial-------#
# sizeGrWindow(12,10)
# #will display correlations and pvalue
# 
# textmatrix=paste(signif(moduleTraitCor,2),"\n(",
#                  signif(moduleTraitPvalue,1),")",sep = "")
# dim(textmatrix)=dim(moduleTraitCor)
# par(mar=c(5,9,5,1))
# 
# #Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = textmatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# 
# 
# 
# #----------------------Enrichment analysis: NOT PERFORMED---------------------------#
# # annot=read.csv(file = "GeneAnnotation.csv")
# # dim(annot)
# # names(annot)
# # probes=names(datExpr)
# # probes2annot=match(probes,annot$substanceBXH)
# 
# #fetching locus link IDs for the probes from annot file
# allLLIDs=annot$LocusLinkID[probes2annot]
# 
# #selecting modules of interest
# intModules=c("brown","red","salmon")
# for (module in intModules)
# {
#   #sselect module probes
#   modGenes=(moduleColors==module)
#   #fetch Entrezcodes(locus link Ids)
#   modLLIDs=allLLIDs[modGenes]
#   #writing to file
#   fileName=paste("LocusLinkIDs-",module,".txt",sep="")
#   write.table(as.data.frame(modLLIDs),file=fileName,
#               row.names = FALSE,col.names = FALSE)
# }
# 
# #background for enrichment analysis, use all probes in analysis
# fileName=paste("LocusLinkIDs-all.txt",sep="")
# write.table(as.data.frame(allLLIDs),file=fileName,
#             row.names = FALSE,col.names = FALSE)
# 
# 
# GOenr=GOenrichmentAnalysis(moduleColors,allLLIDs,organism = "mouse",
#                            nBestP = 10)
# tab=GOenr$bestPTerms[[4]]$enrichment
# names(tab)
# #?GOenrichmentAnalysis
# 
# write.table(tab,file="GOEnrichment.csv",sep=",",quote = TRUE,row.names = FALSE)
# 
# keepCols=c(1,2,5,6,7,12,13)
# screenTab=tab[,keepCols]
# 
# #rounding numbers off by 2 places
# numCols=c(3,4)
# screenTab[,numCols]=signif(apply(screenTab[,numCols],2,as.numeric),2)
# 
# #truncate term name to max  40characters
# screenTab[,7]=substring(screenTab[,7],1,40)
# #shorten column names
# colnames(screenTab)=c("module","size","pval","Bonf","nInTerm","ont","term_name")
# rownames(screenTab)=NULL
# 
# 
# #Set the width of R’s output. The reader should play with this number to obtain satisfactory output.
# options(width=95)
# screenTab
# 
# 
# #--------------------------visualizing gene networks: NOT PERFORMED------#
# #(optional) calculate topological overlkao matrix anew
# #dissimilarity matrix
# dissTOM=1-TOMsimilarityFromExpr(datExpr,power=18)
# #raise to a power for better visibility of strong connections
# plotTOM=dissTOM^7
# 
# #setting diagonals to NA
# diag(plotTOM)=NA
# #plotting
# sizeGrWindow(9,9)
# TOMplot(plotTOM,geneTree,moduleColors,main="Network heatmap plot, all genes")
# 
# #subseting no.of genes forplotting
# nSelect=400
# set.seed(10)#for reproducibility
# select=sample(nGenes,size=nSelect)
# selectTOM=dissTOM[select,select]
# 
# # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
# selectTree = hclust(as.dist(selectTOM), method = "average")
# selectColors = moduleColors[select]
# 
# 
# # Open a graphical window
# sizeGrWindow(9,9)
# # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# # the color palette; setting the diagonal to NA also improves the clarity of the plot
# plotDiss = selectTOM^7;
# diag(plotDiss) = NA;
# myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
# TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=myheatcol)
# 
# #GYAAN
# # #Visualizing the network of eigengenes
# # It is often interesting to study the relationships among the found modules. One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation.
# #The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network.
# #It is usually informative to add a clinical trait (or multiple traits) to the eigengenes to see how the traits fit into the eigengene network
# 
# # Recalculate module eigengenes
# # MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# # 
# # #isolate wt. frrom clinical data
# # weight=as.data.frame(datTraits$weight_g)
# # names(weight)="weight"
# # #adding wt to existing module eigengenes
# # MET = orderMEs(cbind(MEs, fasting_glucose))
# # # Plot the relationships among the eigengenes and the trait
# # sizeGrWindow(5,7.5);
# # par(cex = 0.9)
# # plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
# #                       = 90)
# 
# # Plot the dendrogram
# sizeGrWindow(12,10);
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
#                       plotHeatmaps = FALSE)
# 
# # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,5,5),
#                       plotDendrograms = FALSE, xLabelsAngle = 0)
# 
# #6B. Intramodular analysis: Identifying driver genes 
# 
# # Calculate the module membership and the associated p-values
# 
# # The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# # This quantifies the similarity of all genes on the array to every module.
# 
# module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
# module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
# 
# 
# module.membership.measure.pvals[1:10,1:10]
# 
# 
# # Calculate the gene significance and associated p-values
# 
# gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
# gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
# 
# 
# gene.signf.corr.pvals %>% 
#   as.data.frame() %>% 
#   arrange(V1) %>% 
#   head(25)
# 
# 
# # Using the gene significance you can identify genes that have a high significance for trait of interest 
# # Using the module membership measures you can identify genes with high module membership in interesting modules.

