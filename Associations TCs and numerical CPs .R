# Load required libraries
library(foreign)
library(nnet)
library(ggplot2)
library(reshape2)
library(AER)
library(afex)
library(car)
library(dplyr)

# Load the dataset and create an empty association matrix
# Select the cell types and TCs of interest

# Replace 'CES_with_clinical_pathological_variables' with 'metadata'
# Include the CPs colnames of interest in the rownames
# Include the TCs of interest in the colnames

rownames(association_matrix) = c(colnames(metadata)[c(...)]) 
colnames(association_matrix) = c(colnames(metadata)[c(...)]) 

for(i in rownames(association_matrix)[c(...)]) 
{
  for(j in N:dim(metadata)[2]) # Replace N with column in which TCs start
  {
    metadata_subset = metadata[which(!is.na(metadata[,i])),]
    metadata_subset[,i] = as.numeric(as.character(metadata_subset[,i]))
    metadata_subset = metadata_subset[which(!is.na(metadata_subset[,i])),]
    
    test <- cor.test(metadata_subset[,i], metadata_subset[,j], method = "spearman")
    association_matrix[i,j-N-1] = -log10(test$p.value)*sign(test$estimate) # Replace N with column in which TCs start minus 1
  }
  print(i)
}

View(association_matrix) 

# Write association matrix to file
write.table(association_matrix, file = ".../.txt", sep = "\t", quote = FALSE) 
