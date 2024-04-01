# Load required libraries
library(foreign)
library(nnet)
library(ggplot2)
library(reshape2)
library(AER)
library(afex)
library(car)

# Load metadata and create association matrix
metadata = read.table("...", sep = "\t", header = TRUE)

association_matrix = as.data.frame(matrix(NA, N, M)) # N is the number of clinicopathological paramters (CPs) and M the number of TCs
colnames(association_matrix) = c(colnames(metadata)[c(...)])
rownames(association_matrix) = c(colnames(metadata)[c(...)])

for(i in c(1:dim(association_matrix)[1])[c(1)])
{
  for(j in N:dim(metadata)[2]) # n is the colnumber in which the first TC is
  {
    metadata_subset = metadata[which(!is.na(metadata[,rownames(association_matrix)[i]])),]
    metadata_subset[, rownames(association_matrix)[i]] = as.character(metadata_subset[, rownames(association_matrix)[i]])
    metadata_subset = metadata_subset[which(metadata_subset[, rownames(association_matrix)[i]] != "unknown"),]
    test <- multinom(metadata_subset[,rownames(association_matrix)[i]] ~ metadata_subset[,j], data = metadata_subset)
    association_matrix[i,j-N-1] = -log10(Anova(test, type="III")$`Pr(>Chisq)`) # here the first TC colnumber minus 1
  }
  print(i)
}

# Extract association for specific clinical parameter
assoc = (association_matrix["HERE CP COLNAME",]) 

# Calculate mean of clinical parameter activity scores
mean_of_CP_activity_scores = matrix(NA, length(assoc), 2)
rownames(mean_of_CP_activity_scores) = paste(colnames(assoc), sep = "")
colnames(mean_of_CP_activity_scores) = c("outcome1", "outcome2")
mean_of_CP_activity_scores = as.data.frame(mean_of_CP_activity_scores)

for(CONSENSUS_IC_ in 1:dim(mean_of_CP_activity_scores)[1])
{
  mean_of_CP_activity_scores[CONSENSUS_IC_,"outcome1"] = mean(metadata[which(metadata$clinical_parameters == "outcome1"),rownames(mean_of_CP_activity_scores)[CONSENSUS_IC_]])
  mean_of_CP_activity_scores[CONSENSUS_IC_,"outcome2"] = mean(metadata[which(metadata$clinical_parameters == "outcome2"),rownames(mean_of_CP_activity_scores)[CONSENSUS_IC_]])
} 

# View mean clinical parameter activity scores
View(t(mean_of_CP_activity_scores))
View(t(association_matrix))

# Write association matrix and mean clinical parameter activity scores to file
write.table(association_matrix, file = "...", sep = "\t", quote = FALSE) 
write.table(mean_of_GENDER_INFERRED_activity_scores, file = "...", sep = "\t", quote = FALSE)
