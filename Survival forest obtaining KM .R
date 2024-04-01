# Load required libraries
library(data.table)
library(dendextend)
require("survival")
library(survminer)
library(dplyr)

# Clear workspace
rm(list = ls())

# Read similarity matrix data
mydata = read.table(".../similarity_matrix_between_samples_combined_metric.txt", sep = "\t", header = TRUE)
mydata = as.data.frame(mydata)
rownames(mydata) = mydata$samples
mydata = mydata[,-1] 
rownames(mydata) = colnames(mydata)
mydata = as.matrix(mydata)

# Check for symmetry in the matrix
sum(rownames(mydata)==colnames(mydata))

# Calculate Spearman correlation between samples
spearman_cor_samples = cor(t(mydata))

# Check symmetry after correlation calculation
sum(rownames(mydata) == colnames(spearman_cor_samples))

# Perform hierarchical clustering
hclust_sample_types = hclust(as.dist(1 - spearman_cor_samples), method = "ward.D2")

# Reorder data based on clustering
mydata = mydata[hclust_sample_types$order,]

# Determine clusters
clusters = cutree(hclust_sample_types, k = #number of clusters)

# Read metadata
inputData = read.table("HERE YOUR METADATA")
colnames(inputData) = gsub("\\.", "_", colnames(inputData))
inputData <- inputData[c(..)] # indicate here the columns of interest
names(inputData) <- as.matrix(inputData[1, ])
inputData <- inputData[-1, ]
inputData[] <- lapply(inputData, function(x) type.convert(as.character(x)))
rownames(inputData) = inputData$#here your sample ID
  inputData = inputData[,-1]

# Filter metadata to match with sample data
inputData <- subset(inputData, rownames(inputData) %in% rownames(mydata))

# Add cluster information to metadata
inputData = as.data.frame(cbind(inputData, clusters))

# Generate survival plots for each cluster
pdf("here your file directory.pdf")
for(i in #number of clusters)
    {
      inputData_subset = inputData[,c("survival outcome in months/years", "survival event", paste("cluster_",i, sep = ""))]
      colnames(inputData_subset)[3] =  "Clusters"
      fit <- survfit(Surv(#survival outcome months/years, #survival event) ~ Clusters, data = inputData_subset)
        print(ggsurvplot(fit
                         ,pval = TRUE # Add p-value
                         ,risk.table = TRUE
                         ,tables.height = 0.4
                         ,surv.plot.height = 3
        ))
    }
    dev.off()
    
    # Generate contingency tables
    table(inputData$#lastcluster, rownames(inputData))
          table(rownames(inputData), inputData$lastcluster)
          table(inputData$#secondcluster, rownames(inputData))
                
                # Calculate pairwise log-transformed p-values for survival differences between clusters
                for(clusters in #number of cluster)
                    {
                      survdiff_matrix = matrix(NA, clusters, clusters)
                      for(i in 1:clusters)
                      {
                        for(j in 1:i)
                        {
                          if(i != j)
                          {
                            inputData_subset = inputData[,c("survival outcome in months/years", "survival event", paste("cluster_", clusters, sep = ""))]
                            inputData_subset = inputData_subset[which(inputData_subset[,paste("cluster_", clusters, sep = "" )] %in% c(i, j)),]
                            colnames(inputData_subset)[3] =  "Clusters"
                            
                            fit <- survdiff(Surv(#survival outcome months/years, #survival event) ~ Clusters, data = inputData_subset)
                              survdiff_matrix[i,j] <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
                              survdiff_matrix[j,i] = survdiff_matrix[i,j]
                          }
                        }
                      }
                      
                      write.table(survdiff_matrix, file = paste(".\forest_score_metric_survdiff_between_clusters_",clusters,"_on_projected_mixing_matrix.txt", sep = ""), sep = "\t", quote = FALSE)
                    }
                    
                    # Loop through clusters to generate density plots of log-transformed p-values
                    setwd(".\.")
                    list_of_files = list.files(".\.")
                    for(clusters in #numberofclusters)
                        {
                          forest_score_metric_survdiff_between_clusters_projected = read.table(paste("forest_score_metric_survdiff_between_clusters_",clusters,"_on_projected_mixing_matrix", sep = ""), sep = "\t", header = TRUE)
                          
                          forest_score_metric_survdiff_between_clusters_projected = t(t(-log10(forest_score_metric_survdiff_between_clusters_projected)))
                          
                          forest_score_metric_survdiff_between_clusters_projected = ifelse(is.infinite(forest_score_metric_survdiff_between_clusters_projected), 317, forest_score_metric_survdiff_between_clusters_projected)
                          
                          pdf(paste(".\Density_of_pairwise_survdiff_logtransformed_pvalue_per_method_red_mRNA_green_WMM_black_DFSMM_blue_forest_clusers_",clusters,"_with_projected.pdf", sep = ""))
                          plot(density(as.matrix(forest_score_metric_survdiff_between_clusters_projected), na.rm = TRUE), col="blue")
                          dev.off()
                        }
                        
                        # Calculate and print summary statistics for log-transformed p-values
                        for(clusters in #number of clusters)
                            {
                              forest_score_metric_survdiff_between_clusters_projected = read.table(paste("forest_score_metric_survdiff_between_clusters_",clusters,"_on_projected_mixing_matrix.txt", sep = ""), sep = "\t", header = TRUE)
                              
                              print("**************************")
                              
                              print("forest_score_metric_survdiff_between_clusters_projected")
                              print(length(which(forest_score_metric_survdiff_between_clusters_projected < 0.01/clusters/(clusters-1)*2)))
                              
                              print(clusters)
                              print("clusters")
                              print("**************************")
                            }
                            