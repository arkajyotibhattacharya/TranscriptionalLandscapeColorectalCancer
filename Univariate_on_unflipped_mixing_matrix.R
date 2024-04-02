# Loading required packages
library(survival) # For survival analysis
library(parallel) # For parallel processing

# Read data and adjust column names
data = read.table("/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/Annotations/Annotations_DFS_SURVIVAL_24.03.txt", sep = "\t", header = TRUE)
colnames(data) = gsub("\\.","-", colnames(data))

# Subset data to relevant columns
data = as.data.frame(data)
data = data[,c(1:46)]

# Read consensus mixing matrix
cmm = read.table("/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/ICA_Analysis/Consensus\ mixing\ matrix.txt", sep = "\t", header = TRUE, row.names = "X")

# Adjust row names and subset common samples
rownames(data) = as.character(data$SAMPLE_ID)
common_samples = intersect(rownames(cmm), rownames(data))
data = data[common_samples,]
cmm = cmm[common_samples,]
colnames(cmm) = toupper(gsub("\\.", "_", colnames(cmm)))

# Combine data and mixing matrix
data = as.data.frame(data)
cmm = as.data.frame(cmm)
data = as.data.frame(cbind(data,cmm))

# Prepare matrix for analysis
follow_up= data$DFS_YRS
event = as.numeric(data$DFS_EVENT)
outcome_variable_range=c(47:266)
list_of_predictors=names(data[outcome_variable_range])
newdata = data[,outcome_variable_range]						

# Coxph loop
result <- lapply(newdata, function(x){
  coxph(as.formula(paste("Surv(follow_up,event)~",paste("x")))) # Fit Cox proportional hazards model
})

# Create matrix for results  
full_coef_matrix = NULL
for(i in names(result)){
  each_result_summary = summary(result[[i]])
  full_coef_matrix = rbind(full_coef_matrix, cbind(each_result_summary$coefficients,each_result_summary$conf.int)[1,])
  rownames(full_coef_matrix)[dim(full_coef_matrix)[1]] = i
}

# Write results to file
write.table(full_coef_matrix, file = "/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/Results_04092020/univariate_survival_analysis_on_unflipped_mixing_matrix/Univariate_coxph_on_unflipped_mixing_matrix_weights.txt", sep = "\t")
