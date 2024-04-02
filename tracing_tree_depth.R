# Load necessary libraries
library(MST) # For survival tree modeling
library(stringr) # For string manipulation
library(sqldf) # For SQL-like queries on data frames
library(partykit) # For decision tree visualization

# Load data from specified file path
inputData = read.table("/Users/arkajyotibhattacharya/Work\ outside\ projects/CRC_Dienstman_data/Validation_approach_in_R/Data/Annotations_DFS_SURVIVAL_24_03_with_flipped_consensus_mixing_matrix.txt", sep = "\t", header = TRUE)

# Standardize column names by replacing dots with underscores
colnames(inputData) = gsub("\\.","_", colnames(inputData))

# Filter data based on specific conditions
inputData = inputData[which(inputData$TNM_STAGE==2),]
inputData = inputData[which(inputData$SAMPLE_FROM_COLON_OR_RECTUM=="colon"),]

# Define significant variables
significant_TCs = c("CONSENSUS_IC_208",
                    "CONSENSUS_IC_117",
                    # Other significant variables listed here...
                    "CONSENSUS_IC_166")

# Initialize a similarity matrix with zeros
similarity_matrix_between_samples_combined = matrix(0, dim(inputData)[1], dim(inputData)[1])
colnames(similarity_matrix_between_samples_combined) = rownames(inputData)
rownames(similarity_matrix_between_samples_combined) = rownames(inputData)

# Start time measurement for the loop
time2 = proc.time()[3]

# Initialize a list to store fit results
fit_results=list()

# Initialize a variable to store tree summary
tree_summary_combined = NULL

# Loop for tree fitting (1000 iterations)
for(iteration in 1:1000) {
  
  # Start time measurement for each iteration
  time1 = proc.time()[3]
  
  # Randomly select 5 significant variables
  set.seed(iteration)
  tcs_to_include = sample(significant_TCs, 5)
  
  # Prepare data subset including selected variables
  inputData_subset = inputData[,c("DFS_EVENT", "DFS_YRS", tcs_to_include,colnames(inputData)[c(1:28,31:46)])]
  
  # Handle missing values and convert columns to factors
  for(i in 1:dim(inputData_subset)[2]) {
    length_nas = length(which(is.na(inputData_subset[,i])))
    if(length_nas > 0) {
      inputData_subset[,i] = as.character(inputData_subset[,i])
      inputData_subset[which(is.na(inputData_subset[,i])),i] = "NA"
      inputData_subset[,i] = as.factor(inputData_subset[,i])
    }
  }
  
  # Convert TNM_STAGE column to factor
  inputData_subset$TNM_STAGE = as.factor(inputData_subset$TNM_STAGE)
  
  # Convert character columns to factors
  class_col = sapply(inputData_subset, class)
  for(i in which(class_col=="character")) {
    inputData_subset[,i] = as.factor(inputData_subset[,i])
  }
  
  # Fit a survival tree model
  minsplits = 50
  tree_MST = MST(formula = as.formula(paste("Surv(DFS_YRS, DFS_EVENT) ~ ",paste(tcs_to_include, collapse = "+"),"| SAMPLE_ID")),
                 data = inputData_subset,
                 test = inputData_subset,
                 method = "independence",
                 minsplit = minsplits,
                 minevents = ceiling(minsplits/2),
                 minbucket = ceiling(minsplits/3),
                 selection.method = "test.sample",
                 plot.Ga = FALSE,
                 sortTrees = TRUE,
                 details = FALSE)
  
  # Extract the tree structure
  tree <- getTree(tree_MST , "0")
  allNodes <- nodeids(tree)
  
  # Store samples in each node
  samples_in_allNodes <- lapply(allNodes, function(id, tree) {
    samples_in_each_node <- rownames(data_party(tree, id= id))
    return(samples_in_each_node)
  }, tree)
  
  # Initialize similarity matrix for each node
  similarity_matrix_between_samples_in_each_node = matrix(0, dim(inputData)[1], dim(inputData)[1])
  colnames(similarity_matrix_between_samples_in_each_node) = rownames(inputData)
  rownames(similarity_matrix_between_samples_in_each_node) = rownames(inputData)
  
  # Populate similarity matrix between samples in each node
  for(node in allNodes[-1]) {
    similarity_matrix_between_samples_in_each_node[samples_in_allNodes[[node]],samples_in_allNodes[[node]]] = similarity_matrix_between_samples_in_each_node[samples_in_allNodes[[node]],samples_in_allNodes[[node]]] + 1
  }
  
  # Calculate shared proportion of edges
  shared_proportion_of_edges = similarity_matrix_between_samples_in_each_node
  for(i in 1:dim(shared_proportion_of_edges)[1]) {
    for(j in 1:i) {
      shared_proportion_of_edges[i,j] = 2*similarity_matrix_between_samples_in_each_node[i,j]/(similarity_matrix_between_samples_in_each_node[i,i]+ similarity_matrix_between_samples_in_each_node[j,j])
      shared_proportion_of_edges[j,i] =  shared_proportion_of_edges[i,j]
    }
  }
  
  # Update combined similarity matrix
  similarity_matrix_between_samples_combined = similarity_matrix_between_samples_combined + shared_proportion_of_edges
  
  # Store fit results
  fit_results[[iteration]] = tree
  
  # Find terminal nodes for each sample
  totalNodes <- nodeids(tree, terminal = T)
  samples_in_node <- lapply(totalNodes, function(id, tree) {
    samples_in_each_node <- rownames(data_party(tree, id= id))
    return(samples_in_each_node)
  }, tree)
  
  # Assign terminal node numbers to samples
  for(i in 1:length(totalNodes)) {
    inputData[samples_in_node[[i]],paste("terminal_node", iteration)] = totalNodes[i]
  }
  
  # Extract tree summary
  varid_locations = str_locate_all(pattern ='varid',  as.character(tree)[1])[[1]]
  mean_location = (varid_locations[,1]+varid_locations[,2])/2
  varid_locations = cbind(varid_locations,mean_location)
  id_locations = str_locate_all(pattern ='\\(id = ',  as.character(tree)[1])[[1]]
  mean_location = (id_locations[,1]+id_locations[,2])/2
  id_locations = cbind(id_locations,mean_location)
  
  text_to_search_in = as.character(tree)[1]
  tree_summary = as.data.frame(matrix(NA, dim(varid_locations)[1],3))
  colnames(tree_summary) = c("ID_number", "Var_number", "Var_name")
  for(i in 1:dim(varid_locations)[1]) {
    tree_summary$ID_number[i] = which(id_locations[,3]>varid_locations[i,3])[1]-1
    tree_summary$Var_number[i] = substr(text_to_search_in, str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+1
                                        ,str_locate_all(pattern ='varid = ',  text_to_search_in)[[1]][i,2]+2)
    tree_summary$Var_number[i] = as.numeric(as.character(gsub(",","",tree_summary$Var_number[i])))
  }
  
  tree_summary$Var_name = tcs_to_include[as.numeric(tree_summary$Var_number)-1]
  tab <- sapply(tree_summary$ID_number, function(id) {
    y <- dim(data_party(tree[id]))[1]
  })
  tree_summary$data_size = tab
  
  tree_summary = sqldf("select Var_name, max(data_size) as data_size from tree_summary group by Var_name order by data_size DESC")
  rownames(tree_summary) = tree_summary$Var_name
  tree_summary$iteration = iteration
  tree_summary_combined = rbind(tree_summary_combined,tree_summary)
  
  # Print time taken for each iteration
  print((proc.time()[3]-time1)/60)
}

# Print total time taken
print((proc.time()[3]-time2)/60)

# Convert similarity matrix to data frame
similarity_matrix_between_samples_combined = as.data.frame(similarity_matrix_between_samples_combined)
similarity_matrix_between_samples_combined$samples = rownames(similarity_matrix_between_samples_combined)

# Reorder columns
similarity_matrix_between_samples_combined = similarity_matrix_between_samples_combined[,c(447,1:446)]

# Assign column names
colnames(similarity_matrix_between_samples_combined)[2:447] = inputData$SAMPLE_ID

# Write similarity matrix to file
write.table(similarity_matrix_between_samples_combined, file = "Stage_2_original_similarity_matrix_between_samples_combined.txt", sep = "\t",quote = FALSE)

# Create a copy of similarity matrix and round off values
similarity_matrix_between_samples_combined_v1 = similarity_matrix_between_samples_combined
similarity_matrix_between_samples_combined_v1[,2:447] = round(similarity_matrix_between_samples_combined_v1[,2:447],4)

# Write rounded similarity matrix to file
write.table(similarity_matrix_between_samples_combined_v1, file = "Stage_2_original_similarity_matrix_between_samples_combined_v1.txt", sep = "\t",quote = FALSE)

# Write tree summary to file
write.table(tree_summary, file = "Stage_2_original_tree_summary.txt", sep = "\t",quote = FALSE)

# Save fit results to RData file
save(fit_results, file = "Stage_2_original_fit_results.RData")
