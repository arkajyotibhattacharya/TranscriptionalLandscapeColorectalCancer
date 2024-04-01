library(MST)
library(stringr)
library(sqldf)
library(partykit)
inputData = read.table("/Users/arkajyotibhattacharya/Work\ outside\ projects/CRC_Dienstman_data/Validation_approach_in_R/Data/Annotations_DFS_SURVIVAL_24_03_with_flipped_consensus_mixing_matrix.txt", sep = "\t", header = TRUE)
colnames(inputData) = gsub("\\.","_", colnames(inputData))
inputData = inputData[which(inputData$TNM_STAGE==2),]
inputData = inputData[which(inputData$SAMPLE_FROM_COLON_OR_RECTUM=="colon"),]

significant_TCs = c("CONSENSUS_IC_208",
                    "CONSENSUS_IC_117",
                    "CONSENSUS_IC_193",
                    "CONSENSUS_IC_202",
                    "CONSENSUS_IC_149",
                    "CONSENSUS_IC_21",
                    "CONSENSUS_IC_189",
                    "CONSENSUS_IC_148",
                    "CONSENSUS_IC_80",
                    "CONSENSUS_IC_55",
                    "CONSENSUS_IC_61",
                    "CONSENSUS_IC_111",
                    "CONSENSUS_IC_18",
                    "CONSENSUS_IC_116",
                    "CONSENSUS_IC_153",
                    "CONSENSUS_IC_131",
                    "CONSENSUS_IC_113",
                    "CONSENSUS_IC_123",
                    "CONSENSUS_IC_58",
                    "CONSENSUS_IC_136",
                    "CONSENSUS_IC_17",
                    "CONSENSUS_IC_22",
                    "CONSENSUS_IC_77",
                    "CONSENSUS_IC_44",
                    "CONSENSUS_IC_56",
                    "CONSENSUS_IC_65",
                    "CONSENSUS_IC_126",
                    "CONSENSUS_IC_171",
                    "CONSENSUS_IC_31",
                    "CONSENSUS_IC_69",
                    "CONSENSUS_IC_135",
                    "CONSENSUS_IC_103",
                    "CONSENSUS_IC_157",
                    "CONSENSUS_IC_35",
                    "CONSENSUS_IC_216",
                    "CONSENSUS_IC_74",
                    "CONSENSUS_IC_25",
                    "CONSENSUS_IC_59",
                    "CONSENSUS_IC_125",
                    "CONSENSUS_IC_121",
                    "CONSENSUS_IC_212",
                    "CONSENSUS_IC_5",
                    "CONSENSUS_IC_137",
                    "CONSENSUS_IC_88",
                    "CONSENSUS_IC_52",
                    "CONSENSUS_IC_106",
                    "CONSENSUS_IC_27",
                    "CONSENSUS_IC_38",
                    "CONSENSUS_IC_180",
                    "CONSENSUS_IC_147",
                    "CONSENSUS_IC_115",
                    "CONSENSUS_IC_220",
                    "CONSENSUS_IC_166")

similarity_matrix_between_samples_combined = matrix(0, dim(inputData)[1], dim(inputData)[1])
colnames(similarity_matrix_between_samples_combined) = rownames(inputData)
rownames(similarity_matrix_between_samples_combined) = rownames(inputData)

time2 = proc.time()[3]
fit_results=list()
tree_summary_combined = NULL
# Terminal_node_info_for_each_sample_in_each_tree = NULL
for(iteration in 1:1000)
{
  time1 = proc.time()[3]
  
  set.seed(iteration)
  tcs_to_include = sample(significant_TCs, 5)
  
  inputData_subset = inputData[,c("DFS_EVENT", "DFS_YRS", tcs_to_include,colnames(inputData)[c(1:28,31:46)])]
  # rownames(inputData) = inputData_subset$SAMPLE_ID
  
  for(i in 1:dim(inputData_subset)[2])
  {
    length_nas = length(which(is.na(inputData_subset[,i])))
    if(length_nas>0)
    {
      inputData_subset[,i] = as.character(inputData_subset[,i])
      inputData_subset[which(is.na(inputData_subset[,i])),i] = "NA"
      inputData_subset[,i] = as.factor(inputData_subset[,i])
      
    }
    
  }
  inputData_subset$TNM_STAGE = as.factor(inputData_subset$TNM_STAGE)
  class_col = sapply(inputData_subset, class)
  for(i in which(class_col=="character"))
  {
    inputData_subset[,i] = as.factor(inputData_subset[,i])
  }
  # inputData_subset$SAMPLE_ID = as.character(inputData_subset$SAMPLE_ID)
  minsplits = 50
  tree_MST = MST(formula = as.formula(paste("Surv(DFS_YRS, DFS_EVENT) ~ ",paste(tcs_to_include, collapse = "+"),"| SAMPLE_ID"))
                                        
                                        ,
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
  
  tree <- getTree(tree_MST , "0")
  allNodes <- nodeids(tree)
  samples_in_allNodes <- lapply(allNodes, function(id, tree){
    samples_in_each_node <- rownames(data_party(tree, id= id))
    return(samples_in_each_node)
  },tree)
  
  # similarity_matrix_between_samples_in_each_node = list()
  # 
  # for(node in allNodes[-1])
  # {
  #   similarity_matrix_between_samples_in_each_node[[node]] = matrix(0, dim(inputData)[1], dim(inputData)[1])
  #   colnames(similarity_matrix_between_samples_in_each_node[[node]]) = rownames(inputData)
  #   rownames(similarity_matrix_between_samples_in_each_node[[node]]) = rownames(inputData)
  #   
  #   similarity_matrix_between_samples_in_each_node[[node]][samples_in_allNodes[[node]],samples_in_allNodes[[node]]] = 1
  #   
  # }
  
  similarity_matrix_between_samples_in_each_node = matrix(0, dim(inputData)[1], dim(inputData)[1])
  colnames(similarity_matrix_between_samples_in_each_node) = rownames(inputData)
  rownames(similarity_matrix_between_samples_in_each_node) = rownames(inputData)
  
  for(node in allNodes[-1])
  {
    
    similarity_matrix_between_samples_in_each_node[samples_in_allNodes[[node]],samples_in_allNodes[[node]]] = similarity_matrix_between_samples_in_each_node[samples_in_allNodes[[node]],samples_in_allNodes[[node]]] + 1
    
  }
  
  shared_proportion_of_edges = similarity_matrix_between_samples_in_each_node
  
  for(i in 1:dim(shared_proportion_of_edges)[1])
  {
    for(j in 1:i)
    {
      shared_proportion_of_edges[i,j] = 2*similarity_matrix_between_samples_in_each_node[i,j]/(similarity_matrix_between_samples_in_each_node[i,i]+ similarity_matrix_between_samples_in_each_node[j,j])
      shared_proportion_of_edges[j,i] =  shared_proportion_of_edges[i,j]
    }
  }
  
  similarity_matrix_between_samples_combined = similarity_matrix_between_samples_combined+shared_proportion_of_edges
  fit_results[[iteration]] = tree
  
  #finding terminal nodes for each sample
  totalNodes <- nodeids(tree, terminal = T)
  
  samples_in_node <- lapply(totalNodes, function(id, tree){
    samples_in_each_node <- rownames(data_party(tree, id= id))
    return(samples_in_each_node)
  },tree)
  
  for(i in 1:length(totalNodes))
  {
    inputData[samples_in_node[[i]],paste("terminal_node", iteration)] = totalNodes[i]
  }
  
  
  varid_locations = str_locate_all(pattern ='varid',  as.character(tree)[1])[[1]]
  mean_location = (varid_locations[,1]+varid_locations[,2])/2
  varid_locations = cbind(varid_locations,mean_location)
  id_locations = str_locate_all(pattern ='\\(id = ',  as.character(tree)[1])[[1]]
  mean_location = (id_locations[,1]+id_locations[,2])/2
  id_locations = cbind(id_locations,mean_location)
  
  text_to_search_in = as.character(tree)[1]
  tree_summary = as.data.frame(matrix(NA, dim(varid_locations)[1],3))
  colnames(tree_summary) = c("ID_number", "Var_number", "Var_name")
  for(i in 1:dim(varid_locations)[1])
  {
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
  print((proc.time()[3]-time1)/60)
  
}
print((proc.time()[3]-time2)/60)
similarity_matrix_between_samples_combined = as.data.frame(similarity_matrix_between_samples_combined)
similarity_matrix_between_samples_combined$samples = rownames(similarity_matrix_between_samples_combined)

similarity_matrix_between_samples_combined = similarity_matrix_between_samples_combined[,c(447,1:446)]

colnames(similarity_matrix_between_samples_combined)[2:447] = inputData$SAMPLE_ID
write.table(similarity_matrix_between_samples_combined, file = "Stage_2_original_similarity_matrix_between_samples_combined.txt", sep = "\t",quote = FALSE)
similarity_matrix_between_samples_combined_v1 = similarity_matrix_between_samples_combined
similarity_matrix_between_samples_combined_v1[,2:447] = round(similarity_matrix_between_samples_combined_v1[,2:447],4)
write.table(similarity_matrix_between_samples_combined_v1, file = "Stage_2_original_similarity_matrix_between_samples_combined_v1.txt", sep = "\t",quote = FALSE)

write.table(tree_summary, file = "Stage_2_original_tree_summary.txt", sep = "\t",quote = FALSE)
save(fit_results, file = "Stage_2_original_fit_results.RData")

