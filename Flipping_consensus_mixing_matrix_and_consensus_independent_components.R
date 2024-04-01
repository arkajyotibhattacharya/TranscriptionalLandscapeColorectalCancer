univariate_summary = read.table("/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/Results_04092020/univariate_survival_analysis_on_unflipped_mixing_matrix/Univariate_coxph_on_unflipped_mixing_matrix_weights.txt", sep = "\t", header = TRUE)

flip_or_not = ifelse(univariate_summary$exp.coef.>1, 1,-1)
table(flip_or_not)

consensus_mixing_matrix = read.table("/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/ICA_Analysis/Consensus\ mixing\ matrix.txt", sep = "\t", header = TRUE, row.names = "X")
colnames(consensus_mixing_matrix) = toupper(gsub("\\.", "_", colnames(consensus_mixing_matrix)))
flipped_consensus_mixing_matrix = t(t(consensus_mixing_matrix)*flip_or_not)

table(flipped_consensus_mixing_matrix[,1]/consensus_mixing_matrix[,1])

consensus_independent_components = read.table("/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/ICA_Analysis/Consensus independent components.txt", sep = "\t", header = TRUE, row.names = "X")
colnames(consensus_independent_components) = toupper(gsub("\\.", "_", colnames(consensus_independent_components)))
flipped_consensus_independent_components = t(t(consensus_independent_components)*flip_or_not)

table(flipped_consensus_independent_components[,1]/consensus_independent_components[,1])

write.table(flipped_consensus_mixing_matrix, file = "/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/Results_04092020/flipped_data_on_univariate_survival_analysis/flipped_consensus_mixing_matrix.txt", sep = "\t")
write.table(flipped_consensus_independent_components, file = "/Volumes/Citrix\ Files/Folders/Data/TranscriptionalLandscapeLargeIntestine/Results_04092020/flipped_data_on_univariate_survival_analysis/flipped_consensus_independent_components.txt", sep = "\t")

