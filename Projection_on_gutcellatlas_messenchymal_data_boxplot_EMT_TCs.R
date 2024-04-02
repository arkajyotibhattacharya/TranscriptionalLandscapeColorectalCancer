# Load necessary packages
library(anndata) # For reading AnnData objects
reticulate::install_python() # Install Python environment if not installed already

# Read AnnData object
ad <- anndata::read_h5ad('/Users/arkajyotibhattacharya/Downloads/Mesenchyme_log_counts02_v2.h5ad')

# Load data from h5ad file
library(data.table) # For reading data from h5ad file
data <- read_h5ad(file = "/Users/rudolffehrmann/Downloads/Mesenchyme_log_counts02_v2.h5ad")

# Check class of 'data'
class(data)

# Subset 'data' to get expression matrix
sub_data = data$X
sub_data_matrix = as.matrix(sub_data)

# Set seed for reproducibility
set.seed(1234)

# Randomly select 10% of rows from the expression matrix
sub_data_matrix_10_percent = sub_data_matrix[sample(c(1:dim(sub_data_matrix)[1]),round(0.1*dim(sub_data_matrix)[1])),]

# Get dimensions of 'sub_data'
dim(sub_data)

# Get annotation information
annotation = data$obs

# Write annotation to file
write.table(annotation, file = "/Users/rudolffehrmann/Downloads/annotation.txt", sep = "\t", quote = FALSE)

# Write 'sub_data_matrix_10_percent' to file
write.table(sub_data_matrix_10_percent, file = "/Users/rudolffehrmann/Downloads/sub_data_matrix_10_percent.txt", sep = "\t", quote = FALSE)

# Write 'sub_data_matrix_10_percent' to file using 'fwrite'
fwrite(sub_data_matrix_10_percent, file = "/Users/rudolffehrmann/Downloads/sub_data_matrix.txt", sep = "\t", quote = FALSE)

# Read 'sub_data_matrix_10_percent' from file
sub_data_matrix_10_percent = data.frame(fread("/Users/arkajyotibhattacharya/Downloads/sub_data_matrix_10_percent.txt"), row.names = 1)

# Transpose 'sub_data_matrix_10_percent'
sub_data_matrix_10_percent = t(sub_data_matrix_10_percent)

# Read genomic mapping file
genomic_mapping_file = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Databases/Genomic mapping files/Entrezid_mapping_using_org_Hs_eg_db_21072022.txt"))

# Get common genes between genomic mapping file and 'sub_data_matrix_10_percent'
common_genes = intersect(genomic_mapping_file$SYMBOL, rownames(sub_data_matrix_10_percent))

# Set row names of 'genomic_mapping_file'
rownames(genomic_mapping_file) = genomic_mapping_file$SYMBOL

# Subset 'sub_data_matrix_10_percent' and 'genomic_mapping_file' using common genes
sub_data_matrix_10_percent = sub_data_matrix_10_percent[common_genes,]
genomic_mapping_file = genomic_mapping_file[common_genes,]

# Set row names of 'sub_data_matrix_10_percent' using mapped entrez IDs
rownames(sub_data_matrix_10_percent) = genomic_mapping_file$mapped_entrez_v1

# Read annotation from file
annotation = data.frame(fread(file = "/Users/arkajyotibhattacharya/Downloads/annotation.txt"), row.names = 1)

# Get common samples between annotation and 'sub_data_matrix_10_percent'
common_samples =intersect(rownames(annotation), colnames(sub_data_matrix_10_percent))

# Subset 'annotation' and 'sub_data_matrix_10_percent' using common samples
annotation = annotation[common_samples,]
sub_data_matrix_10_percent = sub_data_matrix_10_percent[,common_samples]

# Write 'sub_data_matrix_10_percent' to file
write.table(sub_data_matrix_10_percent, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/CRC_Dienstman_data/Rebuttal/Results/gutcellatlas_sub_data_matrix_10_percent.txt", sep = "\t", quote = FALSE)

# Write 'annotation' to file
write.table(annotation, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/CRC_Dienstman_data/Rebuttal/Results/gutcellatlas_annotation_10_percent.txt", sep = "\t", quote = FALSE)

# Read 'components' from file
components = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/CRC_Dienstman_data/Data/flipped_consensus_independent_components.txt"), row.names = 1)

# Read genomic mapping file
genomic_mapping_file = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Projects_27092022/Databases/Genomic mapping files/Genomic_Mapping_hgu133plus2_using_jetscore_3003201_v1.txt"))

# Get common genes between genomic mapping file and 'components'
common_genes = intersect(genomic_mapping_file$PROBESET, rownames(components))

# Set row names of 'genomic_mapping_file'
rownames(genomic_mapping_file) = genomic_mapping_file$PROBESET

# Subset 'components' and 'genomic_mapping_file' using common genes
components = components[common_genes,]
genomic_mapping_file = genomic_mapping_file[common_genes,]

# Set row names of 'components' using entrez IDs
rownames(components) = genomic_mapping_file$ENTREZID

# Get indices of EMT_TCs
EMT_TCs_index = paste("CONSENSUS_IC_", c(208, 117, 193, 202, 116,153, 58, 136, 77), sep = "")

# Subset 'components' to get EMT components
components_EMT = components[,EMT_TCs_index]

# Write 'components_EMT' to file
write.table(components_EMT, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/CRC_Dienstman_data/Rebuttal/Results/components_EMT.txt", sep = "\t", quote = FALSE)

# Read corrected mixing matrix from file
corrected_mixing_matrix = data.frame(fread("/Users/arkajyotibhattacharya/Downloads/Project_data_on_independent_components_f98425a4-efeb-4ba9-9413-8bed2f6d9b21/mixing_matrix.tsv"), row.names = 1)

# Get common samples between corrected mixing matrix and annotation
common_samples = intersect(rownames(corrected_mixing_matrix), rownames(annotation))

# Subset corrected mixing matrix and annotation using common samples
corrected_mixing_matrix = corrected_mixing_matrix[common_samples,]
annotation = annotation[common_samples,]

# Combine corrected mixing matrix and annotation into a single dataframe
corrected_mixing_matrix = as.data.frame(cbind(corrected_mixing_matrix,annotation ))

# Order corrected mixing matrix based on annotation
corrected_mixing_matrix = corrected_mixing_matrix[order(corrected_mixing_matrix$annotation),]

# Extract annotation row data
annotation_row_data = as.data.frame(corrected_mixing_matrix[,c(11,12,29)])

# Convert columns of annotation row data to factors
for(i in 1:3) {
  annotation_row_data[,i] = as.factor(annotation_row_data[,i])
}

# Set row names of annotation row data
rownames(annotation_row_data) = rownames(corrected_mixing_matrix)

# Extract relevant data columns
data = corrected_mixing_matrix[,c(1:9)]


# Generate boxplot for CONSENSUS_IC_208
pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/CRC_Dienstman_data/Rebuttal/Plots/corrected_mixing_matrix_208_per_celltype.pdf", width = 12,height = 10)
par(mar = c(12, 5, 4, 2) )
boxplot(CONSENSUS_IC_208 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_208", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19)
dev.off()

# Generate boxplot for EMT components per cell type
pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My Drive/Laptop backup/Work outside projects 09102022/Work outside projects/CRC_Dienstman_data/Rebuttal/Plots/corrected_mixing_matrix_EMT_per_celltype.pdf", width = 36,height = 33)
par(mar = c(23, 5, 4, 2) )
par(mfrow = c(3,3))
for(i in EMT_TCs_index) {
  boxplot(get(i) ~ annotation, data = corrected_mixing_matrix, main = i, ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
}
dev.off()
