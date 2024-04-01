library(anndata)
reticulate::install_python()
ad <- anndata::read_h5ad('/Users/arkajyotibhattacharya/Downloads/Mesenchyme_log_counts02_v2.h5ad')
library(data.table)
data <- read_h5ad(file = "/Users/rudolffehrmann/Downloads/Mesenchyme_log_counts02_v2.h5ad")


class(data)

sub_data = data$X
sub_data_matrix = as.matrix(sub_data)
set.seed(1234)
sub_data_matrix_10_percent = sub_data_matrix[sample(c(1:dim(sub_data_matrix)[1]),round(0.1*dim(sub_data_matrix)[1])),]
dim(sub_data)
annotation = data$obs
write.table(annotation, file = "/Users/rudolffehrmann/Downloads/annotation.txt", sep = "\t", quote = FALSE)
write.table(sub_data_matrix_10_percent, file = "/Users/rudolffehrmann/Downloads/sub_data_matrix_10_percent.txt", sep = "\t", quote = FALSE)
fwrite(sub_data_matrix_10_percent, file = "/Users/rudolffehrmann/Downloads/sub_data_matrix.txt", sep = "\t", quote = FALSE)

sub_data_matrix_10_percent = data.frame(fread("/Users/arkajyotibhattacharya/Downloads/sub_data_matrix_10_percent.txt"), row.names = 1)
sub_data_matrix_10_percent = t(sub_data_matrix_10_percent)
genomic_mapping_file = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Databases/Genomic\ mapping\ files/Entrezid_mapping_using_org_Hs_eg_db_21072022.txt"))
common_genes = intersect(genomic_mapping_file$SYMBOL, rownames(sub_data_matrix_10_percent))
rownames(genomic_mapping_file) = genomic_mapping_file$SYMBOL
sub_data_matrix_10_percent = sub_data_matrix_10_percent[common_genes,]
genomic_mapping_file = genomic_mapping_file[common_genes,]
rownames(sub_data_matrix_10_percent) = genomic_mapping_file$mapped_entrez_v1

annotation = data.frame(fread(file = "/Users/arkajyotibhattacharya/Downloads/annotation.txt"), row.names = 1)
common_samples =intersect( rownames(annotation), colnames(sub_data_matrix_10_percent))

annotation = annotation[common_samples,]

sub_data_matrix_10_percent = sub_data_matrix_10_percent[,common_samples]
write.table(sub_data_matrix_10_percent, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Rebuttal/Results/gutcellatlas_sub_data_matrix_10_percent.txt", sep = "\t", quote = FALSE)
write.table(annotation, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Rebuttal/Results/gutcellatlas_annotation_10_percent.txt", sep = "\t", quote = FALSE)

components = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Data/flipped_consensus_independent_components.txt"), row.names = 1)
genomic_mapping_file = data.frame(fread("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Projects_27092022/Databases/Genomic\ mapping\ files/Genomic_Mapping_hgu133plus2_using_jetscore_3003201_v1.txt"))
common_genes = intersect(genomic_mapping_file$PROBESET, rownames(components))
rownames(genomic_mapping_file) = genomic_mapping_file$PROBESET
components = components[common_genes,]
genomic_mapping_file = genomic_mapping_file[common_genes,]
rownames(components) = genomic_mapping_file$ENTREZID

EMT_TCs_index = paste("CONSENSUS_IC_", c(208, 117, 193, 202, 116,153, 58, 136, 77), sep = "")
components_EMT = components[,EMT_TCs_index]
write.table(components_EMT, file = "/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Rebuttal/Results/components_EMT.txt", sep = "\t", quote = FALSE)


corrected_mixing_matrix = data.frame(fread("/Users/arkajyotibhattacharya/Downloads/Project_data_on_independent_components_f98425a4-efeb-4ba9-9413-8bed2f6d9b21/mixing_matrix.tsv"), row.names = 1)
common_samples = intersect(rownames(corrected_mixing_matrix), rownames(annotation))

corrected_mixing_matrix = corrected_mixing_matrix[common_samples,]
annotation = annotation[common_samples,]
corrected_mixing_matrix = as.data.frame(cbind(corrected_mixing_matrix,annotation ))

corrected_mixing_matrix = corrected_mixing_matrix[order(corrected_mixing_matrix$annotation),]

annotation_row_data = as.data.frame(corrected_mixing_matrix[,c(11,12,29)])
for(i in 1:3)
{
  annotation_row_data[,i] = as.factor(annotation_row_data[,i])
}
# annotation_row_data = as.data.frame(as.factor(data_v1[,c(10,11)]))
rownames(annotation_row_data) = rownames(corrected_mixing_matrix)
data = corrected_mixing_matrix[,c(1:9)]
pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Rebuttal/plots/Heatmap_of_EMT_TCs_in_all_samples_to_show_heterogeneity_with_MSI_CMS_for_annotation_available_samples.pdf", width = 9, height = 8)
pheatmap::pheatmap(
  mat               = data,
  color             = colorRampPalette(c("navy", "white", "red"))(100),
  # breaks            = seq(-0.4, 0.4, length.out = 100),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  annotation_row    = annotation_row_data,
  # annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  cluster_rows=FALSE, cluster_cols=FALSE
)
dev.off()

corrected_mixing_matrix = corrected_mixing_matrix[order(corrected_mixing_matrix$annotation),]

pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Rebuttal/Plots/corrected_mixing_matrix_208_per_celltype.pdf", width = 12,height = 10)
par(mar = c(12, 5, 4, 2) )
boxplot(CONSENSUS_IC_208 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_208", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19)
# axis(1, at = 1:nrow(beta_CRC1), labels = corrected_mixing_matrix, las = 2, srt = 30)
dev.off()



pdf("/Users/arkajyotibhattacharya/Library/CloudStorage/GoogleDrive-a.bhattacharya@rug.nl/My\ Drive/Laptop\ backup/Work\ outside\ projects\ 09102022/Work\ outside\ projects/CRC_Dienstman_data/Rebuttal/Plots/corrected_mixing_matrix_EMT_per_celltype.pdf", width = 36,height = 33)
par(mar = c(23, 5, 4, 2) )
par(mfrow = c(3,3))
boxplot(CONSENSUS_IC_208 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_208", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_117 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_117", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_193 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_193", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_202 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_202", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_116 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_116", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_153 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_153", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_58 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_58", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_136 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_136", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
boxplot(CONSENSUS_IC_77 ~ annotation, data = corrected_mixing_matrix, main = "CONSENSUS_IC_77", ylab = "Activity scores",xlab = "", las = 2, col = "lightblue", notch = TRUE, pch = 19, cex.axis = 2, cex.lab = 2, cex.main = 2)
# axis(1, at = 1:nrow(beta_CRC1), labels = corrected_mixing_matrix, las = 2, srt = 30)
dev.off()




