# Load necessary libraries
library(data.table) # For efficient data manipulation
library(Seurat) # For single-cell RNA-seq analysis
library(ggplot2) # For data visualization
library(RColorBrewer) # For color palettes

# Load corrected mixing matrices for CRC1
corrected_mix_mat = data.frame(fread("Projects/Spatial_transcriptomics_data_Marco/Results/Projection_CRC_components/Results/Subset_3/johnson_p_value_matrix.tsv"), row.names = 1)
corrected_mix_mat = -log10(corrected_mix_mat)
rownames(corrected_mix_mat) = gsub("\\.", "-", rownames(corrected_mix_mat))

corrected_mix_mat_1 = data.frame(fread("Projects/Spatial_transcriptomics_data_Marco/Results/Projection_CRC_components/Results/Subset_4/johnson_p_value_matrix.tsv"), row.names = 1)
corrected_mix_mat_1 = -log10(corrected_mix_mat_1)
rownames(corrected_mix_mat_1) = gsub("\\.", "-", rownames(corrected_mix_mat_1))

corrected_mix_mat = as.data.frame(rbind(corrected_mix_mat, corrected_mix_mat_1))

# Load Visium data
visium <- readRDS("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Data/ALL_VISIUM_FILES_MERGED.rds")

# Define indices for significant TCs
EMT_TCs_index = paste("CONSENSUS_IC_", c(208, 117, 193, 202, 116,153, 58, 136, 77), sep = "")
CN_TCs_index = paste("CONSENSUS_IC_", c(149, 189, 148, 55, 131,56, 126, 171, 69, 157, 35, 216, 25, 52, 38, 180, 115, 166), sep = "")

# Function to plot and save RDS files per image
plot_and_save_rds_files_per_image = function(image_id = "CRC1"
                                             ,TCs_to_plot = EMT_TCs_index
                                             ,to_save = TRUE
                                             ,to_plot = TRUE
                                             ,Title_with_location = paste("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_CRC_components/Results/CRC_components_on", sep = "")
){
  # Get mapping file for the image
  mapping_file <- visium@images[[image_id]]@coordinates
  
  # Check for common samples
  if(length(which(rownames(mapping_file)%in%rownames(corrected_mix_mat))) > 100)
  {
    common_samples = intersect(rownames(mapping_file),rownames(corrected_mix_mat) )
    mapping_file = mapping_file[common_samples,]
    corrected_mix_mat_current = corrected_mix_mat[common_samples,]
    
    # Transform corrected mixing matrix
    corrected_mix_mat_transposed = t(corrected_mix_mat_current)
    corrected_mix_mat_transposed_TCs_only = corrected_mix_mat_transposed[TCs_to_plot,]
    
    # Create Seurat object
    new.seurat.object = CreateSeuratObject(counts = corrected_mix_mat_transposed,
                                           assay = "Spatial")
    
    # Add Visium data to Seurat object
    new.seurat.object@images$image = new(
      Class = 'VisiumV1'
      ,
      assay = "spatial"
      ,key = "image_"
      ,coordinates = visium@images[[image_id]]@coordinates
      ,image = visium@images[[image_id]]@image
      ,scale.factors = visium@images[[image_id]]@scale.factors
      ,spot.radius = visium@images[[image_id]]@spot.radius
    )
    
    # Save files if specified
    if(to_save)
    {
      write.table(corrected_mix_mat_transposed_TCs_only, file = paste(paste(Title_with_location, image_id, collapse = "_"), ".txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
      write.table(mapping_file, file = paste(paste(Title_with_location, image_id, collapse = "_"), "mapping_file.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
      # saveRDS(new.seurat.object, file = paste(paste(Title_with_location, image_id, collapse = "_"), ".rds", sep = ""))
    }
    
    # Plot if specified
    if(to_plot)
    {
      saveRDS(new.seurat.object, file = paste(paste(Title_with_location, image_id, collapse = "_"), ".rds", sep = ""))
      pdf(paste(paste(Title_with_location, image_id, collapse = "_"), "_TCs_to_plot.pdf", sep = ""))
      print(SpatialFeaturePlot(new.seurat.object
                               , features =TCs_to_plot
                               , stroke = 0
                               ,image.alpha = 0.5
                               , pt.size.factor = 1.8
                               , alpha = 0.8
                               ,interactive = FALSE)&
        scale_fill_gradient(limits = c(0,max(corrected_mix_mat_current[,gsub("-", "_", TCs_to_plot)])), low = "white", high = "red" ))
      
      dev.off()
    }
  }
}

# Plot and save for CRC1
plot_and_save_rds_files_per_image(
  image_id = "CRC1"
  ,TCs_to_plot = EMT_TCs_index
  ,to_save = TRUE
  ,to_plot = FALSE
  ,Title_with_location = paste("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_CRC_components/Results/CRC_components_on", sep = "")
)

# Plot and save for CRC2
plot_and_save_rds_files_per_image(
  image_id = "CRC2"
  ,TCs_to_plot = EMT_TCs_index
  ,to_save = TRUE
  ,to_plot = FALSE
  ,Title_with_location = paste("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_CRC_components/Results/CRC_components_on", sep = "")
)

# Plot and save for CRC3
plot_and_save_rds_files_per_image(
  image_id = "CRC3"
  ,TCs_to_plot = EMT_TCs_index
  ,to_save = TRUE
  ,to_plot = FALSE
  ,Title_with_location = paste("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Results/Projection_CRC_components/Results/CRC_components_on", sep = "")
)

# Additional plot for CRC1
pdf(paste(paste(Title_with_location, image_id, collapse = "_"), "_TCs_to_plot.pdf", sep = ""), height = 18, width = 18)
SpatialFeaturePlot(new.seurat.object
                   , features =TCs_to_plot
                   , stroke = 0
                   ,image.alpha = 0.5
                   , pt.size.factor = 1.8
                   , alpha = 0.8
                   ,interactive = FALSE)&
  scale_fill_gradient(limits = c(0,max(corrected_mix_mat_current[,gsub("-", "_", TCs_to_plot)])), low = "white", high = "red" )
dev.off()

# Plot heatmaps for CRC1 CN TCs
for(TC in CN_TCs_index)
{
  time1 = proc.time()[3]
  heatmap_data = matrix(NA, max(mapping_file$row)+1, max(mapping_file$col)+1)
  rownames_corrected_mix_mat = rownames(corrected_mix_mat)
  for(i in 1:dim(corrected_mix_mat)[1])
  {
    mapping_index = which(rownames(mapping_file)==rownames_corrected_mix_mat[i])
    
    heatmap_data[mapping_file$row[mapping_index]+1, mapping_file$col[mapping_index]+1] = corrected_mix_mat[i,TC]
    
  }
  
  pdf(paste0("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Plots/CRC_samples/CRC_1_CN_TCs_",TC,".pdf"))
  image(heatmap_data,useRaster=TRUE)
  dev.off()
  
  print(TC)
  print((proc.time()[3] - time1)/60)
}

# Plot heatmaps for CRC1 EMT TCs
for(TC in EMT_TCs_index)
{
  time1 = proc.time()[3]
  heatmap_data = matrix(NA, max(mapping_file$row)+1, max(mapping_file$col)+1)
  rownames_corrected_mix_mat = rownames(corrected_mix_mat)
  for(i in 1:dim(corrected_mix_mat)[1])
  {
    mapping_index = which(rownames(mapping_file)==rownames_corrected_mix_mat[i])
    
    heatmap_data[mapping_file$row[mapping_index]+1, mapping_file$col[mapping_index]+1] = corrected_mix_mat[i,TC]
    
  }
  
  pdf(paste0("/home/arkajyotibhattacharya/Projects/Spatial_transcriptomics_data_Marco/Plots/CRC_samples/CRC_1_EMT_TCs_",TC,".pdf"))
  image(heatmap_data,useRaster=TRUE)
  dev.off()
  
  print(TC)
  print((proc.time()[3] - time1)/60)
}

# Continue for CRC2 and CRC3...
