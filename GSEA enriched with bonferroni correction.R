# Load required library
library(data.table)

# First, import the enrichment score matrix of interest and TCs of interest

# Calculate Spearman correlation for genesets and samples
spearman_cor_genesets = cor(t(enrichment_score_matrix), method = "spearman")
spearman_cor_samples = cor(enrichment_score_matrix, method = "spearman")

# Check for row and column consistency
sum(rownames(enrichment_score_matrix) == rownames(spearman_cor_genesets))
sum(colnames(enrichment_score_matrix) == colnames(spearman_cor_samples))

# Perform hierarchical clustering for genesets and samples
hclust_genesets = hclust(as.dist(1 - spearman_cor_genesets), method = "ward.D2")
hclust_sample_types = hclust(as.dist(1 - spearman_cor_samples), method = "ward.D2")

# Reorder enrichment score matrix based on geneset clustering
enrichment_score_matrix = enrichment_score_matrix[hclust_genesets$order, ]

# If Bonferroni correction is wanted, continue to this part of the script:
high_enriched_genesets = NULL
for (i in 1:ncol(enrichment_score_matrix)) {
  high_enriched_genesets = union(high_enriched_genesets, rownames(enrichment_score_matrix)[which(abs(enrichment_score_matrix[,i]) > -log10(0.01/N))]) # Adjust N to data
}
enrichment_score_matrix = enrichment_score_matrix[high_enriched_genesets, ]
