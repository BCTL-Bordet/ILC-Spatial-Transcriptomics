#### Load Required Libraries ####
library(STutility)
library(sjmisc)

#### Load Merged Object ####
merged_object <- readRDS("object_merged_clusters_different_res.RDS")

#### Function to Compute Co-occurrence Metrics ####
compute_co_occurrence <- function(data) {
  
  # Create an empty data frame to store results
  co_occurrence_results <- data.frame(matrix(ncol = 13))
  colnames(co_occurrence_results) <- c("Sample", "Fat_tissue", "High_TILs_stroma", "Cellular_stroma", 
                                       "Acellular_stroma", "Vessels", "Canal_galactophore", "In_situ", 
                                       "Nerve", "Lymphocyte", "Microcalcification", "Apocrine metaplasia", "Immune")
  
  # Iterate over each sample
  for (sample in unique(data$orig.ident)) {
    
    sample_data <- data[data$orig.ident == sample, ]
    
    # Ensure there are tumor-containing spots
    if (sum(sample_data$Tumor) > 0) {
      
      sample_data <- sample_data[sample_data$Tumor > 0, ]  # Keep only tumor-containing spots
      total_tumor_spots <- nrow(sample_data)  # Total tumor spots in the sample
      
      # Extract relevant tissue types for co-occurrence analysis
      sample_annotations <- sample_data[, c("Fat_tissue", "High_TILs_stroma", "Cellular_stroma", 
                                            "Acellular_stroma", "Vessels", "Canal_galactophore", "In_situ", 
                                            "Nerve", "Lymphocyte", "Microcalcification", "Apocrine metaplasia", "Immune")]
      
      # Convert presence to binary format
      sample_binary <- dicho(sample_annotations, dich.by = 0, append = FALSE, as.num = TRUE)
      
      # Compute normalized co-occurrence score
      co_occurrence_scores <- colSums(sample_binary) / total_tumor_spots
      
      # Store results
      co_occurrence_results <- rbind(co_occurrence_results, c(sample, co_occurrence_scores))
    }
  }
  
  # Remove NA values and set row names
  co_occurrence_results <- na.omit(co_occurrence_results)
  rownames(co_occurrence_results) <- co_occurrence_results$Sample
  co_occurrence_results$Sample <- NULL
  
  return(co_occurrence_results)
}

#### Compute Co-occurrence Metrics ####
co_occurrence_scores <- compute_co_occurrence(merged_object@meta.data)

