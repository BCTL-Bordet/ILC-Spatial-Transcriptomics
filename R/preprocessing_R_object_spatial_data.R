#### Load Required Libraries ####
library(Seurat)
library(data.table)
library(STutility)
library(dplyr)
library(stringr)

#### Define Sample IDs ####
samples_ids <- c("1", "3", "5", "7", seq(33,48), seq(50,72))
list_samples <- vector("list", length(samples_ids))

#### Process Each Sample ####
for (i in seq_along(samples_ids)) {
  sample_id <- samples_ids[i]
  
  # Define file names. All files are outputs of spaceranger
  samples <- paste0("ST", sample_id, "/outs/filtered_feature_bc_matrix.h5")
  spotfiles <- paste0("ST", sample_id, "/outs/spatial/tissue_positions_list.csv")
  imgs <- paste0("ST", sample_id, "/outs/spatial/tissue_hires_image.png")
  json <- paste0("ST", sample_id, "/outs/spatial/scalefactors_json.json")
  
  # Create an info table
  infoTable <- data.frame(samples, spotfiles, imgs, json, stringsAsFactors = FALSE)
  
  # Load data into STutility object
  object_st <- InputFromTable(
    infotable = infoTable, 
    platform = "Visium", 
    min.spot.feature.count = 200, 
    min.gene.count = 1, 
    type = "2"
  )
  
  object_st@meta.data$orig.ident <- sample_id
  
  # Load spatial annotations. tissue_positions_list_annotation.csv was obtained using the python script "spatial_annotation.py", where the morphological annotation was transferred to each ST spot
  spatial_locations_annotated <- fread(paste0("ST", sample_id, "/outs/spatial/tissue_positions_list_annotation.csv"))
  colnames(spatial_locations_annotated) <- c('barcode', 'in_tissue', 'array_row', 'array_col', 'row_pxl', 'col_pxl', 
                                             "Tumor", "Necrosis", "Fat_tissue", "High_TILs_stroma", "Cellular_stroma", 
                                             "Acellular_stroma", "Vessels", "Artefact", "Canal_galactophore", "Nodule_lymphoid", 
                                             "In_situ", "Nerve", "Lymphocyte", "Hole", "Microcalcification", "Out", "Apocrine_metaplasia")
  
  # Match annotations with barcodes
  annotations <- spatial_locations_annotated[spatial_locations_annotated$barcode %in% str_split_fixed(rownames(object_st@meta.data), "_", 2)[,1], 
                                             c("Tumor", "Necrosis", "Fat_tissue", "High_TILs_stroma", "Cellular_stroma", 
                                               "Acellular_stroma", "Vessels", "Artefact", "Canal_galactophore", "Nodule_lymphoid", 
                                               "In_situ", "Nerve", "Lymphocyte", "Hole", "Microcalcification", "Out", "Apocrine_metaplasia")]
  
  object_st@meta.data <- cbind(object_st@meta.data, annotations)
  list_samples[[i]] <- object_st
}

#### Preprocessing and Feature Selection ####
list_var_features <- list()
for (i in seq_along(list_samples)) {
  # Remove artefact spots
  list_samples[[i]] <- SubsetSTData(list_samples[[i]], spots = rownames(list_samples[[i]]@meta.data[list_samples[[i]]@meta.data$Hole < 0.3, ]))
  list_samples[[i]] <- SubsetSTData(list_samples[[i]], spots = rownames(list_samples[[i]]@meta.data[list_samples[[i]]@meta.data$Artefact < 0.3, ]))
  list_samples[[i]] <- SubsetSTData(list_samples[[i]], spots = rownames(list_samples[[i]]@meta.data[list_samples[[i]]@meta.data$Out < 0.3, ]))
  
  # Remove unwanted genes
  RPL.genes <- rownames(list_samples[[i]])[grepl("RPL", rownames(list_samples[[i]]))]
  RPS.genes <- rownames(list_samples[[i]])[grepl("RPS", rownames(list_samples[[i]]))]
  MT.genes <- rownames(list_samples[[i]])[grepl("MT-", rownames(list_samples[[i]]))]
  MTRNR.genes <- rownames(list_samples[[i]])[grepl("MTRNR", rownames(list_samples[[i]]))]
  metabolic.genes <- unique(c(RPL.genes, RPS.genes, MT.genes, MTRNR.genes))
  list_samples[[i]] <- list_samples[[i]][!(rownames(list_samples[[i]]) %in% metabolic.genes), ]
  
  # Load and process images
  list_samples[[i]] <- LoadImages(list_samples[[i]], verbose = TRUE, time.resolve = FALSE)
  ImagePlot(list_samples[[i]], method = "raster", darken = TRUE, type = "raw")
  list_samples[[i]] <- MaskImages(list_samples[[i]])
  
  # Normalize and select variable features
  list_samples[[i]] <- SCTransform(list_samples[[i]])
  list_samples[[i]] <- FindVariableFeatures(list_samples[[i]], selection.method = "vst", nfeatures = 2000)
  list_var_features[[i]] <- as.vector(list_samples[[i]]@assays$SCT@var.features)
}

#### Merge All Processed Objects ####
merged_object <- MergeSTData(
  x = list_samples[[1]], 
  y = list_samples[-1], 
  merge.data = TRUE, 
  add.spot.ids = paste0("sample", samples_ids)
)
