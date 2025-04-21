options(repos = c(CRAN = "https://cloud.r-project.org"))
options(stringsAsFactors = FALSE)

install.packages(c("hdf5r", "openxlsx", "devtools", "remotes"), dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(
  c("DirichletMultinomial", "TFBSTools", "JASPAR2020"),
  ask = FALSE, force = TRUE
)

remotes::install_github("satijalab/seurat@seurat5")
remotes::install_github("satijalab/seurat-data")
remotes::install_github("satijalab/azimuth")
remotes::install_github("mojaveazure/seurat-disk")



library(Matrix)
library(Seurat)
library(SeuratData)
library(Azimuth)
library(SeuratDisk)
library(hdf5r)
options(timeout = 600)

bm <- RunAzimuth(query = "human_cd34_bone_marrow.h5ad", reference = "bonemarrowref")

bm <- NormalizeData(
  object = bm,
  assay  = "RNA"
)

counts_mat <- GetAssayData(bm, assay="RNA", layer="counts")     
norm_mat   <- GetAssayData(bm, assay="RNA", layer="data")      


meta_df <- bm@meta.data


clean <- CreateSeuratObject(
  counts    = counts_mat,
  meta.data = meta_df,
  project   = bm@project.name
)


clean <- SetAssayData(
  object   = clean,
  assay    = "RNA",
  slot     = "data",
  new.data = norm_mat
)


for(a in setdiff(Assays(clean), "RNA")) clean[[a]] <- NULL

#-------Data Export-------
# 1) raw counts
counts_mat <- GetAssayData(clean, assay = "RNA", layer = "counts")
writeMM(counts_mat, file = "counts.mtx")

# 2) normalized data
norm_mat   <- GetAssayData(clean, assay = "RNA", layer = "data")
writeMM(norm_mat, file = "data.mtx")

# 3) barcodes / cells
write.table(
  colnames(counts_mat),
  file      = "barcodes.tsv",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 4) genes / features
write.table(
  rownames(counts_mat),
  file      = "genes.tsv",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# 5) metadata
write.csv(
  clean@meta.data,
  file      = "meta.csv",
  quote     = FALSE,
  row.names = TRUE
)
#-------Data Export--------