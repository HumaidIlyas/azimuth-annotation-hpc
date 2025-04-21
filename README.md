# Assessment - Azimuth Mapping of Human CD34 $+$  Bone Marrow

## 1. Introduction

This document describes the pipeline used to map a human CD34 $+$  bone marrow single-cell RNA-seq dataset onto the Azimuth bone-marrow reference, and to reassemble the results into an AnnData (.h5ad) fle for downstream analysis.

## 2. Steps Taken

The entire project was done on the BigRed200 HPC Firstly, Data is downloaded from the HuBMAP. After that, in the r file, the h5ad file is directly read when Azimuth is run.

### 2.1 Piece-wise Export & Reassembly

After Azimuth is done running, the rds object, which we are calling our rds file (bm) is converted into pieces and is exported to Python, where these pieces are reassembled

#### 2.1.1 Export from Seurat
The rds file is broken down into 5 pieces 
![Figure 1: Export from R](r_out.jpg)

#### 2.1.2 Reassemble in Python
Which is reassembled in python for analysis
![Figure 2: Combine in Python](py_in.jpg)

### 2.2 Analyzing the Data

After all the processes are done, the data is analysed by having the number of cells per cell type and visualizing the distribution of cells per cell type.

## 3. Issues Encountered

### 3.1 Insufficient Memory on Local Machine

Problem: Running RunAzimuth() on the local laptop resulted in out-of-memory errors due to the large size of the bone marrow dataset and only 16 GB of RAM available.
![Figure 3: Combine in Python](py_in.jpg)
Solution: Switched to the BigRed200 HPC cluster, providing sufficient memory to successfully complete the Azimuth run.

### 3.2 Persistent HDF5 Conversion Error

Problem: After obtaining the annotated Seurat object, every attempt to convert it back to .h5ad with Convert() failed with the following HDF5 API error:
![Figure 4: Error Exporting h5ad Directly](err.jpg)

#### Attempts to Diagnose:

• Injected metadata via Excel export/import.
• Tried switching between slots (data vs. counts) and specifying layers explicitly.
• Downgraded to Seurat v4 (Assay class) and forced conversion of Assay5 to Assay.
• Experimented with writing out Matrix Market (.mtx) dumps for counts and normalized data.

None of these approaches removed the HDF5 copy error.

#### Final Workaround

##### 1) Dump the two key matrices and metadata to disk in portable formats:

• Raw counts: counts.mtx
• Normalized data: data.mtx
• Gene list: genes.tsv
• Cell barcodes: barcodes.tsv
• Metadata: meta.csv

##### 2) Reconstruct the annotated AnnData object entirely in Python:

The manual reassembly bypassed the faulty HDF5 copy step and produced a valid .h5ad

## 3. Results

The final annotated h5ad contains 20 distinct cell types. The most abundant populations were GMP (approximately 1065 cells), LMPP (1138), and HSC (1028), while rare populations such as MAIT and Stromal were each represented by a single cell.


![Figure 5: Distribution of Cells per Cell Type](Distribution.jpg)


• Initial Azimuth mapping succeeded on HPC after local memory limits.
• Direct SeuratDisk conversion to .h5ad was blocked by an HDF5 copy bug.
• Workaround: export assay matrices and metadata separately, then reassemble in Python.
• The resulting AnnData file is fully compatible with Scanpy or Anndata for downstream analyses.

