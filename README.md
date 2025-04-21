# Assessment - Azimuth Mapping of Human CD34 $+$  Bone Marrow

## Humaid Ilyas

Monday 21 $st$  April, 2025

## 1. Introduction

This document describes the pipeline used to map a human CD34 $+$  bone marrow single-cell RNA-seq dataset onto the Azimuth bone-marrow reference, and to reassemble the results into an AnnData (.h5ad) fle for downstream analysis.

## 1. Steps Taken

The entire project was done on the BigRed200 HPC Firstly, Data is downloaded from the HuBMAP After that, in the r fle the h5ad fle is directly read when Azimuth is ran.

### 1.1 Piece-wise Export & Reassembly

After Azimuth is done running, the rds object which we are calling bm is converted into pieces and is exported to python in the python frontm these peices are reassembled

#### 1.1.1 Export from Seurat

bm &lt;- readRDS("bm_azimuth_results.rds")

counts_mat $<-$  GetAssayData(bm, assay="RNA", $\text {slot="counts")}$ 

norm_mat $<-$  GetAssayData(bm, $\text {assay="RNA"}$ , $\text {slot="data")}$ 

writeMM(counts_mat, "counts.mtx")

writeMM(norm_mat,"data.mtx")

write.table(rownames(norm_mat), "genes.tsv", $\text {quote=FALSE}$ , $row.names=FALSE,$  col.names $=$ FALSE)

write.table(colnames(norm_mat), "barcodes.tsv" $,quote=FALSE$ , $row.names=FALSE,$  col.names=FALSE)

write.csv(bm@meta.data, "meta.csv", quot $e=FALSE)$ 

#### 1.1.2 Reassemble in Python

from scipy import io

import pandas as pd

import anndata as ad

counts $=$  io.mmread("counts.mtx").T.tocsr()

norm $=$ io.mmread("data.mtx").T.tocsr()

genes $=pd.read_csv("genes.tsv",$  header $=$ None)[0].astype(str)

barcodes $=$  pd.read_csv("barcodes.tsv", header=None)[0].astype(str)

meta $=pd.read_csv("meta.csv"$ , index_col $=0$ )

adata $=$  ad.AnnData(

$\mathrm {X}=\text {norm,}$ 

$obs=meta,$ 

$\text {var=pd.}$ DataFrame(index $=$ genes)

)

adata.raw $=$  ad.AnnData $(\mathrm {X}=\text {counts}$ , $\text {var=adata.var}$ , obs=adata.obs)

adata.layers["counts"] $=$  counts

adata.write("annotated_expr_manual.h5ad")

### 1.2 Anlyzing the Data

After all the processes are done, the data is analysed by having the number of cells per cell type and visualizing the distribution of cells per cell type.

## 2. Issues Encountered

### 2.1 Insufcient Memory on Local Machine

Problem: Running RunAzimuth() on the local laptop resulted in out-of-memory errors due to the large size of the bone marrow dataset and only 16 GB of RAM available.

Solution: Switched to the BigRed200 HPC cluster, which provided sufcient memory to complete the Azimuth run successfully.

### 2.2 Persistent HDF5 Conversion Error

Problem: After obtaining the annotated Seurat object, every attempt to convert it back to .h5ad with Convert() failed with the following HDF5 API error:

Error in assay.group&#36;obj_copy_to(dst_loc $=dfile$ $,dst_name="X"$ $,src_name=x.data):$ 

HDF5-API Errors:

error #000: H5O.c in H5Ocopy(): ... Unable to copy object

...

error #005: Source object not found

Calls: Convert ... Convert.H5Seurat -&gt; H5SeuratToH5AD -&gt; &lt;Anonymous&gt; -&gt; .Call

#### Attempts to Diagnose:

• Injected metadata via Excel export/import.

• Tried switching between slots (data vs. counts) and specifying layers explicitly.

• Downgraded to Seurat v4 (Assay class) and forced conversion of Assay5 to Assay.

• Experimented with writing out Matrix Market (.mtx) dumps for counts and normalized data.

None of these approaches removed the HDF5 copy error.

#### Final Workaround

1. Dump the two key matrices and metadata to disk in portable formats:

• Raw counts: counts.mtx

• Normalized data: data.mtx

• Gene list: genes.tsv

• Cell barcodes: barcodes.tsv

• Metadata: meta.csv

2. Reconstruct the annotated AnnData object entirely in Python:

import scipy.io as io

import pandas as pd

import anndata as ad

counts $=$  io.mmread("counts.mtx").T.tocsr()

norm $.$ 

genes $=$  pd.read_csv("genes.tsv", header $=$ None)[0].astype(str)

barcodes $=$  pd.read_csv("barcodes.tsv", header $=$ None)[0].astype(str)

meta $=$  pd.read_csv("meta.csv", index_col $=0$ )

adata $=$  ad.AnnData $(\mathrm {X}=\text {norm},$  laye $s=\{"$ counts": counts}, obs $=$ meta, $\text {var=pd.}$ DataFrame(index $=$ genes))adata.write("annotated_expr_manual.h5ad")

\end{minted}

\item This manual reassembly bypassed the faulty HDF5 copy step and produced a valid \texttt{.h5ad

## 3. Results

The fnal annotated expr manual.h5ad contains 20 distinct cell types. The most abundant popula tions were GMP (approximately 1065 cells), LMPP (1138), and HSC (1028), while rare populations such as MAIT and Stromal were each represented by a single cell.


![Figure 1: Distribution of Cells per Cell Type](Distribution.jpg)



• Initial Azimuth mapping succeeded on HPC after local memory limits.

• Direct SeuratDisk conversion to .h5ad was blocked by an HDF5 copy bug.

• Workaround: export assay matrices and metadata separately, then reassemble in Python.

• The resulting AnnData fle is fully compatible with Scanpy or Anndata for downstream analyses.

