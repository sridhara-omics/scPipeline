pkgname <- "scPipeline"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('scPipeline')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ConvertGeneIdentifiers")
### * ConvertGeneIdentifiers

flush(stderr()); flush(stdout())

### Name: ConvertGeneIdentifiers
### Title: Convert Gene Identifiers in a Seurat Object
### Aliases: ConvertGeneIdentifiers

### ** Examples

## No test: 
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
# Convert RefSeq IDs to gene symbols
seurat_obj_converted <- ConvertGeneIdentifiers(seurat_obj, id_type = "refseq", \
                                               to_id_type = "symbol")
## End(No test)



cleanEx()
nameEx("ReactomeData")
### * ReactomeData

flush(stderr()); flush(stdout())

### Name: ReactomeData
### Title: Reactome Data Analysis for Seurat Object
### Aliases: ReactomeData

### ** Examples

## No test: 
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
# Reactome Analysis
seurat_reactome <- ReactomeData(seurat_obj)
## End(No test)



cleanEx()
nameEx("SeuratLowDim")
### * SeuratLowDim

flush(stderr()); flush(stdout())

### Name: SeuratLowDim
### Title: Create a Low dimensional Seurat object from scaled seurat object
### Aliases: SeuratLowDim

### ** Examples

## No test: 
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
## End(No test)



cleanEx()
nameEx("SeuratMarkers")
### * SeuratMarkers

flush(stderr()); flush(stdout())

### Name: SeuratMarkers
### Title: A thresholded markers list for better calculation of DE genes
### Aliases: SeuratMarkers

### ** Examples

## No test: 
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
# Create Markers list
seurat_markers <- SeuratMarkers(seurat_obj)
## End(No test)



cleanEx()
nameEx("SeuratPreprocess")
### * SeuratPreprocess

flush(stderr()); flush(stdout())

### Name: SeuratPreprocess
### Title: Preprocess count data and create a Seurat object
### Aliases: SeuratPreprocess

### ** Examples

## No test: 
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
## End(No test)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
