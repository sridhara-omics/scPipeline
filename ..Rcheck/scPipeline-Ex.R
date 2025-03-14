pkgname <- "scPipeline"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "scPipeline-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('scPipeline')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ConvertGeneIdentifiers")
### * ConvertGeneIdentifiers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ConvertGeneIdentifiers
### Title: Convert Gene Identifiers in a Seurat Object
### Aliases: ConvertGeneIdentifiers

### ** Examples

## No test: 
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
library(Seurat)
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
# Convert RefSeq IDs to gene symbols
seurat_obj_converted <- ConvertGeneIdentifiers(
  seurat_obj,
  id_type = "refseq",
  to_id_type = "symbol"
)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ConvertGeneIdentifiers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReactomeData")
### * ReactomeData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReactomeData
### Title: Reactome Data Analysis for Seurat Object
### Aliases: ReactomeData

### ** Examples

## No test: 
library(Seurat)
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
# Reactome Analysis
seurat_reactome <- ReactomeData(seurat_obj)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReactomeData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SeuratLowDim")
### * SeuratLowDim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SeuratLowDim
### Title: Create a Low dimensional Seurat object from scaled seurat object
### Aliases: SeuratLowDim

### ** Examples

## No test: 
library(Seurat)
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SeuratLowDim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SeuratMarkers")
### * SeuratMarkers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SeuratMarkers
### Title: A thresholded markers list for better calculation of DE genes
### Aliases: SeuratMarkers

### ** Examples

## No test: 
library(Seurat)
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
seurat_obj <- SeuratLowDim(seurat_obj)
# Create Markers list
seurat_markers <- SeuratMarkers(seurat_obj)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SeuratMarkers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SeuratPreprocess")
### * SeuratPreprocess

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SeuratPreprocess
### Title: Preprocess count data and create a Seurat object
### Aliases: SeuratPreprocess

### ** Examples

## No test: 
library(Seurat)
# Read 10X counts data from matrix.mtx, barcodes.tsv and genes.tsv
counts <- Read10X(data.dir = "../inst/extdata", gene.column = 1)

# Create Seurat object without batch correction
seurat_obj <- SeuratPreprocess(counts)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SeuratPreprocess", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
