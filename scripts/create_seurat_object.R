#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
})

####################
#   arg parsing.   #
# ##################

#Need to disect the argumments passed from Snakemake: 
    #[1] --snake_outs=/path/to/outs/fastqs_A10/outs
    #[2] --output_rds=/path/to/snake_outs/seurat_objects/fastqs_A10.rds


args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL, required = FALSE) { #function looks in args for an argument /
                                                              #that matches the pattern --<name>=<value> and extracts <value>
  pat <- paste0("^--", name, "=") 
  hit <- grep(pat, args, value = TRUE) #`grep` returns the elements of args that match the pattern
  if (length(hit) == 0) {
    if (required) stop(sprintf("Missing required argument --%s", name), call. = FALSE)
    return(default)
  }
  sub(pat, "", hit[1]) #return the value after the '='
}

snake_outs    <- get_arg("snake_outs",    required = TRUE) #call function 2x
output_rds <- get_arg("output_rds", required = TRUE)

message("snake_outs: ", snake_outs)
message("output_rds: ", output_rds)


#########################
# Paths & sanity checks #
#########################

#Define paths to required input files and check they exist.


h5_path   <- file.path(snake_outs, "filtered_peak_bc_matrix.h5")
meta_path <- file.path(snake_outs, "singlecell.csv")
frag_path <- file.path(snake_outs, "fragments.tsv.gz")
frag_tbi  <- paste0(frag_path, ".tbi")

required <- c(h5_path, meta_path, frag_path)
missing  <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop(sprintf("Missing required file(s):\n- %s", paste(missing, collapse = "\n- ")), call. = FALSE)
}

###########################
# Load matrix & metadata #
###########################

# load the peak-by-cell count matrix and cell-level metadata


message("Reading peak-by-cell matrix (H5)…")
counts <- Read10X_h5(filename = h5_path)

#h5 file can sometimes return a list whose elements are matrices. One of those elements is typically named "Peaks" 
#and that’s the peak-by-cell matrix we want. Downstream in the script we do things like colnames(counts) and pass 
#counts into CreateChromatinAssay which expects counts to be a matrix, not a list. 

if (is.list(counts)) {
  if ("Peaks" %in% names(counts)) {
    counts <- counts[["Peaks"]]
  } else {
    # Fallback: take the first element if there isn't a "Peaks" name
    counts <- counts[[1]]
  }
}

message("Reading cell-level metadata (singlecell.csv)…")
cell_metadata <- read.csv(
  file = meta_path,
  header = TRUE,
  row.names = 1
)

# ensure colnames(counts) and rownames(cell_metadata) are the same set and in the same order
barcodes_mat  <- colnames(counts)
barcodes_meta <- rownames(cell_metadata)

common <- intersect(barcodes_mat, barcodes_meta)
if (length(common) == 0) {
  stop("No overlapping barcodes between matrix and metadata. Check inputs.", call. = FALSE)
}
if (length(common) < length(barcodes_mat)) {
  warning(sprintf("Subsetting to %d barcodes present in both matrix and metadata (of %d in matrix).",
                  length(common), length(barcodes_mat)))
  counts        <- counts[, common, drop = FALSE]
}
if (length(common) < length(barcodes_meta)) {
  cell_metadata <- cell_metadata[common, , drop = FALSE]
}

###########################
# Fragments index (tabix) #
# #########################

#    Creates a .tbi index for the fragments file if missing.
#    Tabix index accelerates random access used by Signac (e.g., TSS enrichment, coverage tracks).

if (!file.exists(frag_tbi)) {
  message("Fragments index (.tbi) not found; creating with Signac::IndexFragments()…")
  # This will create a .tbi in-place; doesn't require object
  Signac::IndexFragments(fragments = frag_path)
}


# ####################################
# Build ChromatinAssay & Seurat object #
# ####################################

#    Create a Signac ChromatinAssay object to store the ATAC-seq data, then wrap it in a Seurat object.  


message("Creating ChromatinAssay…")
chrom_assay <- CreateChromatinAssay(
  counts     = counts,
  sep        = c(":", "-"),
  fragments  = frag_path,
  min.cells  = 10,
  min.features = 200
)

message("Wrapping into Seurat object and attaching metadata…")
seurat_obj <- CreateSeuratObject(
  counts   = chrom_assay,
  assay    = "peaks",
  meta.data = cell_metadata
)
DefaultAssay(seurat_obj) <- "peaks"

# ########################
# Gene annotations (hg38) #
# #########################


 #   Fetch gene annotation data from Ensembl (EnsDb.Hsapiens.v86) and add to the Seurat object.


message("Adding hg38 gene annotations (EnsDb.Hsapiens.v86)…")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"  # e.g., 'chr1'
genome(annotations)          <- "hg38"
Annotation(seurat_obj)       <- annotations

# ##########
# Save Obj #
# ##########


#    Save the Seurat object as an RDS file.


dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
message("Saving Seurat object to: ", output_rds)
saveRDS(seurat_obj, file = output_rds)

message("Done.")
