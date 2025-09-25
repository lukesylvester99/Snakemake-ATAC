#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(harmony)
})

# -----------------
# Arg parsing
# -----------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL, required = FALSE) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0) {
    if (required) stop(sprintf("Missing required argument --%s", name), call. = FALSE)
    return(default)
  }
  sub(pat, "", hit[1])
}

# Prefer many --input_rds=... flags; fall back to single --inputs="a,b c"
input_flags <- grep("^--input_rds=", args, value = TRUE)
if (length(input_flags) > 0) {
  inputs_vec <- sub("^--input_rds=", "", input_flags)
} else {
  input <- get_arg("inputs", required = TRUE)  # only required if no flags present
  inputs_vec <- unlist(strsplit(input, "[,\\s]+"))
  inputs_vec <- inputs_vec[nzchar(inputs_vec)]
}

output_rds <- get_arg("output_rds", required = TRUE)
message("inputs: ", paste(inputs_vec, collapse = " "))
message("output_rds: ", output_rds)

# Basic validation
if (length(inputs_vec) < 2) {
  stop("Need at least two input RDS files for merging. Got: ", length(inputs_vec), call. = FALSE)
}
missing <- inputs_vec[!file.exists(inputs_vec)]
if (length(missing)) {
  stop("Missing input files:\n  - ", paste(missing, collapse = "\n  - "), call. = FALSE)
}

# Load objects
seurat_list <- lapply(inputs_vec, readRDS)
names(seurat_list) <- sub("(_filtered_cells)?\\.rds$", "", basename(inputs_vec))
message("Loaded ", length(seurat_list), " Seurat objects.")


#stamp the replicate name into the cells

rep_labels <- names(seurat_list)

# helper: derive timepoint (e.g., "S1_1" -> "S1")
tp_from <- function(x) sub("^([A-Za-z0-9]+)_.*$", "\\1", x)

for (lab in rep_labels) {
  o <- seurat_list[[lab]]
  o$replicate <- lab
  o$timepoint <- tp_from(lab)
  o$source_rds <- basename(inputs_vec[which(names(seurat_list) == lab)])
  seurat_list[[lab]] <- o
}

message("Stamped meta: 'replicate' and 'timepoint' (and 'source_rds') on each object.") 


# ---------------------------
# create shared peak sets across replicates
# ---------------------------

# Extract per-replicate peak sets
peaks_list <- lapply(seurat_list, function(o) {
  gr <- granges(o) #pull peak coordinates
  seqlevelsStyle(gr) <- "UCSC"                 
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")  # drops decoys/alt contigs
  gr <- gr[width(gr) > 0] # drop any zero/negative widths
  gr <- unique(gr) #de-duplicate
  sort(gr)
})

stopifnot(all(vapply(peaks_list, inherits, logical(1), what = "GRanges"))) #must be granges for reduce() below
message("Extracted peak sets from each object: ",
        paste(vapply(peaks_list, length, integer(1)), collapse = ", "),
        " peaks per replicate.")

# Concatenate as GRangesList, then flatten to a single GRanges
all_peaks <- unlist(GenomicRanges::GRangesList(peaks_list), use.names = FALSE)

#find shared peaks across all replicates (within 50bp)
shared_peaks <- GenomicRanges::reduce(
  all_peaks,
  min.gapwidth  = 50,  # merge ranges that are ≤50 bp apart, bridging tiny gaps
  ignore.strand = TRUE # strand-agnostic merging
)

seqlevelsStyle(shared_peaks) <- "UCSC"
shared_peaks <- keepStandardChromosomes(shared_peaks, pruning.mode = "coarse") #drop any nonstandard chromosomes
shared_peaks <- sort(shared_peaks)

message("Shared peak set built: ", length(shared_peaks), " peaks.")

#_---------------------------
# Re-quantify counts in shared peaks
# ---------------------------   

recounted <- vector("list", length(seurat_list)) #pre-allocate a list to hold the recounted objects
names(recounted) <- names(seurat_list)

for (i in seq_along(seurat_list)) {
  #pull the object into o and its label into lab
  o   <- seurat_list[[i]]
  lab <- names(seurat_list)[i]

  frags <- Fragments(o)
  stopifnot(length(frags) >= 1)
  frag_obj <- frags[[1]]
  stopifnot(inherits(frag_obj, "Fragment"))

  # collect the cell barcodes present
  cells_vec <- colnames(o)
  stopifnot(length(cells_vec) > 0)

  # build new peak x cell matrix on shared_peaks
  mtx <- FeatureMatrix(
    fragments = list(frag_obj),
    features  = shared_peaks,
    cells     = cells_vec
  )

  # reuse annotation if present; it’s fine if NULL
  anno <- tryCatch(Annotation(o), error = function(e) NULL)

  assay <- CreateChromatinAssay(
    counts    = mtx,
    fragments = frag_obj,
    annotation = anno
  )

  # carry over original metadata for those cells
  md <- o@meta.data[match(colnames(mtx), rownames(o@meta.data)), , drop = FALSE]

  #build a fresh Seurat object from the new chromatin assay
  so <- CreateSeuratObject(
    counts    = assay,
    assay     = "peaks",
    meta.data = md,
    project   = paste0(o$timepoint[1], "_recount")
  )
  DefaultAssay(so) <- "peaks"

  # (re)stamp to be safe
  so$replicate <- o$replicate
  so$timepoint <- o$timepoint

  recounted[[i]] <- so #store the recount result for this replicate
}
message("Recounted all replicates against the shared peak set.")

# ---------------------------
# Merge all replicates into one object
# ---------------------------   

# merge any number of replicates
if (length(recounted) == 1L) {
  merged <- recounted[[1]] #If there’s only one replicate, nothing to merge
} else {
  merged <- merge(
    x = recounted[[1]],
    y = recounted[-1], #Everything except the first becomes the y list
    add.cell.ids = names(recounted) #prepend IDs to every cell barcode (T1_rep1_AAAC...)
  )
}
DefaultAssay(merged) <- "peaks"

# TF-IDF + LSI on merged counts
merged <- FindTopFeatures(merged, min.cutoff = "q0")
merged <- RunTFIDF(merged)
merged <- RunSVD(merged)

# LSI matrix (cells x dims)
lsi <- Embeddings(merged, "lsi")
md  <- merged@meta.data

# Run Harmony on LSI
H <- harmony::HarmonyMatrix(
  data_mat          = lsi,
  meta_data         = md,
  vars_use          = "replicate",
  do_pca            = FALSE,
  max.iter.harmony  = 20
)

# Attach as Seurat reduction and continue
merged[["harmony"]] <- CreateDimReducObject(
  embeddings = H,
  key        = "harmony_",
  assay      = DefaultAssay(merged)
)

set.seed(1)
# UMAP on the Harmony space
merged <- RunUMAP(merged, reduction = "harmony", dims = 2:ncol(H))

#save the merged object
saveRDS(merged, file = output_rds)
message("Saved merged object to: ", output_rds)
