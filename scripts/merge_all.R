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
# Arg parsing (same style as merge_timepoints)
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
  input <- get_arg("inputs", required = TRUE)
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

# -----------------
# Load objects & stamp replicate/timepoint from filenames
# (identical pattern to merge_timepoints)
# -----------------
seurat_list <- lapply(inputs_vec, readRDS)
names(seurat_list) <- sub("(_filtered_cells)?\\.rds$", "", basename(inputs_vec))
message("Loaded ", length(seurat_list), " Seurat objects.")

rep_labels <- names(seurat_list)
tp_from <- function(x) sub("^([A-Za-z0-9]+)_.*$", "\\1", x)

for (lab in rep_labels) {
  o <- seurat_list[[lab]]
  o$replicate  <- lab
  o$timepoint  <- tp_from(lab)
  o$source_rds <- basename(inputs_vec[which(names(seurat_list) == lab)])
  seurat_list[[lab]] <- o
}
message("Stamped meta: 'replicate', 'timepoint', 'source_rds' on each object.")

# ---------------------------
# Build a shared peak set across ALL inputs (same style as merge_timepoints)
# ---------------------------
peaks_list <- lapply(seurat_list, function(o) {
  gr <- granges(o)
  seqlevelsStyle(gr) <- "UCSC"
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- gr[width(gr) > 0]
  gr <- unique(gr)
  sort(gr)
})

stopifnot(all(vapply(peaks_list, inherits, logical(1), what = "GRanges")))
message("Extracted peak sets: ",
        paste(vapply(peaks_list, length, integer(1)), collapse = ", "),
        " peaks per replicate.")

all_peaks <- unlist(GenomicRanges::GRangesList(peaks_list), use.names = FALSE)

shared_peaks <- GenomicRanges::reduce(
  all_peaks,
  min.gapwidth  = 50,
  ignore.strand = TRUE
)
seqlevelsStyle(shared_peaks) <- "UCSC"
shared_peaks <- keepStandardChromosomes(shared_peaks, pruning.mode = "coarse")
shared_peaks <- sort(shared_peaks)
message("Shared peak set built: ", length(shared_peaks), " peaks.")

# ---------------------------
# Re-quantify counts in shared peaks (same logic as merge_timepoints)
# ---------------------------
recounted <- vector("list", length(seurat_list))
names(recounted) <- names(seurat_list)

for (i in seq_along(seurat_list)) {
  o   <- seurat_list[[i]]
  lab <- names(seurat_list)[i]

  frags <- Fragments(o)
  stopifnot(length(frags) >= 1)
  frag_obj <- frags[[1]]
  stopifnot(inherits(frag_obj, "Fragment"))

  cells_vec <- colnames(o)
  stopifnot(length(cells_vec) > 0)

  # Build new peak x cell matrix on shared peaks
  mtx <- FeatureMatrix(
    fragments = list(frag_obj),
    features  = shared_peaks,
    cells     = cells_vec
  )

  # Carry over metadata for those cells, matched to mtx columns
  md <- o@meta.data[match(colnames(mtx), rownames(o@meta.data)), , drop = FALSE]

  # ---- tiny guard to keep Seurat happy ----
  # Ensure rownames(meta) == colnames(counts) exactly (same order & values)
  if (nrow(md) != ncol(mtx)) {
    # If match introduced NAs (cells dropped), subset both by non-NA rows
    ok <- !is.na(rownames(md))
    md <- md[ok, , drop = FALSE]
    mtx <- mtx[, ok, drop = FALSE]
  }
  rownames(md) <- colnames(mtx)
  stopifnot(identical(rownames(md), colnames(mtx)))
  # -----------------------------------------

  anno <- tryCatch(Annotation(o), error = function(e) NULL)

  assay <- CreateChromatinAssay(
    counts     = mtx,
    fragments  = frag_obj,
    annotation = anno
  )

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

  # quick diag
  message(sprintf("[recount] %s: %d cells after recount", lab, ncol(so)))

  recounted[[i]] <- so
}
message("Recounted all replicates against the shared peak set.")

# ---------------------------
# Merge all replicates into one object (same as merge_timepoints)
# ---------------------------
if (length(recounted) == 1L) {
  merged <- recounted[[1]]
} else {
  merged <- merge(
    x = recounted[[1]],
    y = recounted[-1],
    add.cell.ids = names(recounted)
  )
}
DefaultAssay(merged) <- "peaks"

# TF-IDF + LSI + Harmony + UMAP (same as merge_timepoints)
merged <- FindTopFeatures(merged, min.cutoff = "q0")
merged <- RunTFIDF(merged)
merged <- RunSVD(merged)

lsi <- Embeddings(merged, "lsi")
md  <- merged@meta.data

H <- harmony::HarmonyMatrix(
  data_mat          = lsi,
  meta_data         = md,
  vars_use          = "replicate",
  do_pca            = FALSE,
  max.iter.harmony  = 20
)

merged[["harmony"]] <- CreateDimReducObject(
  embeddings = H,
  key        = "harmony_",
  assay      = DefaultAssay(merged)
)

set.seed(1)
merged <- RunUMAP(merged, reduction = "harmony", dims = 2:ncol(H))

saveRDS(merged, file = output_rds)
message("Saved merged object to: ", output_rds)
