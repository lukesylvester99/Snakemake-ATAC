#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  library(grid)   # for optional table/text pages if needed later
})

# -----------------------------
# Global plot theme & settings
# -----------------------------
theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
      axis.title = element_text(size = 12),
      plot.margin = margin(10, 12, 10, 12),
      panel.grid.minor = element_blank()
    )
)

# QC threshold constants (easy to tweak)
CUT_NCOUNT_MIN <- 1000
CUT_NCOUNT_MAX <- 100000
CUT_TSS_MIN    <- 4
CUT_NS_MAX     <- 4
CUT_FRIP_MIN   <- 20       # percent
CUT_BL_MAX     <- 0.05     # fraction

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

input_rds  <- get_arg("input_rds",  required = TRUE)
output_pdf <- get_arg("output_pdf", required = TRUE)
output_rds <- get_arg("output_rds", required = TRUE)

message("input_rds: ", input_rds)
message("output_pdf: ", output_pdf)
message("output_rds: ", output_rds)

# -----------------------------------
# Load object & basic sanity checks
# -----------------------------------
seurat_obj <- readRDS(input_rds)
if (!"peaks" %in% Assays(seurat_obj)) {
  stop("Expected a 'peaks' ChromatinAssay in the Seurat object.", call. = FALSE)
}
DefaultAssay(seurat_obj) <- "peaks"

frags <- Fragments(seurat_obj, assay = "peaks")
if (length(frags) == 0) {
  stop("No fragment objects attached to assay 'peaks'.", call. = FALSE)
}
frag_path <- frags[[1]]@path  # Path() accessor not exported in Signac 1.15.0
if (is.null(frag_path) || !nzchar(frag_path) || !file.exists(frag_path)) {
  stop(sprintf("Fragments file missing or not found at: %s", as.character(frag_path)), call. = FALSE)
}

# -----------------------------
# Open PDF (safe close on exit)
# -----------------------------
dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)
pdf(file = output_pdf, width = 8.5, height = 11, pointsize = 12)
on.exit({
  try(dev.off(), silent = TRUE)
}, add = TRUE)

# ---------------------------------------------------
# Helper: text page with sane margins (no overlap)
# ---------------------------------------------------
print_text_page <- function(title, lines) {
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(1.5, 1.5, 3.5, 1.5))  # extra top margin
  plot.new()
  title(main = title, cex.main = 1.2, line = 1)
  # Add leading newlines so stats don’t overlap the title
  txt <- paste(c("", "", lines), collapse = "\n")
  mtext(txt, side = 3, adj = 0, line = -0.5, cex = 0.8, family = "mono")
}

# ----------------------
# Object summary page
# ----------------------
print_text_page(
  "Object summary",
  c(
    paste("Cells:", ncol(seurat_obj)),
    paste("Peaks:", nrow(seurat_obj)),
    paste("Fragments file:", basename(frag_path))
  )
)

# -----------------------
# TSS enrichment
# -----------------------
message("Computing TSS enrichment…")
seurat_obj <- TSSEnrichment(seurat_obj)

p_tss <- DensityScatter(
  object = seurat_obj,
  x = "nCount_peaks",
  y = "TSS.enrichment",
  log_x = TRUE,
  quantiles = TRUE
) +
  ggtitle("Library Size vs TSS Enrichment") +
  xlab("nCount_peaks (log10)") +
  ylab("TSS.enrichment") +
  geom_hline(yintercept = CUT_TSS_MIN, linetype = "dashed") +
  geom_vline(xintercept = c(CUT_NCOUNT_MIN, CUT_NCOUNT_MAX), linetype = "dashed") +
  coord_cartesian(
    ylim = c(0, max(CUT_TSS_MIN * 1.5, quantile(seurat_obj$TSS.enrichment, 0.99, na.rm = TRUE)))
  )

print(p_tss)

tss_vals <- seurat_obj$TSS.enrichment
ncount   <- seurat_obj$nCount_peaks
tss_stats <- c(
    sprintf("QC cutoffs: TSS.enrichment > %d", CUT_TSS_MIN),
   "",
   "",
  "TSS.enrichment (summary):",
  capture.output(summary(tss_vals)),
  "",
"",
  "TSS.enrichment (deciles):",
  capture.output(quantile(tss_vals, probs = seq(0, 1, 0.1), na.rm = TRUE))
)
ncount_stats <- c(
    sprintf("QC cutoffs: %d < nCount_peaks < %d", CUT_NCOUNT_MIN, CUT_NCOUNT_MAX),
  "",
  "",
  "nCount_peaks (summary):",
  capture.output(summary(ncount)),
  "",
  "",
  "nCount_peaks (deciles):",
  capture.output(quantile(ncount, probs = seq(0, 1, 0.1), na.rm = TRUE))
)
print_text_page("Summary: TSS & nCount_peaks", c(tss_stats, "", ncount_stats))

# -----------------------
# Nucleosome signal
# -----------------------
message("Computing nucleosome signal…")
seurat_obj <- NucleosomeSignal(seurat_obj)

seurat_obj$nucleosome_group <- ifelse(seurat_obj$nucleosome_signal > CUT_NS_MAX, "NS > cutoff", "NS ≤ cutoff")
p_frag_hist <- FragmentHistogram(object = seurat_obj, group.by = "nucleosome_group") +
  ggtitle(sprintf("Fragment Size Distribution by Nucleosome Signal (cutoff = %g)", CUT_NS_MAX)) +
  theme(legend.position = "top")

print(p_frag_hist)

ns_vals <- seurat_obj$nucleosome_signal
ns_stats <- c(
    sprintf("QC cutoff: nucleosome_signal < %g", CUT_NS_MAX),
    "",
    "",
  "nucleosome_signal (summary):",
  capture.output(summary(ns_vals)),
  "",
  "",
  "nucleosome_signal (deciles):",
  capture.output(quantile(ns_vals, probs = seq(0, 1, 0.1), na.rm = TRUE))
)
print_text_page("Summary: Nucleosome signal", ns_stats)

# -------------------------------
# Percent reads in peaks (FRIP)
# -------------------------------
message("Computing pct_reads_in_peaks…")
if (!all(c("peak_region_fragments", "passed_filters") %in% colnames(seurat_obj[[]]))) {
  warning("Metadata missing 'peak_region_fragments' or 'passed_filters'. Percent reads in peaks will be NA where missing.")
}
seurat_obj$pct_reads_in_peaks <- with(
  seurat_obj@meta.data,
  100 * (peak_region_fragments / pmax(passed_filters, 1))
)

hist(
  seurat_obj$pct_reads_in_peaks,
  breaks = 60,
  main = sprintf("Percent Reads in Peaks per Cell (cutoff = %d%%)", CUT_FRIP_MIN),
  xlab  = "Percent"
)
abline(v = CUT_FRIP_MIN, lty = 2)

frip_vals <- seurat_obj$pct_reads_in_peaks
frip_stats <- c(
   sprintf("QC cutoff: pct_reads_in_peaks > %d%%", CUT_FRIP_MIN),
  "",
  "",
  "pct_reads_in_peaks (summary):",
  capture.output(summary(frip_vals)),
  "",
  "",
  "pct_reads_in_peaks (deciles):",
  capture.output(quantile(frip_vals, probs = seq(0, 1, 0.1), na.rm = TRUE))
)
print_text_page("Summary: FRIP (pct_reads_in_peaks)", frip_stats)

# -----------------------
# Blacklist fraction
# -----------------------
message("Computing blacklist fraction…")
data("blacklist_hg38_unified", package = "Signac", envir = environment())
seurat_obj$blacklist_ratio <- FractionCountsInRegion(
  object  = seurat_obj,
  assay   = "peaks",
  regions = blacklist_hg38_unified
)

hist(
  seurat_obj$blacklist_ratio,
  breaks = 60,
  main = sprintf("Blacklist Ratio per Cell (cutoff = %.2f)", CUT_BL_MAX),
  xlab  = "Fraction of fragments in blacklist regions"
)
abline(v = CUT_BL_MAX, lty = 2)

bl_vals <- seurat_obj$blacklist_ratio
bl_stats <- c(
    sprintf("QC cutoff: blacklist_ratio < %.2f", CUT_BL_MAX),
  "",
  "",
  "blacklist_ratio (summary):",
  capture.output(summary(bl_vals)),
  "",
"",
  "blacklist_ratio (deciles):",
  capture.output(quantile(bl_vals, probs = seq(0, 1, 0.1), na.rm = TRUE))
)
print_text_page("Summary: Blacklist fraction", bl_stats)

message("QC report written: ", output_pdf)

# -----------------------------------------
# Filtering: apply thresholds, save object
# -----------------------------------------
filtered_cells <- subset(
  x = seurat_obj,
  subset =
    nCount_peaks > CUT_NCOUNT_MIN &
    nCount_peaks < CUT_NCOUNT_MAX &
    pct_reads_in_peaks > CUT_FRIP_MIN &
    blacklist_ratio < CUT_BL_MAX &
    nucleosome_signal < CUT_NS_MAX &
    TSS.enrichment > CUT_TSS_MIN
)

# Record the sample name inside the object (from the filename)
sample_id <- sub("_filtered_cells\\.rds$", "", basename(output_rds))
filtered_cells@misc$sample_id <- sample_id

# Save cleaned object
dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(filtered_cells, file = output_rds)
message("Cleaned Seurat saved: ", output_rds)
