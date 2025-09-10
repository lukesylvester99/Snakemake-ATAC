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
#  arg parsing #
# ##################

#Need to disect the argumments passed from Snakemake: 
#        --input_rds snake_outs/seurat_objects/sample.rds 
#       --output_pdf snake_outs/qc_reports/sample_qc_report.pdf


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

input_rds    <- get_arg("input_rds",    required = TRUE) #call function 3x
output_pdf <- get_arg("output_pdf", required = TRUE)
output_rds  <- get_arg("output_rds", required = TRUE)

message("input_rds: ", input_rds)
message("output_pdf: ", output_pdf)
message("output_rds: ", output_rds)

##############
# QC Metrics #
##############

#Load the Seurat object and generate QC metrics and plots.


# ---------- Load object ----------
seurat_obj <- readRDS(input_rds) #load seurat object
if (!"peaks" %in% Assays(seurat_obj)) {
  stop("Expected a 'peaks' ChromatinAssay in the Seurat object.", call. = FALSE)
}
DefaultAssay(seurat_obj) <- "peaks"

frag_path <- Fragments(seurat_obj, assay = "peaks") #check that fragments file is attached to object
if (length(frag_path) == 0 || !file.exists(Signac::Path(frags[[1]]))) {
  stop("Fragments file not attached or missing. Ensure create_seurat_object.R set 'fragments='.", call. = FALSE)
}

# ---------- Open a PDF device ----------
dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)
pdf(file = output_pdf, width = 8.5, height = 11)
on.exit(dev.off(), add = TRUE)

# A helper to print text blocks as a page
print_text_page <- function(title, lines) {
  plot.new(); title(main = title, cex.main = 1.2)
  txt <- paste(lines, collapse = "\n")
  mtext(txt, side = 3, adj = 0, line = -2, cex = 0.75, family = "mono")
}

# ---------- TSS enrichment ----------
message("Computing TSS enrichment…")
seurat_obj <- TSSEnrichment(seurat_obj)

# Scatter: counts vs TSS
p_tss <- DensityScatter(
  object = seurat_obj,
  x = "nCount_peaks",
  y = "TSS.enrichment",
  log_x = TRUE,
  quantiles = TRUE
) + ggtitle("Library Size vs TSS Enrichment") +
    xlab("nCount_peaks (log)") + ylab("TSS.enrichment")

print(p_tss)

tss_vals <- seurat_obj$TSS.enrichment #extract TSS and nCount values from metadata
ncount   <- seurat_obj$nCount_peaks
tss_stats <- c(capture.output(summary(tss_vals)),
               "",
               paste("Quantiles (TSS.enrichment):"),
               capture.output(quantile(tss_vals, probs = seq(0, 1, 0.1), na.rm = TRUE)))
ncount_stats <- c(capture.output(summary(ncount)),
                  "",
                  paste("Quantiles (nCount_peaks):"),
                  capture.output(quantile(ncount, probs = seq(0, 1, 0.1), na.rm = TRUE)))
print_text_page("Summary: TSS & nCount_peaks", c(tss_stats, "", ncount_stats))

# ---------- Nucleosome signal ----------
message("Computing nucleosome signal…")
seurat_obj <- NucleosomeSignal(seurat_obj)

seurat_obj$nucleosome_group <- ifelse(seurat_obj$nucleosome_signal > 4, "NS > 4", "NS <= 4")
p_frag_hist <- FragmentHistogram(object = seurat_obj, group.by = "nucleosome_group") +
  ggtitle("Fragment size distribution by nucleosome signal")
print(p_frag_hist)

ns_vals <- seurat_obj$nucleosome_signal
ns_stats <- c(capture.output(summary(ns_vals)),
              "",
              "Quantiles (nucleosome_signal):",
              capture.output(quantile(ns_vals, probs = seq(0, 1, 0.1), na.rm = TRUE)))
print_text_page("Summary: Nucleosome signal", ns_stats)

# ---------- percent reads in peaks  ----------
message("Computing pct_reads_in_peaks…")
if (!all(c("peak_region_fragments","passed_filters") %in% colnames(seurat_obj[[]]))) {
  warning("Metadata missing 'peak_region_fragments' or 'passed_filters'. Percent reads in peaks will be NA where missing.")
}
seurat_obj$pct_reads_in_peaks <- with(seurat_obj@meta.data,
  100 * (peak_region_fragments / pmax(passed_filters, 1))
)

hist(seurat_obj$pct_reads_in_peaks,
     main = "Percent Reads in Peaks per Cell", xlab = "Percent",
     breaks = 50, col = "steelblue")
frip_vals <- seurat_obj$pct_reads_in_peaks
frip_stats <- c(capture.output(summary(frip_vals)),
                "",
                "Quantiles (pct_reads_in_peaks):",
                capture.output(quantile(frip_vals, probs = seq(0, 1, 0.1), na.rm = TRUE)))
print_text_page("Summary: FRIP (pct_reads_in_peaks)", frip_stats)

# ---------- Blacklist fraction ----------
message("Computing blacklist fraction…")
data("blacklist_hg38_unified", package = "Signac", envir = environment())
seurat_obj$blacklist_ratio <- FractionCountsInRegion(
  object  = seurat_obj,
  assay   = "peaks",
  regions = blacklist_hg38_unified
)

hist(seurat_obj$blacklist_ratio,
     breaks = 50,
     main = "Blacklist Ratio per Cell",
     xlab = "Fraction of fragments in blacklist regions",
     col = "skyblue", border = "white")
bl_vals <- seurat_obj$blacklist_ratio
bl_stats <- c(capture.output(summary(bl_vals)),
              "",
              "Quantiles (blacklist_ratio):",
              capture.output(quantile(bl_vals, probs = seq(0, 1, 0.1), na.rm = TRUE)))
print_text_page("Summary: Blacklist fraction", bl_stats)

# ---------- Close PDF ----------
dev.off()

message("QC report written: ", output_pdf)


################################
# Filter the cells based on QC #
################################

#filter the Seurat object to remove low-quality cells based on QC metrics.
#Will return a new Seurat object with only high-quality cells, 'filtered_cells'

#QC metrics used:
#    - nCount_peaks > 1000 & < 100000
#    - pct_reads_in_peaks > 20
#   - blacklist_ratio < 0.05
#    - nucleosome_signal < 4
#    - TSS.enrichment > 4


filtered_cells <- subset(
  x = seurat_obj,
  subset = nCount_peaks > 1000 & 
    nCount_peaks < 100000 & 
    pct_reads_in_peaks > 20 & 
    blacklist_ratio < 0.05 & 
    nucleosome_signal < 4 & 
    TSS.enrichment > 4 
  ) 

# record the sample name inside the object (from the filename)
sample_id <- sub("_filtered_cells\\.rds$", "", basename(output_rds))
filtered_cells@misc$sample_id <- sample_id

# ensure the output dir exists and save the cleaned object
dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(filtered_cells, file = output_rds)
message("Cleaned Seurat saved: ", output_rds)
