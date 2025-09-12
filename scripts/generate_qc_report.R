#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  library(grid)
})

# -----------------------------
# Global plot theme & settings
# -----------------------------
theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
      axis.title = element_text(size = 12),
      plot.margin = margin(30, 12, 10, 12),
      panel.grid.minor = element_blank()
    )
)

# QC thresholds 
CUT_NCOUNT_MIN <- 1000
CUT_NCOUNT_MAX <- 100000
CUT_TSS_MIN    <- 4
CUT_NS_MAX     <- 4
CUT_PTC_MIN   <- 20       # percent
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
frag_path <- frags[[1]]@path  # Path() not exported in Signac 1.15.0
if (is.null(frag_path) || !nzchar(frag_path) || !file.exists(frag_path)) {
  stop(sprintf("Fragments file missing or not found at: %s", as.character(frag_path)), call. = FALSE)
}

# ----------
# Open PDF 
# ----------
dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)
pdf(file = output_pdf, width = 8.5, height = 11, pointsize = 12)
on.exit({ try(dev.off(), silent = TRUE) }, add = TRUE)

# -------------------------------------------------------------
# Helper: settings for text page with summary statistics of QC
# -------------------------------------------------------------
print_text_page <- function(title, lines, lines_per_page = 40, body_top = 0.82, left = 0.04) {
  n <- length(lines)
  if (n == 0) lines <- ""
  pages <- ceiling(n / lines_per_page)
  for (i in seq_len(pages)) {
    i_start <- (i - 1) * lines_per_page + 1
    i_end   <- min(i * lines_per_page, n)
    chunk   <- lines[i_start:i_end]

    grid::grid.newpage()
    # Title (append page x/y if multiple pages)
    page_title <- if (pages > 1) sprintf("%s  —  page %d/%d", title, i, pages) else title
    grid::grid.text(page_title, x = 0.5, y = 0.95,
                    gp = grid::gpar(fontface = "bold", cex = 1.2))
    # Body
    grid::grid.text(paste(chunk, collapse = "\n"),
                    x = left, y = body_top, just = c("left","top"),
                    gp = grid::gpar(cex = 0.8, family = "mono"))
  }
}

# -------------------------------------------------------------
# Helper: Percent pass/fail stats
# -------------------------------------------------------------

pass_stats <- function(values, pass_logical) {
  n_total <- sum(!is.na(values))
  n_pass  <- sum(pass_logical, na.rm = TRUE)
  n_fail  <- n_total - n_pass
  n_na    <- sum(is.na(values))
  pct     <- if (n_total > 0) 100 * n_pass / n_total else NA_real_
  list(n_total = n_total, n_pass = n_pass, n_fail = n_fail, n_na = n_na, pct = pct)
}

fmt_pass_line <- function(label, s, detail = TRUE) {
  base <- sprintf("%-22s: %5.1f%% pass", label, s$pct)
  if (detail) {
    paste0(base, sprintf("  (n=%d / N=%d; NA=%d)", s$n_pass, s$n_total, s$n_na))
  } else {
    base
  }
}

# ----------------------
# Object summary page
# ----------------------
print_text_page(
  "Object summary",
  c(
    paste("Cells:", ncol(seurat_obj)),
    paste("Peaks:", nrow(seurat_obj)),
    paste("Fragments file:", basename(frag_path)),
    "",
    sprintf("QC thresholds: nCount_peaks in (%d, %d), TSS > %g, NS < %g, PTC > %d%%, Blacklist < %.2f",
            CUT_NCOUNT_MIN, CUT_NCOUNT_MAX, CUT_TSS_MIN, CUT_NS_MAX, CUT_PTC_MIN, CUT_BL_MAX)
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

# capture vectors for later summary
tss_vals <- seurat_obj$TSS.enrichment
ncount   <- seurat_obj$nCount_peaks

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

# capture vectors for later summary
ns_vals <- seurat_obj$nucleosome_signal

# -----------------------
# Percent reads in peaks 
# ------------------------
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
  main = sprintf("Percent Reads in Peaks per Cell (cutoff = %d%%)", CUT_PTC_MIN),
  xlab  = "Percent"
)
abline(v = CUT_PTC_MIN, lty = 2)

ptc_vals <- seurat_obj$pct_reads_in_peaks

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


# ----- Pass/fail per metric -----
pass_ncount <- (ncount > CUT_NCOUNT_MIN) & (ncount < CUT_NCOUNT_MAX)
pass_tss    <- (tss_vals > CUT_TSS_MIN)
pass_ns     <- (ns_vals  < CUT_NS_MAX)
pass_ptc   <- (ptc_vals > CUT_PTC_MIN)
pass_bl     <- (bl_vals  < CUT_BL_MAX)

s_ncount <- pass_stats(ncount,   pass_ncount)
s_tss    <- pass_stats(tss_vals, pass_tss)
s_ns     <- pass_stats(ns_vals,  pass_ns)
s_ptc   <- pass_stats(ptc_vals,pass_ptc)
s_bl     <- pass_stats(bl_vals,  pass_bl)

# Overall: all criteria simultaneously (NA-safe “and”)
all_mat <- cbind(pass_ncount, pass_tss, pass_ns, pass_ptc, pass_bl)
row_ok  <- apply(all_mat, 1, function(x) all(x %in% TRUE))  # FALSE if any criterion FALSE; NA if any NA
overall_values <- ifelse(rowSums(!is.na(all_mat)) == ncol(all_mat), TRUE, NA)  # marks rows with complete data
s_overall <- pass_stats(overall_values, row_ok)


# -----------------------------------------
# ONE summary page with all percentiles
# -----------------------------------------
fmt_quantiles <- function(x, label, probs = seq(0, 1, 0.1)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  c(
    label,
    paste(names(q), sprintf("= %.3f", as.numeric(q)))
  )
}

summary_lines <- c(
  sprintf("QC thresholds:"),
  sprintf("  - nCount_peaks: %d < x < %d", CUT_NCOUNT_MIN, CUT_NCOUNT_MAX),
  sprintf("  - TSS.enrichment: x > %g", CUT_TSS_MIN),
  sprintf("  - nucleosome_signal: x < %g", CUT_NS_MAX),
  sprintf("  - pct_reads_in_peaks: x > %d%%", CUT_PTC_MIN),
  sprintf("  - blacklist_ratio: x < %.2f", CUT_BL_MAX)
)

summary_lines <- c(
  summary_lines,
  fmt_quantiles(ncount,   "nCount_peaks:"),
  "",
  fmt_quantiles(tss_vals, "TSS.enrichment:"),
  "",
  fmt_quantiles(ns_vals,  "nucleosome_signal:"),
  "",
  fmt_quantiles(ptc_vals,"pct_reads_in_peaks:"),
  "",
  fmt_quantiles(bl_vals,  "blacklist_ratio:")
)

print_text_page("QC Summary", summary_lines, lines_per_page = 38)

message("QC report written: ", output_pdf)


print_text_page(
  "QC Threshold Pass Rates",
  c(
    sprintf("Cutoffs:  nCount_peaks (%d, %d),  TSS > %g,  NS < %g,  PTC > %d%%,  Blacklist < %.2f",
            CUT_NCOUNT_MIN, CUT_NCOUNT_MAX, CUT_TSS_MIN, CUT_NS_MAX, CUT_PTC_MIN, CUT_BL_MAX),
    "",
    fmt_pass_line("nCount_peaks",       s_ncount),
    fmt_pass_line("TSS.enrichment",     s_tss),
    fmt_pass_line("nucleosome_signal",  s_ns),
    fmt_pass_line("pct_reads_in_peaks", s_ptc),
    fmt_pass_line("blacklist_ratio",    s_bl),
    "",
    fmt_pass_line("ALL criteria",       s_overall)
  )
)

# -----------------------------------------
# Filtering: apply thresholds, save object
# -----------------------------------------
filtered_cells <- subset(
  x = seurat_obj,
  subset =
    nCount_peaks > CUT_NCOUNT_MIN &
    nCount_peaks < CUT_NCOUNT_MAX &
    pct_reads_in_peaks > CUT_PTC_MIN &
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
