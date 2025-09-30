#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np

# -----------------
# Arg parsing
# -----------------
parser = argparse.ArgumentParser(description="Process mgatk outputs")
parser.add_argument(
    "--mgatk_dir",
    required=True,
    help="Path to mgatk output directory containing the final files "
         "(e.g., .../mgatk/<sample>/final)"
)
parser.add_argument(
    "--out_dir",
    required=True,
    help="Directory to write processed mgatk summaries (four TSVs)"
)
args = parser.parse_args()

# Normalize paths (useful for logs)
args.mgatk_dir = os.path.abspath(args.mgatk_dir)
args.out_dir   = os.path.abspath(args.out_dir)

print(f"mgatk_dir: {args.mgatk_dir}")
print(f"out_dir: {args.out_dir}")

# Derive sample_id robustly: .../mgatk/<sample>/final -> <sample>
sample_id = os.path.basename(os.path.dirname(args.mgatk_dir))

# Ensure base output directory exists
out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

#-----------------
# Process mgatk output
#-----------------

############################
# Filter cells by coverage #
############################

#thresholds
MIN_CELL_MEAN_DEPTH = 10 #caleb used 10 in paper

# Paths from args
path_depth_table = os.path.join(args.mgatk_dir, f"{sample_id}.depthTable.txt")
out_dir = args.out_dir


# Load depth table (two columns: barcode and mean depth)
depth_df = pd.read_csv(
    path_depth_table,
    sep="\t",
    header=None,
    names=["barcode", "mean_depth"] #explicitly name columns
)

# Drop incomplete rows and force mean_depth to numeric
depth_df = depth_df.dropna()
depth_df["mean_depth"] = pd.to_numeric(depth_df["mean_depth"], errors="coerce") #coerce -> errors to NaN
depth_df = depth_df.dropna(subset=["mean_depth"]) #drop rows where mean_depth is NaN

# Filter cells by coverage threshold (good_cells is a one-column Series of barcodes from cells with enough coverage)
good_cells = depth_df.loc[depth_df["mean_depth"] >= MIN_CELL_MEAN_DEPTH, "barcode"].astype(str)
print(f"Cells passing depth ≥{MIN_CELL_MEAN_DEPTH}×: {len(good_cells)} / {len(depth_df)}")

# Converts the Series into a one-column DataFrame with column name "barcode" and saves as TSV
os.makedirs(out_dir, exist_ok=True)
good_cells.to_frame(name="barcode").to_csv(
    os.path.join(out_dir, "cells_passing_depth.tsv"),
    sep="\t",
    index=False
)

############################
# Inspect & filter variants #
############################

path_var_stats = os.path.join(args.mgatk_dir, f"{sample_id}.variant_stats.tsv.gz")

# Load variant stats
var_stats = pd.read_csv(path_var_stats, sep="\t", compression="infer")
print("Variant stats columns:", list(var_stats.columns))

# thresholds 
MIN_CELLS_WITH_VARIANT = 5
MIN_STRAND_CORR = 0.5 #caleb used 0.65 in paper
MIN_LOG10_VMR  = -2.0

# Expected columns in variant_stats
expected_cols = [
    "position", "nucleotide", "variant", "vmr", "mean", "variance",
    "n_cells_conf_detected", "n_cells_over_5", "n_cells_over_10",
    "n_cells_over_20", "n_cells_over_95", "max_heteroplasmy",
    "strand_correlation", "mean_coverage"
]

# Fail fast if schema is not as expected
missing = [c for c in expected_cols if c not in var_stats.columns]
if missing:
    raise ValueError(f"variant_stats missing expected columns: {missing}")

# Ensure numeric dtype where appropriate
numeric_cols = [
    "position", "vmr", "mean", "variance", "n_cells_conf_detected",
    "n_cells_over_5", "n_cells_over_10", "n_cells_over_20", "n_cells_over_95",
    "max_heteroplasmy", "strand_correlation", "mean_coverage"
]
for c in numeric_cols:
    var_stats[c] = pd.to_numeric(var_stats[c], errors="coerce")

# Compute log10(vmr) safely (vmr <= 0 -> NaN -> will fail the > -2 filter)
log10_vmr = np.log10(var_stats["vmr"].where(var_stats["vmr"] > 0))

#build boolean mask of variants to keep
keep = (
    (var_stats["n_cells_conf_detected"] >= MIN_CELLS_WITH_VARIANT) &
    (var_stats["strand_correlation"] >= MIN_STRAND_CORR) &
    (log10_vmr > MIN_LOG10_VMR)
)

filtered_variants = var_stats.loc[keep].copy()
print(f"Variants passing filters: {len(filtered_variants)} / {len(var_stats)}")

# Save filtered table
os.makedirs(out_dir, exist_ok=True)
filtered_variants.to_csv(os.path.join(out_dir, "variant_stats_filtered.tsv"), sep="\t", index=False)


############################
# per-cell heteroplasmy table #
############################

path_cell_heteroplasmic_df = os.path.join(args.mgatk_dir, f"{sample_id}.cell_heteroplasmic_df.tsv.gz")

# Load per-cell heteroplasmy calls
cell_het = pd.read_csv(
    path_cell_heteroplasmic_df,
    sep="\t",
    compression="infer",
    header=0,
    index_col=0
)

# force rows (barcodes) and columns (variants) to str
cell_het.index = cell_het.index.astype(str)
cell_het.columns = cell_het.columns.astype(str)

# Keep only depth-QC’d cells (rows)
cell_het = cell_het.loc[cell_het.index.isin(set(good_cells))]

# Keep only passing variants (columns)
pass_variants = set(filtered_variants["variant"].astype(str).unique())
cols_to_keep = [c for c in cell_het.columns if c in pass_variants]
cell_het_filt = cell_het.loc[:, cols_to_keep].copy()

# Coerce all entries to numeric (in case of any stray non-numeric values)
cell_het_filt = cell_het_filt.apply(pd.to_numeric, errors="coerce").fillna(0.0)

print(f"Per-cell heteroplasmy matrix after filtering: {cell_het_filt.shape[0]} cells × {cell_het_filt.shape[1]} variants")

# Save filtered VAF matrix: rows=cells, columns=variants, values=VAF
cell_het_filt.to_csv(os.path.join(out_dir, "cell_heteroplasmy_filtered.tsv"), sep="\t")

############################
# Variant-level summaries  #
############################

#for each variant, how many unique cells is it found in,
# and then order the variants from most common to least common
var_counts = (
    cell_het_filt[cell_het_filt > 0]
    .count(axis=0)
    .sort_values(ascending=False)
    .rename("n_cells_with_variant")
)

#allele fraction distribution for each variant across all cells
#tells you # cells, avg and median heteroplasmy, and maximum heteroplasmy
var_vaf_stats = (
    cell_het_filt
    .agg(["count","mean","median","max"])
    .T
    .rename(columns={"count":"n_obs"})
)

#combine the two summaries into one table
#Each row = one variant, with both prevalence (# cells) and VAF stats
summary = var_counts.to_frame().join(var_vaf_stats, how="outer").fillna(0)
summary = summary.sort_values("n_cells_with_variant", ascending=False)
summary.to_csv(os.path.join(out_dir, "variant_summary.tsv"), sep="\t")

############
# READ ME #
############

#outputs a set of files to the output directory:
# 1. cells_passing_depth.tsv: list of cell barcodes passing depth filter

# 2. variant_stats_filtered.tsv: filtered variant stats after applying strand cor & VMR thresholds

# 3. cell_heteroplasmy_filtered.tsv: A wide VAF matrix after two filters:
            # Rows (cells): only barcodes in cells_passing_depth.tsv
            # Columns (variants): only variants in variant_stats_filtered.tsv

# 4. variant_summary.tsv: A per-variant summary joining:
            # n_cells_with_variant (how many unique cells carry the variant)
            # VAF stats across cells: n_obs, mean, median, max
