configfile: "config.yaml"
workdir: "../.."
from collections import defaultdict

#-------- Load config variables and helper functions --------#

SAMPLES     = config["samples"]            
FASTQ_ROOT = config["path_to_fastqs"]         
OUT_ROOT    = config["out_root"]           
REFERENCE   = config["reference"]     
SEURAT_INPUTS = config["seurat_inputs"]
LOG_ROOT = config["log_dir"]

#helper function to get full paths to seurat input files. Called in cellranger_count and create_seurat_object rules
def seurat_input_paths(sample): 
    return [f"{OUT_ROOT}/{sample}/outs/{fname}" for fname in SEURAT_INPUTS]

""" 
    Following functions extract timepoint from sample name.
    Assumes sample names are formatted as 'S1_1', 'S2_2', etc.
    Returns the timepoint part, e.g., 'S1', 'S2'.
    Called in the merge_timepoint rule.
"""

def timepoint_of(sample: str) -> str:
    # Expects names like S1_1, S2_2, etc. Returns 'S1', 'S2', ...
    return sample.split("_")[0]

# Group replicates by timepoint
TP_REPS = defaultdict(list)
for s in SAMPLES:
    TP_REPS[timepoint_of(s)].append(s)

TIMEPOINTS = sorted(TP_REPS.keys())

#-------- Rules --------#
rule all:
    input:
         expand("{out_root}/seurat_objects_merged/{tp}_merged.rds",
               out_root=OUT_ROOT, tp=TIMEPOINTS)


rule cellranger_count:
    """
        Rule runs all fastq files through cellranger-atac count. It assumes a directory structure like:

        └── parent directory (patient ID)/
            ├── Raw fastqs/
            │   ├── time point 1
            │   ├── tiepoint 2
            │   └── timepoint 3
            ├── outs/
            │   ├── time point 1/
            │   │   └── cellranger files
            │   ├── timepoint 2/
            │   │   └── cellranger files
            │   └── timepoint 3/
            │       └── cellranger files
            └── workflows/
                ├── snakefiles/
                │   ├── snakefile
                │   └── config.yaml
                └── scripts/
                    ├── Rscript1
                    └── Rscript2
        """

    input:
        fastq_dir=lambda wc: f"{FASTQ_ROOT}/{wc.sample}"
    output:
        done = touch(f"{OUT_ROOT}/{{sample}}/.cellranger_done")
    params:
        outdir=OUT_ROOT,
        reference=REFERENCE,
        run_id=lambda wc: wc.sample
    threads: 16
    resources:
        mem_mb=64000
    log:
        run = f"{LOG_ROOT}" + "/cellranger_count/{sample}.log"
    shell:
          r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname "{log.run}")"

        (
          echo "==== cellranger_count START $(date) ===="
          echo "HOST: $(hostname)"
          echo "SAMPLE: {wildcards.sample}"
          echo "FASTQ_DIR: {input.fastq_dir}"
          echo "OUT_ROOT: {params.outdir}"
          echo "REFERENCE: {params.reference}"
          echo "THREADS: {threads}"
          echo

          test -d "{input.fastq_dir}" || (echo "FASTQ folder not found: {input.fastq_dir}" >&2; exit 1)

          cd "{params.outdir}"

          cellranger-atac count \
            --id="{params.run_id}" \
            --reference="{params.reference}" \
            --fastqs="{input.fastq_dir}" \
            --localcores {threads}
         
          echo
          echo "Marking completion: {wildcards.sample}/.cellranger_done"
          touch "{wildcards.sample}/.cellranger_done"

          echo "==== cellranger_count END $(date) ===="
        ) &> "{log.run}"
        """


rule create_seurat_object:
    """
        Rule to create a Seurat object from Cell Ranger ATAC output.
        Object will be loaded with all necessary matrices, metadata, and annotations.
        Object will be saved as an RDS file.
    """
    conda: "../envs/seurat.yaml"
    input:
        done = f"{OUT_ROOT}/{{sample}}/.cellranger_done"
    output:
        rds=f"{OUT_ROOT}" + "/seurat_objects/{sample}.rds"
    params:
        snake_outs=lambda wc: f"{OUT_ROOT}/{wc.sample}/outs",
        anno_cache = "workflows/cache/hg38_annotations_ucsc.rds"
    threads: 4  
    log:
        run = f"{LOG_ROOT}" + "/create_seurat_obj/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUT_ROOT}/seurat_objects" "$(dirname "{log.run}")"

        (
        echo "==== create_seurat_object START $(date) ===="
        echo "SAMPLE: {wildcards.sample}"
        echo "snake_outs: {params.snake_outs}"
        echo "OUTPUT_RDS: {output.rds}"
        echo

        # Show the conda rule's R
        which Rscript
        Rscript -e 'cat("R.home(): ", R.home(), "\n"); cat(".libPaths():\n"); print(.libPaths())'
        Rscript -e 'library(Signac); cat("Signac version: ", as.character(packageVersion("Signac")), "\n")'

        Rscript workflows/scripts/create_seurat_object.R \
            --snake_outs="{params.snake_outs}" \
            --output_rds="{output.rds}"\
            --anno_cache="{params.anno_cache}"

        echo "==== create_seurat_object END $(date) ===="
        ) &> "{log.run}"
        """


rule qc_metrics:
    """
        Rule to generate QC metrics and plots from Seurat object.
        Generates a PDF report with various QC plots and statistics, as well
        as a cleaned seurat object with QC metrics added.  
    """
    conda: "../envs/seurat.yaml"
    input:
        rds = f"{OUT_ROOT}" + "/seurat_objects/{sample}.rds"
    output:
        pdf     = f"{OUT_ROOT}" + "/qc_reports/{sample}_qc_report.pdf",
        cleaned = f"{OUT_ROOT}" + "/seurat_objects_clean/{sample}_filtered_cells.rds"
    threads: 6
    log:
        run = f"{LOG_ROOT}" + "/qc_metrics/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUT_ROOT}/qc_reports" "{OUT_ROOT}/seurat_objects_clean" "$(dirname "{log.run}")"

        (
          echo "==== qc_metrics START $(date) ===="
          echo "SAMPLE: {wildcards.sample}"
          echo "INPUT_RDS: {input.rds}"
          echo "OUTPUT_PDF: {output.pdf}"
          echo "OUTPUT_RDS (clean): {output.cleaned}"
          echo

          Rscript workflows/scripts/generate_qc_report.R \
            --input_rds="{input.rds}" \
            --output_pdf="{output.pdf}" \
            --output_rds="{output.cleaned}"

          echo "==== qc_metrics END $(date) ===="
        ) &> "{log.run}"
        """


rule merge_timepoint:
    """
    Merge cleaned Seurat objects for replicates within each timepoint.
    Produces one merged object per timepoint (e.g., S1_merged.rds).
    """
    conda: "../envs/merge.yaml"
    input:
        lambda wc: expand(
            f"{OUT_ROOT}/seurat_objects_clean/{{s}}_filtered_cells.rds",
            s=TP_REPS[wc.timepoint]
        )
    params:
        # Build one --input_rds="..." flag per file (robust, no giant argv token)
        input_flags=lambda wc: ' '.join(
            f'--input_rds="{p}"'
            for p in expand(
                f"{OUT_ROOT}/seurat_objects_clean/{{s}}_filtered_cells.rds",
                s=TP_REPS[wc.timepoint]
            )
        )
    output:
        merged = f"{OUT_ROOT}/seurat_objects_merged/{{timepoint}}_merged.rds"
    threads: 6
    log:
        run = f"{LOG_ROOT}/merge/{{timepoint}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUT_ROOT}/seurat_objects_merged" "$(dirname "{log.run}")"

        (
          echo "==== merge_timepoint START $(date) ===="
          echo "TIMEPOINT: {wildcards.timepoint}"
          echo "N_INPUTS: $(echo {input} | wc -w)"
          for f in {input}; do echo "  - $f"; done
          echo "OUTPUT_RDS (merged): {output.merged}"
          echo

          Rscript workflows/scripts/merge.R \
            {params.input_flags} \
            --output_rds="{output.merged}"

          echo "==== merge_timepoint END $(date) ===="
        ) &> "{log.run}"
        """
