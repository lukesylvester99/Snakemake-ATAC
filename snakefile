configfile: "config.yaml"
workdir: "../.."

#-------- Load config variables --------#

SAMPLES     = config["samples"]            
FASTQ_ROOT = config["path_to_fastqs"]         
OUT_ROOT    = config["out_root"]           
REFERENCE   = config["reference"]     
SEURAT_INPUTS = config["seurat_inputs"]
LOG_ROOT = config["log_dir"]

def seurat_input_paths(sample): #helper function to get full paths to seurat input files. Called in cellranger_count and create_seurat_object rules
    return [f"{OUT_ROOT}/{sample}/outs/{fname}" for fname in SEURAT_INPUTS]


#-------- Rules --------#
rule all:
    input:
        expand(f"{OUT_ROOT}/seurat_objects/{{sample}}.rds", sample=SAMPLES)


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
    conda: "workflows/envs/seurat.yaml"
    input:
       expand("{out_root}/{sample}/.cellranger_done",
           out_root=OUT_ROOT, sample=SAMPLES)
           

    output:
        rds=f"{OUT_ROOT}" + "/seurat_objects/{sample}.rds"
    params:
        snake_outs=lambda wc: f"{OUT_ROOT}/{wc.sample}/outs"
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

          Rscript workflows/scripts/create_seurat_object.R \
            --snake_outs="{params.snake_outs}" \
            --output_rds="{output.rds}"

          echo "==== create_seurat_object END $(date) ===="
        ) &> "{log.run}"
        """
