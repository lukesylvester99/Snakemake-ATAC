configfile: "config.yaml"
workdir: "../.."

SAMPLES     = config["samples"]            
FASTQ_ROOT = config["path_to_fastqs"]         
OUT_ROOT    = config["out_root"]           
REFERENCE   = config["reference"]     
SEURAT_INPUTS = config["seurat_inputs"]      

rule all:
    input:
        # final QC reports
        expand("{out_root}/qc_reports/{sample}_qc_report.pdf",
               out_root=OUT_ROOT, sample=SAMPLES),
        # final cleaned Seurat objects
        expand("{out_root}/seurat_objects_clean/{sample}_filtered_cells.rds",
               out_root=OUT_ROOT, sample=SAMPLES)



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
        snake_outs_dir = directory(f"{OUT_ROOT}/{{sample}}/outs"),
        done = touch(f"{OUT_ROOT}/{{sample}}/.cellranger_done")
    params:
        outdir=OUT_ROOT,
        reference=REFERENCE,
        run_id=lambda wc: wc.sample
    threads: 16
    resources:
        mem_mb=64000

    shell:
        r"""
        set -euo pipefail

        test -d "{input.fastq_dir}" || (echo "FASTQ folder not found: {input.fastq_dir}" && exit 1)

        mkdir -p "{params.outdir}"
        cd "{params.outdir}"

        cellranger-atac count \
          --id="{params.run_id}" \
          --reference="{params.reference}" \
          --fastqs="{input.fastq_dir}" \
          --localcores {threads}
         
        touch "{wildcards.sample}/.cellranger_done"
        """

rule create_seurat_object:
    """
        Rule to create a Seurat object from Cell Ranger ATAC output.
        Object will be loaded with all necessary matrices, metadata, and annotations.
        Object will be saved as an RDS file.
    """
    input:
        # ensure Cell Ranger finished
        snake_outs_dir = lambda wc: f"{OUT_ROOT}/{wc.sample}/outs/fragments.tsv.gz",
        done = lambda wc: f"{OUT_ROOT}/{wc.sample}/.cellranger_done"

    output:
        rds=f"{OUT_ROOT}" + "/seurat_objects/{sample}.rds"
    params:
        snake_outs=lambda wc: f"{OUT_ROOT}/{wc.sample}/outs"
    threads: 4  
    shell:
        r"""
        mkdir -p "{OUT_ROOT}/seurat_objects"
        Rscript workflows/scripts/create_seurat_object.R \
          --snake_outs "{params.snake_outs}" \
          --output_rds "{output.rds}"
        """

