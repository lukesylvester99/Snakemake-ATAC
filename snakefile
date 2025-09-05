configfile: "config.yaml"

SAMPLES     = config["samples"]            
FASTQ_ROOT  = config["fastq_root"]          
OUT_ROOT    = config["out_root"]           
REFERENCE   = config["reference"]           

rule all:
    input:
        expand("results/{samples}/qc/filtered.{samples}.seurat.rds", samples=SAMPLES)

rule cellranger_count:
    """
        Rule runs all fastq files through cellranger-atac count. It assumes a directory structure like:

        └── parent directory (patient ID)/
            ├── Raw fastqs/
            │   ├── time point 1
            │   ├── tiepoint 2
            │   └── timepoint 3
            └── outs/
                ├── time point 1/
                │   └── cellranger files
                ├── timepoint 2/
                │   └── cellranger files
                └── timepoint 3/
                    └── cellranger files
        """

    input:
        fastq_dir=lambda wc: f"{FASTQ_ROOT}/{wc.sample}"
    output:
        html=f"{OUT_ROOT}" + "/{sample}/outs/web_summary.html"
    params:
        outdir=OUT_ROOT,
        reference=REFERENCE,
        run_id=lambda wc: wc.sample
    threads: 8
    resources:
        mem_mb=64000
    conda:
        "envs/cellranger-atac.yaml"   # optional; ensure cellranger-atac is on PATH
    shell:
        r"""
        test -d "{input.fastq_dir}" || (echo "FASTQ folder not found: {input.fastq_dir}" && exit 1)

        mkdir -p "{params.outdir}"
        cd "{params.outdir}"

        cellranger-atac count \
          --id="{params.run_id}" \
          --reference="{params.reference}" \
          --fastqs="{input.fastq_dir}" \
          --localcores {threads} \
          --localmem {int(resources.mem_mb/1024)}
        """