rule all:
    input:
        expand("results/{samples}/qc/filtered.{samples}.seurat.rds", samples=SAMPLES)