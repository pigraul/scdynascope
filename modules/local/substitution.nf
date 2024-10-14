process SUBSTITUTION {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1'
    container "raulee/sgr-python-samtools"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.tsv"),  emit: substitution_stat
    tuple val(meta), path("*.json"), emit: json

    script:

    """
    substitution.py \\
        --bam ${bam} \\
        --outdir ./ \\
        --sample ${meta.id}
    """
}