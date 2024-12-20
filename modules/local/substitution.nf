process SUBSTITUTION {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1'
    container "raulee/sgr-python-samtools"

    input:
    tuple val(meta), path(bam), path(snp)
    val bgfiles

    output:
    tuple val(meta), path("*.TC_substitution.tsv"),  emit: substitution_stat
    tuple val(meta), path("*.json"), emit: json

    script:
    def args = task.ext.args ?: ''
    def extra_bg = bgfiles.join(' ')

    """
    substitution.py \\
        --bam ${bam} \\
        --outdir ./ \\
        --sample ${meta.id} \\
        --bg_snp ${snp} $extra_bg  $args 
    """
}