process CONVERSION_MERGE {
    tag "$meta.id"
    label 'process_median'

    conda 'conda-forge::pandas==2.2.1 bioconda::samtools=1.20'
    container "raulee/sgr-python-samtools"

    input:
    tuple val(meta), path(input_bams, stageAs: "?/*")
    tuple val(meta), path(input_csvs, stageAs: "?/*")
    val snp_threshold
    val snp_min_depth

    output:
    tuple val(meta), path("${meta.id}.PosTag.bam") , emit: conversion_bam
    tuple val(meta), path("${meta.id}.PosTag.csv") , emit: conversion_csv
    tuple val(meta), path("${meta.id}.snp.csv")    , emit: conversion_snp
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def csvList = input_csvs.join(',')

    """

    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        ${meta.id}.PosTag.bam \\
        ${input_bams}
    

    conversion_merge.py \\
        --csvlist ${input_csvs} \\
        --sample ${meta.id} \\
        --outdir ./ \\
        --snp_threshold ${snp_threshold} \\
        --snp_min_depth ${snp_min_depth}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS


    """

}