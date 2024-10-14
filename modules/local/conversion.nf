process CONVERSION {
    tag "$meta.id"
    label 'process_high'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1'
    container "raulee/sgr-python-samtools"

    input:
    tuple val(meta), path(bam), path(bclist)
    path gtf
    val conversion_type
    val basequalilty


    output:
    tuple val(meta), path("${meta.id}_split/*.bam"), emit: bam_chunks
    tuple val(meta), path("${meta.id}_split/*.csv"), emit: csv_chunks
    //path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bam_chunks_dir = "${meta.id}_split"

    """
    mkdir ${bam_chunks_dir}

    conversion.py \\
        --bam ${bam} \\
        --gtf ${gtf} \\
        --bclist ${bclist}  \\
        --outdir ${bam_chunks_dir} \\
        --conversion_type ${conversion_type} \\
        --basequalilty ${basequalilty}

    """

}