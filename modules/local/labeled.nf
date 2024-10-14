process LABELED {
    tag "$meta.id"
    label 'process_high_memory'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 conda-forge::scipy==1.14.0'
    container "raulee/sgr-python-samtools"

    input:
    tuple val(meta), path(bam), path(snp), path(matrix_dir)
    val bgfiles


    output:
    tuple val(meta), path("${meta.id}.matrix/"),               emit: matrix2
    tuple val(meta), path("${meta.id}.matrix/labeled"),        emit: labeled_matrix
    tuple val(meta), path("${meta.id}.matrix/unlabeled"),      emit: unlabeled_matrix
    tuple val(meta), path("${meta.id}.labeled_detail.txt.gz"), emit: labeled_detail
    tuple val(meta), path("*.json"),                           emit: json
    //path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extra_bg = bgfiles.join(' ')
    
    """
    labeled.py \\
        --bam ${bam} \\
        --outdir ./ \\
        --sample ${meta.id} \\
        --filtered_matrix ${matrix_dir} \\
        --bg_snp ${snp} $extra_bg
    """

}