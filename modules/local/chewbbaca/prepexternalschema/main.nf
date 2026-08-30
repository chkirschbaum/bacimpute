process CHEWBBACA_PREPEXTERNALSCHEMA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.5.4--pyh106432d_0':
        'quay.io/biocontainers/chewbbaca:3.5.4--pyh106432d_0' }"

    input:
    tuple val(meta), path(schema_dir)
    val adapted_schema_dir
    path training_file

    output:
    path "*_invalid_alleles.txt"        , emit: invalid_alleles
    path "*_invalid_genes.txt"          , emit: invalid_genes
    path "*_summary_stats.tsv"          , emit: summary_stats
    tuple val(meta), path(adapted_schema_dir), emit: schema_dir
    tuple val("${task.process}"), val("chewbbaca"), eval("chewie --version 2>&1 | sed 's/chewBBACA version: //'"), topic: versions, emit: versions_chewbbaca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    chewbbaca.py \\
        PrepExternalSchema \\
        --cpu ${task.cpus} \\
        $args \\
        --schema-directory ${schema_dir} \\
        --output-directory ${adapted_schema_dir} \\
        --ptf ${training_file}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args
    
    touch ${adapted_schema_dir}_invalid_alleles.txt
    touch ${adapted_schema_dir}_invalid_genes.txt
    touch ${adapted_schema_dir}_summary_stats.tsv
    mkdir -p ${adapted_schema_dir}
    touch ${adapted_schema_dir}/.schema_config
    """
}
