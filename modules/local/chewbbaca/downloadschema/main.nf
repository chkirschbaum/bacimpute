process CHEWBBACA_DOWNLOADSCHEMA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.5.4--pyh106432d_0':
        'quay.io/biocontainers/chewbbaca:3.5.4--pyh106432d_0' }"

    input:
    val meta
    val species
    val schema
    val schema_dir

    output:
    tuple val(meta), path(schema_dir), emit: schema_dir
    tuple val("${task.process}"), val("chewbbaca"), eval("chewie --version 2>&1 | sed 's/chewBBACA version: //'"), topic: versions, emit: versions_chewbbaca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    chewBBACA.py \\
        DownloadSchema \\
        --cpu ${task.cpus} \\
        $args \\
        --species-id ${species} \\
        --schema-id ${schema} \\
        --download-folder ${schema_dir}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args

    mkdir -p ${schema_dir}
    touch ${schema_dir}/.schema_config
    """
}
