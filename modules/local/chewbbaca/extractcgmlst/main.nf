process CHEWBBACA_EXTRACTCGMLST {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.5.4--pyh106432d_0':
        'quay.io/biocontainers/chewbbaca:3.5.4--pyh106432d_0' }"

    input:
    path results_alleles

    output:
    path "Presence_Absence.tsv"         , emit: presence_absence
    path "mdata_stats.tsv"              , emit: mdata_stats
    path "cgMLST95.tsv"                 , emit: cgmlst95
    path "cgMLSTschema95.txt"           , emit: cgmlstschema95
    path "cgMLST99.tsv"                 , emit: cgmlst99
    path "cgMLSTschema99.txt"           , emit: cgmlstschema99
    path "cgMLST100.tsv"                , emit: cgmlst100
    path "cgMLSTschema100.txt"          , emit: cgmlstschema100
    path "cgMLST.html"                  , emit: cgmlst_report
    tuple val("${task.process}"), val("chewbbaca"), eval("chewie --version 2>&1 | sed 's/chewBBACA version: //'"), topic: versions, emit: versions_chewbbaca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    chewBBACA.py \\
        ExtractCgMLST \\
        $args \\
        --input-file ${results_alleles} \\
        --output-directory results \\
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args
    
    touch Presence_Absence.tsv
    touch mdata_stats.tsv
    touch cgMLST95.tsv
    touch cgMLSTschema95.txt
    touch cgMLST99.tsv
    touch cgMLSTschema99.txt
    touch cgMLST100.tsv
    touch cgMLSTschema100.txt
    touch cgMLST.html
    """
}
