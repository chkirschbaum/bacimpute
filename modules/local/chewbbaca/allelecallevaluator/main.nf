process CHEWBBACA_ALLELECALLEVALUATOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.5.4--pyh106432d_0':
        'quay.io/biocontainers/chewbbaca:3.5.4--pyh106432d_0' }"

    input:
    path allelecall_dir
    tuple val(meta), path(scheme)

    output:
    path "allelecall_report.html"           , emit: allelecall_report
    path "masked_profiles.tsv"              , emit: masked_profiles
    path "cgMLST_profiles.tsv"              , emit: cgMLST_profiles
    path "report_bundle.js"                 , emit: report_bundle
    path "presence_absence.tsv"             , emit: presence_absence    , optional:true
    path "distance_matrix_symmetric.tsv"    , emit: distance_matrix     , optional:true
    path "protein_msa.fasta"                , emit: protein_msa         , optional:true
    path "cgMLST.tree"                      , emit: tree                , optional:true
    tuple val("${task.process}"), val("chewbbaca"), eval("chewie --version 2>&1 | sed 's/chewBBACA version: //'"), topic: versions, emit: versions_chewbbaca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    chewBBACA.py \\
        AlleleCallEvaluator \\
        --cpu ${task.cpus} \\
        $args \\
        --input-files ${allelecall_dir} \\
        --schema-directory ${scheme} \\
        --output-directory results

    mv results/allelecall_report.html allelecall_report.html 
    mv results/masked_profiles.tsv masked_profiles.tsv
    mv results/cgMLST_profiles.tsv cgMLST_profiles.tsv
    mv results/report_bundle.js report_bundle.js

    # Handle optional output files
    [ -f results/presence_absence.tsv ] && mv results/presence_absence.tsv presence_absence.tsv
    [ -f results/distance_matrix_symmetric.tsv ] && mv results/distance_matrix_symmetric.tsv distance_matrix_symmetric.tsv
    [ -f results/protein_msa.fasta ] && mv results/protein_msa.fasta protein_msa.fasta
    [ -f results/cgMLST.tree ] && mv results/cgMLST.tree cgMLST.tree
    """

    stub:
    """
    echo $args

    touch allelecall_report.html
    touch masked_profiles.tsv
    touch cgMLST_profiles.tsv
    touch report_bundle.js

    # Optional files
    touch presence_absence.tsv
    touch distance_matrix_symmetric.tsv
    touch protein_msa.fasta
    touch cgMLST.tree
    """
}
