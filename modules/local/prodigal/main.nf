process PRODIGAL_TRAIN {
    label 'process_single'

    input:
    path genome
    val species

    output:
    path "${species}.trn", emit: training_file

    script:
    """
    prodigal \
        -i "${genome}" \
        -t "${species}.trn" \
        -p single \
        -o /dev/null
    """
}