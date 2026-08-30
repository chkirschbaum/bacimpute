process MASKING_PROFILES {
    label 'process_single'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.output}/${params.masking_outdir}", mode: params.publish_dir_mode

    input:
        path train_profiles
        path test_profiles
        val model
        path model_file
        path report_template

    output:
        path "test_report.html",   emit: test_report

    script:
    """
    echo "IMPUTING MISSING ALLELES"
    
    masking.py \
        --ri ${train_profiles} \
        --ti ${test_profiles} \
        -m ${model} \
        -f ${model_file} \
        --cpu ${task.cpus} \
        --template ${report_template} \
        -r "test_report.html"
    """

    stub:
    """
    touch test_report.html
    """
}