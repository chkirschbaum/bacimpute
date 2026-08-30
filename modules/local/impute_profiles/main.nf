process IMPUTE_PROFILES {
    label 'process_single'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.output}/${params.impute_outdir}", mode: params.publish_dir_mode

    input:
        path profiles
        val model
        path model_file

    output:
        path "imputed.tsv",   emit: imputed

    script:
    """
    echo "IMPUTING MISSING ALLELES"
    
    impute.py \
    -i ${profiles} \
    -m ${model} \
    -f ${model_file} \
    --cpu ${task.cpus} \
    -o "imputed.tsv"
    """

    stub:
    """
    touch imputed.tsv
    """
}