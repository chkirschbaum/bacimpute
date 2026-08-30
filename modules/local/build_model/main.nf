process BUILD_MODEL {
    label 'process_single'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.output}/${params.model_outdir}", mode: params.publish_dir_mode

    input:
        path profiles
        path contigs
        val ref
        val model
        val mode

    output:
        path "model.pkl",   emit: model

    script:
    """
    echo "BUILDING A NEW MODEL"
    
    build_model.py \
        -i ${profiles} \
        -c ${contigs} \
        -r ${ref} \
        -m ${model} \
        --mode ${mode} \
        -o "model.pkl"
    """

    stub:
    """
    touch model.pkl
    """
}