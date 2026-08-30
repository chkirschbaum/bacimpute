#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PRODIGAL_TRAIN }                  from './modules/local/prodigal'

include { CHEWBBACA_CREATESCHEMA }          from './modules/nf-core/chewbbaca/createschema'
include { 
    CHEWBBACA_ALLELECALL as CHEWBBACA_ALLELECALL_CREATE;
    CHEWBBACA_ALLELECALL as CHEWBBACA_ALLELECALL_BUILD;  
    CHEWBBACA_ALLELECALL as CHEWBBACA_ALLELECALL_IMPUTE 
}                                           from './modules/nf-core/chewbbaca/allelecall'

include { CHEWBBACA_DOWNLOADSCHEMA }        from './modules/local/chewbbaca/downloadschema'
include { CHEWBBACA_PREPEXTERNALSCHEMA }    from './modules/local/chewbbaca/prepexternalschema'
include { CHEWBBACA_EXTRACTCGMLST }         from './modules/local/chewbbaca/extractcgmlst'
include { CHEWBBACA_ALLELECALLEVALUATOR }   from './modules/local/chewbbaca/allelecallevaluator'

include { BUILD_MODEL }                     from './modules/local/build_model'
include { IMPUTE_PROFILES }                 from './modules/local/impute_profiles'
include { MASKING_PROFILES }                from './modules/local/masking_profiles'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // Help message
    if (params.help) { exit 0, helpMSG() }

    // Validate params
    validateParameters()

    if (params.training_file.endsWith('.fna')){
        ch_training_genome = channel.fromPath( file("${params.training_file}", checkIfExists: true) )

        PRODIGAL_TRAIN(
            ch_training_genome, 
            params.species
        )

        ch_training_file = PRODIGAL_TRAIN.out.training_file
    } else if (params.training_file.endsWith('.trn')) {
        ch_training_file = channel.fromPath( file("${params.training_file}", checkIfExists: true) )
    }

    if ((params.new_schema && params.schema_source == 'create') || params.new_model) {
        def ref_path = file(params.ref_dir)
        def ref_pattern = ref_path.isDirectory()
            ? "${ref_path}/*.{fa,fna,fasta,fa.gz,fna.gz,fasta.gz}"
            : params.ref_dir

        ch_ref_fastas = channel
            .fromPath(ref_pattern, checkIfExists: true)
            .collect()
            .map { fastas -> tuple([id: 'reference_set'], fastas) }
    }

    if (params.new_schema) {
        if (params.schema_source == 'chewie-ns') {
            //
            // Download schema from Chewie-NS
            //
            CHEWBBACA_DOWNLOADSCHEMA(
                [id: 'new_schema'],
                params.species,
                1,
                'new_schema'
            )
            ch_schema = CHEWBBACA_DOWNLOADSCHEMA.out.schema_dir
        }

        if (params.schema_source == 'external') {
            //
            // Import schema into chewBBACA
            //
            ch_external_schema = channel
                .fromPath(params.ext_schema_dir, checkIfExists: true)
                .map { schema -> tuple([id: 'external_schema'], schema) }

            CHEWBBACA_PREPEXTERNALSCHEMA(
                ch_external_schema,
                'new_schema',
                ch_training_file
            )
            ch_schema = CHEWBBACA_PREPEXTERNALSCHEMA.out.schema_dir
        }

        if (params.schema_source == 'create') {
            //
            // Create a new schema with chewBBACA
            //
            ch_empty = channel.value([])

            CHEWBBACA_CREATESCHEMA(
                ch_ref_fastas,
                ch_training_file,
                ch_empty
            )
            ch_schema = CHEWBBACA_CREATESCHEMA.out.schema
            
            CHEWBBACA_ALLELECALL_CREATE(
                ch_ref_fastas,
                ch_schema
            )

            ch_create_alleles = CHEWBBACA_ALLELECALL_CREATE.out.alleles.map {
                meta, alleles -> alleles
            }

            CHEWBBACA_EXTRACTCGMLST(
                ch_create_alleles
            )
            ch_loci_95    = CHEWBBACA_EXTRACTCGMLST.out.cgmlstschema95
            ch_loci_100   = CHEWBBACA_EXTRACTCGMLST.out.cgmlstschema100
        } 
    } else {
        ch_schema = channel
            .fromPath(params.schema_dir, checkIfExists: true)
            .map { schema -> tuple([id: 'schema'], schema) }
    }

    //
    // Build a new model
    //
    if (params.new_model) {
        CHEWBBACA_ALLELECALL_BUILD(
            ch_ref_fastas,
            ch_schema
        )

        ch_contigs = CHEWBBACA_ALLELECALL_BUILD.out.contigs_info.map {
            meta, contigs -> contigs
        }

        ch_build_profiles = CHEWBBACA_ALLELECALL_BUILD.out.alleles.map {
            meta, alleles -> alleles
        }

        ch_allelecall_dir = CHEWBBACA_ALLELECALL_BUILD.out.stats
        .mix(
            CHEWBBACA_ALLELECALL_BUILD.out.contigs_info,
            CHEWBBACA_ALLELECALL_BUILD.out.alleles,
            CHEWBBACA_ALLELECALL_BUILD.out.log,
            CHEWBBACA_ALLELECALL_BUILD.out.paralogous_counts,
            CHEWBBACA_ALLELECALL_BUILD.out.paralogous_loci,
            CHEWBBACA_ALLELECALL_BUILD.out.cds_coordinates,
            CHEWBBACA_ALLELECALL_BUILD.out.invalid_cds,
            CHEWBBACA_ALLELECALL_BUILD.out.loci_summary_stats
        )
        .map { meta, result -> result }
        .collect()

        if (params.mode == 'clonal'){
            //
            // Allele Call Evaluation - if working on 100% core
            //
            CHEWBBACA_ALLELECALLEVALUATOR(
                ch_allelecall_dir,
                ch_schema
            )
            ch_profiles = CHEWBBACA_ALLELECALLEVALUATOR.out.cgMLST_profiles
        } else {
            ch_profiles = ch_build_profiles
        }

        BUILD_MODEL(
            ch_profiles, 
            ch_contigs, 
            params.ref,
            params.model,
            params.mode
        )
        ch_model = BUILD_MODEL.out.model
    } else {
        ch_model = channel.fromPath( file("${params.model_file}", checkIfExists: true) )
    }

    //
    // Testing of a new model / species
    //
    if (params.test){
        ch_test_profiles = channel.fromPath( "${params.test_dir}", checkIfExists: true )
        ch_report_template = channel.value( file("${projectDir}/bin/utils/test_report_template.html", checkIfExists: true) )
        MASKING_PROFILES(
            ch_profiles, 
            ch_test_profiles,
            params.model,
            ch_model,
            ch_report_template
        )
    }

    //
    // Imputation
    //
    if (params.input.endsWith('.tsv')) {
        ch_impute_profiles = channel.fromPath(params.input, checkIfExists: true)
    } else {
        def input_path = file(params.input)
        def input_pattern = input_path.isDirectory()
            ? "${input_path}/*.{fa,fna,fasta,fa.gz,fna.gz,fasta.gz}"
            : params.input

        ch_input_fastas = channel
            .fromPath(input_pattern, checkIfExists: true)
            .collect()
            .map { fastas -> tuple([id: 'imputation'], fastas) }

        CHEWBBACA_ALLELECALL_IMPUTE(
            ch_input_fastas,
            ch_schema
        )

        ch_impute_profiles = CHEWBBACA_ALLELECALL_IMPUTE.out.alleles.map {
            meta, profiles -> profiles
        }
    }

    IMPUTE_PROFILES (
        ch_impute_profiles,
        params.model,
        ch_model
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateParameters() {
    Set valid_params = [
        'cores', 'max_cores', 'memory', 'version', 'help', 'profile',  
        'input', 'mode', 'new_model', 'model', 'model_file', 'ref',
        'new_schema', 'schema_source', 'schema_dir', 
        'species', 'ext_schema_dir', 'ref_dir', 'training_file', 
        'test', 'test_dir',
        'output', 'model_outdir', 'schema_outdir', 
        'allelecall_outdir', 'impute_outdir', 'masking_outdir',
        'runinfo_dir', 'publish_dir_mode', 
        'conda_cache_dir', 'singularity_cache_dir',
        'cloudProcess', 'cloud-process', 'trace_timestamp'
    ]

    def parameter_diff = params.keySet() - valid_params
    if (parameter_diff.size() != 0) {
        error("Parameter(s) $parameter_diff is/are not valid in the pipeline!\n")
    }

    if (params.profile) { 
        error("--profile is wrong please use -profile!\n")
    }

    if (!params.input) {
        error("Missing required parameter: --input")
    }

    if (!(params.mode in ['clonal', 'variable'])) {
        error("Invalid mode: '${params.mode}'. Choose 'clonal' or 'variable'.")
    }

    if (!(params.model in ['imc', 'maxst'])) {
        error("Invalid model: '${params.model}'. Choose 'imc' or 'maxst'.")
    }

    if (params.new_schema) {
        if (!(params.schema_source in ['chewie-ns', 'external', 'create'])) {
            error("--new_schema requires --schema_source chewie-ns, external, or create.")
        }

        if (params.schema_source == 'chewie-ns' && !params.species) {
            error("--schema_source chewie-ns requires --species.")
        }

        if (params.schema_source == 'external' && !params.ext_schema_dir) {
            error("--schema_source external requires --ext_schema_dir.")
        }

        if (params.schema_source == 'create' && !params.ref_dir) {
            error("--schema_source create requires --ref_dir.")
        }

        if (
            params.schema_source in ['external', 'create'] &&
            !params.training_file
        ) {
            error("--schema_source ${params.schema_source} requires --training_file.")
        }
    }
    else if (!params.schema_dir) {
        error("An existing --schema_dir is required unless --new_schema is enabled.")
    }

    if (params.new_model) {
        if (!params.ref_dir) {
            error("--new_model requires --ref_dir.")
        }

        if (!params.ref) {
            error("--new_model requires --ref.")
        }
    }
    else if (!params.model_file) {
        error(
            "An existing --model_file is required unless --new_model is enabled."
        )
    }

    if (
        params.training_file &&
        !(params.training_file.endsWith('.trn') ||
        params.training_file.endsWith('.fna'))
    ) {
        error("--training_file must end in .trn or .fna.")
    }
    if (params.training_file.endsWith('.fna') && !params.species) {
        error("Creating a Prodigal training file from .fna requires --species.")
    }

    if (params.test && !params.new_model) {
        error("--test currently requires --new_model.")
    }
    if (params.test && !params.test_dir) {
        error("--test requires --test_dir.")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMSG() {
    def c_green = "\033[0;32m"
    def c_reset = "\033[0m"
    def c_yellow = "\033[0;33m"
    def c_blue = "\033[0;34m"
    def c_dim = "\033[2m"
    log.info """
    ____________________________________________________________________________________________
    
    ${c_blue}Robert Koch Institute, MF1 Bioinformatics${c_reset}

    Workflow: BACimpute

    Imputation of Bacterial Genomes based on their cgMLST profiles

    ${c_yellow}Usage examples:${c_reset}
    nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local --input 'test/profiles.tsv'

    ${c_yellow}Input options for runs:${c_reset}
    ${c_green} --input ${c_reset}           REQUIRED! Path to the input directory or file.
                        [ default: $params.input ]
    ${c_green} --mode ${c_reset}            Define if you are working on a clonal or variable species.
                        [ default: $params.mode ]

    ${c_green} --schema_dir ${c_reset}      Path to the directory of the existing schema.
                        [ default: $params.schema_dir ]
    ${c_green} --new_schema ${c_reset}      Set to true if a new schema should be generated.
                        [ default: $params.new_schema ]
    ${c_green} --schema_source ${c_reset}   Define the source of your schema. Choose 'chewie-ns' to download a schema
                        from there, 'external' to import a schema from Ridom or Enterobase, and 'create' to build a new
                        individual schema from your reference set.
                        [ default: $params.schema_source ]
    ${c_green} --species ${c_reset}         Species name or ID to download schema for.
                        Required for --schema_source 'chewie-ns'.
                        More information available in the README.
                        [ default: $params.species ]
    ${c_green} --ext_schema_dir ${c_reset}  Path to the directory of the external schema to import.
                        Required for --schema_source 'external'.
                        [ default: $params.ext_schema_dir ]
    ${c_green} --ref_dir ${c_reset}         Path to the directory with the sequences which should be used to create 
                        the new schema.
                        Required for --schema_source 'create'.
                        [ default: $params.ref_dir ]
    ${c_green} --training_file ${c_reset}   Either the GCF file to build a training file on or a training file for
                        the species.
                        Required for --schema_source 'create'.
                        [ default: $params.training_file ]

    ${c_green} --model ${c_reset}           Define if you are using a 'imc' or a 'maxst' model.
                        [ default: $params.model ]
    ${c_green} --model_file ${c_reset}      Path to the model file.
                        [ default: $params.model_file ]
    ${c_green} --new_model ${c_reset}       Set to true if a new model should be generated.
                        [ default: $params.new_model ]
    ${c_green} --ref ${c_reset}             Name of the sequence that should be used as reference.
                        Required for --new_model true.
                        [ default: $params.ref ]

    ${c_green} --test ${c_reset}          Species to download schema for. Required for --schema_source 'chewie-ns'.
                        More information available in the README.
                        [ default: $params.test ]

    ${c_yellow}Computing options:${c_reset}
    ${c_green} --cores ${c_reset}           Max cores per process for local use 
                        [ default: $params.cores ]
    ${c_green} --max_cores ${c_reset}       Max cores used on the machine for local use 
                        [ default: $params.max_cores ]
    ${c_green} --memory ${c_reset}          Max memory in GB for local use
                        [ default: $params.memory ]
    
    ${c_yellow}Output options:${c_reset}
    ${c_green} --output ${c_reset}                  Name of the result folder 
                                [ default: $params.output ]
    ${c_green} --publish_dir_mode ${c_reset}        Mode of output publishing: 'copy', 'symlink' 
                                [ default: $params.publish_dir_mode ]
                                ${c_dim}With 'symlink' results are lost when removing the work directory.${c_reset}
    
    ${c_yellow}Caching:${c_reset}
    ${c_green} --conda_cache_dir ${c_reset}         Location for storing the conda environments 
                                [ default: $params.conda_cache_dir ]
    
    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executors${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executor${c_reset} (choose one):
        local
        slurm
    
    ${c_blue}Engines${c_reset} (choose one):
        conda
        mamba
    """
}
