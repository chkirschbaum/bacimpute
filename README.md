<div id="top"></div>

<div align="center">
<h1 align="center"> BACimpute </h1>
<h3 align="center"> cgMLST Based Imputation for Bacterial Genes </h3>
</div>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

BACimpute predicts missing allele calls in bacterial core-genome multilocus sequence typing (cgMLST) profiles. It combines [chewBBACA](https://github.com/B-UMMI/chewBBACA) schema and allele-calling workflows with statistical models trained on a reference collection of bacterial genomes.

**! IMPORTANT !**
**BACimpute is under active development and is not yet a production release.**

- [How it works](#how-it-works)
- [Getting Started](#getting-started)
  - [Quick Installation](#quick-installation)
  - [Get / Update BACimpute](#get--update-bacimpute)
  - [Call help](#call-help)
- [Running BACimpute](#running-bacimpute)
  - [Get a cgMLST Scheme](#get-a-cgmlst-scheme)
    - [Download from Chewie-NS](#download-from-chewie-ns)
    - [Import from Ridom or Enterobase](#import-from-ridom-or-enterobase)
    - [Create a new Schema](#create-a-new-schema)
  - [Get an new Imputation Model](#get-an-new-imputation-model)
    - [Test your new Model](#test-your-new-model)
- [Parameter List](#parameter-list)
- [Citations](#citations)
- [Contact](#contact)


# How it works

The workflow has the following main functionalities:

1. Obtain a cgMLST schema by downloading one from Chewie-NS, adapting an external schema, creating a schema from a reference genome collection, or supplying an existing chewBBACA schema.
2. Train and optionally evaluate an imputation model from the resulting allelic profiles.
3. Call alleles for new assemblies and impute loci reported as missing.

BACimpute provides two models:

- **Inhomogeneous Markov chain (IMC):** uses the physical order of loci in a reference genome. A missing locus is inferred from neighbouring loci; consecutive missing loci are handled with a forward-backward procedure.
- **Maximum spanning tree (MaxST):** connects loci according to normalized mutual information. Missing values are inferred by passing information between parent and child loci in the rooted tree.

Before model fitting, allele identifiers are converted to frequency ranks learned from the reference collection. Alleles not observed during training are represented as new alleles, while chewBBACA `LNF` calls are treated as missing values.


# Getting Started

## Quick Installation

To run the pipeline, you need to have `Nextflow` and either `conda`, `Docker` or `Singularity`.

<details><summary><strong>Click!</strong> If you want to install <code>Nextflow</code> directly, you can use the following one-liner. </summary>

```bash
wget -qO- https://get.nextflow.io | bash
```
</details>

<details><summary><strong>Click!</strong> If you want to set up <code>conda</code> to run the pipeline and install all other dependencies through it, you can use the following steps. </summary>

Use the following bash commands if you are working on **Linux**:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Use the following bash commands if you are working on **Mac**:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

Then, `Nextflow` an be installed over `conda`:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```
</details>

## Get / Update BACimpute

```bash
nextflow pull chkirschbaum/bacimpute
```

## Call help

```bash
nextflow run chkirschbaum/bacimpute -r <version> --help
```

# Running BACimpute

BACimpute takes as input either a directory with fasta sequences or a `results_alleles.tsv` from `chewBBACA AlleleCall`.

```bash
nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local \
    --input 'profiles.tsv'
```

## Get a cgMLST Scheme

### Download from Chewie-NS

```bash
nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local \
    --input 'profiles.tsv' \
    --new_schema \
    --schema_source 'chewie-ns' \
    --species 'Listeria monocytogenes' \
    --model_file 'path/to/model'
```

<details><summary><strong>Click!</strong> If you want to see which schemes you can download from Chewie-NS. </summary>

The following schemes are available from Chewie-NS:

| Species                    | Species ID (`-sp`) | Schema ID (`-sc`) |
| -------------------------- | -----------------: | ----------------: |
| *Streptococcus pyogenes*   |                  1 |                 1 |
| *Acinobacter baumanii*     |                  2 |                 1 |
| *Aliarcobacter butzleri*   |                  3 |                 1 |
| *Bacillus anthracis*       |                  4 |                 1 |
| *Brucella*                 |                  5 |                 1 |
| *Campylobacter jejuni*     |                  6 |                 1 |
| *Clostridium chauvoei*     |                  7 |                 1 |
| *Clostridium neonatale*    |                  8 |                 1 |
| *Clostridium perfringens*  |                  9 |                 1 |
| *Escherichia coli*         |                 10 |                 1 |
| *Klebsiella oxytoca*       |                 11 |                 1 |
| *Legionella longbeachae*   |                 12 |                 1 |
| *Neisseria meningitidis*   |                 13 |                 1 |
| *Salmonella enterica*      |                 14 |                 1 |
| *Shewanella*               |                 15 |                 1 |
| *Yersinia enterocolitica*  |                 16 |                 1 |
| *Chlamydia trachomatis*    |                 17 |                 1 |
| *Listeria monocytogenes*   |                 18 |                 1 |
| *Legionella pneumophila*   |                 19 |                 1 |
| *Yersinia pestis*          |                 20 |                 1 |
</details>

### Import from Ridom or Enterobase

```bash
nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local \
    --input 'profiles.tsv' \
    --new_schema \
    --schema_source 'external' \
    --ext_schema_dir 'path/to/external_schema' \
    --model_file 'path/to/model'
```

### Create a new Schema

Creating a schema requires either:

- a Prodigal training file ending in `.trn`; or
- a representative genome (GCF from NCBI RefSeq) ending in `.fna`, from which BACimpute generates a training file.

Several training files are provided in [`data/prodigal_training_files`](data/prodigal_training_files) or via `chewBBACA`.

```bash
nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local \
    --input 'profiles.tsv' \
    --new_schema \
    --schema_source 'create' \
    --ref_dir 'path/to/reference_set' \
    --training_file 'path/to/training_or_GCF_file' \
    --model_file 'path/to/model'
```

## Get an new Imputation Model

```bash
nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local \
    --input 'profiles.tsv' \
    --schema_dir 'path/to/schema' \
    --model 'imc' \
    --new_model \
    --ref 'refseq_name'
```

### Test your new Model

You can test your new built model on the reference set you build the model on and a test set and get an interactive HTML report with the results.

```bash
nextflow run chkirschbaum/bacimpute -r <version> -profile conda,local \
    --input 'profiles.tsv' \
    --schema_dir 'path/to/schema' \
    --model 'imc' \
    --new_model \
    --ref 'refseq_name' \
    --ref_dir 'path/to/reference_set' \
    --test_dir 'path/to/test_set'
```


# Parameter List

| Parameter | Default | Description |
| --- | --- | --- |
| `--input` | `''` | FASTA assembly or quoted glob of assemblies to impute. |
| `--mode` | `clonal` | Species handling: `clonal` or `variable`. |
| `--schema_dir` | `''` | Existing chewBBACA schema directory. Used when `--new_schema` is false. |
| `--new_schema` | `false` | Obtain or create a schema during the run. |
| `--schema_source` | `''` | Schema route: `chewie-ns`, `external`, or `create`. |
| `--species` | `''` | Species identifier for Chewie-NS and name used for a generated Prodigal training file. |
| `--external_schema_dir` | `''` | External schema to adapt with chewBBACA. |
| `--ref_dir` | `''` | Reference genome collection used for schema creation and/or model training. |
| `--training_file` | `''` | Prodigal `.trn` file or `.fna` genome used to produce one. |
| `--model` | `maxst` | Model type: `maxst` or `imc`. |
| `--model_file` | `''` | Existing `.pkl` model. Used when `--new_model` is false. |
| `--new_model` | `false` | Train a model from the reference collection. |
| `--ref` | `''` | Reference sample identifier used to order loci. Required for training. |
| `--test` | `false` | Run masking-based model evaluation. |
| `--output` | `results` | Top-level result directory. |
| `--help` | `false` | Print the command-line help message. |

# Citations

COMING SOON

Additionally, an extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

# Contact

Did you find a bug? 🐛 Suggestion / Feedback / Feature request? 👨‍💻
Please visit [GitHub Issues](https://github.com/rki-mf1/viruswarn-sc2/issues)

For business inquiries or professional support requests, please feel free to contact us!