#!/usr/bin/env python

import argparse

from jinja2 import Template

from utils.utils import load_imc
from utils.utils import load_maxst

from utils.data import read_chewie
from utils.data import freq_rank_alleles
from utils.data import get_missing_idx
from utils.data import mask_per_sample
from utils.data import mask_per_locus

from utils.imc import impute_missing as imc_impute

from utils.maxst import impute_missing as maxst_impute

from utils.report import accuracy
from utils.report import plot_summary


def main():
    """
    Main function

    Returns:
        None
    """
    parser = argparse.ArgumentParser(
        description='Mask and impute on profiles',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--ri', '--train_input',
        required=True,
        dest="train_input",
        help='Path to results_alleles.tsv from chewBBACA'
    )
    parser.add_argument(
        '--ti', '--test_input',
        required=True,
        dest="test_input",
        help='Path to results_alleles.tsv from chewBBACA'
    )
    parser.add_argument(
        '-m', '--model',
        required=True,
        help='Use an IMC or a MaxST model'
    )
    parser.add_argument(
        '-f', '--file',
        required=True,
        help='Path to the model file (.pkl)'
    )
    parser.add_argument(
        '--mask',
        default='sample',
        help='Masking per sample or per locus'
    )
    parser.add_argument(
        '--cpu',
        required=True,
        help='CPU cores to run the process'
    )
    parser.add_argument(
        '--template',
        required=True,
        help='Path to the HTML report template'
    )
    parser.add_argument(
        '-r', '--report',
        default='imputation.html',
        help='HTML report for imputation'
    )
    args = parser.parse_args()

    profiles_df = read_chewie(
        path=args.train_input
    )

    if args.model == "imc":
        species, freq_rank_dict, transition_array, _, cg_list, ref = load_imc(
            path=args.file,
        )
    elif args.model == "maxst":
        species, freq_rank_dict, rooted_mst, n_maxalleles, node_to_genome_idx, tree_idx_to_node, genome_idx_to_node, cg_list, ref = load_maxst(
            path=args.file,
        )

    allele_counts = [len(freq_rank_dict[locus]) for locus in cg_list]
    model_summary = {
        "model": {"imc": "IMC", "maxst": "MaxST"}[args.model],
        "n_loci": len(cg_list),
        "min_alleles": min(allele_counts),
        "max_alleles": max(allele_counts),
        "mean_alleles": sum(allele_counts) / len(allele_counts),
        "reference": ref,
        # Placeholder until species is stored in and loaded from the model file.
        "species": "Not available",
    }

    profiles_df = profiles_df[
        [loci for loci in cg_list if loci in profiles_df.columns]
    ]

    freq_df, _ = freq_rank_alleles(
        profiles=profiles_df, 
        ref_rank=freq_rank_dict, 
        mode="run"
    )
    freq_array = freq_df.to_numpy()

    test_df = read_chewie(
        args.test_input
    )

    test_df = test_df[
        test_df.columns.intersection(freq_df.columns)
    ]

    test_df, _ = freq_rank_alleles(
        profiles=test_df, 
        ref_rank=freq_rank_dict, 
        mode="run"
    )

    test_array = test_df.to_numpy()

    thresh = [0.05, 0.1, 0.15, 0.2, 0.25]

    train_model_all = []
    train_baseline_all = []
    test_model_all = []
    test_baseline_all = []
    train_accuracy = []
    test_accuracy = []

    for t in thresh:
        if args.mask == "sample":
            masked_df = mask_per_sample(freq_df, percentage=t)
        elif args.mask == "locus":
            masked_df = mask_per_locus(freq_df, percentage=t)
        masked_array = masked_df.to_numpy()

        baseline_df = masked_df.copy()
        baseline_df.fillna(1, inplace=True)
        baseline_array = baseline_df.to_numpy()

        missing_idx = get_missing_idx(
            profiles=masked_df,
        )

        if args.model == "imc":
            imputed_array = imc_impute(
                transition=transition_array, 
                samples=masked_array, 
                missing_idx=missing_idx,
                n_jobs=args.cpu,
            )
        elif args.model == "maxst":
            imputed_array = maxst_impute(
                rooted_mst=rooted_mst,
                samples=masked_array,
                missing_idx=missing_idx,
                n_maxalleles=n_maxalleles,
                node_to_genome_idx=node_to_genome_idx,
                tree_idx_to_node=tree_idx_to_node,
                genome_idx_to_node=genome_idx_to_node,
                n_jobs=args.cpu, 
            )

        acc_baseline, acc_per_col_baseline, acc_model, acc_per_col_model = accuracy(
            original=freq_array,
            imputed_baseline=baseline_array,
            imputed_model=imputed_array,
            masked=masked_array
        )
        train_baseline_all.extend(acc_per_col_baseline)
        train_model_all.extend(acc_per_col_model)
        train_accuracy.append({
            "threshold": f"{t:.0%}",
            "baseline": acc_baseline,
            "model": acc_model,
        })

    fig_train = plot_summary(
        n_loci=len(freq_df.columns),
        imputed_baseline_all=train_baseline_all,
        imputed_model_all=train_model_all,
        model=args.model,
    )

    for t in thresh:
        if args.mask == "sample":
            masked_df = mask_per_sample(test_df, percentage=t)
        elif args.mask == "locus":
            masked_df = mask_per_locus(test_df, percentage=t)
        masked_array = masked_df.to_numpy()

        baseline_df = masked_df.copy()
        baseline_df.fillna(1, inplace=True)
        baseline_array = baseline_df.to_numpy()

        missing_idx = get_missing_idx(
            profiles=masked_df,
        )

        if args.model == "imc":
            imputed_array = imc_impute(
                transition=transition_array, 
                samples=masked_array, 
                missing_idx=missing_idx,
                n_jobs=args.cpu,
            )
        elif args.model == "maxst":
            imputed_array = maxst_impute(
                rooted_mst=rooted_mst,
                samples=masked_array,
                missing_idx=missing_idx,
                n_maxalleles=n_maxalleles,
                node_to_genome_idx=node_to_genome_idx,
                tree_idx_to_node=tree_idx_to_node,
                genome_idx_to_node=genome_idx_to_node,
                n_jobs=args.cpu, 
            )

        acc_baseline, acc_per_col_baseline, acc_model, acc_per_col_model = accuracy(
            original=test_array,
            imputed_baseline=baseline_array,
            imputed_model=imputed_array,
            masked=masked_array
        )
        test_baseline_all.extend(acc_per_col_baseline)
        test_model_all.extend(acc_per_col_model)
        test_accuracy.append({
            "threshold": f"{t:.0%}",
            "baseline": acc_baseline,
            "model": acc_model,
        })

    fig_test = plot_summary(
        n_loci=len(test_df.columns),
        imputed_baseline_all=test_baseline_all,
        imputed_model_all=test_model_all,
        model=args.model,
    )

    plotly_jinja_data = {
        "violin_train":fig_train.to_html(full_html=False),
        "violin_test":fig_test.to_html(full_html=False),
        "train_accuracy": train_accuracy,
        "test_accuracy": test_accuracy,
        "train_n_sequences": len(freq_df),
        "test_n_sequences": len(test_df),
        "model_summary": model_summary,
    }

    with open(args.report, "w", encoding="utf-8") as output_file:
        with open(args.template, encoding="utf-8") as template_file:
            template = Template(template_file.read())
            output_file.write(template.render(plotly_jinja_data))


if __name__ == "__main__":
    main()
