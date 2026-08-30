#!/usr/bin/env python

import argparse

import pandas as pd

from utils.utils import load_imc
from utils.utils import load_maxst

from utils.data import read_chewie
from utils.data import freq_rank_alleles
from utils.data import get_missing_idx

from utils.imc import impute_missing as imc_impute

from utils.maxst import impute_missing as maxst_impute


def main():
    """
    Main function

    Returns:
        None
    """
    parser = argparse.ArgumentParser(
        description='Impute on a profile with missing alleles',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
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
        '--cpu',
        required=True,
        help='CPU cores to run the process'
    )
    parser.add_argument(
        '-o', '--output',
        default='imputed.tsv',
        help='TSV file with the imputed profiles'
    )
    parser.add_argument(
        '-r', '--report',
        default='imputation.html',
        help='HTML report for imputation'
    )
    args = parser.parse_args()

    profiles_df = read_chewie(
        path=args.input
    )

    if args.model == "imc":
        species, freq_rank_dict, transition_array, _, cg_list, ref = load_imc(
            path=args.file,
        )
    elif args.model == "maxst":
        species, freq_rank_dict, rooted_mst, n_maxalleles, node_to_genome_idx, tree_idx_to_node, genome_idx_to_node, cg_list, ref = load_maxst(
            path=args.file,
        )

    profiles_df = profiles_df[
        [loci for loci in cg_list if loci in profiles_df.columns]
    ]

    freq_df, _ = freq_rank_alleles(
        profiles=profiles_df, 
        ref_rank=freq_rank_dict, 
        mode="run"
    )
    freq_array = freq_df.to_numpy()

    missing_idx = get_missing_idx(
        profiles=freq_df,
    )

    if args.model == "imc":
        imputed_array = imc_impute(
            transition=transition_array, 
            samples=freq_array, 
            missing_idx=missing_idx,
            n_jobs=args.cpu,
        )
    elif args.model == "maxst":
        imputed_array = maxst_impute(
            rooted_mst=rooted_mst,
            samples=freq_array,
            missing_idx=missing_idx,
            n_maxalleles=n_maxalleles,
            node_to_genome_idx=node_to_genome_idx,
            tree_idx_to_node=tree_idx_to_node,
            genome_idx_to_node=genome_idx_to_node,
            n_jobs=args.cpu, 
        )

    imputed_df = pd.DataFrame(
        imputed_array,
        index=freq_df.index,
        columns=freq_df.columns,
    )

    imputed_df.to_csv(
        args.output,
        sep='\t'
    )


if __name__ == "__main__":
    main()