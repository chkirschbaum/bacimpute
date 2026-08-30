#!/usr/bin/env python

import argparse

from utils.preprocess import mutual_info_matrix

from utils.data import read_chewie
from utils.data import order_to_ref
from utils.data import freq_rank_alleles

from utils.imc import build_markov_chain

from utils.maxst import build_rooted_mst
from utils.maxst import get_transitions

from utils.utils import save_imc
from utils.utils import save_maxst
from utils.utils import save_cytoscape


def main():
    """
    Main function

    Returns:
        None
    """
    parser = argparse.ArgumentParser(
        description='Build a new Markov Chain or Maximum Spanning Tree for Impuation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to cgMLST_profiles.tsv / results_alleles.tsv from chewBBACA'
    )
    parser.add_argument(
        '-c', '--contigs',
        required=True,
        help='Path to results_contigsInfo.tsv from chewBBACA'
    )
    parser.add_argument(
        '-r', '--ref',
        required=True,
        help='Name of the reference sequence in the dataset'
    )
    parser.add_argument(
        '-m', '--model',
        required=True,
        help='Build an IMC or a MaxST model'
    )
    parser.add_argument(
        '--mode',
        default='clonal',
        help='Running on a clonal or variable species'
    )
    parser.add_argument(
        "--species",
        default="Not available",
    )
    parser.add_argument(
        '--cytoscape',
        default=False,
        help='Generate a JSON file to visualize in cytoscape'
    )
    parser.add_argument(
        '-o', '--output',
        default='model.pkl',
        help='Pickle file with the new Maximum Spanning Tree model'
    )
    args = parser.parse_args()

    if args.mode == "variable":
        profiles_df = read_chewie(
            path=args.input,
            species_mode="variable",
            inf_mode="keep"
        )
        profiles_df = order_to_ref(
            profiles_df=profiles_df, 
            ref_id=args.ref, 
            contigs_path=args.contigs
        )

        freq_df, freq_rank_dict = freq_rank_alleles(
            profiles=profiles_df,
            cg=95
        )
    
    else:
        profiles_df = read_chewie(
            path=args.input,
            inf_mode="keep"
        )
        profiles_df = order_to_ref(
            profiles_df=profiles_df, 
            ref_id=args.ref, 
            contigs_path=args.contigs
        )

        freq_df, freq_rank_dict = freq_rank_alleles(
            profiles=profiles_df
        )

    loci_list = profiles_df.columns

    if args.model == "imc":
        transition, max_alleles_array = build_markov_chain(
            freq_rank_profiles=freq_df,
        )

        save_imc(
            species=args.species,
            freq_rank_dict=freq_rank_dict,
            transition_array=transition,
            max_alleles_array=max_alleles_array,
            cg_list=profiles_df.columns,
            ref=args.ref,
            path=args.output,
        )

    elif args.model == "maxst":
        mi_df = mutual_info_matrix(
            freq_profile_df=freq_df,
        )

        rooted_mst, node_to_genome_idx, tree_idx_to_node, genome_idx_to_node = build_rooted_mst(
            mi_df=mi_df,
            loci_order=loci_list,
            root=loci_list[0],
        )
        rooted_mst, n_maxalleles = get_transitions(
            freq_rank_profiles=freq_df, 
            rooted_mst=rooted_mst,
        )

        save_maxst(
            species=args.species,
            freq_rank_dict=freq_rank_dict,
            rooted_mst=rooted_mst,
            n_maxalleles=n_maxalleles,
            node_to_genome_idx=node_to_genome_idx,
            tree_idx_to_node=tree_idx_to_node,
            genome_idx_to_node=genome_idx_to_node,
            cg_list=loci_list,
            ref=args.ref,
            path=args.output,
        )

        if args.cytoscape:
            save_cytoscape(
                G=rooted_mst,
            )



if __name__ == "__main__":
    main()