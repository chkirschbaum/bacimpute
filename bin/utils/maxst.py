#!/usr/bin/env python

import numpy as np
import pandas as pd

import networkx as nx

from joblib import Parallel, delayed


def build_rooted_mst(
        mi_df: pd.DataFrame,
        loci_order: list,
        root: str,
        algo: str = 'prim',
) -> tuple[nx.DiGraph, dict, dict, dict, dict]:
    """
    Build a rooted maximum spanning tree from a given normalized mutual information matrix.

    Args:
        path (str): Path to the normalized mutual information matrix file.
        loci_order (list): List of loci in the order they appear in the genome.
        root (str): Name of the locus to use as the root node of the tree. Normally the first locus in the genome.
        algo (str): The algorithm to use for building the maximum spanning tree. Defaults to 'prim'. 
            Other options are 'kruskal' and 'boruvka'.
    
    Returns:
        tuple[nx.DiGraph, dict, dict, dict, dict]: A tuple containing the rooted maximum spanning tree and the mapping dictionaries.
    """
    mask = np.triu(np.ones(mi_df.shape, dtype=bool), k=1)

    edge_df = (
        mi_df.where(mask)
        .stack()
        .rename("weight")
        .reset_index()
        .dropna(subset=["weight"])
    )
    edge_df.columns = ["source", "target", "weight"]

    graph = nx.from_pandas_edgelist(
        edge_df,
        source='source',
        target='target',
        edge_attr='weight',
        create_using=nx.Graph()
    )

    mst = nx.maximum_spanning_tree(graph, algorithm=algo)

    rooted_mst = nx.bfs_tree(mst, source=root)

    genome_idx = {node: i for i, node in enumerate(loci_order)}
    nx.set_node_attributes(rooted_mst, genome_idx, "genome_index")
    node_to_genome_idx = {node: att["genome_index"] for node, att in rooted_mst.nodes(data=True)}
    genome_idx_to_node = {att["genome_index"]: node for node, att in rooted_mst.nodes(data=True)}

    tree_idx = {node: i for i, node in enumerate(list(rooted_mst.nodes()))}
    nx.set_node_attributes(rooted_mst, tree_idx, "tree_index")
    tree_idx_to_node = {att["tree_index"]: node for node, att in rooted_mst.nodes(data=True)}

    return rooted_mst, node_to_genome_idx, tree_idx_to_node, genome_idx_to_node


def get_transitions(
        freq_rank_profiles: pd.DataFrame,
        rooted_mst: nx.Graph,
        epsilon: float = 0.02,
) -> nx.DiGraph:
    """ 
    Building a start probability and transition matrix from locus to locus based on the given reference profile.  

    Args:
        freq_rank_profiles (pd.Dataframe): Dataframe of the reference profiles of shape (n_samples, n_loci).
        rooted_mst (nx.Graph): The rooted maximum spanning tree as defined by `build_rooted_mst()`.
        epsilon (float): Pseudo count value to avoid `0.` values in the transitions.

    Returns:
        nx.DiGraph: A directed graph representing the rooted maximum spanning tree with transition probabilities.
    """
    max_alleles_array = freq_rank_profiles.max()
    max_alleles = max(max_alleles_array)

    if max_alleles < 255:
        freq_rank_profiles = freq_rank_profiles.astype(np.uint8)
    elif max_alleles < 65535:
        freq_rank_profiles = freq_rank_profiles.astype(np.uint16)

    mst_with_transitions = rooted_mst.copy()

    for locus in rooted_mst.nodes():
        ### Parent node
        if rooted_mst.in_degree(locus) != 0:
            parent = list(rooted_mst.predecessors(locus))
            mst_with_transitions.nodes[locus]["parent"] = parent
        ### Root
        else:
            mst_with_transitions.nodes[locus]["parent"] = []

        ### Children nodes
        if rooted_mst.out_degree(locus) != 0:
            children = list(rooted_mst.successors(locus))
            mst_with_transitions.nodes[locus]["children"] = children

            for c in children:
                col1 = freq_rank_profiles[locus]
                col1_max = int(max_alleles_array[locus])
                col2 = freq_rank_profiles[c]
                col2_max = int(max_alleles_array[c])
        
                ### Get allele ID combinations for two consecutive columns
                transition_temp = pd.crosstab(col2, col1)
        
                ### Include 0 and all valid alleles
                transition_temp = transition_temp.reindex(
                    index=range(col2_max + 1),
                    columns=range(col1_max + 1),
                    fill_value=0,
                )

                ### Apply pseudocount
                transition_temp = transition_temp.astype(float) + epsilon
        
                ### Normalize
                transition_temp = transition_temp.div(
                    transition_temp.sum(axis=0),
                    axis=1,
                )
        
                ### Uniform length for all transition matrices
                transition_temp = transition_temp.reindex(
                    index=range(max_alleles + 1),
                    columns=range(max_alleles + 1),
                    fill_value=0,
                )
                mst_with_transitions.nodes[locus][c] = transition_temp.T.to_numpy(dtype=np.float32)
        ### Leaves
        else:
            mst_with_transitions.nodes[locus]["children"] = []

    return mst_with_transitions, max_alleles


def impute_single(
        sample: np.ndarray,
        locus: int,
        rooted_mst: nx.DiGraph,
) -> int:
    """
    Imputes a single missing allele.

    Args:
        sample (np.ndarray): The array (1D) of shape (n_loci) with the missing allele.
        locus (int): Position of the locus with the missing allele.
        rooted_mst (nx.DiGraph): The rooted maximum spanning tree as defined by `build_rooted_mst()`.

    Returns:
        int: The imputed allele for the missing locus.
    """
    parent = rooted_mst.nodes[locus]["parent"]
    children = rooted_mst.nodes[locus]["children"]

    ### Missing at root
    if rooted_mst.in_degree(locus) == 0:
        imputed_allele = np.prod(
            [
                rooted_mst.nodes[locus][c][
                    :, sample[rooted_mst.nodes[c]["genome_index"]]
                ]
                for c in children
            ]
        )                
        return int(imputed_allele)

    ### Missing at leave
    elif rooted_mst.out_degree(locus) == 0:
        imputed_allele = np.nanargmax(
            rooted_mst.nodes[parent[0]][locus][
                sample[rooted_mst.nodes[parent[0]]["genome_index"]], :
            ]
        )
        return int(imputed_allele)

    ### Searching MAP for missing node
    else:
        p_impute = rooted_mst.nodes[parent[0]][locus][
            sample[rooted_mst.nodes[parent[0]]["genome_index"]], :
        ]
        c_impute = np.prod(
            [
                rooted_mst.nodes[locus][c][
                    :, sample[rooted_mst.nodes[c]["genome_index"]]
                ]
                for c in children
            ]
        )   

        posterior = (p_impute * c_impute) / np.nansum(p_impute * c_impute)

        if np.nansum(posterior) == 0:
            return int(np.nanargmax(p_impute))

        return int(np.nanargmax(posterior))


def impute_stretches(
        missing_nodes_list: list,
        sample: np.ndarray,
        rooted_mst: nx.DiGraph,
        n_maxalleles: int
) -> list:
    """
    Imputes a missing stretch of alleles.

    Args:
        missing_stretch_idx (list): List of tree indices of the missing stretch.
        sample (np.ndarray): The array (1D) of shape (n_loci) with the missing alleles.
        rooted_mst (nx.DiGraph): The rooted maximum spanning tree as defined by `build_rooted_mst()`.

    Returns:
        list: The predicted alleles for the missing stretch.
    """
    postorder = missing_nodes_list[::-1]
    ### Sets for faster lookup
    missing = set(missing_nodes_list)

    ### Upward Messaging
    up_dict = {}
    for node in postorder:
        children = rooted_mst.nodes[node]["children"]
        child_messages = []

        for child in children:
            if child in missing:
                message = (
                    rooted_mst.nodes[node][child] @ up_dict[child]
                )

            else:
                ### Observed allele
                child_allele = sample[
                    rooted_mst.nodes[child]["genome_index"]
                ]
                message = rooted_mst.nodes[node][child][
                    :, child_allele
                ]

            child_messages.append(message)

        if child_messages:
            up_dict[node] = np.prod(child_messages, axis=0)
        else:
            ### Leaf node 
            up_dict[node] = np.ones(n_maxalleles)

    ### Downward Messaging
    down_dict = {}
    for node in missing_nodes_list:
        parent = rooted_mst.nodes[node]["parent"]

        if not parent:
            ### Root node
            down_dict[node] = np.ones(n_maxalleles)

        elif parent[0] not in missing:
            ### Observed allele
            parent_allele = sample[
                rooted_mst.nodes[parent[0]]["genome_index"]
            ]
            down_dict[node] = rooted_mst.nodes[parent[0]][node][
                int(parent_allele), :
            ]

        ### Integrgrate all information of subtree besides node itself
        children = rooted_mst.nodes[node]["children"]
        for child in children:
            if child not in missing:
                continue

            parent_messages = []
            parent_messages.append(down_dict[node])

            for sibling in children:
                if sibling == child:
                    continue

                elif sibling in missing:
                    message = (
                        rooted_mst.nodes[node][sibling] @ up_dict[sibling]
                    )

                else:
                    sibling_allele = sample[
                        rooted_mst.nodes[sibling]["genome_index"]
                    ]
                    message = rooted_mst.nodes[node][sibling][
                        :, sibling_allele
                    ]
                parent_messages.append(message)

            product = np.prod(parent_messages, axis=0)
            down_dict[child] = product @ rooted_mst.nodes[node][child]
        

    predictions = {}
    for node in missing_nodes_list:

        posterior = up_dict[node] * down_dict[node]

        if np.nansum(posterior) == 0:
            posterior = down_dict[node]

        predictions[node] = np.nanargmax(posterior)

    return predictions


def split_missing(
    rooted_mst: nx.DiGraph,
    missing_nodes: list,
) -> tuple[list, list]:
    """
    Split the missing tree indices into single missing alleles and stretches of missing alleles.

    The important part is that even though we have the order of the missing alleles in the tree index, 
    consecutive values in the tree index do not necessarily mean that they are in a parent-child relationship in the tree.

    Args:
        rooted_mst (nx.DiGraph): The rooted maximum spanning tree as defined by `build_rooted_mst()`.
        missing_nodes (list): The sorted list of tree indices of the missing alleles.

    Returns:
        tuple[list, list]: The first list contains the tree indices of the single missing alleles, 
            and the second list contains the lists of tree indices of the stretches of missing alleles.
    """
    missing_single = []
    missing_stretches = []

    ### Sets for faster lookup
    missing = set(missing_nodes)

    for node in missing:
        parent = rooted_mst.nodes[node]["parent"]
        ### Parent node also missing - subtree gets added in iteration for parent node
        if len(parent) != 0 and parent[0] in missing:
            continue
        

        nodes = [node]
        stretch = []

        while nodes:
            current = nodes.pop(0)

            stretch.append(rooted_mst.nodes[current]["tree_index"])

            ### Consider branching nodes
            ### Add all missing children
            nodes.extend(
                child
                for child in rooted_mst.nodes[current]["children"]
                if child in missing
            )

        if len(stretch) == 1:
            missing_single.append(stretch[0])
        else:
            missing_stretches.append(sorted(stretch))

    return missing_single, missing_stretches


def _impute_sample(
        rooted_mst: nx.DiGraph,
        sample: np.ndarray,
        missing_idx: np.ndarray,
        n_maxalleles: int,
        node_to_genome_idx: dict,
        tree_idx_to_node: dict,
        genome_idx_to_node: dict,
) -> np.ndarray:
    """ 
    Predicting missing samples with a Maximum Spanning Tree. Missing samples are marked as NA.

    Args:
        rooted_mst (nx.DiGraph): The rooted maximum spanning tree as defined by `build_rooted_mst()`.
        sample (np.ndarray): The array (1D) of shape (n_loci,) of the sample to predict.
        missing_idx (np.ndarray): The array (1D) of shape (n_missingidx,) of the missing alleles
            generated from `get_missing_idx()`.
        node_to_genome_idx (dict): A dictionary mapping the nodes of the tree to the genome indices.
        tree_idx_to_node (dict): A dictionary mapping the tree indices to the nodes of the tree.
        genome_idx_to_node (dict): A dictionary mapping the genome indices to the nodes of the tree.

    Returns:
        np.ndarray: Array (2D) of shape (n_loci,) with the predicted sample.
    """
    imputed_sample = sample.copy()

    if len(missing_idx) == 0:
        return imputed_sample

    missing_nodes = [genome_idx_to_node[idx] for idx in missing_idx]

    missing_single, missing_stretches = split_missing(rooted_mst, missing_nodes)

    for t_i in missing_single:
        ### Impute missing single alleles
        locus = tree_idx_to_node[t_i]
        g_i = node_to_genome_idx[locus]

        imputed_sample[g_i] = impute_single(
            sample=sample,
            locus=locus, 
            rooted_mst=rooted_mst, 
        )

    for stretch in missing_stretches:
        stretch = [tree_idx_to_node[idx] for idx in stretch]
        ### Impute missing stretches of alleles
        imputed_stretch = impute_stretches(
            missing_nodes_list=stretch,
            sample=sample,
            rooted_mst=rooted_mst,
            n_maxalleles=n_maxalleles+1,
        )

        for node, imputed_allele in imputed_stretch.items():
            g_i = node_to_genome_idx[node]
            imputed_sample[g_i] = imputed_allele
                
    return imputed_sample

def impute_missing(
    rooted_mst: nx.DiGraph,
    samples: np.ndarray,
    missing_idx: np.ndarray,
    n_maxalleles: int,
    node_to_genome_idx: dict,
    tree_idx_to_node: dict,
    genome_idx_to_node: dict,
    n_jobs: int = 1,
) -> np.ndarray:
    """
    Predicting missing samples with a Maximum Spanning Tree. Missing samples are marked as NA.

    Args:
        samples (np.ndarray): The array (2D) of shape (n_samples, n_loci) of the samples to predict.
        missing_idx (np.ndarray): The array (2D) of shape (n_samples, n_missingidx) of the missing alleles
            generated from `get_missing_idx()`.
        rooted_mst (nx.DiGraph): The rooted maximum spanning tree as defined by `build_rooted_mst()`.
        node_to_genome_idx (dict): A dictionary mapping the nodes of the tree to the genome indices.
        tree_idx_to_node (dict): A dictionary mapping the tree indices to the nodes of the tree.
        genome_idx_to_tree_idx (dict): A dictionary mapping the genome indices to the tree indices.
        n_jobs (int): The number of jobs to run in parallel. Defaults to 1.

    Returns:
        np.ndarray: Array (2D) of shape (n_samples, n_loci) with the predicted samples.
    """
    imputed = Parallel(
        ### Do not use `prefer="threads"` as it changes the results!
        n_jobs=n_jobs,
        backend="loky",
    )(
        delayed(_impute_sample)(
            rooted_mst,
            sample,
            idx,
            n_maxalleles,
            node_to_genome_idx,
            tree_idx_to_node,
            genome_idx_to_node,
        )
        for sample, idx in zip(samples, missing_idx)
    )

    return np.asarray(imputed)