#!/usr/bin/env python

import numpy as np

import networkx as nx

import pickle
import json


def save_imc(
        species: str,
        freq_rank_dict: dict,
        transition_array: np.ndarray,
        max_alleles_array: np.ndarray,
        cg_list: list,
        ref: str, 
        path: str = "model.pkl",
) -> None:
    """
    Saves the IMC model to a file using pickle.

    Args:
        freq_rank_dict (dict): The frequency mask build on the reference set of the used model. 
        transition_array (np.ndarray): The transition matrix (3D) of shape (n_loci, n_maxalleles, n_maxalleles).
        max_alleles_array (np.ndarray): The maximum number of alleles for each locus.
        cg_list (list): The ordered loci of the core genome.
        ref (str): Name of the reference sequence used.
        path (str): The path to the file to save the model to. Defaults to "model.pkl".
        
    Returns:
        None
    """
    model = {
        "species": species,
        "freq_rank_dict": freq_rank_dict,
        "transition_array": transition_array,
        "max_alleles_array": max_alleles_array,
        "cg_list": cg_list,
        "ref": ref
    }

    with open(path, "wb") as f:
        pickle.dump(model, f)

    print(f"The IMC model was saved successfully in {path}")


def load_imc(
        path: str = "",
) -> tuple[dict, np.ndarray, np.ndarray, list, str]:
    """
    Loads the IMC model from a file using pickle.

    Args:
        path (str): The path to the file to load the model from.

    Returns:
        tuple[dict, np.ndarray, np.ndarray, list, str]: The loaded model components.
    """
    if path == "":
        raise ValueError(
            "Please add the path to your input file to the parameter `path`."
        )

    with open(path, "rb") as f:
        model = pickle.load(f)

    species = model["species"]
    freq_rank_dict = model["freq_rank_dict"]
    transition_array = model["transition_array"]
    max_alleles_array = model["max_alleles_array"]
    cg_list = model["cg_list"]
    ref = model["ref"]

    return species, freq_rank_dict, transition_array, max_alleles_array, cg_list, ref


def save_maxst(
        species: str,
        freq_rank_dict: dict,
        rooted_mst: nx.DiGraph,
        n_maxalleles: int,
        node_to_genome_idx: dict,
        tree_idx_to_node: dict,
        genome_idx_to_node: dict,
        cg_list: list,
        ref: str,
        path: str = "model.pkl",
) -> None:
    """
    Saves the MaxST model to a file using pickle.

    Args:
        freq_rank_dict (dict): The frequency mask build on the reference set of the used model.
        rooted_mst (nx.DiGraph): The rooted minimum spanning tree.
        node_to_genome_idx (dict): Mapping from node to genome index.
        tree_idx_to_node (dict): Mapping from tree index to node.
        cg_list (list): The ordered loci of the core genome.
        ref (str): Name of the reference sequence used.
        path (str): The path to the file to save the model to. Defaults to "model.pkl".

    Returns:
        None
    """
    model = {
        "species": species,
        "freq_rank_dict": freq_rank_dict,
        "rooted_mst": rooted_mst,
        "n_maxalleles": n_maxalleles,
        "node_to_genome_idx": node_to_genome_idx,
        "tree_idx_to_node": tree_idx_to_node,
        "genome_idx_to_node": genome_idx_to_node,
        "cg_list": cg_list,
        "ref": ref
    }

    with open(path, "wb") as f:
        pickle.dump(model, f)

    print(f"The MaxST model was saved successfully in {path}")

    return None


def load_maxst(
        path: str = "",
) -> tuple[dict, nx.DiGraph, dict, dict, dict, list, str]:
    """
    Loads the MaxST model from a file using pickle.
    
    Args:
        path (str): The path to the file to load the model from.

    Returns:
        tuple[dict, nx.DiGraph, dict, dict, dict, list, str]: The loaded model components.
    """
    if path == "":
        raise ValueError(
            "Please add the path to your input file to the parameter `path`."
        )

    with open(path, "rb") as f:
        model = pickle.load(f)

    species = model["species"]
    freq_rank_dict = model["freq_rank_dict"]
    rooted_mst = model["rooted_mst"]
    n_maxalleles = model["n_maxalleles"]
    node_to_genome_idx = model["node_to_genome_idx"]
    tree_idx_to_node = model["tree_idx_to_node"]
    genome_idx_to_node = model["genome_idx_to_node"]
    cg_list = model["cg_list"]
    ref = model["ref"]

    return species, freq_rank_dict, rooted_mst, n_maxalleles, node_to_genome_idx, tree_idx_to_node, genome_idx_to_node, cg_list, ref


def save_cytoscape(
        G: nx.Graph,
        path: str = "maxst.json",
) -> None:
    """
    Save the networkx graph in JSON file for cytoscape.

    Args:
        G (nx.Graph): The nextworkx graph
        path (str): The path for the JSON file

    Returns:
        None
    """
    cyto = nx.cytoscape_data(G)

    with open(path, 'w') as f:
        json.dump(cyto, f, indent=2)

    print(f"The MaxST Tree was saved successfully in {path} for visualization in Cytoscape.")

    return None