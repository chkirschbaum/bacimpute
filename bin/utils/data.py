#!/usr/bin/env python

import numpy as np
import pandas as pd


def _convert_int(
        x: str,
) -> int | str:
    """
    Helper function to convert a string to an integer if possible, otherwise return the original string.

    Args:
        x (str): The input string to convert.

    Returns:
        int or str: The converted integer if possible, otherwise the original string.
    """
    try:
        return int(x)
    
    except (ValueError):
        return x
    

def read_chewie(
    path: str,
    species_mode: str = "clonal", 
    inf_mode: str = "drop",
) -> pd.DataFrame:
    """
    Reads in an output tsv fromn chewBBACA and prepares it for further use.

    Args:
        path (str): Path to the input file from chewBBACA (cgMLST_profiles.tsv, results_alleles.tsv, results_contigsInfo.tsv, ...).
        species_mode (str, optional): Defines if the species is clonal or variable. Defaults to "clonal".
        inf_mode (str, optional): Defines if inferred alleles should be kept or dropped. Defaults to "drop".

    Returns:
        pd.DataFrame: Prepared input dataframe for following tasks.
    """
    chewie_df = pd.read_csv(path, sep='\t')
    chewie_df.rename(columns={"FILE":""}, inplace=True)
    chewie_df.set_index(chewie_df.columns[0], inplace=True)

    if inf_mode == "keep":
        ### results_alleles.tsv can hold INF alleles
        ### We can just include them into the model when they are in the Reference Set
        chewie_df = chewie_df.astype(str).replace(r"INF-", "", regex=True)
        
    if species_mode == "variable":
        ### Calculate percentage of special chewBBACA categories in every locus
        percent_series = chewie_df.astype(str).apply(
            lambda col: ~col.str.contains(r"LNF|PLOT3|PLOT5|PAMA|LOTSC|NIPH|NIPHEM|ASM|ALM", regex=True)
        ).mean() * 100
        ### Get loci which appear in less than 95% of the sequences
        sparse_loci_list = percent_series[percent_series < 95].index.tolist()
        ### Remove the loci from the dataframe
        chewie_df = chewie_df.drop(columns=sparse_loci_list)
        chewie_df = chewie_df.map(_convert_int)

    return chewie_df


def get_pos(
        contig: pd.DataFrame,
) -> tuple[int, int]:
    """
    Get start and end position of the loci in a sequence from chewBBACA results_contigsInfo.tsv.

    Args:
        contig (pd.DataFrame): Dataframe from chewBBACA results_contigsInfo.tsv prepared with `read_chewie()`.

    Returns:
        tuple[int, int]: Start and end position of the loci in the sequence.
    """
    pos = contig.split("&")[1]
    start, end = pos.split("-")

    return int(start), int(end)


def order_to_ref(
        profiles_df: pd.DataFrame,
        ref_id: str,
        contigs_path: str,
) -> pd.DataFrame:
    """
    Orders the loci of the given profiles dataframe to the order of a reference profile.

    Args:
        profiles_df (pd.DataFrame): Dataframe of the sample profiles to order of shape (n_samples, n_loci).
        ref_id (str): The ID of the reference profile in the chewBBACA results_contigsInfo.tsv.
        contigs_path (str): Path to the chewBBACA results_contigsInfo.tsv.

    Returns:
        pd.DataFrame: The ordered dataframe of the sample profiles of shape (n_samples, n_loci).
    """
    contigs = read_chewie(contigs_path)
    contigs = contigs[contigs.columns.intersection(profiles_df.columns)]

    ref = contigs.loc[ref_id]

    loci_order = sorted(ref.index, key=lambda loci: get_pos(ref[loci]))

    profiles_df_ordered = profiles_df[loci_order]    

    return profiles_df_ordered


def recategorize_chewie(
        profiles: pd.DataFrame,
        cg: int = 100,
) -> pd.DataFrame:
    """
    Recategorizes the special chewBBACA categories to new alleles (0) for INF, PLOT3 & PLOT5, PAMA, LOTSC, NIPH & NIPHEM, ASM & ALM and 
    missing alleles (NA) for LNF for further processing.

    In a special case where a model is built for a variable species like Neisseria gonorrhoeae and we use 95% core and the chewBBACA results_alleles.tsv,
    INF alleles are kept as they are in the reference set and therefore can be used for imputation.

    Args:
        profiles (pd.DataFrame): Dataframe of shape (n_samples, n_loci) of the sample profiles to recategorize.
        cg (int, optional): Defining if we are using 100% or 95% core genome. Defaults to 100.

    Returns:
        pd.DataFrame: The recategorized dataframe of shape (n_samples, n_loci) of the sample profiles.
    """
    recategorize = profiles.copy()

    if cg == 100:
        # Locus not found means missing
        recategorize = recategorize.replace("LNF", pd.NA)

        # Inferred alleles are new alleles
        # Set everything with INF-* and *  in front to 0
        recategorize = recategorize.replace(r".*\*.*", 0, regex=True)

        # Define special cases as new alleles
        chewie = ["PLOT3", "PLOT5", "PAMA", "LOTSC", "NIPH", "NIPHEM", "ASM", "ALM"]
        recategorize = recategorize.replace(chewie, 0)
    else:
        # Define special cases as new alleles
        chewie = ["LNF", "PLOT3", "PLOT5", "PAMA", "LOTSC", "NIPH", "NIPHEM", "ASM", "ALM"]
        recategorize = recategorize.replace(chewie, 0)

    # Int columns instead of object
    non_num_cols = recategorize.select_dtypes(include=["object"]).columns
    recategorize[non_num_cols] = recategorize[non_num_cols].apply(
        lambda col: pd.to_numeric(col, errors="coerce")
    ).astype("Int64") 

    return recategorize


def freq_rank_alleles(
        profiles: pd.DataFrame,
        ref_rank: dict = {},
        mode: str = "ref",
        cg: int = 100,
) -> tuple[pd.DataFrame, dict]:
    """
    Mask the cgMLST profiles by frequency of the appearing IDs in the reference set and therefore also reducing the number of observations.
    If you apply a mask to a new dataset in mode `run`, new alleles (not in the reference set used to build the model) are set to 0.

    Args:
        refs (pd.DataFrame): Dataframe of the reference profiles of shape (n_samples, n_loci).
        mode (str): Choose between "ref" if you prepare a reference set to build a model and "run" if you prepare a dataset for imputation. 
            Defaults to "ref".
        ref_rank (dict, optional): Only needed for mode "run". The dict with the frequency mask build on the reference set of the used model. 
            Defaults to {}.

    Returns:
        tuple[pd.DataFrame, dict]: Tuple containing the masked dataframe and the frequency ranks dict.
    """
    freq_rank = profiles.copy()
    rank = {}

    if mode == "run" and ref_rank == {}:
        raise ValueError(
            "The parameter ref_rank can not be an empty dict in mode `run`. Please provide the mask of the model you use."
        )

    elif mode == "ref":
        freq_rank = recategorize_chewie(freq_rank, cg=cg)

        for col in profiles.columns:
            freq = (freq_rank[col].value_counts().rank(method='first',ascending=False).astype("Int64"))

            rank_temp = freq.to_dict()

            # Freq adjust col of reference set
            freq_rank[col] = freq_rank[col].map(rank_temp)
            # Add rank for col to dict
            rank[col] = rank_temp

        return freq_rank, rank
    
    elif mode == "run":
        freq_rank = recategorize_chewie(freq_rank)

        for col in profiles.columns:
            ranked = freq_rank[col].map(ref_rank[col])
            freq_rank[col] = ranked.where(freq_rank[col].notna(), pd.NA).fillna(0).astype("Int64")
        
        return freq_rank, ref_rank
    
    else:
        raise ValueError(
            f"Invalid mode: {mode}. Please choose between `ref` and `run`."
        )

def mask_per_sample(
        profiles: pd.DataFrame,
        percentage: float = 0.05,
) -> pd.DataFrame:
    """
    Masks a dataframe of test profiles per row to check the performance of the model or approaches chosen for comparison.
    Missing alleles are defined by NA.

    Args:
        profiles (pd.DataFrame): The dataframe of sample profiles to mask of shape (n_samples, n_loci).
        percentage (float, optional): The percentage of every sample profile that should be masked. Defaults to 0.05.

    Returns:
        pd.DataFrame: The masked dataframe for the given percentage of samples in every row of shape (n_samples, n_loci).
    """
    masked = profiles.astype("Int64").copy()
    n_cols = masked.shape[1]
    n_mask = int(np.floor(n_cols * percentage))

    for i in range(masked.shape[0]):
        mask_indices = np.random.choice(n_cols, size=n_mask, replace=False)
        masked.iloc[i, mask_indices] = pd.NA
            
    return masked


def mask_per_locus(
        profiles: pd.DataFrame,
        percentage: float = 0.05,
) -> pd.DataFrame:
    """
    Masks a dataframe of test profiles per column to check the performance of the model or approaches chosen for comparison.
    Missing alleles are defined by NA.

    Args:
        profiles (pd.DataFrame): The dataframe of sample profiles to mask of shape (n_samples, n_loci).
        percentage (float, optional): The percentage of every sample profile that should be masked. Defaults to 0.05.

    Returns:
        pd.DataFrame: The masked dataframe for the given percentage of samples in every row of shape (n_samples, n_loci).
    """
    masked = profiles.astype("Int64").copy()
    n_rows = masked.shape[0]
    n_mask = int(np.floor(n_rows * percentage))

    for j in range(masked.shape[1]):
        mask_indices = np.random.choice(n_rows, size=n_mask, replace=False)
        masked.iloc[mask_indices, j] = pd.NA
            
    return masked


def get_missing_idx(
        profiles: pd.DataFrame,
) -> np.ndarray:
    """ 
    Defines an array with the positions of the missing alleles from the given sample profiles.
    Missing alleles are defined by NA.
    
    Args:
        profiles (pd.DataFrame): Dataframe of shape (n_samples, n_loci) of sample profiles with missing alleles.

    Returns:
        np.ndarray: Array (2D) of shape (n_samples, n_missingidx) with the positions of the missing alleles for every sample.
    """
    na_mask = profiles.isna().to_numpy()
    missing_idx = [list(np.flatnonzero(sample)) for sample in na_mask]

    return missing_idx