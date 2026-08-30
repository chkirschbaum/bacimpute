#!/usr/bin/env python

import numpy as np
import pandas as pd

from sklearn.metrics import normalized_mutual_info_score

import warnings


def mutual_info_matrix(
        freq_profile_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate the normalized mutual information between all columns of a dataframe.

    Args:
        freq_profile_df (pd.DataFrame): With `freq_rank_alleles` prepared dataframe of cgMLST allele profile from chewBBACA allele calling

    Returns:
        pd.DataFrame: Mutual information between all columns
    """
    allele_array = freq_profile_df.to_numpy()

    ### Initalize MI matrix
    col_names = freq_profile_df.columns
    col_len = len(col_names)
    mi_matrix = np.zeros((col_len, col_len), dtype=np.float32)

    ### Loop through cols
    for i in range(col_len):
        for j in range(i, col_len):
            true = allele_array[:, i]
            pred = allele_array[:, j]
            ### Mask to consider only complete pairs like correlation function
            mask = ~(np.isnan(true) | np.isnan(pred))
            mi = normalized_mutual_info_score(true[mask], pred[mask])
            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi

    mi_matrix_df = pd.DataFrame(mi_matrix, index=col_names, columns=col_names)

    return mi_matrix_df
