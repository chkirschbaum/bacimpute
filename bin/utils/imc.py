#!/usr/bin/env python

import numpy as np
import pandas as pd

from joblib import Parallel, delayed


def build_markov_chain(
        freq_rank_profiles: pd.DataFrame,
        epsilon: float = 0.02,
) -> tuple[np.ndarray, np.ndarray]:
    """ 
    Building a start probability and transition matrix from locus to locus based on the given reference profile.  

    Args:
        freq_rank_profiles (pd.Dataframe): Dataframe of the reference profiles of shape (n_samples, n_loci).
        epsilon (float): Pseudo count value to avoid `0.` values in the transitions.

    Returns:
        tuple[np.ndarray, np.ndarray]: The tuple contains the arrays for the transition matrix (3D) 
            of shape (n_loci - 1, n_maxalleles, n_maxalleles) and max alleles per locus (1D) of shape (n_loci).
    """
    states = freq_rank_profiles.columns 

    ### Get maximum number of alleles per locus
    max_alleles_array = freq_rank_profiles.max()
    max_alleles = max(max_alleles_array)
    
    if max_alleles < 255:
        freq_rank_profiles = freq_rank_profiles.astype(np.uint8)
    elif max_alleles < 65535:
        freq_rank_profiles = freq_rank_profiles.astype(np.uint16)

    ### For every state, every allele has an own transition array
    transition = np.zeros(
        (len(states) - 1, max_alleles+1, max_alleles+1), 
        dtype=np.float32
    )

    for locus in range(len(states) - 1):
        col1 = freq_rank_profiles[states[locus]]
        col1_max = int(max_alleles_array.iloc[locus])
        col2 = freq_rank_profiles[states[locus + 1]]
        col2_max = int(max_alleles_array.iloc[locus + 1])

        ### Get allele ID combinations for two consecutive columns
        ### col1 rows, col2 columns
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

        transition[locus] = transition_temp.T.to_numpy(dtype=np.float32)

    return transition, max_alleles_array


def impute_single(
        sample: np.ndarray,
        locus: int, 
        transition: np.ndarray,
) -> int:
    """
    Imputes a single missing allele.
    
    Args:
        sample (np.ndarray): The array (1D) of shape (n_loci) with the missing allele.
        locus (int): Position of the locus with the missing allele.
        transition (np.ndarray): The transition matrix (3D) of shape (n_loci, n_maxalleles, n_maxalleles)
            as defined in `build_markov_chain()`.

    Returns:
        int: The predicted allele at the given locus.
    """
    ### Missing at first locus
    if locus == 0:
        return int(
            np.nanargmax(transition[locus][:, sample[locus + 1]])
        )

    ### Missing at last locus
    elif locus == (len(sample) - 1):
        return int(
            np.nanargmax(transition[locus - 1][sample[locus - 1], :])
        )

    ### Searching MAP for missing x2
    else:
        ### i is observed state at x1, j at x3
        i = sample[locus - 1]
        j = sample[locus + 1]

        ### p12 is transitions from x1 to x2, p23 from x2 to x3
        p12 = transition[locus - 1][i, :]
        p23 = transition[locus][:, j]

        posterior = p12 * p23

        if np.nansum(posterior) == 0:
            return int(np.nanargmax(p12))

        return int(np.nanargmax(posterior))
                

def impute_stretches(
        missing_stretch_idx: list,
        sample: np.ndarray,
        transition: np.ndarray,
) -> list:
    """
    Imputes a missing stretch of alleles.
    
    Args:
        missing_stretch_idx (list): List of the indices of the missing alleles on the stretch of missing alleles.
        sample (np.ndarray): The array (1D) of shape (n_loci) with the missing allele.
        transition (np.ndarray): The transition matrix (3D) of shape (n_loci, n_maxalleles, n_maxalleles)
            as defined in `build_markov_chain()`.

    Returns:
        list: List with the imputed missing alleles of the stretch.
    """
    n_maxalleles = transition.shape[1]

    start = missing_stretch_idx[0]
    end = missing_stretch_idx[-1]

    n_missing = end - start + 1

    if 0 in missing_stretch_idx:
        ### Missing stretch is at the start of the sample
        missing_stretch_idx = np.append(missing_stretch_idx, end + 1)
        j = sample[end + 1]

        backward = np.empty((n_missing, n_maxalleles), dtype=float)

        ### Only Backward pass
        backward[-1] = transition[end][:, j]

        for k in range(n_missing - 2, -1, -1):
            locus = start + k
            backward[k] = transition[locus] @ backward[k + 1]

        return np.nanargmax(backward, axis=1)

    elif (len(sample) - 1) in missing_stretch_idx:
        ### Missing stretch is at the end of the sample
        missing_stretch_idx = np.insert(missing_stretch_idx, 0, start - 1)
        i = sample[start - 1]

        forward = np.empty((n_missing, n_maxalleles), dtype=float)

        ### Only Forward pass
        forward[0] = transition[start - 1][i, :]

        for k in range(1, n_missing):
            locus = start + k - 1
            forward[k] = forward[k - 1] @ transition[locus]

        return np.nanargmax(forward, axis=1)

    else:
        missing_stretch_idx = np.insert(missing_stretch_idx, 0, start - 1)
        missing_stretch_idx = np.append(missing_stretch_idx, end + 1)

        i = sample[start - 1]
        j = sample[end + 1]

        forward = np.empty((n_missing, n_maxalleles), dtype=float)
        backward = np.empty((n_missing, n_maxalleles), dtype=float)

        ### Forward pass
        forward[0] = transition[start - 1][i, :]

        for k in range(1, n_missing):
            locus = start + k - 1
            forward[k] = forward[k - 1] @ transition[locus]

        ### Backward pass
        backward[-1] = transition[end][:, j]

        for k in range(n_missing - 2, -1, -1):
            locus = start + k
            backward[k] = transition[locus] @ backward[k + 1]

        # Get MAP estimate for each missing locus
        posterior = forward * backward

        return np.nanargmax(posterior, axis=1)


def _impute_sample(
    sample: np.ndarray,
    missing_idx: np.ndarray,
    transition: np.ndarray,
) -> np.ndarray:
    """
    Imputes a single sample with missing alleles.

    Args:
        sample (np.ndarray): The array (1D) of shape (n_loci,) with the missing alleles.
        missing_idx (np.ndarray): The array (1D) of shape (n_missingidx,) of the missing alleles generated 
            from `get_missing_idx()`.
        transition (np.ndarray): The transition matrix (3D) of shape (n_loci, n_maxalleles, n_maxalleles) 
            as defined in `build_markov_chain()`.

    Returns:
        np.ndarray: Array (1D) of shape (n_loci,) with the predicted sample.
    """
    imputed_sample = sample.copy()

    ### Sample with nothing to impute
    if len(missing_idx) == 0:
        return imputed_sample

    ### Split for missing single alleles and missing stretches of alleles
    breaks = np.where(np.diff(missing_idx) != 1)[0] + 1

    for loci in np.split(missing_idx, breaks):
        if len(loci) == 1:
            ### Impute missing single alleles
            locus = loci[0]
            imputed_sample[locus] = impute_single(
                sample=sample,
                locus=locus,
                transition=transition,
            )

        else:
            ### Impute missing stretches of alleles
            imputed_sample[loci] = impute_stretches(
                missing_stretch_idx=loci,
                sample=sample,
                transition=transition,
            )

    return imputed_sample


def impute_missing(
    transition: np.ndarray,
    samples: np.ndarray,
    missing_idx: np.ndarray,
    n_jobs: int = 1,
) -> np.ndarray:
    """
    Imputes missing alleles in the given samples using a Markov chain. Missing alleles are marked as NA.

    Args:
        transition (np.ndarray): The transition matrix (3D) of shape (n_loci, n_maxalleles, n_maxalleles) as defined in `build_markov_chain()`.
        samples (np.ndarray): The array (2D) of shape (n_samples, n_loci) of samples to predict.
        missing_idx (np.ndarray): The array (2D) of shape (n_samples, n_missingidx) of missing allele positions generated from `get_missing_idx()`.
        n_jobs (int): The number of jobs to run in parallel. Default is 1.

    Returns:
        np.ndarray: Array (2D) of shape (n_samples, n_loci) with the predicted samples.
    """
    imputed = Parallel(
        ### Do not use `prefer="threads"` as it changes the results!
        n_jobs=n_jobs,
        backend="loky",
    )(
        delayed(_impute_sample)(
            sample,
            idx,
            transition,
        )
        for sample, idx in zip(samples, missing_idx)
    )

    return np.asarray(imputed)