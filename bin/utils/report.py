#!/usr/bin/env python

import numpy as np
import pandas as pd

import plotly.graph_objects as go


def get_acc(
        original: np.ndarray,
        imputed: np.ndarray,
        missing: np.ndarray,
        missing_per_col: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """
    Calculates overall accuracy and accuracy per locus.

    Args:
        original (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the original profiles.
        imputed (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the imputed profiles.
        masked (np.ndarray): Array (2D) of shape (n_samples, n_loci) where there were a certain percentage of 
            every original profiles was masked as missing.

    Returns:
        tuple[float, np.ndarray, np.ndarray]: The overall accuracy, the two arrays (1D) of shape (n_loci) with the 
            accuracy per column and the count of correct imputed alleles per column.
    """
    ### Overall accuracy
    correct_count = np.sum(original[missing] == imputed[missing])
    total_count = np.sum(missing)
    acc = (correct_count / total_count) * 100

    ### Accuracy per loci 
    correct_imputations = (imputed == original) & missing
    correct_per_col = correct_imputations.sum(axis=0)
    acc_per_col = np.where(missing_per_col > 0, (100 * correct_per_col / missing_per_col), 0)

    return acc, acc_per_col, correct_per_col


def accuracy(
        original: np.ndarray,
        imputed_baseline: np.ndarray, 
        imputed_model: np.ndarray,
        masked: np.ndarray,
) -> tuple[float, np.ndarray, float, np.ndarray]:
    """
    Compares the results of an imputed profile to the original profile array of the Baseline and the model.
    Calculates overall accuracy, plots the accuracy per locus and a histogram with the accuracy fractions.

    Args:
        original (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the original profiles.
        imputed_baseline (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the imputed profiles from the baseline.
        imputed_model (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the imputed profiles from the model.
        masked (np.ndarray): Array (2D) of shape (n_samples, n_loci) where there were a certain percentage of 
            every original profiles was masked as missing.

    Returns:
        tuple[np.ndarray, np.ndarray]: Two arrays (1D) of shape (n_loci) with tyhe accuracy per locus for the baseline and the model. 
    """
    missing = pd.isna(masked)
    missing_per_col = missing.sum(axis=0)

    acc_baseline, acc_per_col_baseline, _ = get_acc(original, imputed_baseline, missing, missing_per_col)

    acc_model, acc_per_col_model, _ = get_acc(original, imputed_model, missing, missing_per_col)

    return acc_baseline, acc_per_col_baseline, acc_model, acc_per_col_model


def plot_summary(
        n_loci: int,
        imputed_baseline_all: np.ndarray, 
        imputed_model_all: np.ndarray,
        model: str,
        legend_width: int = 80,
        title: str = "",
) -> None:
    """
    Args:
        n_loci: Number of loci per sample.
        imputed_baseline (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the imputed profiles from the baseline.
        imputed_model (np.ndarray): Array (2D) of shape (n_samples, n_loci) with the imputed profiles from the Markov chain.

    Returns:
        None
    """
    fig_vio = go.Figure()

    fig_vio.add_trace(
        go.Violin(
            y=imputed_model_all,
            x=(["5% Masked"]*n_loci)+(["10% Masked"]*n_loci)+(["15% Masked"]*n_loci+["20% Masked"]*n_loci)+(["25% Masked"]*n_loci),
            marker_color="#5BBCD6",
            name="Imputation Model" + model.upper(),
            scalegroup="Imputation Model",
            legendgroup="Imputation Model",
            showlegend=True,
            side="positive",
            line_color="#5BBCD6",
        )
    )

    fig_vio.add_trace(
        go.Violin(
            y=imputed_baseline_all,
            x=(["5% Masked"]*n_loci)+(["10% Masked"]*n_loci)+(["15% Masked"]*n_loci+["20% Masked"]*n_loci)+(["25% Masked"]*n_loci),
            marker_color="#FF0000",
            name="Baseline",
            scalegroup="Imputation Model",
            legendgroup="Baseline",
            showlegend=True,
            side="negative",
            line_color="#FF0000",
        )
    )

    fig_vio.update_layout(
        boxmode="group",
        title=title,
        yaxis_title="Accuracy per Locus (%)",
        plot_bgcolor='#FFFFFF',
        hovermode="x unified",
        legend=dict(
            orientation="h",
            entrywidth=legend_width,
            yanchor="bottom",
            y=1.05,
            xanchor="right",
            x=1,
            font=dict(size=14),
        ),
        xaxis=dict(title_font=dict(size=16), tickfont=dict(size=14)),
        yaxis=dict(title_font=dict(size=16), tickfont=dict(size=14)),
    )

    return fig_vio