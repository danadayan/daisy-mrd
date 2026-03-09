"""
daisy_mrd.lspv.identify
=======================
Final LSPV extraction and summary visualisation.

After clonality labelling, LSPVs are defined as variants that are:

1. **Clonal** — present in the entire leukemic cell population
   (``clonality == "clonal"``).
2. **Non-coding** — no amino-acid-level annotation in ``hgvsp_short``
   (coding driver mutations are excluded so that only passenger
   variants remain).

A pie chart summarising clonal vs sub-clonal variant counts is also
produced here.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.figure
import pandas as pd


# ---------------------------------------------------------------------------
# LSPV extraction
# ---------------------------------------------------------------------------

def extract_lspvs(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract leukemia-specific passenger variants (LSPVs).

    Keeps only variants that are both clonal and non-coding.

    Parameters
    ----------
    vcf_df : pd.DataFrame
        DataFrame after :func:`~daisy_mrd.lspv.gmm.label_clonality`.
        Must have ``clonality`` and ``hgvsp_short`` columns.

    Returns
    -------
    pd.DataFrame
        LSPV table, reset index.
    """
    # Step 1: clonal only
    df = vcf_df[vcf_df["clonality"] == "clonal"].copy()

    # Step 2: non-coding only (hgvsp_short is NaN / None for non-coding)
    df = df[df["hgvsp_short"].isna()].copy()

    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def lspv_summary(
    vcf_df: pd.DataFrame,
    lspv_df: pd.DataFrame,
    patient_id: str = "",
) -> pd.DataFrame:
    """
    Build a one-row summary DataFrame for a single patient / sample.

    Parameters
    ----------
    vcf_df : pd.DataFrame
        Full filtered variant table (before LSPV extraction).
    lspv_df : pd.DataFrame
        LSPV table from :func:`extract_lspvs`.
    patient_id : str

    Returns
    -------
    pd.DataFrame
        Columns: ``patient``, ``total_variants``, ``clonal``,
        ``subclonal``, ``n_lspvs``.
    """
    clonal_counts = vcf_df["clonality"].value_counts()
    return pd.DataFrame(
        {
            "patient": [patient_id],
            "total_variants": [len(vcf_df)],
            "clonal": [clonal_counts.get("clonal", 0)],
            "subclonal": [clonal_counts.get("subclonal", 0)],
            "n_lspvs": [len(lspv_df)],
        }
    )


# ---------------------------------------------------------------------------
# Pie chart
# ---------------------------------------------------------------------------

def plot_clonality_pie(
    vcf_df: pd.DataFrame,
    patient_id: str = "",
) -> matplotlib.figure.Figure:
    """
    Plot a pie chart of clonal vs sub-clonal variant counts.

    Parameters
    ----------
    vcf_df : pd.DataFrame
        Must have a ``clonality`` column.
    patient_id : str
        Used in the plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    color_map = {
        "clonal": "#2563eb",
        "subclonal": "#7dd3fc",
    }

    counts = vcf_df["clonality"].value_counts()

    # Keep a stable ordering: clonal first, then subclonal
    labels = [k for k in ["clonal", "subclonal"] if k in counts.index]
    sizes = [counts[k] for k in labels]
    colors = [color_map[k] for k in labels]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(
        sizes,
        labels=labels,
        autopct="%1.1f%%",
        colors=colors,
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.2},
    )
    title = "Clonal vs Sub-clonal variants"
    if patient_id:
        title += f"\nPatient {patient_id}"
    ax.set_title(title, fontsize=13)

    fig.tight_layout()
    return fig
