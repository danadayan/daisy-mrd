"""
daisy_mrd.lspv.pon
==================
Panel of Normals (PoN) filtering.

Variants that appear recurrently in a panel of normal samples are likely
sequencing artefacts. For each variant position, the PoN stores the
background noise rate ``P_N`` (observed in normal samples). A binomial
test is used to decide whether the observed alt-allele count in the
tumour sample is significantly above this background.

A variant is **retained** if the binomial p-value ≤ 0.05 (i.e. the
observed count is significantly higher than expected noise).

The built-in PoN was built from the same Ultima Genomics platform used
in the original study. Users can supply their own PoN CSV via the
``pon_path`` argument of the pipeline.

PoN CSV format
--------------
Must contain at minimum these columns::

    CHROM   POS   P_N

Where ``P_N`` is the per-position noise rate (float in [0, 1]).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import binomtest


# ---------------------------------------------------------------------------
# PoN loading
# ---------------------------------------------------------------------------

def load_pon(pon_path: str | Path) -> pd.DataFrame:
    """
    Load a Panel of Normals CSV file.

    Adds a composite ``pos`` key (``CHROM_POS``) used for fast merging.

    Parameters
    ----------
    pon_path : str or Path

    Returns
    -------
    pd.DataFrame
        PoN table with columns ``CHROM``, ``POS``, ``P_N``, ``pos``.

    Raises
    ------
    ValueError
        If required columns are missing.
    """
    pon = pd.read_csv(pon_path)
    required = {"CHROM", "POS", "P_N"}
    missing = required - set(pon.columns)
    if missing:
        raise ValueError(f"PoN file is missing required columns: {missing}")

    pon["pos"] = pon["CHROM"].astype(str) + "_" + pon["POS"].astype(str)
    return pon


# ---------------------------------------------------------------------------
# Binomial noise test
# ---------------------------------------------------------------------------

def _binomial_pvalues(df: pd.DataFrame) -> list[float]:
    """
    Compute one-sided (greater) binomial p-value for each variant.

    Tests whether ``mut_reads`` out of ``tot_reads`` is significantly
    greater than the background noise rate ``P_N``.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns ``mut_reads``, ``tot_reads``, ``P_N``.

    Returns
    -------
    list[float]
    """
    pvalues: list[float] = []
    for _, row in df.iterrows():
        result = binomtest(
            k=int(row["mut_reads"]),
            n=int(row["tot_reads"]),
            p=float(row["P_N"]),
            alternative="greater",
        )
        pvalues.append(result.pvalue)
    return pvalues


# ---------------------------------------------------------------------------
# Main filter function
# ---------------------------------------------------------------------------

def filter_pon(
    vcf_df: pd.DataFrame,
    pon: pd.DataFrame,
    pvalue_threshold: float = 0.05,
) -> pd.DataFrame:
    """
    Filter variants using the Panel of Normals.

    Steps:

    1. Add a composite position key and merge with the PoN to get
       per-position noise rates.
    2. Remove rows with zero total reads (cannot be tested).
    3. Run a one-sided binomial test: keep variants whose alt-allele
       count is significantly above the PoN noise rate.

    Parameters
    ----------
    vcf_df : pd.DataFrame
        DataFrame after :func:`~daisy_mrd.lspv.reads.get_reads`.
        Must have ``CHROM``, ``POS``, ``tot_reads``, ``mut_reads``.
    pon : pd.DataFrame
        PoN table loaded by :func:`load_pon`.
    pvalue_threshold : float
        Variants with p-value ≤ this threshold are kept.
        Default: 0.05.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame. The ``pos`` and ``P_N`` helper columns are
        dropped before returning.
    """
    df = vcf_df.copy()
    df["pos"] = df["CHROM"].astype(str) + "_" + df["POS"].astype(str)

    # Merge to get P_N for each position
    merged = pd.merge(df, pon[["pos", "P_N"]], on="pos", how="inner")

    # Remove rows with no coverage
    merged = merged[merged["tot_reads"] > 0].copy()
    merged["tot_reads"] = merged["tot_reads"].astype(int)
    merged["mut_reads"] = merged["mut_reads"].astype(int)

    # Binomial test
    merged["noise_pvalue"] = _binomial_pvalues(merged)

    # Keep only variants that pass
    filtered = merged[merged["noise_pvalue"] <= pvalue_threshold].copy()
    filtered = filtered.drop(columns=["pos", "P_N"], errors="ignore")

    return filtered.reset_index(drop=True)
