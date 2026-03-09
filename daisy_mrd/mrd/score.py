"""
daisy_mrd.mrd.score
===================
Compute the DAISY-MRD score and compare it against background noise.

The DAISY-MRD score for a patient is:

    score = Σ(read_num_match) / Σ(Read_Depth)

where the sums run over all LSPV positions that pass the noise filters.

Background noise estimation
---------------------------
The score is evaluated against two noise distributions:

* **Local noise** — the same patient's remission reads piled up against
  *other* patients' LSPV positions.  This tests how much signal the
  patient's remission sample shows at positions that are *not* their own
  leukemia markers.
* **Global noise** — all patients' self-comparison scores evaluated at
  *this* patient's LSPV positions. This tests how noisy these particular
  positions are across the cohort.

Both distributions are illustrated in
:func:`plot_noise_distributions`.

Single-patient vs multi-patient
--------------------------------
For a single patient (no cross-comparisons):

    >>> result = compute_mrd_score(filtered_df)
    >>> print(result.score)

For a cohort (cross-comparison noise estimation):

    >>> scores = compute_cohort_scores(patient_filtered_dfs)
    >>> results_df, figures = compute_cohort_scores(patient_filtered_dfs)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Single-patient score
# ---------------------------------------------------------------------------

@dataclass
class MrdScore:
    """
    DAISY-MRD score for a single patient at a single time-point.

    Attributes
    ----------
    patient_id : str
    total_read_depth : int
        Sum of ``Read_Depth`` across all LSPV positions after filtering.
    total_alt_reads : int
        Sum of ``read_num_match`` across all LSPV positions.
    score : float
        DAISY-MRD score = ``total_alt_reads / total_read_depth``.
        ``nan`` if ``total_read_depth`` == 0.
    n_lspv_positions : int
        Number of LSPV positions that contributed to the score.
    """

    patient_id: str
    total_read_depth: int
    total_alt_reads: int
    score: float
    n_lspv_positions: int


def compute_mrd_score(
    filtered_df: pd.DataFrame,
    patient_id: str = "",
    use_noise_pvalue_filter: bool = True,
) -> MrdScore:
    """
    Compute the DAISY-MRD score from a filtered pileup DataFrame.

    If ``use_noise_pvalue_filter`` is ``True`` and a ``noise_pvalue``
    column is present, the numerator (alt reads) is restricted to
    positions with ``noise_pvalue < 0.05``, while the denominator
    (total depth) uses all positions. This matches the behaviour of the
    ``no_xy_no_germline_no_pon`` filter layer.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        Output of one of the :mod:`~daisy_mrd.mrd.filters` layers.
        Must have ``Read_Depth`` and ``read_num_match``.
    patient_id : str
    use_noise_pvalue_filter : bool
        Default: ``True``.

    Returns
    -------
    MrdScore
    """
    df = filtered_df.copy()
    total_depth = int(df["Read_Depth"].sum())

    if use_noise_pvalue_filter and "noise_pvalue" in df.columns:
        numerator_df = df[df["noise_pvalue"] < 0.05]
    else:
        numerator_df = df

    total_alt = int(numerator_df["read_num_match"].sum())
    score = total_alt / total_depth if total_depth > 0 else float("nan")

    return MrdScore(
        patient_id=patient_id,
        total_read_depth=total_depth,
        total_alt_reads=total_alt,
        score=score,
        n_lspv_positions=len(df),
    )


# ---------------------------------------------------------------------------
# Cohort cross-comparison
# ---------------------------------------------------------------------------

def compute_cohort_scores(
    patient_filtered_dfs: dict[str, pd.DataFrame],
    use_noise_pvalue_filter: bool = True,
) -> pd.DataFrame:
    """
    Compute a cross-comparison MRD score matrix for a cohort of patients.

    For each patient ``P`` and each other patient ``Q``, the score is:

        score(P, Q) = Σ(read_num_match of P's reads at Q's LSPV positions)
                      / Σ(Read_Depth of P's reads at Q's LSPV positions)

    The **diagonal** ``score(P, P)`` is the true DAISY-MRD score.
    Off-diagonal values form the noise background.

    Parameters
    ----------
    patient_filtered_dfs : dict[str, pd.DataFrame]
        Mapping of ``patient_id`` → filtered pileup DataFrame.
        Each DataFrame must contain the columns from
        :func:`~daisy_mrd.mrd.readcount.apply_read_counts` and must
        have been produced by piling up **that patient's remission BAM**
        against **that patient's LSPVs**.

        .. note::
            For full cross-comparison noise estimation, you need pileup
            data for every (remission sample, LSPV set) combination.
            This function handles the simpler case where each DataFrame
            already represents one (remission, LSPV) pairing.

    use_noise_pvalue_filter : bool

    Returns
    -------
    pd.DataFrame
        Matrix with ``Patient`` column and one ``Compared_to_{id}``
        column per patient. Diagonal = true MRD scores.
    """
    patient_ids = list(patient_filtered_dfs.keys())
    columns = ["Patient"] + [f"Compared_to_{p}" for p in patient_ids]
    results: list[dict] = []

    for pid, df in patient_filtered_dfs.items():
        row: dict = {"Patient": pid}

        # Self-comparison: true MRD score
        mrd = compute_mrd_score(df, patient_id=pid,
                                use_noise_pvalue_filter=use_noise_pvalue_filter)
        row[f"Compared_to_{pid}"] = mrd.score

        # Cross-comparisons are only possible if we have per-(remission, LSPV)
        # pileup files for all pairs. In the single-pair case we fill with NaN.
        for other_id in patient_ids:
            if other_id != pid:
                row.setdefault(f"Compared_to_{other_id}", float("nan"))

        results.append(row)

    return pd.DataFrame(results, columns=columns)


def compute_scores_from_pileup_matrix(
    pileup_matrix: dict[str, dict[str, pd.DataFrame]],
    use_noise_pvalue_filter: bool = True,
) -> pd.DataFrame:
    """
    Compute a full cross-comparison MRD score matrix from a nested
    dict of pileup DataFrames.

    This function handles the **full** cross-comparison case where you
    have piled up every remission BAM against every patient's LSPV set.

    Parameters
    ----------
    pileup_matrix : dict[str, dict[str, pd.DataFrame]]
        ``pileup_matrix[remission_id][lspv_patient_id]`` = filtered
        pileup DataFrame for that (remission, LSPV) pair.
    use_noise_pvalue_filter : bool

    Returns
    -------
    pd.DataFrame
        Rows = remission samples, columns = LSPV patient sets.
    """
    remission_ids = list(pileup_matrix.keys())
    lspv_ids = list(next(iter(pileup_matrix.values())).keys())

    columns = ["Patient"] + [f"Compared_to_{lid}" for lid in lspv_ids]
    results: list[dict] = []

    for rid in remission_ids:
        row: dict = {"Patient": rid}
        for lid in lspv_ids:
            df = pileup_matrix[rid].get(lid)
            if df is None or df.empty:
                row[f"Compared_to_{lid}"] = float("nan")
                continue

            total_depth = int(df["Read_Depth"].sum())
            if use_noise_pvalue_filter and "noise_pvalue" in df.columns:
                numerator_df = df[df["noise_pvalue"] < 0.05]
            else:
                numerator_df = df
            total_alt = int(numerator_df["read_num_match"].sum())
            row[f"Compared_to_{lid}"] = (
                total_alt / total_depth if total_depth > 0 else float("nan")
            )
        results.append(row)

    return pd.DataFrame(results, columns=columns)


# ---------------------------------------------------------------------------
# Noise distribution plots
# ---------------------------------------------------------------------------

def plot_noise_distributions(
    scores_df: pd.DataFrame,
    patient_id: str,
    output_dir: str | Path | None = None,
) -> matplotlib.figure.Figure:
    """
    Plot the MRD score for one patient against local and global noise
    distributions.

    * **Local noise** (plum histogram): scores of this patient's remission
      reads against all *other* patients' LSPV sets.
    * **Global noise** (orange histogram): scores of all patients' remission
      reads against *this* patient's LSPV set.
    * **MRD value** (blue vertical line): the self-comparison score.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Cross-comparison matrix from :func:`compute_cohort_scores` or
        :func:`compute_scores_from_pileup_matrix`.
    patient_id : str
        The patient whose MRD value to highlight.
    output_dir : str, Path, or None
        If provided, the figure is saved as a PDF in this directory.

    Returns
    -------
    matplotlib.figure.Figure
    """
    self_col = f"Compared_to_{patient_id}"

    if patient_id not in scores_df["Patient"].values:
        raise ValueError(f"Patient '{patient_id}' not found in scores_df['Patient'].")
    if self_col not in scores_df.columns:
        raise ValueError(f"Column '{self_col}' not found in scores_df.")

    patient_row = scores_df[scores_df["Patient"] == patient_id].iloc[0]
    self_score = patient_row[self_col]

    # Local noise: this patient's row, excluding self-comparison
    local_noise = patient_row.drop(["Patient", self_col]).values.astype(float)
    local_noise = local_noise[~np.isnan(local_noise)]
    local_noise = local_noise[local_noise != self_score]

    # Global noise: all rows' score in self_col, excluding self
    global_noise = scores_df.loc[
        scores_df["Patient"] != patient_id, self_col
    ].values.astype(float)
    global_noise = global_noise[~np.isnan(global_noise)]

    n_bins = max(10, min(30, int(np.sqrt(len(local_noise) + 1))))

    fig, ax = plt.subplots(figsize=(7, 5))

    if len(local_noise) > 0:
        ax.hist(local_noise, bins=n_bins, alpha=0.5, color="plum",
                label="Local noise", density=True)
    if len(global_noise) > 0:
        ax.hist(global_noise, bins=n_bins, alpha=0.5, color="orange",
                label="Global noise", density=True)

    ax.axvline(
        x=self_score,
        color="cornflowerblue",
        linewidth=3,
        label=f"MRD score ({self_score:.4f})",
    )

    ax.set_title(f"Patient {patient_id} — DAISY-MRD score", fontsize=13)
    ax.set_xlabel("VAF (DAISY-MRD score)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.legend(bbox_to_anchor=(1, 1.01), loc="upper left", fontsize=10)
    ax.grid(False)
    fig.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            output_dir / f"mrd_noise_{patient_id}.pdf",
            format="pdf",
            dpi=300,
            bbox_inches="tight",
        )

    return fig
