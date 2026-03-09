"""
daisy_mrd.mrd.pileup
====================
Run ``samtools mpileup`` over LSPV positions in a remission BAM/CRAM,
then parse the resulting pileup file into a pandas DataFrame.

Two public functions:

* :func:`run_mpileup`  — calls samtools and writes a ``.pileup`` file
* :func:`pileup_to_df` — parses a ``.pileup`` file into a DataFrame

The pileup format produced by ``samtools mpileup -r chrom:start-end``
emits *two* lines per region when the ``-r`` flag targets a single
position: a header line and a data line. :func:`pileup_to_df` keeps
only the data lines (every second row, index 1, 3, 5 …).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

PILEUP_COLUMNS = ["CHROM", "POS", "REF", "Read_Depth", "Read_Bases", "Base_Qualities"]


# ---------------------------------------------------------------------------
# Position formatting
# ---------------------------------------------------------------------------

def format_positions(lspv_df: pd.DataFrame) -> list[str]:
    """
    Convert an LSPV DataFrame into a list of samtools region strings.

    samtools ``-r`` expects ``chrom:start-end`` where coordinates are
    **1-based** and the region is half-open on the left, so to target
    a single SNV at position ``POS`` we use ``chrom:POS-1-POS``.

    Parameters
    ----------
    lspv_df : pd.DataFrame
        Must contain ``CHROM`` and ``POS`` columns.

    Returns
    -------
    list[str]
        e.g. ``["chr1:99-100", "chr3:2049-2050", …]``
    """
    positions: list[str] = []
    for _, row in lspv_df.iterrows():
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        positions.append(f"{chrom}:{pos - 1}-{pos}")
    return positions


# ---------------------------------------------------------------------------
# samtools mpileup runner
# ---------------------------------------------------------------------------

def run_mpileup(
    cram_or_bam: str | Path,
    reference: str | Path,
    lspv_df: pd.DataFrame,
    output_pileup: str | Path,
    samtools_path: str = "samtools",
    extra_flags: list[str] | None = None,
) -> Path:
    """
    Run ``samtools mpileup`` over all LSPV positions in a remission
    BAM or CRAM file.

    Each LSPV position is queried individually with ``-r chrom:start-end``
    and all output is appended to a single ``.pileup`` file. This mirrors
    the original approach and gives one pileup entry per LSPV position.

    Parameters
    ----------
    cram_or_bam : str or Path
        Path to the remission alignment file (``.cram`` or ``.bam``).
        An index file (``.crai`` / ``.bai``) must be present alongside it.
    reference : str or Path
        Path to the reference FASTA used for CRAM decoding (``-f`` flag).
    lspv_df : pd.DataFrame
        LSPV table produced by Step 1.  Must contain ``CHROM`` and
        ``POS`` columns.
    output_pileup : str or Path
        Path where the output ``.pileup`` file will be written.
    samtools_path : str
        Path to the ``samtools`` binary. Defaults to ``"samtools"``
        (assumes it is on ``$PATH``).
    extra_flags : list[str] or None
        Any additional flags to pass to ``samtools mpileup``, e.g.
        ``["--min-BQ", "20"]``.

    Returns
    -------
    Path
        Path to the written pileup file.

    Raises
    ------
    subprocess.CalledProcessError
        If samtools exits with a non-zero return code for any position.
    FileNotFoundError
        If ``cram_or_bam`` or ``reference`` do not exist.
    """
    cram_or_bam = Path(cram_or_bam)
    reference = Path(reference)
    output_pileup = Path(output_pileup)

    if not cram_or_bam.exists():
        raise FileNotFoundError(f"BAM/CRAM not found: {cram_or_bam}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")

    positions = format_positions(lspv_df)
    output_pileup.parent.mkdir(parents=True, exist_ok=True)

    base_cmd = [samtools_path, "mpileup", "-f", str(reference)]
    if extra_flags:
        base_cmd.extend(extra_flags)

    with open(output_pileup, "w") as out_fh:
        for position in positions:
            cmd = base_cmd + ["-r", position, str(cram_or_bam)]
            subprocess.run(cmd, shell=False, check=True, stdout=out_fh)

    return output_pileup


# ---------------------------------------------------------------------------
# Pileup parser
# ---------------------------------------------------------------------------

def pileup_to_df(
    pileup_file: str | Path,
    output_csv: str | Path | None = None,
) -> pd.DataFrame:
    """
    Parse a samtools pileup file into a pandas DataFrame.

    When samtools mpileup is called with ``-r`` per position, it emits
    two lines per region: a header line and a data line.  This function
    keeps only the **data lines** (every second row starting at index 1).

    Parameters
    ----------
    pileup_file : str or Path
    output_csv : str, Path, or None
        If provided, the DataFrame is saved as CSV to this path.

    Returns
    -------
    pd.DataFrame
        Columns: ``CHROM``, ``POS``, ``REF``, ``Read_Depth``,
        ``Read_Bases``, ``Base_Qualities``.
    """
    pileup_file = Path(pileup_file)

    df = pd.read_csv(
        pileup_file,
        sep="\t",
        header=None,
        names=PILEUP_COLUMNS,
    )

    # Keep only data rows (every second row, 0-indexed: 1, 3, 5, …)
    df = df.iloc[1::2].reset_index(drop=True)

    # Coerce types
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["Read_Depth"] = pd.to_numeric(df["Read_Depth"], errors="coerce").fillna(0).astype(int)

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_csv, index=False)

    return df


# ---------------------------------------------------------------------------
# Merge pileup with LSPV original ALT
# ---------------------------------------------------------------------------

def merge_lspv_alts(
    pileup_df: pd.DataFrame,
    lspv_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge the pileup DataFrame with the LSPV table to attach the
    original ALT allele (``OG_ALT``) to each pileup row.

    This is needed so downstream read counting knows which base to look
    for in the pileup ``Read_Bases`` string.

    Parameters
    ----------
    pileup_df : pd.DataFrame
        Output of :func:`pileup_to_df`.
    lspv_df : pd.DataFrame
        LSPV table (must have ``CHROM``, ``POS``, ``ALT``).

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with ``OG_ALT`` column added.
    """
    alt_df = lspv_df[["CHROM", "POS", "ALT"]].copy()
    alt_df = alt_df.rename(columns={"ALT": "OG_ALT"})
    alt_df["POS"] = pd.to_numeric(alt_df["POS"], errors="coerce")

    merged = pd.merge(pileup_df, alt_df, on=["CHROM", "POS"], how="inner")
    return merged
