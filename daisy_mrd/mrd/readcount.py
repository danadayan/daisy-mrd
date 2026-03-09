"""
daisy_mrd.mrd.readcount
=======================
Count ALT-supporting reads from the samtools pileup ``Read_Bases``
column, then attach per-variant counts to the pileup DataFrame.

The pileup ``Read_Bases`` field uses a compact encoding:

* ``.`` / ``,``  — reference base match (forward / reverse strand)
* ``A C G T`` (upper) — mismatch on forward strand
* ``a c g t`` (lower) — mismatch on reverse strand
* ``^<qual>``    — read start (followed by mapping-quality character)
* ``$``          — read end
* ``+Nseq``      — insertion of *N* bases *seq*
* ``-Nseq``      — deletion of *N* bases *seq*
* ``*``          — deletion placeholder

This module provides:

* :func:`extract_specific_indel` — count one specific indel pattern
* :func:`count_variants`         — count all base-level mismatches for one row
* :func:`match_alts`             — select the right count column per variant
* :func:`apply_read_counts`      — apply the above to an entire DataFrame
"""

from __future__ import annotations

import re
from collections import Counter

import pandas as pd


# ---------------------------------------------------------------------------
# Indel counting
# ---------------------------------------------------------------------------

def extract_specific_indel(bases_str: str, target_indel: str) -> int:
    """
    Count occurrences of a specific indel in a pileup Read_Bases string.

    Parameters
    ----------
    bases_str : str
        Raw ``Read_Bases`` string from samtools pileup.
    target_indel : str
        The ALT allele. Insertions are longer than 1 character (e.g.
        ``"AT"`` means the inserted sequence is ``"T"``). Deletions are
        represented by a ``"-"`` prefix in your ALT field.

    Returns
    -------
    int
        Number of reads supporting this exact indel.
    """
    # In VCF, an insertion ALT like "AT" means REF="A", ALT="AT" → inserted bases = "T"
    # The pileup encodes this as "+1T".
    is_insertion = len(target_indel) > 1 and not target_indel.startswith("-")
    target_seq = target_indel[1:] if is_insertion else target_indel.lstrip("-")

    if is_insertion:
        pattern = rf"\+{len(target_seq)}({re.escape(target_seq)})"
    else:
        pattern = rf"-{len(target_seq)}({re.escape(target_seq)})"

    matches = re.finditer(pattern, bases_str, re.IGNORECASE)
    return sum(1 for _ in matches)


# ---------------------------------------------------------------------------
# Per-row variant counting
# ---------------------------------------------------------------------------

def count_variants(
    bases_str: str,
    ref_base: str,
    og_alt: str,
) -> tuple[int, str, int, int, int, int, int]:
    """
    Count non-reference bases in a pileup Read_Bases string.

    For **indels** (``len(og_alt) > 1``), the function delegates to
    :func:`extract_specific_indel` and returns zero for all SNV counts.

    For **SNVs**, the function cleans the string of pileup encoding
    characters and counts occurrences of each non-reference base.

    Parameters
    ----------
    bases_str : str
        Raw pileup ``Read_Bases`` field.
    ref_base : str
        Reference base (``REF`` column).
    og_alt : str
        Expected ALT allele (``OG_ALT`` column).

    Returns
    -------
    tuple of 7 values:
        ``(total_mutations, mutated_bases_str, A_count, T_count,
          C_count, G_count, indel_count)``
    """
    bases_str = str(bases_str)
    og_alt = str(og_alt).upper()

    # --- Indel ---
    if len(og_alt) > 1:
        indel_count = extract_specific_indel(bases_str, og_alt)
        return (0, "", 0, 0, 0, 0, indel_count)

    # --- SNV: clean pileup encoding ---
    s = bases_str
    s = re.sub(r"-\d+[a-zA-Z]+", "", s)   # remove deletions
    s = re.sub(r"\*\w*", "", s)            # remove deletion placeholders
    s = re.sub(r"\+\d+[A-Za-z]+", "", s)  # remove insertions
    s = re.sub(r"\^.", "", s)              # remove read-start markers
    s = re.sub(r"[$<>]", "", s)            # remove read-end and other markers

    # Reference bases are '.' (forward) and ',' (reverse); everything
    # else is a mismatch.
    mutations = [b.upper() for b in s if b not in ".,"]
    base_counts = Counter(mutations)
    total = sum(base_counts.values())

    return (
        total,
        ",".join(sorted(base_counts.keys())),
        base_counts.get("A", 0),
        base_counts.get("T", 0),
        base_counts.get("C", 0),
        base_counts.get("G", 0),
        0,  # indel_count = 0 for SNVs
    )


# ---------------------------------------------------------------------------
# Match ALT to the right count column
# ---------------------------------------------------------------------------

def match_alts(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each row, look up the read count for the expected ALT allele
    and write it to a new ``read_num_match`` column.

    * For SNVs: reads the pre-computed ``{ALT}_count`` column
      (e.g. ``A_count``, ``T_count``).
    * For indels: uses ``indel_count``.

    Parameters
    ----------
    df : pd.DataFrame
        Must have ``OG_ALT``, ``A_count``, ``T_count``, ``C_count``,
        ``G_count``, ``indel_count`` columns.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with ``read_num_match`` column added.
    """
    df = df.copy()
    df["read_num_match"] = 0

    for idx, row in df.iterrows():
        og_alt = str(row["OG_ALT"]).upper()
        if len(og_alt) == 1:
            col = f"{og_alt}_count"
            if col in df.columns:
                df.at[idx, "read_num_match"] = row[col]
        else:
            df.at[idx, "read_num_match"] = row["indel_count"]

    return df


# ---------------------------------------------------------------------------
# Apply to full DataFrame
# ---------------------------------------------------------------------------

def apply_read_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply :func:`count_variants` and :func:`match_alts` to an entire
    pileup DataFrame.

    Adds columns: ``read_num``, ``ALT``, ``A_count``, ``T_count``,
    ``C_count``, ``G_count``, ``indel_count``, ``read_num_match``.

    Parameters
    ----------
    df : pd.DataFrame
        Merged pileup DataFrame (output of
        :func:`~daisy_mrd.mrd.pileup.merge_lspv_alts`).  Must have
        ``Read_Bases``, ``REF``, ``OG_ALT``.

    Returns
    -------
    pd.DataFrame
    """
    results = df.apply(
        lambda row: count_variants(
            row["Read_Bases"], row["REF"], row["OG_ALT"]
        ),
        axis=1,
    )

    (
        df["read_num"],
        df["ALT"],
        df["A_count"],
        df["T_count"],
        df["C_count"],
        df["G_count"],
        df["indel_count"],
    ) = zip(*results)

    df = match_alts(df)
    return df
