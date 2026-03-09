"""
daisy_mrd.lspv.filter
=====================
VCF-level hard filters applied before clonality analysis.

Functions
---------
filter_vcf_pass         : Keep only PASS variants (file → file)
remove_germline_variants: Drop variants present in gnomAD (AF ≥ 1e-3)
remove_rs               : Drop variants with dbSNP rs IDs
remove_indels           : Keep SNVs only (drop indels via VEP SO term)
apply_hard_filters      : Convenience wrapper that applies all three
                          DataFrame-level filters in sequence
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# File-level filter (PASS)
# ---------------------------------------------------------------------------

def filter_vcf_pass(input_vcf: str | Path, output_vcf: str | Path) -> Path:
    """
    Write a new VCF containing only variants whose FILTER column is ``PASS``.

    All header lines (starting with ``#``) are preserved unchanged.

    Parameters
    ----------
    input_vcf : str or Path
        Path to the input VCF (plain or ``.gz``).
    output_vcf : str or Path
        Path where the filtered VCF will be written.

    Returns
    -------
    Path
        Path to the written output file.
    """
    input_vcf = Path(input_vcf)
    output_vcf = Path(output_vcf)

    written = 0
    with open(input_vcf) as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            fields = line.strip().split("\t")
            if len(fields) > 6 and fields[6] == "PASS":
                outfile.write(line)
                written += 1

    return output_vcf


# ---------------------------------------------------------------------------
# DataFrame-level filters
# ---------------------------------------------------------------------------

def remove_germline_variants(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove variants that are likely germline based on gnomAD allele frequency.

    A variant is kept if **either**:

    * Its gnomAD allele frequency (``AF`` column, added by ``annotate_vcf``)
      is below 1 × 10⁻³, **or**
    * It was not found in gnomAD at all (``GNOMAD`` column == ``"NO"``).

    Also restricts to biallelic SNVs by requiring ``ALT`` to be a single
    character (removes multi-allelic sites before the indel filter).

    Parameters
    ----------
    vcf_df : pd.DataFrame
        DataFrame produced by :func:`~daisy_mrd.utils.read_vcf` after
        gnomAD annotation.

    Returns
    -------
    pd.DataFrame
        Filtered copy of the input DataFrame.
    """
    df = vcf_df.copy()

    # Biallelic SNVs only (single-character ALT)
    df = df[df["ALT"].str.len() == 1]

    df["AF"] = pd.to_numeric(df.get("AF", pd.Series(dtype=float)), errors="coerce")

    gnomad_col = df.get("GNOMAD", pd.Series("NO", index=df.index))
    keep = (df["AF"] < 1e-3) | (gnomad_col == "NO")

    return df[keep].copy()


def remove_rs(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove variants that carry a dbSNP ``rs`` identifier in the ID column.

    These are known germline variants regardless of their gnomAD frequency.

    Parameters
    ----------
    vcf_df : pd.DataFrame

    Returns
    -------
    pd.DataFrame
    """
    mask = ~vcf_df["ID"].astype(str).str.startswith("rs")
    return vcf_df[mask].copy()


def remove_indels(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove indel variants detected via the ``VARIANT_TYPE`` tag in the INFO
    field (set by VEP or the variant caller).

    If ``VARIANT_TYPE`` is absent from INFO, the variant is kept (it is
    assumed to be an SNV).

    Parameters
    ----------
    vcf_df : pd.DataFrame

    Returns
    -------
    pd.DataFrame
    """
    def _is_indel(info_str: str) -> bool:
        for part in str(info_str).split(";"):
            if part.startswith("VARIANT_TYPE="):
                return "indel" in part.split("=", 1)[1].lower()
        return False

    mask = ~vcf_df["INFO"].apply(_is_indel)
    return vcf_df[mask].copy()


def apply_hard_filters(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply all three DataFrame-level hard filters in sequence:

    1. :func:`remove_germline_variants`
    2. :func:`remove_rs`
    3. :func:`remove_indels`

    Parameters
    ----------
    vcf_df : pd.DataFrame
        Annotated VCF DataFrame (must have ``GNOMAD`` and ``AF`` columns
        if gnomAD annotation was run; otherwise germline removal falls
        back to gnomAD-absent logic only).

    Returns
    -------
    pd.DataFrame
        Hard-filtered DataFrame.
    """
    df = remove_germline_variants(vcf_df)
    df = remove_rs(df)
    df = remove_indels(df)
    return df
