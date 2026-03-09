"""
daisy_mrd.utils
===============
Shared utility functions used across the pipeline.
"""

from __future__ import annotations

import importlib.resources
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# VCF reader
# ---------------------------------------------------------------------------

def read_vcf(vcf_path: str | Path) -> pd.DataFrame:
    """
    Read a VCF file into a pandas DataFrame.

    Header lines (``##``) are skipped. The column-header line (``#CHROM``)
    is used to name columns, with the leading ``#`` stripped so the first
    column is named ``CHROM``.

    Extra tab-separated columns appended after standard VCF columns
    (e.g. by ``annotate_vcf``) are preserved.

    Parameters
    ----------
    vcf_path : str or Path
        Path to the VCF file. Plain ``.vcf`` and gzip-compressed
        ``.vcf.gz`` files are both supported.

    Returns
    -------
    pd.DataFrame
        One row per variant. All values are strings; callers are
        responsible for type conversion.

    Raises
    ------
    ValueError
        If the file contains no ``#CHROM`` header line.
    """
    vcf_path = Path(vcf_path)
    open_fn = _gzip_open if vcf_path.suffix == ".gz" else open

    header: list[str] = []
    rows: list[list[str]] = []

    with open_fn(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                continue
            rows.append(line.strip().split("\t"))

    if not header:
        raise ValueError(f"No #CHROM header line found in {vcf_path}")

    # Pad short rows so DataFrame construction never fails
    n_cols = len(header)
    padded = [row + [""] * (n_cols - len(row)) for row in rows]

    return pd.DataFrame(padded, columns=header)


def _gzip_open(path: Path, mode: str = "rt"):
    import gzip
    return gzip.open(path, mode)


# ---------------------------------------------------------------------------
# Bundled data helpers
# ---------------------------------------------------------------------------

def get_builtin_pon_path() -> Path:
    """
    Return the path to the built-in Panel of Normals (PoN) CSV bundled
    with the package (``daisy_mrd/data/pon_default.csv``).

    Returns
    -------
    Path

    Raises
    ------
    FileNotFoundError
        If the bundled PoN file is missing from the installation.
    """
    try:
        # Python 3.9+ importlib.resources API
        ref = importlib.resources.files("daisy_mrd.data").joinpath("pon_default.csv")
        path = Path(str(ref))
    except AttributeError:
        # Fallback for Python 3.8
        import pkg_resources
        path = Path(pkg_resources.resource_filename("daisy_mrd", "data/pon_default.csv"))

    if not path.exists():
        raise FileNotFoundError(
            f"Built-in PoN not found at {path}. "
            "If you installed from source, make sure 'data/pon_default.csv' exists."
        )
    return path


def resolve_pon_path(pon_path: str | Path | None) -> Path:
    """
    Resolve which PoN to use.

    Parameters
    ----------
    pon_path : str, Path, or None
        User-supplied PoN path. Pass ``None`` to use the built-in PoN.

    Returns
    -------
    Path
        Resolved path to the PoN CSV.
    """
    if pon_path is None:
        return get_builtin_pon_path()
    return Path(pon_path)


# ---------------------------------------------------------------------------
# Output directory helper
# ---------------------------------------------------------------------------

def ensure_output_dir(output_dir: str | Path) -> Path:
    """
    Create *output_dir* (and parents) if it does not exist, then return it
    as a ``Path``.
    """
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path
