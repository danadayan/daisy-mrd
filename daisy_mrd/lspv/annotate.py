"""
daisy_mrd.lspv.annotate
=======================
Variant annotation: VEP (optional) and gnomAD germline frequency lookup.

VEP
---
VEP annotation is **optional**. Users who have already run VEP externally
(or who use a different annotation tool) can skip this step entirely and
pass their annotated VCF directly to the pipeline.

When ``run_vep`` is called it invokes VEP via **Docker**, using the
``ensemblorg/ensembl-vep`` image. A local VEP cache directory can be
supplied; otherwise VEP queries its online database (slower).

gnomAD
------
gnomAD annotation adds two columns to the VCF:

* ``GNOMAD``   — ``"YES"`` if the position is found in gnomAD, else ``"NO"``
* ``gnomad_AF``— allele frequency reported by gnomAD (``"."`` if not found)

Lookup can be performed against:

1. A **local gnomAD index** (pre-built ``chrom_pos.txt`` files; fastest).
2. The **gnomAD GraphQL API** (no local files needed; slower, rate-limited).

Both modes can be combined: local lookup first, API as fallback.
"""

from __future__ import annotations

import glob
import os
import subprocess
import time
from pathlib import Path

import requests

GNOMAD_API = "https://gnomad.broadinstitute.org/api"


# ---------------------------------------------------------------------------
# VEP (optional)
# ---------------------------------------------------------------------------

def run_vep(
    input_vcf: str | Path,
    output_vcf: str | Path,
    cache_dir: str | Path = "~/.vep",
    assembly: str = "GRCh38",
    use_cache: bool = False,
) -> subprocess.CompletedProcess:
    """
    Annotate a VCF with Ensembl VEP via Docker.

    This step is **entirely optional**. If you have already annotated your
    VCF with VEP (or a different tool), skip this function and pass the
    annotated VCF directly to :func:`annotate_vcf` or the pipeline.

    The Docker image ``ensemblorg/ensembl-vep`` must be available. The
    input VCF and output VCF must reside in the same directory, because
    only that directory is mounted into the container.

    Parameters
    ----------
    input_vcf : str or Path
    output_vcf : str or Path
    cache_dir : str or Path
        Path to a local VEP cache directory. Only used when
        ``use_cache=True``. Defaults to ``~/.vep``.
    assembly : str
        Genome assembly, e.g. ``"GRCh38"`` or ``"GRCh37"``.
    use_cache : bool
        If ``True``, pass ``--cache --offline`` to VEP (requires a local
        cache download). If ``False``, VEP queries its online database.

    Returns
    -------
    subprocess.CompletedProcess
        Contains ``returncode``, ``stdout`` and ``stderr``.

    Raises
    ------
    RuntimeError
        If VEP exits with a non-zero return code.
    """
    input_vcf = Path(input_vcf).resolve()
    output_vcf = Path(output_vcf).resolve()
    workdir = input_vcf.parent
    cache_dir = Path(cache_dir).expanduser().resolve()

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{cache_dir}:/opt/vep/.vep",
        "-v", f"{workdir}:/data",
        "ensemblorg/ensembl-vep",
        "vep",
        "--assembly", assembly,
        "--input_file", f"/data/{input_vcf.name}",
        "--output_file", f"/data/{output_vcf.name}",
        "--vcf",
        "--symbol",
        "--hgvs",
        "--canonical",
        "--pick",
        "--everything",
        "--no_stats",
    ]

    if use_cache:
        cmd += ["--cache", "--offline"]
    else:
        cmd += ["--database"]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"VEP failed (return code {result.returncode}).\n"
            f"STDERR:\n{result.stderr}"
        )

    return result


# ---------------------------------------------------------------------------
# gnomAD annotation
# ---------------------------------------------------------------------------

def _query_gnomad_api(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    dataset: str = "gnomad_r3",
) -> tuple[bool, float | None]:
    """Query the gnomAD GraphQL API for a single variant."""
    query = """
    query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
      variant(variantId: $variantId, dataset: $datasetId) {
        genome { af }
      }
    }
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    try:
        r = requests.post(
            GNOMAD_API,
            json={"query": query, "variables": {"variantId": variant_id, "datasetId": dataset}},
            timeout=10,
        )
        data = r.json()
        variant = data.get("data", {}).get("variant")
        if variant and variant.get("genome"):
            return True, variant["genome"]["af"]
    except Exception:
        pass
    return False, None


def _load_local_gnomad(gnomad_path: str | Path) -> dict[str, set[int]]:
    """
    Load a pre-built gnomAD position index into memory.

    Expects files matching the pattern::

        gnomad.genomes.v3.1.2.sites.chr*.vcf.bgz.chrom_pos.txt

    Each file contains whitespace-separated ``(chrom, pos)`` pairs,
    one per line.

    Parameters
    ----------
    gnomad_path : str or Path
        Directory containing the ``chrom_pos.txt`` files.

    Returns
    -------
    dict[str, set[int]]
        Mapping of chromosome name → set of positions.
    """
    chrom_dict: dict[str, set[int]] = {}
    pattern = os.path.join(
        str(gnomad_path),
        "gnomad.genomes.v3.1.2.sites.chr*.vcf.bgz.chrom_pos.txt",
    )
    for txt_file in glob.glob(pattern):
        # Extract chromosome from filename, e.g. "chr1"
        chrom = os.path.basename(txt_file).split(".")[6]
        positions: set[int] = set()
        with open(txt_file) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    positions.add(int(parts[1]))
        chrom_dict[chrom] = positions
    return chrom_dict


def annotate_vcf(
    input_vcf: str | Path,
    output_vcf: str | Path,
    gnomad_path: str | Path | None = None,
    use_api: bool = True,
    api_dataset: str = "gnomad_r3",
    api_sleep: float = 0.05,
) -> Path:
    """
    Annotate a VCF with gnomAD presence and allele frequency.

    Adds two tab-separated columns to each variant line:

    * ``GNOMAD``    — ``"YES"`` / ``"NO"``
    * ``gnomad_AF`` — numeric AF or ``"."``

    Parameters
    ----------
    input_vcf : str or Path
    output_vcf : str or Path
    gnomad_path : str, Path, or None
        Directory containing pre-built ``chrom_pos.txt`` index files
        (see :func:`_load_local_gnomad`). When ``None``, local lookup
        is skipped.
    use_api : bool
        Whether to fall back to the gnomAD GraphQL API for variants
        not found in the local index (or when no local index is given).
        Disable to avoid network calls.
    api_dataset : str
        gnomAD dataset identifier passed to the API, e.g. ``"gnomad_r3"``
        (GRCh38) or ``"gnomad_r2_1"`` (GRCh37).
    api_sleep : float
        Seconds to sleep between API requests to avoid rate-limiting.

    Returns
    -------
    Path
        Path to the written output file.
    """
    input_vcf = Path(input_vcf)
    output_vcf = Path(output_vcf)

    local_db: dict[str, set[int]] | None = None
    if gnomad_path is not None:
        local_db = _load_local_gnomad(gnomad_path)

    with open(input_vcf) as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            # Preserve meta-info lines
            if line.startswith("##"):
                outfile.write(line)
                continue

            # Augment column header
            if line.startswith("#CHROM"):
                outfile.write(line.rstrip() + "\tGNOMAD\tgnomad_AF\n")
                continue

            fields = line.strip().split("\t")
            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            gnomad_found = False
            gnomad_af: float | None = None

            # --- Local lookup ---
            if local_db is not None:
                if chrom in local_db and pos in local_db[chrom]:
                    gnomad_found = True

            # --- API fallback ---
            if not gnomad_found and use_api:
                gnomad_found, gnomad_af = _query_gnomad_api(
                    chrom, pos, ref, alt, dataset=api_dataset
                )
                time.sleep(api_sleep)

            gnomad_flag = "YES" if gnomad_found else "NO"
            af_str = str(gnomad_af) if gnomad_af is not None else "."

            outfile.write(line.rstrip() + f"\t{gnomad_flag}\t{af_str}\n")

    return output_vcf
