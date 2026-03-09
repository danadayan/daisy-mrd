"""
daisy_mrd.lspv
==============
LSPV identification pipeline: takes a diagnosis VCF and returns a
list of leukemia-specific passenger variants (LSPVs).

Modules
-------
filter      : VCF-level hard filters (PASS, SNV-only, germline removal)
annotate    : VEP (optional) and gnomAD annotation
reads       : Extract read depth and VAF from FORMAT fields
pon         : Panel of Normals binomial filter
gmm         : Gaussian Mixture Model clonality classification
identify    : Final LSPV extraction (clonal + non-coding)
pipeline    : End-to-end orchestrator
"""

from daisy_mrd.lspv.pipeline import run_lspv_pipeline  # noqa: F401

__all__ = ["run_lspv_pipeline"]
