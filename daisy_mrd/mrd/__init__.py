"""
daisy_mrd.mrd
=============
DAISY-MRD scoring pipeline (Step 2).

Takes the LSPV list from Step 1 and a remission BAM/CRAM file,
piles up reads at LSPV positions using ``samtools mpileup``,
applies multi-layer noise filtering, and computes the DAISY-MRD score.

Modules
-------
pileup      : Run samtools mpileup + parse pileup files
readcount   : Count ALT-supporting reads from Read_Bases strings
filters     : Multi-layer noise filtering (context, sex chr, depth,
              germline, PoN, high-VAF)
score       : DAISY-MRD score calculation + noise distribution plots
pipeline    : End-to-end orchestrators (single patient + cohort)
"""

from daisy_mrd.mrd.pipeline import run_mrd_cohort, run_mrd_single  # noqa: F401

__all__ = ["run_mrd_single", "run_mrd_cohort"]
