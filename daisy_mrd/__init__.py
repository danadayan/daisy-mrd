"""
DAISY-MRD
=========
Distributed Analysis of Integrated Sites for Yielding MRD.

A WGS-based pipeline for measurable residual disease (MRD) monitoring
in pediatric AML using leukemia-specific passenger variants (LSPVs).

Pipeline overview
-----------------
Step 1  LSPV identification  (daisy_mrd.lspv)
    • Filter VCF (PASS variants only)
    • Optional VEP annotation
    • gnomAD germline removal
    • Panel of Normals (PoN) filtering
    • GMM-based clonality classification
    • LSPV extraction (pan-clonal, non-coding variants)

Step 2  MRD scoring  (daisy_mrd.mrd)
    • samtools mpileup over LSPV positions in remission BAM/CRAM
    • Multi-layer noise filtering (context, sex chr, depth, germline, PoN, VAF)
    • DAISY-MRD score = Σ(alt reads) / Σ(total depth) across all LSPV positions

Basic usage
-----------
>>> from daisy_mrd import run_lspv_pipeline, run_mrd_single
>>>
>>> # Step 1: identify LSPVs from diagnosis VCF
>>> lspv_result = run_lspv_pipeline(
...     vcf_path="patient_001_diagnosis.vcf",
...     output_dir="results/patient_001/",
...     patient_id="001",
... )
>>>
>>> # Step 2: compute DAISY-MRD score from remission BAM/CRAM
>>> mrd_result = run_mrd_single(
...     lspv_csv="results/patient_001/001_lspvs.csv",
...     remission_bam="patient_001_remission.cram",
...     reference="GRCh38.fa",
...     output_dir="results/patient_001/mrd/",
...     patient_id="001",
... )
>>> print(mrd_result.score.score)
"""

__version__ = "0.1.0"
__author__ = "Dana Dayan, Yosef Maruvka"

from daisy_mrd.lspv.pipeline import run_lspv_pipeline  # noqa: F401
from daisy_mrd.mrd.pipeline import run_mrd_cohort, run_mrd_single  # noqa: F401

__all__ = ["run_lspv_pipeline", "run_mrd_single", "run_mrd_cohort", "__version__"]
