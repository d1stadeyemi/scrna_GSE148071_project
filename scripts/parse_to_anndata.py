"""
Stage 2 — Parse GEO expression matrices into a merged AnnData object.

Input:  42 per-sample .txt.gz files in data/raw/extracted/
Output: data/processed/merged_raw.h5ad

Assumptions (verified pre-run):
- All samples share identical gene sets (29,528 genes)
- Matrix orientation: genes x cells (rows x cols)
- First row: cell barcodes (with sample-specific numeric prefix)
- First column: gene names
- Values: raw integer counts
- Barcode format: {sample_index}_{barcode} — globally unique across samples
"""

import os
import re
import glob
import argparse
import logging

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def parse_filename(filepath):
    """Extract GSM accession and patient ID from filename.
    
    Expected format: GSM4453576_P1_exp.txt.gz
    """
    basename = os.path.basename(filepath)
    match = re.match(r"(GSM\d+)_(P\d+)_exp\.txt\.gz", basename)
    if not match:
        raise ValueError(f"Unexpected filename format: {basename}")
    return match.group(1), match.group(2)   # gsm_id, patient_id


def read_sample(filepath):
    """Read a single sample matrix into AnnData.
    
    Returns AnnData of shape (n_cells, n_genes) with sparse X matrix.
    """
    gsm_id, patient_id = parse_filename(filepath)

    df = pd.read_csv(filepath, sep="\t", index_col=0, compression="gzip")
    # df shape: (n_genes, n_cells) — transpose for AnnData convention
    
    adata = ad.AnnData(
        X        = csr_matrix(df.values.T.astype(np.float32)),
        obs      = pd.DataFrame(index=df.columns.tolist()),
        var      = pd.DataFrame(index=df.index.tolist()),
    )

    adata.obs["patient_id"]  = patient_id
    adata.obs["gsm_id"]      = gsm_id
    adata.obs["sample_id"]   = f"{gsm_id}_{patient_id}"

    log.info(f"  {patient_id} ({gsm_id}): {adata.n_obs:>5} cells x {adata.n_vars} genes")
    return adata


def main(extracted_dir, output_path):
    files = sorted(glob.glob(os.path.join(extracted_dir, "*.txt.gz")))
    
    if len(files) == 0:
        raise FileNotFoundError(f"No .txt.gz files found in {extracted_dir}")
    
    log.info(f"Found {len(files)} sample files")

    # ── Parse all samples ────────────────────────────────────────────────────
    adatas = []
    for f in files:
        adata = read_sample(f)
        adatas.append(adata)

    total_cells = sum(a.n_obs for a in adatas)
    log.info(f"Total cells before merge: {total_cells:,}")

    # ── Concatenate ──────────────────────────────────────────────────────────
    # inner join is safe here — confirmed identical gene sets across all samples
    log.info("Concatenating all samples (inner join)...")
    merged = ad.concat(
        adatas,
        join    = "inner",
        merge   = "same",
    )

    # ── Sanity checks ────────────────────────────────────────────────────────
    assert merged.n_obs == total_cells, \
        f"Cell count mismatch after concat: {merged.n_obs} vs {total_cells}"
    assert merged.n_vars == adatas[0].n_vars, \
        f"Gene count mismatch after concat: {merged.n_vars} vs {adatas[0].n_vars}"
    assert merged.obs["patient_id"].nunique() == len(files), \
        f"Expected {len(files)} patients, got {merged.obs['patient_id'].nunique()}"

    log.info(f"Merged AnnData: {merged.n_obs:,} cells x {merged.n_vars:,} genes")
    log.info(f"Samples represented: {merged.obs['patient_id'].nunique()}")

    # ── Cell counts per patient ──────────────────────────────────────────────
    counts = merged.obs["patient_id"].value_counts().sort_index()
    log.info("Cell counts per patient:\n" + counts.to_string())

    # ── Save ─────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    merged.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved to {output_path}")

    # ── Cleanup extracted files to save HPC storage ──────────────────────────
    log.info("Cleaning up extracted .txt.gz files...")
    for f in files:
        os.remove(f)
    
    extracted_dir_path = extracted_dir.rstrip("/")
    if os.path.isdir(extracted_dir_path) and not os.listdir(extracted_dir_path):
        os.rmdir(extracted_dir_path)
        log.info(f"Removed empty directory: {extracted_dir_path}")

    log.info("Stage 2 complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse GEO scRNA-seq matrices to AnnData"
    )
    parser.add_argument(
        "--extracted_dir",
        default = "data/raw/extracted",
        help    = "Directory containing per-sample .txt.gz files"
    )
    parser.add_argument(
        "--output",
        default = "data/processed/merged_raw.h5ad",
        help    = "Output path for merged AnnData (.h5ad)"
    )
    args = parser.parse_args()
    main(args.extracted_dir, args.output)