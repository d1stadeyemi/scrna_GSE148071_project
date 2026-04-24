"""
Stage 2 — Parse Per-Sample Expression Matrices to AnnData

GSE148071 provides 42 per-sample expression matrices:
  GSM4453576_P1_exp.txt.gz → Patient 1
  GSM4453577_P2_exp.txt.gz → Patient 2
  ...
  GSM4453617_P42_exp.txt.gz → Patient 42

Each file:
  - Rows = genes (29,528 genes)
  - Columns = cell barcodes
  - Values = UMI counts (integers)

Barcode uniqueness:
  Each sample's barcodes are prefixed with a sample-specific number
  (1_, 3_, 4_...) already embedded in the GEO files. These are globally
  unique within the merged matrix.

Expected output:
  ~89,887 cells × 29,527 genes before QC
  (paper reports 90,406 — minor GEO update since 2021)
"""

import os
import re
import gzip
import argparse
import logging

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def parse_sample_matrix(filepath, patient_id, gsm_id):
    """
    Parse one per-sample expression matrix using pandas.
    Format varies slightly across GEO submissions — pandas handles edge cases
    more robustly than manual line parsing.
    
    Matrix orientation: genes x cells (rows=genes, cols=cells)
    Output: cells x genes AnnData
    """
    import pandas as pd
    
    log.info(f"    Reading {os.path.basename(filepath)}...")
    
    # Read with pandas — handles variable whitespace, empty first fields, etc.
    df = pd.read_csv(
        filepath,
        sep        = "\t",
        index_col  = 0,    # first column = gene names
        compression = "gzip",
    )
    
    # df is now genes x cells — transpose to cells x genes
    genes    = df.index.tolist()
    barcodes = df.columns.tolist()
    
    log.info(f"    {len(genes):,} genes x {len(barcodes):,} cells")
    
    X = csr_matrix(df.values.T.astype(np.int32))
    
    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = genes
    
    adata.obs["patient_id"] = patient_id
    adata.obs["gsm_id"]     = gsm_id
    adata.obs["sample_id"]  = f"{gsm_id}_{patient_id}"
    
    return adata


def extract_patient_info(filename):
    """
    Extract patient ID and GSM ID from filename.
    Pattern: GSM4453576_P1_exp.txt.gz → ("GSM4453576", "P1")
    """
    basename = os.path.basename(filename)
    m = re.match(r"(GSM\d+)_(P\d+)_exp\.txt\.gz", basename)
    if m:
        return m.group(1), m.group(2)
    raise ValueError(f"Cannot parse patient info from filename: {basename}")


def main(input_dir, metadata_path, output_path):

    # Find all expression files
    files = sorted([
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith("_exp.txt.gz")
    ])
    log.info(f"Found {len(files)} expression files in {input_dir}")

    if len(files) == 0:
        raise FileNotFoundError(
            f"No *_exp.txt.gz files found in {input_dir}"
        )

    # Load metadata
    metadata = pd.read_csv(metadata_path, sep="\t")
    meta_lookup = dict(zip(metadata["gsm_id"], metadata.to_dict("records")))

    # Parse each sample
    adatas = []
    for i, filepath in enumerate(files, 1):
        gsm_id, patient_id = extract_patient_info(filepath)
        log.info(f"  [{i:>2}/{len(files)}] {patient_id} ({gsm_id})...")

        adata = parse_sample_matrix(filepath, patient_id, gsm_id)
        log.info(f"    {adata.n_obs:>5} cells x {adata.n_vars:,} genes")

        # Attach any metadata available
        if gsm_id in meta_lookup:
            for col, val in meta_lookup[gsm_id].items():
                if col not in ["patient_id", "gsm_id"]:
                    adata.obs[col] = val

        adatas.append(adata)

    # Merge all samples
    log.info("Merging all samples...")
    merged = ad.concat(adatas, join="outer", fill_value=0)
    merged.var_names_make_unique()

    log.info(f"Merged AnnData: {merged.n_obs:,} cells x {merged.n_vars:,} genes")
    log.info(f"Patients: {merged.obs['patient_id'].nunique()}")

    # Cell counts per patient
    log.info("Cells per patient:")
    for pid, count in merged.obs["patient_id"].value_counts().sort_index().items():
        log.info(f"  {pid}: {count:,}")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    merged.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved to {output_path}")
    log.info("Stage 2 complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing *_exp.txt.gz files")
    parser.add_argument("--metadata",  required=True,
                        help="Patient metadata TSV from parse_metadata.py")
    parser.add_argument("--output",    required=True,
                        help="Output h5ad path")
    args = parser.parse_args()
    main(args.input_dir, args.metadata, args.output)
