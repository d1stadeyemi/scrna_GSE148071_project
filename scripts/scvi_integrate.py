"""
Stage 4b — scVI Integration and Batch Correction

Input:  data/processed/normalized.h5ad
Output: data/processed/scvi_integrated.h5ad
        models/scvi_model/

Steps:
- Set up scVI model using raw counts + patient_id as batch key
- Train on GPU if available, fallback to CPU
- Extract latent representation (adata.obsm['X_scvi'])
- Save trained model for reproducibility
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import scvi

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def main(input_path, output_path, model_dir,
         n_latent, n_epochs, batch_size):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── Subset to HVGs for scVI ───────────────────────────────────────────────
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    log.info(f"Using {adata_hvg.n_vars:,} HVGs for scVI")

    # ── Set up scVI model ─────────────────────────────────────────────────────
    log.info("Setting up scVI model...")
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer      = "counts",       # raw counts required
        batch_key  = "patient_id",   # batch correction across 42 patients
    )

    model = scvi.model.SCVI(
        adata_hvg,
        n_latent = n_latent,
        n_layers = 2,
        n_hidden = 128,
    )

    log.info(f"Model architecture:")
    log.info(f"  n_latent: {n_latent}")
    log.info(f"  n_layers: 2")
    log.info(f"  n_hidden: 128")
    log.info(f"  batch_key: patient_id")

    # ── Train ─────────────────────────────────────────────────────────────────
    log.info("Training scVI model...")
    model.train(
        max_epochs          = n_epochs,
        batch_size          = batch_size,
        early_stopping      = True,
        early_stopping_patience = 10,
        plan_kwargs         = {"lr": 1e-3},
    )

    log.info("Training complete.")
    log.info(f"  Final ELBO: {model.history['elbo_train'].values[-1]}")

    # ── Extract latent representation ─────────────────────────────────────────
    log.info("Extracting scVI latent representation...")
    latent = model.get_latent_representation()
    adata.obsm["X_scvi"] = latent
    log.info(f"  Latent shape: {latent.shape}")

    # ── Save model ────────────────────────────────────────────────────────────
    os.makedirs(model_dir, exist_ok=True)
    model.save(model_dir, overwrite=True)
    log.info(f"Model saved to {model_dir}")

    # ── Save AnnData ──────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved integrated AnnData to {output_path}")
    log.info("Stage 4b complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",      default="data/processed/normalized.h5ad")
    parser.add_argument("--output",     default="data/processed/scvi_integrated.h5ad")
    parser.add_argument("--model_dir",  default="models/scvi_model")
    parser.add_argument("--n_latent",   type=int, default=30)
    parser.add_argument("--n_epochs",   type=int, default=400)
    parser.add_argument("--batch_size", type=int, default=128)
    args = parser.parse_args()
    main(args.input, args.output, args.model_dir,
         args.n_latent, args.n_epochs, args.batch_size)