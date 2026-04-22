"""
Stage 4b — scVI Integration (Alternative Track Only)

scVI (single-cell Variational Inference) models each cell's gene expression
as a function of:
  1. Biological variation (latent variables z)
  2. Technical variation (batch = patient_id)

The latent space z is batch-free and used for UMAP/clustering.

If a trained model already exists at model_dir, it is loaded instead
of retraining. This allows the post-training steps (latent extraction,
h5ad save) to rerun quickly without GPU.

Reference: Lopez et al. 2018 (Nature Methods)
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import scvi
import torch

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def main(input_path, output_path, model_dir,
         n_latent, n_epochs, batch_size):

    log.info(f"[ALTERNATIVE TRACK] Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Device logging
    log.info(f"  PyTorch: {torch.__version__}")
    log.info(f"  CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        log.info(f"  GPU: {torch.cuda.get_device_name(0)}")
    else:
        log.info("  Running on CPU (training will take ~60-90 min)")

    # ── Subset to HVGs for scVI ───────────────────────────────────────────────
    # scVI trains on HVGs only — reduces noise, speeds training
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    log.info(f"  HVG subset: {adata_hvg.n_vars:,} genes")

    # ── Setup scVI ────────────────────────────────────────────────────────────
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer     = "counts",       # raw counts required by scVI's NB model
        batch_key = "patient_id",   # 42 patients as batches
    )

    # ── Load or train ─────────────────────────────────────────────────────────
    model_pt = os.path.join(model_dir, "model.pt")

    if os.path.exists(model_pt):
        log.info(f"  Found existing model at {model_dir} — loading...")
        model = scvi.model.SCVI.load(model_dir, adata=adata_hvg)
        log.info("  Model loaded.")
    else:
        log.info("  No existing model — training from scratch...")
        model = scvi.model.SCVI(
            adata_hvg,
            n_latent = n_latent,
            n_layers = 2,
            n_hidden = 128,
        )
        log.info(f"  Architecture: n_latent={n_latent}, n_layers=2, n_hidden=128")
        log.info(f"  Batch key: patient_id ({adata_hvg.obs['patient_id'].nunique()} batches)")

        model.train(
            max_epochs              = n_epochs,
            batch_size              = batch_size,
            early_stopping          = True,
            early_stopping_patience = 10,
            plan_kwargs             = {"lr": 1e-3},
        )
        log.info("  Training complete.")

        # Log final ELBO
        try:
            elbo = float(model.history["elbo_train"].values[-1].flat[0])
            log.info(f"  Final train ELBO: {elbo:.4f}")
        except Exception:
            log.info("  ELBO not available from history")

        os.makedirs(model_dir, exist_ok=True)
        model.save(model_dir, overwrite=True)
        log.info(f"  Model saved to {model_dir}")

    # ── Extract latent representation ─────────────────────────────────────────
    log.info("  Extracting latent representation...")
    latent = model.get_latent_representation()
    adata.obsm["X_scvi"] = latent
    log.info(f"  X_scvi shape: {latent.shape}")

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"  Saved to {output_path}")
    log.info("Stage 4b [alternative] complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="scVI integration for alternative track"
    )
    parser.add_argument("--input",      required=True)
    parser.add_argument("--output",     required=True)
    parser.add_argument("--model_dir",  required=True)
    parser.add_argument("--n_latent",   type=int, default=30)
    parser.add_argument("--n_epochs",   type=int, default=400)
    parser.add_argument("--batch_size", type=int, default=128)
    args = parser.parse_args()
    main(args.input, args.output, args.model_dir,
         args.n_latent, args.n_epochs, args.batch_size)
