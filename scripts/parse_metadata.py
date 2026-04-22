"""
Parse GEO Series Matrix — Extract Patient Metadata

GEO provides only age and gender for GSE148071.
Clinical metadata (histology, mutation, smoking, stage) is in Supplementary
Table 1 of the paper (PDF only, not machine-readable from GEO).

This script extracts what is available from the series matrix file and
creates a clean TSV for downstream use.
"""

import os
import argparse
import logging
import gzip
import re
import pandas as pd

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def parse_series_matrix(filepath):
    """
    Parse !Sample_characteristics_ch1 rows from GEO series matrix.
    Returns dict of {field: [values per sample]}.
    """
    log.info(f"Parsing {filepath}...")
    characteristics = {}
    sample_ids = []

    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.strip()

            # Sample IDs row
            if line.startswith('"ID_REF"') or line.startswith("!Sample_geo_accession"):
                parts = line.split("\t")
                sample_ids = [p.strip('"') for p in parts[1:]]

            # Characteristics rows
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                values = [p.strip('"') for p in parts[1:]]

                # Determine field from first value
                if values:
                    first = values[0]
                    if ":" in first:
                        field = first.split(":")[0].strip()
                        parsed = [v.split(":", 1)[1].strip()
                                  if ":" in v else ""
                                  for v in values]
                        characteristics[field] = parsed

    log.info(f"  Found {len(sample_ids)} samples")
    log.info(f"  Available fields: {list(characteristics.keys())}")
    return sample_ids, characteristics


def build_metadata_table(sample_ids, characteristics):
    """
    Build patient metadata dataframe from parsed fields.
    Maps GSM IDs to patient IDs (P1-P42 in order).
    """
    # GSM IDs map sequentially to patient IDs
    # GSM4453576=P1, GSM4453577=P2, ..., GSM4453617=P42
    patient_ids = [f"P{i+1}" for i in range(len(sample_ids))]

    df = pd.DataFrame({
        "patient_id": patient_ids,
        "gsm_id":     sample_ids,
    })

    # Add available characteristics
    for field, values in characteristics.items():
        if len(values) == len(sample_ids):
            col_name = field.lower().replace(" ", "_")
            df[col_name] = values
        else:
            log.warning(f"  Field '{field}' has {len(values)} values "
                        f"(expected {len(sample_ids)}) — skipping")

    log.info(f"  Metadata table: {df.shape[0]} rows x {df.shape[1]} cols")
    log.info(f"  Columns: {list(df.columns)}")
    return df


def main(input_path, output_path):
    sample_ids, characteristics = parse_series_matrix(input_path)
    df = build_metadata_table(sample_ids, characteristics)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    log.info(f"Saved metadata to {output_path}")
    log.info("\nMetadata preview:")
    log.info(df.to_string(index=False))

    # Note about missing clinical data
    log.info(
        "\nNOTE: Per-patient clinical metadata (histology, mutation, smoking,"
        "\nstage, treatment) is only available in Supplementary Table 1 of the"
        "\npaper (PDF format). It is not machine-readable from GEO."
        "\nThese fields are required for Fig 2d/2e, 5a/5b analysis."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True,
                        help="Path to GSE148071_series_matrix.txt.gz")
    parser.add_argument("--output", required=True,
                        help="Output TSV path")
    args = parser.parse_args()
    main(args.input, args.output)
