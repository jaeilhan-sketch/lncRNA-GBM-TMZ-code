#!/usr/bin/env python3
"""
Phase 1 (continued): Download IDH mutation and MGMT methylation status.

Sources:
  - cBioPortal API for TCGA-GBM molecular annotations
  - GDC supplementary files (fallback)

Outputs:
  - data/clinical/molecular_annotations.tsv
  - data/clinical/sample_info_annotated.tsv  (merged with molecular data)
"""

import json
import requests
import pandas as pd
from pathlib import Path

PROJECT_DIR = Path("/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ")
CLINICAL_DIR = PROJECT_DIR / "data" / "clinical"

CBIOPORTAL_API = "https://www.cbioportal.org/api"
STUDY_ID = "gbm_tcga"


def get_molecular_annotations():
    """Query cBioPortal for IDH and MGMT status."""
    print("[1/3] Querying cBioPortal for molecular annotations...")

    # Get clinical data from cBioPortal
    url = f"{CBIOPORTAL_API}/studies/{STUDY_ID}/clinical-data"
    params = {
        "clinicalDataType": "SAMPLE",
        "projection": "DETAILED",
    }
    headers = {"Accept": "application/json"}

    response = requests.get(url, params=params, headers=headers)
    response.raise_for_status()
    clinical_data = response.json()

    # Parse into DataFrame
    rows = {}
    for entry in clinical_data:
        sid = entry["sampleId"]
        pid = entry["patientId"]
        attr = entry["clinicalAttributeId"]
        val = entry["value"]
        if sid not in rows:
            rows[sid] = {"sample_id_cbio": sid, "patient_id_cbio": pid}
        rows[sid][attr] = val

    df = pd.DataFrame(rows.values())

    # Select relevant columns (names vary by study)
    col_map = {}
    for col in df.columns:
        col_lower = col.lower()
        if "idh" in col_lower and "status" in col_lower:
            col_map[col] = "IDH_status"
        elif "idh" in col_lower and ("mutation" in col_lower or "mut" in col_lower):
            col_map[col] = "IDH_mutation"
        elif "mgmt" in col_lower and ("methyl" in col_lower or "status" in col_lower):
            col_map[col] = "MGMT_status"
        elif "subtype" in col_lower or "expression_subtype" in col_lower:
            col_map[col] = "expression_subtype"

    df_mol = df[["sample_id_cbio", "patient_id_cbio"] + list(col_map.keys())].copy()
    df_mol.rename(columns=col_map, inplace=True)

    out_path = CLINICAL_DIR / "molecular_annotations.tsv"
    df_mol.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    print(f"  Columns available: {list(df_mol.columns)}")

    return df_mol


def merge_with_sample_info(df_mol):
    """Merge molecular annotations with existing sample_info."""
    print("[2/3] Merging molecular annotations with sample info...")

    sample_info_path = CLINICAL_DIR / "sample_info.tsv"
    if not sample_info_path.exists():
        print("  ERROR: sample_info.tsv not found. Run 01_query_gdc_metadata.py first.")
        return None

    df_sample = pd.read_csv(sample_info_path, sep="\t")

    # Match via submitter_id → patient_id_cbio (TCGA barcode: TCGA-XX-XXXX)
    df_sample["patient_barcode"] = df_sample["submitter_id"].str[:12]
    df_mol["patient_barcode"] = df_mol["patient_id_cbio"].str[:12]

    df_merged = df_sample.merge(
        df_mol.drop(columns=["sample_id_cbio", "patient_id_cbio"], errors="ignore"),
        on="patient_barcode",
        how="left"
    )

    # Standardize IDH status
    if "IDH_status" in df_merged.columns:
        df_merged["IDH_status"] = df_merged["IDH_status"].fillna("Unknown")
    elif "IDH_mutation" in df_merged.columns:
        df_merged["IDH_status"] = df_merged["IDH_mutation"].fillna("Unknown")
    else:
        df_merged["IDH_status"] = "Unknown"

    # Standardize MGMT status
    if "MGMT_status" in df_merged.columns:
        df_merged["MGMT_status"] = df_merged["MGMT_status"].fillna("Unknown")
    else:
        df_merged["MGMT_status"] = "Unknown"

    print(f"\n  IDH status distribution:")
    print(f"  {df_merged['IDH_status'].value_counts().to_string()}")
    print(f"\n  MGMT status distribution:")
    print(f"  {df_merged['MGMT_status'].value_counts().to_string()}")

    return df_merged


def apply_exclusion_criteria(df):
    """Apply final exclusion criteria: remove IDH-mutant cases."""
    print("[3/3] Applying exclusion criteria...")

    n_before = len(df)

    # Exclude IDH-mutant (WHO 2021 classification)
    idh_mut_patterns = ["mutant", "mut", "yes", "positive", "idh1", "r132h"]
    df_filtered = df[~df["IDH_status"].str.lower().isin(idh_mut_patterns)].copy()

    n_excluded_idh = n_before - len(df_filtered)
    print(f"  Excluded IDH-mutant: {n_excluded_idh}")
    print(f"  Remaining samples: {len(df_filtered)}")

    # Final response distribution
    print(f"\n  Final TMZ Response distribution:")
    print(f"  {df_filtered['tmz_response'].value_counts().to_string()}")

    out_path = CLINICAL_DIR / "sample_info_final.tsv"
    df_filtered.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")

    return df_filtered


if __name__ == "__main__":
    df_mol = get_molecular_annotations()
    df_merged = merge_with_sample_info(df_mol)
    if df_merged is not None:
        df_final = apply_exclusion_criteria(df_merged)

        print("\n" + "=" * 60)
        print("Molecular annotation complete.")
        print("Review: data/clinical/sample_info_final.tsv")
        print("Next: Download raw FASTQ files with 03_download_fastq.sh")
        print("=" * 60)
