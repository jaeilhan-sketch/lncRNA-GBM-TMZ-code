#!/usr/bin/env python3
"""
Phase 1: Query GDC API for TCGA-GBM RNA-seq metadata and clinical data.

Outputs:
  - data/clinical/gdc_rnaseq_manifest.tsv     (file manifest)
  - data/clinical/clinical_metadata_raw.tsv    (raw clinical data)
  - data/clinical/sample_info.tsv              (curated sample info)
"""

import json
import os
import sys
import requests
import pandas as pd
from pathlib import Path

# ── Paths ──
PROJECT_DIR = Path("/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ")
CLINICAL_DIR = PROJECT_DIR / "data" / "clinical"
CLINICAL_DIR.mkdir(parents=True, exist_ok=True)

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_CASES_ENDPOINT = "https://api.gdc.cancer.gov/cases"

# ============================================================================
# 1. Query RNA-seq file metadata from GDC
# ============================================================================
def query_rnaseq_files():
    """Query GDC for TCGA-GBM RNA-seq files."""
    print("[1/3] Querying GDC for TCGA-GBM RNA-seq files...")

    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-GBM"}},
            {"op": "=", "content": {"field": "data_category", "value": "Transcriptome Profiling"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "experimental_strategy", "value": "RNA-Seq"}},
            {"op": "=", "content": {"field": "access", "value": "open"}},
        ]
    }

    fields = [
        "file_id", "file_name", "file_size", "data_format",
        "cases.case_id", "cases.submitter_id",
        "cases.samples.sample_id", "cases.samples.sample_type",
        "cases.samples.portions.analytes.aliquots.aliquot_id",
        "analysis.workflow_type",
    ]

    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": 1000,
    }

    response = requests.get(GDC_FILES_ENDPOINT, params=params)
    response.raise_for_status()
    data = response.json()

    hits = data["data"]["hits"]
    print(f"  Found {len(hits)} RNA-seq files.")

    rows = []
    for hit in hits:
        case = hit["cases"][0] if hit.get("cases") else {}
        sample = case.get("samples", [{}])[0] if case.get("samples") else {}
        rows.append({
            "file_id": hit["file_id"],
            "file_name": hit["file_name"],
            "file_size": hit.get("file_size", 0),
            "data_format": hit.get("data_format", ""),
            "workflow_type": hit.get("analysis", {}).get("workflow_type", ""),
            "case_id": case.get("case_id", ""),
            "submitter_id": case.get("submitter_id", ""),
            "sample_id": sample.get("sample_id", ""),
            "sample_type": sample.get("sample_type", ""),
        })

    df_files = pd.DataFrame(rows)
    out_path = CLINICAL_DIR / "gdc_rnaseq_manifest.tsv"
    df_files.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    return df_files


# ============================================================================
# 2. Query clinical data for all TCGA-GBM cases
# ============================================================================
def query_clinical_data():
    """Query comprehensive clinical data from GDC."""
    print("[2/3] Querying GDC for clinical metadata...")

    filters = {
        "op": "=",
        "content": {"field": "project.project_id", "value": "TCGA-GBM"}
    }

    fields = [
        "case_id", "submitter_id",
        "demographic.gender", "demographic.race", "demographic.ethnicity",
        "demographic.vital_status", "demographic.days_to_death",
        "diagnoses.age_at_diagnosis",
        "diagnoses.primary_diagnosis",
        "diagnoses.tumor_grade",
        "diagnoses.tumor_stage",
        "diagnoses.morphology",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.days_to_recurrence",
        "diagnoses.progression_or_recurrence",
        "diagnoses.treatments.treatment_type",
        "diagnoses.treatments.therapeutic_agents",
        "diagnoses.treatments.treatment_outcome_of_treatment",
        "diagnoses.treatments.days_to_treatment_start",
        "diagnoses.treatments.days_to_treatment_end",
    ]

    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": 1000,
    }

    response = requests.get(GDC_CASES_ENDPOINT, params=params)
    response.raise_for_status()
    data = response.json()

    hits = data["data"]["hits"]
    print(f"  Found {len(hits)} cases.")

    rows = []
    for hit in hits:
        demo = hit.get("demographic", {})
        diag = hit.get("diagnoses", [{}])[0] if hit.get("diagnoses") else {}
        treatments = diag.get("treatments", [])

        # Find TMZ treatment
        tmz_treatment = None
        for tx in treatments:
            agents = tx.get("therapeutic_agents", "") or ""
            if "temozolomide" in agents.lower() or "tmz" in agents.lower():
                tmz_treatment = tx
                break

        rows.append({
            "case_id": hit["case_id"],
            "submitter_id": hit.get("submitter_id", ""),
            "gender": demo.get("gender", ""),
            "race": demo.get("race", ""),
            "vital_status": demo.get("vital_status", ""),
            "days_to_death": demo.get("days_to_death", None),
            "age_at_diagnosis": diag.get("age_at_diagnosis", None),
            "primary_diagnosis": diag.get("primary_diagnosis", ""),
            "tumor_grade": diag.get("tumor_grade", ""),
            "days_to_last_follow_up": diag.get("days_to_last_follow_up", None),
            "days_to_recurrence": diag.get("days_to_recurrence", None),
            "progression_or_recurrence": diag.get("progression_or_recurrence", ""),
            "n_treatments": len(treatments),
            "tmz_treated": tmz_treatment is not None,
            "tmz_outcome": tmz_treatment.get("treatment_outcome_of_treatment", "") if tmz_treatment else "",
            "tmz_days_start": tmz_treatment.get("days_to_treatment_start", None) if tmz_treatment else None,
            "tmz_days_end": tmz_treatment.get("days_to_treatment_end", None) if tmz_treatment else None,
            "all_treatments": "; ".join(
                [f"{tx.get('treatment_type', '')}:{tx.get('therapeutic_agents', '')}" for tx in treatments]
            ),
        })

    df_clinical = pd.DataFrame(rows)
    out_path = CLINICAL_DIR / "clinical_metadata_raw.tsv"
    df_clinical.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    return df_clinical


# ============================================================================
# 3. Curate sample info: define TMZ response groups
# ============================================================================
def curate_sample_info(df_files, df_clinical, pfs_cutoff_months=6):
    """Merge file + clinical data and define TMZ response groups."""
    print("[3/3] Curating sample info and defining TMZ response groups...")

    # Merge on case_id
    df = df_files.merge(df_clinical, on=["case_id", "submitter_id"], how="inner")

    # Filter: Primary Tumor only
    df = df[df["sample_type"] == "Primary Tumor"].copy()
    print(f"  Primary Tumor samples: {len(df)}")

    # Filter: TMZ treated
    df_tmz = df[df["tmz_treated"]].copy()
    print(f"  TMZ-treated samples: {len(df_tmz)}")

    # Calculate PFS (days)
    # PFS = days_to_recurrence (if available) or days_to_death or days_to_last_follow_up
    def calc_pfs(row):
        if pd.notna(row["days_to_recurrence"]) and row["days_to_recurrence"] > 0:
            return row["days_to_recurrence"]
        elif row["vital_status"] == "Dead" and pd.notna(row["days_to_death"]):
            return row["days_to_death"]
        elif pd.notna(row["days_to_last_follow_up"]):
            return row["days_to_last_follow_up"]
        return None

    df_tmz["pfs_days"] = df_tmz.apply(calc_pfs, axis=1)
    df_tmz["pfs_months"] = df_tmz["pfs_days"] / 30.44  # approximate

    # Define response groups
    cutoff_days = pfs_cutoff_months * 30.44
    df_tmz["tmz_response"] = df_tmz["pfs_days"].apply(
        lambda x: "Responder" if pd.notna(x) and x >= cutoff_days
        else ("NonResponder" if pd.notna(x) else "Unknown")
    )

    # Calculate OS
    def calc_os(row):
        if row["vital_status"] == "Dead" and pd.notna(row["days_to_death"]):
            return row["days_to_death"], 1
        elif pd.notna(row["days_to_last_follow_up"]):
            return row["days_to_last_follow_up"], 0
        return None, None

    os_data = df_tmz.apply(calc_os, axis=1, result_type="expand")
    df_tmz["os_days"] = os_data[0]
    df_tmz["os_status"] = os_data[1]  # 1=dead, 0=censored

    # Age group
    df_tmz["age_years"] = df_tmz["age_at_diagnosis"] / 365.25
    df_tmz["age_group"] = pd.cut(
        df_tmz["age_years"], bins=[0, 50, 65, 100], labels=["<50", "50-65", ">65"]
    )

    # Summary
    print(f"\n  TMZ Response distribution:")
    print(f"  {df_tmz['tmz_response'].value_counts().to_string()}")
    print(f"\n  Samples with PFS data: {df_tmz['pfs_days'].notna().sum()}")

    # Remove Unknown response
    df_final = df_tmz[df_tmz["tmz_response"] != "Unknown"].copy()
    print(f"\n  Final cohort: {len(df_final)} samples")

    out_path = CLINICAL_DIR / "sample_info.tsv"
    df_final.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")

    # Also save a summary table for Table 1
    summary_path = CLINICAL_DIR / "cohort_summary.tsv"
    summary = df_final.groupby("tmz_response").agg(
        n=("case_id", "count"),
        age_mean=("age_years", "mean"),
        age_sd=("age_years", "std"),
        male_pct=("gender", lambda x: (x == "male").mean() * 100),
        pfs_median=("pfs_months", "median"),
        os_median_days=("os_days", "median"),
    ).round(2)
    summary.to_csv(summary_path, sep="\t")
    print(f"  Saved: {summary_path}")

    return df_final


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    df_files = query_rnaseq_files()
    df_clinical = query_clinical_data()
    df_sample = curate_sample_info(df_files, df_clinical, pfs_cutoff_months=6)

    print("\n" + "=" * 60)
    print("Phase 1 complete.")
    print("Next steps:")
    print("  1. Review data/clinical/sample_info.tsv")
    print("  2. Check MGMT methylation status (cBioPortal or GDC supplements)")
    print("  3. Check IDH mutation status and exclude IDH-mutant cases")
    print("  4. Download raw FASTQ via gdc-client (controlled access)")
    print("=" * 60)
