#!/usr/bin/env python
"""
main.py — Transcriptomic Drug Repurposing Engine (TDRE)

Usage:
    python main.py "Disease Name"

Example:
    python main.py "Friedreich ataxia"
"""

import os
import sys
import argparse

# Load .env BEFORE tool imports so all os.getenv() calls inside functions
# see the values from the .env file.
from dotenv import load_dotenv
load_dotenv()

from tools.disease_lookup  import lookup_disease
from tools.geo_search      import (search_geo_datasets,
                                   score_dataset_suitability,
                                   select_dataset_interactive)
from tools.geo_download    import download_geo_dataset, assign_sample_groups_interactive
from tools.deg_analysis    import run_deseq2
from tools.signature       import build_signature
from tools.drug_repurposing import find_drug_candidates
from utils.provenance      import ProvenanceLog


# ── Helpers ───────────────────────────────────────────────────────────────────

def slugify(name: str) -> str:
    """Convert a disease name to a filesystem-safe directory name."""
    return (name.lower()
                .replace(" ", "_")
                .replace("/", "-")
                .replace("'", "")
                .replace("(", "")
                .replace(")", ""))


def make_narrative(prov: ProvenanceLog, top_candidates: list[dict]) -> str:
    """
    Generate a plain-text summary from the provenance log.
    Inserted into the Markdown report as the Executive Summary section.
    In a future version this will be written by an LLM.
    """
    steps    = prov.data.get("pipeline_steps", {})
    dataset  = prov.data.get("dataset", {})
    deg      = steps.get("deg_analysis", {})
    sig      = steps.get("signature_building", {})
    drugs    = steps.get("drug_repurposing", {})

    lines = [
        f"This analysis identified {drugs.get('n_candidates', 0)} drug repurposing "
        f"candidates for {prov.disease_name} using transcriptomic data from "
        f"{dataset.get('accession', 'GEO')} ({dataset.get('title', '')}).",

        f"Differential expression analysis compared "
        f"{deg.get('n_disease_samples', '?')} disease samples against "
        f"{deg.get('n_control_samples', '?')} controls, identifying "
        f"{deg.get('n_significant_up', 0)} upregulated and "
        f"{deg.get('n_significant_down', 0)} downregulated genes. "
        f"The disease signature comprised {sig.get('n_up', 0)} upregulated and "
        f"{sig.get('n_down', 0)} downregulated genes after filtering and symbol mapping.",
    ]

    high_conf = [d for d in top_candidates if d.get("confidence") == "high"]
    if high_conf:
        names = ", ".join(d["drug_name"] for d in high_conf[:5])
        lines.append(
            f"High-confidence candidates (q < 0.05) include: {names}."
        )
    elif top_candidates:
        names = ", ".join(d["drug_name"] for d in top_candidates[:3])
        lines.append(
            f"Top candidates by signature reversal score include: {names}. "
            f"Note that no candidates reached high-confidence (q < 0.05) — "
            f"results should be interpreted with caution."
        )

    if prov.has_warnings():
        lines.append(
            f"{len(prov.data['warnings'])} quality control warning(s) were raised "
            f"during the analysis — see the QC Warnings section below."
        )

    return " ".join(lines)


def _header(title: str):
    print(f"\n── {title} {'─' * max(1, 60 - len(title))}")


def _fatal(message: str):
    print(f"\n  ✗  {message}")
    print("  Pipeline cannot continue.\n")
    sys.exit(1)


# ── Pipeline ──────────────────────────────────────────────────────────────────

def run_pipeline(disease_name: str):

    slug        = slugify(disease_name)
    data_dir    = os.path.join("data",    slug)
    reports_dir = os.path.join("reports", slug)

    os.makedirs(data_dir,    exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    prov = ProvenanceLog(disease_name, output_dir=reports_dir)

    print(f"\n{'='*65}")
    print(f"  TDRE — Transcriptomic Drug Repurposing Engine")
    print(f"  Disease: {disease_name}")
    print(f"  Data:    {data_dir}")
    print(f"  Reports: {reports_dir}")
    print(f"{'='*65}")

    # ── Step 1: Disease lookup ────────────────────────────────────────────────
    _header("Step 1: Disease Lookup")
    lookup      = lookup_disease(disease_name, data_dir=data_dir, provenance=prov)
    feasibility = lookup["feasibility"]

    print(f"  Feasibility score: {feasibility['feasibility_score']}/100")
    for flag in feasibility["flags"]:
        print(f"  ⚠  {flag}")

    if not feasibility["proceed_recommended"]:
        print("\n  ⚠  This disease did not meet the recommended feasibility threshold.")
        response = input("  Proceed anyway? (y/n): ").strip().lower()
        if response != "y":
            print("  Exiting.")
            return

    # ── Step 2: GEO dataset search and selection ──────────────────────────────
    _header("Step 2: GEO Dataset Search")
    datasets = search_geo_datasets(disease_name, max_results=50, provenance=prov)
    print(f"  Found {len(datasets)} datasets — scoring suitability...")
    scored   = [score_dataset_suitability(d, provenance=prov) for d in datasets]
    selected = select_dataset_interactive(scored, provenance=prov)
    accession = selected["geo_accession"]

    # ── Step 3: Download and sample group assignment ──────────────────────────
    _header(f"Step 3: Download {accession}")
    download = download_geo_dataset(accession, data_dir=data_dir, provenance=prov)

    if download["count_matrix_path"] is None:
        _fatal("Could not obtain a count matrix from this dataset.")

    if download["issues"]:
        for issue in download["issues"]:
            print(f"  ⚠  {issue}")

    groups = assign_sample_groups_interactive(
        download["metadata_path"], provenance=prov
    )

    if not groups["ready_for_deg"]:
        _fatal(
            f"Insufficient samples for DESeq2 — "
            f"disease: {groups['group_counts'].get('disease', 0)}, "
            f"control: {groups['group_counts'].get('control', 0)}. "
            f"Need at least 3 per group."
        )

    # Update provenance dataset record with actual group counts
    prov.record_dataset(
        accession   = accession,
        title       = download["title"],
        n_samples   = download["n_samples"],
        n_disease   = groups["group_counts"].get("disease", 0),
        n_control   = groups["group_counts"].get("control", 0),
        organism    = selected.get("organism", ""),
        data_source = download["data_source"]
    )

    # ── Step 4: Differential expression ──────────────────────────────────────
    _header("Step 4: Differential Expression Analysis")
    # Merge title-keyed and GSM-keyed dicts so run_deseq2 can match
    # count matrix columns regardless of whether they are sample titles or
    # GSM accession IDs (varies by dataset).
    merged_groups = {**groups["sample_groups"], **groups["sample_groups_by_gsm"]}
    deg = run_deseq2(
        count_matrix_path = download["count_matrix_path"],
        sample_groups     = merged_groups,
        output_prefix     = accession,
        data_dir          = data_dir,
        provenance        = prov
    )

    if "error" in deg:
        _fatal(f"DESeq2 failed: {deg['error']}")

    print(f"  Genes tested:        {deg['n_genes_tested']}")
    print(f"  Significant up:      {deg['n_significant_up']}")
    print(f"  Significant down:    {deg['n_significant_down']}")

    if deg["n_significant_up"] + deg["n_significant_down"] == 0:
        _fatal("No significant DEGs found. Cannot build a transcriptomic signature.")

    # ── Step 5: Signature building ────────────────────────────────────────────
    _header("Step 5: Signature Building")
    signature = build_signature(deg["results_path"], provenance=prov)

    print(f"  Signature quality:   {signature['signature_quality']}")
    print(f"  Upregulated genes:   {signature['n_up']}")
    print(f"  Downregulated genes: {signature['n_down']}")

    if signature["n_up"] == 0 and signature["n_down"] == 0:
        _fatal("Empty signature — cannot query drug database.")

    # ── Step 6: Drug repurposing ──────────────────────────────────────────────
    _header("Step 6: Drug Repurposing")
    drug_result = find_drug_candidates(
        signature     = signature,
        output_prefix = f"{accession}_drug_candidates",
        data_dir      = data_dir,
        provenance    = prov
    )

    print(f"  Candidates found:    {drug_result['n_candidates']}")
    if drug_result.get("result_url"):
        print(f"  Interactive results: {drug_result['result_url']}")

    # ── Step 7: Report ────────────────────────────────────────────────────────
    _header("Step 7: Generating Report")
    narrative   = make_narrative(prov, drug_result["top_candidates"])
    report_path = prov.generate_report(drug_result["top_candidates"], narrative)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print(f"  Pipeline complete.")
    print(f"  Report:   {report_path}")
    print(f"  Data:     {data_dir}")
    print(f"{'='*65}")

    if drug_result["top_candidates"]:
        print(f"\n  Top drug candidates:\n")
        print(f"  {'Drug':<32} {'Score':>8}  {'Q-value':>8}  {'Conf':<8}  MoA")
        print(f"  {'-'*75}")
        for drug in drug_result["top_candidates"]:
            name = drug["drug_name"]
            if len(name) > 32:
                name = name[:29] + "..."
            print(
                f"  {name:<32} "
                f"{str(drug.get('score', '')):>8}  "
                f"{str(round(drug.get('qvalue') or 1.0, 4)):>8}  "
                f"{drug.get('confidence', ''):<8}  "
                f"{drug.get('moa', '')}"
            )
        print()


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TDRE — Transcriptomic Drug Repurposing Engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  python main.py "Friedreich ataxia"'
    )
    parser.add_argument(
        "disease",
        type=str,
        help='Name of the rare/orphan disease to analyse'
    )
    args = parser.parse_args()
    run_pipeline(args.disease)
