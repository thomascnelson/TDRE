"""
run_pipeline.py — End-to-end drug repurposing pipeline orchestrator.

Usage:
    python run_pipeline.py
    python run_pipeline.py --accession GSE297335 --disease "Barth Syndrome"

The pipeline runs four tools in sequence:
    1. GEO download        → count matrix + sample metadata
    2. DESeq2 DEG          → differential gene expression results
    3. Signature building  → disease transcriptomic signature
    4. Drug repurposing    → L1000FWD drug candidates

All outputs land in data/{accession}/ and a full provenance log + human-
readable Markdown report are written to reports/.
"""

import argparse
import importlib
import json
import os
import sys

from dotenv import load_dotenv

load_dotenv()

# Ensure project root is on the path so tool imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Tool filenames start with digits, which aren't valid Python identifiers.
# Use importlib to load them by their actual on-disk names.
_tools_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools")
sys.path.insert(0, _tools_dir)

_geo      = importlib.import_module("1_geo_download")
_deg      = importlib.import_module("2_deg_analysis")
_sig      = importlib.import_module("3_signature")
_drug     = importlib.import_module("4_drug_repurposing")

download_geo_dataset             = _geo.download_geo_dataset
assign_sample_groups_interactive = _geo.assign_sample_groups_interactive
run_deseq2                       = _deg.run_deseq2
build_signature                  = _sig.build_signature
find_drug_candidates             = _drug.find_drug_candidates

from utils.provenance import ProvenanceLog


# ── Helpers ───────────────────────────────────────────────────────────────────

def _prompt(label: str, default: str = "") -> str:
    """Prompt the user for input, showing a default value."""
    suffix = f" [{default}]" if default else ""
    raw = input(f"  {label}{suffix}: ").strip()
    return raw if raw else default


def _banner(title: str):
    width = max(65, len(title) + 8)
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


def _step_banner(n: int, title: str):
    print(f"\n{'─'*65}")
    print(f"  STEP {n}: {title}")
    print(f"{'─'*65}")


# ── Main pipeline ─────────────────────────────────────────────────────────────

def run_pipeline(accession: str, disease_name: str):
    """
    Execute the full drug repurposing pipeline for the given GEO accession
    and disease name.
    """

    data_dir = os.path.join("data", accession)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs("reports", exist_ok=True)

    _banner(f"Drug Repurposing Pipeline — {disease_name}  ({accession})")
    print(f"  Data directory : {data_dir}/")
    print(f"  Reports        : reports/")

    # ── Provenance log ────────────────────────────────────────────────────────
    prov = ProvenanceLog.get_or_create(disease_name)

    # ── Step 1: GEO Download ─────────────────────────────────────────────────
    _step_banner(1, "GEO Download")
    try:
        result = download_geo_dataset(accession, data_dir=data_dir, provenance=prov)
    except Exception as exc:
        print(f"\n  ERROR in GEO download: {exc}")
        prov.record_warning("geo_download", f"Fatal error: {exc}")
        sys.exit(1)

    if result.get("count_matrix_path") is None:
        print("\n  Cannot continue — no count matrix obtained.")
        for issue in result.get("issues", []):
            print(f"  {issue}")
        sys.exit(1)

    count_matrix_path = result["count_matrix_path"]
    metadata_path     = result["metadata_path"]

    print(f"\n  Title   : {result['title']}")
    print(f"  Samples : {result['n_samples']}")
    print(f"  Source  : {result['data_source']}")
    if result.get("issues"):
        for issue in result["issues"]:
            print(f"  WARNING : {issue}")

    # ── Step 1b: Sample Group Assignment ─────────────────────────────────────
    _step_banner("1b", "Sample Group Assignment")
    try:
        groups_result = assign_sample_groups_interactive(
            metadata_path, provenance=prov
        )
    except Exception as exc:
        print(f"\n  ERROR in sample assignment: {exc}")
        prov.record_warning("sample_assignment", f"Fatal error: {exc}")
        sys.exit(1)

    if not groups_result.get("ready_for_deg"):
        counts = groups_result.get("group_counts", {})
        print(f"\n  Not enough samples for DESeq2:")
        print(f"    disease: {counts.get('disease', 0)}  "
              f"control: {counts.get('control', 0)}  (need ≥3 each)")
        print("  Pipeline cannot continue.")
        sys.exit(1)

    groups_path = groups_result["groups_file"]
    counts      = groups_result["group_counts"]

    # Update provenance with final sample counts
    prov.update_dataset_counts(
        n_disease = counts.get("disease", 0),
        n_control = counts.get("control", 0)
    )

    # ── Step 2: Differential Gene Expression ─────────────────────────────────
    _step_banner(2, "Differential Gene Expression (DESeq2)")
    try:
        deg_result = run_deseq2(
            count_matrix_path = count_matrix_path,
            groups_path       = groups_path,
            output_prefix     = accession,
            data_dir          = data_dir,
            provenance        = prov
        )
    except Exception as exc:
        print(f"\n  ERROR in DESeq2: {exc}")
        prov.record_warning("deg_analysis", f"Fatal error: {exc}")
        sys.exit(1)

    if "error" in deg_result:
        print(f"\n  DESeq2 failed: {deg_result['error']}")
        sys.exit(1)

    print(f"\n  Genes tested       : {deg_result['n_genes_tested']}")
    print(f"  Significant up     : {deg_result['n_significant_up']}")
    print(f"  Significant down   : {deg_result['n_significant_down']}")

    if deg_result["n_significant_up"] + deg_result["n_significant_down"] == 0:
        print("\n  No significant DEGs found. Cannot build a meaningful signature.")
        print("  Consider relaxing thresholds (edit run_pipeline.py padj/lfc params).")
        sys.exit(1)

    deg_results_path = deg_result["results_path"]

    # ── Step 3: Disease Signature ─────────────────────────────────────────────
    _step_banner(3, "Disease Transcriptomic Signature")
    try:
        signature = build_signature(
            deg_results_path = deg_results_path,
            provenance       = prov
        )
    except Exception as exc:
        print(f"\n  ERROR building signature: {exc}")
        prov.record_warning("signature_building", f"Fatal error: {exc}")
        sys.exit(1)

    print(f"\n  Quality            : {signature['signature_quality']}")
    print(f"  Upregulated genes  : {signature['n_up']}")
    print(f"  Downregulated genes: {signature['n_down']}")

    if signature["n_up"] + signature["n_down"] == 0:
        print("\n  Signature is empty — cannot query L1000FWD.")
        sys.exit(1)

    # ── Step 4: Drug Repurposing ──────────────────────────────────────────────
    _step_banner(4, "Drug Repurposing (L1000FWD)")
    try:
        drug_result = find_drug_candidates(
            signature     = signature,
            output_prefix = f"{accession}_drug_candidates",
            data_dir      = data_dir,
            provenance    = prov
        )
    except Exception as exc:
        print(f"\n  ERROR in drug repurposing: {exc}")
        prov.record_warning("drug_repurposing", f"Fatal error: {exc}")
        # Don't exit — we can still generate a partial report
        drug_result = {"n_candidates": 0, "top_candidates": [], "source_used": "error"}

    if drug_result.get("result_url"):
        print(f"\n  Interactive results: {drug_result['result_url']}")
    print(f"  Candidates found   : {drug_result.get('n_candidates', 0)}")

    # ── Final Report ──────────────────────────────────────────────────────────
    _step_banner("✓", "Generating Report")

    # Load full drug candidate list from saved JSON for the report
    drug_candidates = drug_result.get("top_candidates", [])
    drug_json_path  = drug_result.get("json_path")
    if drug_json_path and os.path.exists(drug_json_path):
        try:
            with open(drug_json_path) as f:
                drug_candidates = json.load(f)
        except Exception:
            pass  # fall back to the in-memory top_candidates

    report_path = prov.generate_report(drug_candidates)

    # ── Summary ───────────────────────────────────────────────────────────────
    _banner("Pipeline Complete")
    print(f"  Provenance log : {prov.log_path}")
    print(f"  Report         : {report_path}")
    print(f"  Data directory : {data_dir}/")
    print()
    print("  Key results:")
    print(f"    DEGs identified : {deg_result['n_significant_up']} up, "
          f"{deg_result['n_significant_down']} down")
    print(f"    Signature genes : {signature['n_up']} up, "
          f"{signature['n_down']} down ({signature['signature_quality']})")
    print(f"    Drug candidates : {drug_result.get('n_candidates', 0)} "
          f"(source: {drug_result.get('source_used', 'unknown')})")
    if drug_result.get("result_url"):
        print(f"    Online results  : {drug_result['result_url']}")
    print()


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Drug repurposing pipeline: GEO → DESeq2 → signature → L1000FWD"
    )
    parser.add_argument("--accession", "-a",
                        help="GEO accession (e.g. GSE297335)")
    parser.add_argument("--disease",   "-d",
                        help="Disease name for provenance log (e.g. 'Barth Syndrome')")
    args = parser.parse_args()

    accession    = args.accession
    disease_name = args.disease

    # Interactive prompts if not provided on command line
    if not accession or not disease_name:
        print("\n  Drug Repurposing Pipeline")
        print("  ─────────────────────────")
        if not accession:
            accession = _prompt("GEO accession (e.g. GSE297335)")
        if not disease_name:
            disease_name = _prompt("Disease name (e.g. Barth Syndrome)")

    if not accession or not disease_name:
        print("  ERROR: accession and disease name are required.")
        sys.exit(1)

    accession = accession.strip().upper()

    run_pipeline(accession, disease_name)


if __name__ == "__main__":
    main()
