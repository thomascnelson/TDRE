# tools/geo_search.py
from Bio import Entrez
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.provenance import ProvenanceLog

def search_geo_datasets(disease_name: str, max_results: int = 50, provenance: ProvenanceLog = None) -> list[dict]:
    """
    Search GEO for RNA-seq datasets relevant to the disease.
    Returns structured metadata for each candidate dataset.
    """
    Entrez.email = os.getenv("EMAIL_ADDRESS")
    Entrez.api_key = os.getenv("NCBI_API_KEY")

    # GEO DataSets database
    search_term = (
        f'("{disease_name}"[title]) OR '
        f'("{disease_name}"[summary])'
        #f'("Homo sapiens"[organism])'
        #f'("expression profiling by high throughput sequencing"[DataSet Type])'
    )

    handle = Entrez.esearch(db="gds", term=search_term, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    datasets = []
    for gds_id in record["IdList"]:
        handle = Entrez.esummary(db="gds", id=gds_id)
        summary = Entrez.read(handle)
        handle.close()

        if summary:
            s = summary[0]
            datasets.append({
                "geo_accession": s.get("Accession", ""),
                "title": s.get("title", ""),
                "summary": s.get("summary", "")[:500],
                "organism": s.get("taxon", ""),
                "n_samples": int(s.get("n_samples", 0)),
                "platform": s.get("GPL", ""),
                "submission_date": s.get("PDAT", ""),
                "gds_id": gds_id
            })

    if provenance:
        provenance.record_step("geo_search", {
            "disease_name":    disease_name,
            "n_datasets_found": len(datasets)
        })

    return datasets

def score_dataset_suitability(dataset: dict, provenance: ProvenanceLog = None) -> dict:
    """
    Heuristic scoring of a GEO dataset for DEG analysis suitability.
    The LLM will re-evaluate this with reasoning, but hard rules
    filter out obviously unsuitable datasets first.
    """
    score = 0
    issues = []
    positives = []

    # Sample size: need at least 3 per group (6 total) for DESeq2
    n = dataset["n_samples"]
    if n >= 20:
        score += 30
        positives.append(f"Good sample size ({n} samples)")
    elif n >= 10:
        score += 20
        positives.append(f"Adequate sample size ({n} samples)")
    elif n >= 6:
        score += 10
        issues.append(f"Borderline sample size ({n} samples) — DESeq2 results may be underpowered")
    else:
        score -= 30
        issues.append(f"Insufficient sample size ({n} samples) — need at least 6")

    # Human data preferred
    if "Homo sapiens" in dataset.get("organism", ""):
        score += 30
        positives.append("Human data")
    elif "Mus musculus" in dataset.get("organism", ""):
        score += 10
        issues.append("Mouse data — gene ortholog mapping will be required")

    # Check summary for case/control language
    summary_lower = dataset["summary"].lower()
    control_keywords = ["control", "healthy", "normal", "wild-type", "versus"]
    if any(k in summary_lower for k in control_keywords):
        score += 20
        positives.append("Appears to include control group")
    else:
        issues.append("No clear mention of control group in summary")

    if provenance:
        recommended = score >= 40
        if not recommended:
            provenance.record_rejected_dataset(
                accession = dataset.get("geo_accession", ""),
                title     = dataset.get("title", ""),
                reason    = "; ".join(issues) if issues else "Low suitability score"
            )
        for issue in issues:
            provenance.record_warning("geo_search",
                f"{dataset.get('geo_accession', '')}: {issue}")

    return {
        **dataset,
        "suitability_score": score,
        "issues": issues,
        "positives": positives,
        "recommended": score >= 40
    }


def select_dataset_interactive(scored_datasets: list[dict],
                               provenance: ProvenanceLog = None) -> dict:
    """
    Present scored datasets to the user and ask them to select one for analysis.
    Shows all metadata for each dataset so the user can make an informed choice.
    Returns the selected dataset dict (including geo_accession for download).
    """
    if not scored_datasets:
        raise ValueError("No datasets to select from")

    # Sort by suitability score descending, show only recommended, cap at 5
    ranked = sorted(scored_datasets, key=lambda d: d["suitability_score"], reverse=True)
    ranked = [d for d in ranked if d["recommended"]][:5]

    if not ranked:
        print("  ⚠  No datasets met the recommendation threshold. Showing top 5 overall.")
        ranked = sorted(scored_datasets, key=lambda d: d["suitability_score"], reverse=True)[:5]

    print("\n" + "="*70)
    print("  GEO DATASET SELECTION")
    print("="*70)
    print("  Review each dataset and enter the number of the one to analyse.")
    print("="*70)

    for i, ds in enumerate(ranked, 1):
        print(f"\n  [{i}]  {ds['geo_accession']}  —  Score: {ds['suitability_score']}"
              f"  {'✓ Recommended' if ds['recommended'] else '⚠  Not recommended'}")
        print(f"  Title:    {ds['title']}")
        print(f"  Organism: {ds['organism']}")
        print(f"  Samples:  {ds['n_samples']}")
        print(f"  Platform: {ds['platform']}")
        print(f"  Submitted:{ds['submission_date']}")
        print(f"  URL:      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={ds['geo_accession']}")
        print(f"  Summary:  {ds['summary']}")
        if ds["positives"]:
            for p in ds["positives"]:
                print(f"    ✓  {p}")
        if ds["issues"]:
            for iss in ds["issues"]:
                print(f"    ⚠  {iss}")
        print()

    print("="*70)

    while True:
        raw = input(
            f"  Select dataset [1-{len(ranked)}] or enter a GEO accession directly: "
        ).strip()

        if raw.upper().startswith("GSE") or raw.upper().startswith("GDS"):
            # User typed an accession directly — build a minimal dict
            selected = next((d for d in ranked if d["geo_accession"] == raw.upper()), None)
            if selected is None:
                print(f"  Accession {raw.upper()} not in the list above — using it anyway.")
                selected = {"geo_accession": raw.upper(), "title": "", "suitability_score": 0,
                            "n_samples": 0, "organism": "", "platform": "",
                            "submission_date": "", "summary": "",
                            "issues": [], "positives": [], "recommended": False}
            break

        try:
            idx = int(raw) - 1
            if 0 <= idx < len(ranked):
                selected = ranked[idx]
                break
            else:
                print(f"  Please enter a number between 1 and {len(ranked)}.")
        except ValueError:
            print(f"  Please enter a number between 1 and {len(ranked)}, or a GEO accession.")

    print(f"\n  Selected: {selected['geo_accession']}  —  {selected['title']}")
    print("="*70 + "\n")

    if provenance:
        provenance.record_step("dataset_selection", {
            "selected_accession": selected["geo_accession"],
            "selected_title":     selected["title"],
            "suitability_score":  selected["suitability_score"],
            "n_candidates_shown": len(ranked)
        })
        provenance.record_dataset(
            accession   = selected["geo_accession"],
            title       = selected["title"],
            n_samples   = selected["n_samples"],
            n_disease   = 0,   # not yet known — updated after sample assignment
            n_control   = 0,
            organism    = selected["organism"],
            data_source = "pending"
        )

    return selected


### TEST ROUTINE
if __name__ == "__main__":
    import json
    from dotenv import load_dotenv
    load_dotenv()

    # Create the provenance log
    prov = ProvenanceLog.get_or_create("Friedreich ataxia")

    geo_search_results = search_geo_datasets("Friedreich ataxia", max_results=20, provenance=prov)
    scored_results = [score_dataset_suitability(d, provenance=prov) for d in geo_search_results]

    selected = select_dataset_interactive(scored_results, provenance=prov)
    print(f"Proceeding with: {selected['geo_accession']}")



