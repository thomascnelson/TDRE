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
    Entrez.email = "EMAIL_ADDRESS"
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
        if not result["recommended"]:
            provenance.record_rejected_dataset(
                accession = result.get("geo_accession", ""),
                title     = result.get("title", ""),
                reason    = "; ".join(result["issues"]) if result["issues"]
                            else "Low suitability score"
            )
        for issue in result["issues"]:
            provenance.record_warning("geo_search",
                f"{result.get('geo_accession', '')}: {issue}")

    return {
        **dataset,
        "suitability_score": score,
        "issues": issues,
        "positives": positives,
        "recommended": score >= 40
    }


### TEST ROUTINE
if __name__ == "__main__":
    import json
    from dotenv import load_dotenv
    load_dotenv()

    # Create the provenance log
    prov = ProvenanceLog.get_or_create("Friedreich ataxia")

    geo_search_results = search_geo_datasets("Friedreich ataxia", max_results=20, provenance = prov)
    scored_results = []
    
    for item in geo_search_results:
        scored_results.append(score_dataset_suitability(item))
    sorted_scored_results = sorted(scored_results, key=lambda d: d['suitability_score'], reverse=True)
    print(json.dumps(sorted_scored_results, indent=4))



