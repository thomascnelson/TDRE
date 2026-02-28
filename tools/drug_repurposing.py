# tools/drug_repurposing.py
import requests
import pandas as pd
import json
import os
import re
import time
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.provenance import ProvenanceLog

L1000FWD_URL = "https://maayanlab.cloud/l1000fwd/"

RESOURCES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources")

# Module-level cache — loaded once on first call
_hub_cache = None

def _load_repurposing_hub() -> dict:
    """
    Load the Broad Drug Repurposing Hub annotation files from resources/.
    Returns a dict with:
      brd_to_name  : broad_id → pert_iname
      name_to_info : pert_iname → {moa, target, clinical_phase}
    Results are cached after the first load.
    """
    global _hub_cache
    if _hub_cache is not None:
        return _hub_cache

    brd_to_name  = {}
    name_to_info = {}

    # Find annotation files by prefix — tolerates different date suffixes
    sample_file = None
    drug_file   = None
    if os.path.isdir(RESOURCES_DIR):
        for fname in os.listdir(RESOURCES_DIR):
            if fname.startswith("repo-sample-annotation"):
                sample_file = os.path.join(RESOURCES_DIR, fname)
            elif fname.startswith("repo-drug-annotation"):
                drug_file = os.path.join(RESOURCES_DIR, fname)

    # Sample annotation: broad_id → pert_iname
    if sample_file:
        try:
            df = pd.read_csv(sample_file, sep="\t", comment="!", dtype=str).fillna("")
            for _, row in df.iterrows():
                brd  = row.get("broad_id", "").strip()
                name = row.get("pert_iname", "").strip()
                if brd and name:
                    brd_to_name[brd] = name
                    # Also index by the short BRD prefix (BRD-XXXXXXXX)
                    # in case L1000FWD uses a different batch suffix
                    short = "-".join(brd.split("-")[:2])
                    if short not in brd_to_name:
                        brd_to_name[short] = name
            print(f"  Repurposing Hub: loaded {len(brd_to_name)} BRD ID entries")
        except Exception as e:
            print(f"  Warning: could not load sample annotation file: {e}")

    # Drug annotation: pert_iname → moa, target, clinical_phase
    if drug_file:
        try:
            df = pd.read_csv(drug_file, sep="\t", comment="!", dtype=str).fillna("")
            for _, row in df.iterrows():
                name = row.get("pert_iname", "").strip()
                if name:
                    name_to_info[name.lower()] = {
                        "moa":            row.get("moa",            ""),
                        "target":         row.get("target",         ""),
                        "clinical_phase": row.get("clinical_phase", "")
                    }
            print(f"  Repurposing Hub: loaded {len(name_to_info)} drug annotation entries")
        except Exception as e:
            print(f"  Warning: could not load drug annotation file: {e}")

    _hub_cache = {"brd_to_name": brd_to_name, "name_to_info": name_to_info}
    return _hub_cache


def _pubchem_lookup(brd_id: str) -> str:
    """
    Fallback: look up a BRD ID in PubChem.
    Returns the common name, or the BRD ID itself if not found.
    """
    try:
        url  = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                f"{brd_id}/property/IUPACName,Title/JSON")
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            props = resp.json().get("PropertyTable", {}).get("Properties", [])
            if props:
                return props[0].get("Title") or props[0].get("IUPACName") or brd_id
    except Exception:
        pass
    return brd_id


def resolve_drug_names(sig_ids: list[str]) -> dict:
    """
    Resolve BRD compound IDs to drug names and annotations.

    Lookup order:
      1. Broad Drug Repurposing Hub (local files in resources/) — best coverage,
         also provides moa, target, and clinical_phase
      2. PubChem API — fallback for compounds not in the Hub
      3. Raw BRD ID — last resort
    """
    hub = _load_repurposing_hub()
    brd_to_name  = hub["brd_to_name"]
    name_to_info = hub["name_to_info"]

    parsed        = {}
    brd_to_sigids = {}

    for sig_id in sig_ids:
        parts     = sig_id.split(":")
        brd_id    = parts[1] if len(parts) > 1 else ""
        prefix    = parts[0].split("_")
        cell_line = prefix[1] if len(prefix) > 1 else ""
        timepoint = prefix[2] if len(prefix) > 2 else ""
        dose      = parts[2] if len(parts) > 2 else ""

        parsed[sig_id] = {
            "sig_id":         sig_id,
            "brd_id":         brd_id,
            "cell_id":        cell_line,
            "pert_time":      timepoint,
            "pert_dose":      dose,
            "pert_iname":     brd_id,   # overwritten below if resolved
            "moa":            "",
            "target":         "",
            "clinical_phase": ""
        }

        if brd_id:
            brd_to_sigids.setdefault(brd_id, []).append(sig_id)

    unique_brds    = list(brd_to_sigids.keys())
    hub_resolved   = 0
    pubchem_resolved = 0
    unresolved     = 0

    for brd_id in unique_brds:
        # Try exact match, then short prefix match
        short = "-".join(brd_id.split("-")[:2])
        name  = brd_to_name.get(brd_id) or brd_to_name.get(short)

        if name:
            hub_resolved += 1
        else:
            # Fall back to PubChem
            name = _pubchem_lookup(brd_id)
            time.sleep(0.1)   # respect PubChem rate limit
            if not name.startswith("BRD-"):
                pubchem_resolved += 1
            else:
                unresolved += 1

        # Look up rich metadata by drug name
        info = name_to_info.get(name.lower(), {})

        # Apply to all sig_ids sharing this BRD ID
        for sig_id in brd_to_sigids[brd_id]:
            parsed[sig_id]["pert_iname"]     = name
            parsed[sig_id]["moa"]            = info.get("moa",            "")
            parsed[sig_id]["target"]         = info.get("target",         "")
            parsed[sig_id]["clinical_phase"] = info.get("clinical_phase", "")

    print(f"  Resolved {hub_resolved} via Repurposing Hub, "
          f"{pubchem_resolved} via PubChem, "
          f"{unresolved} unresolved out of {len(unique_brds)} unique BRD IDs")

    return parsed

def parse_l1000fwd_results(raw: dict, top_n: int = 50) -> list[dict]:
    """
    Parse L1000FWD topn results into a clean ranked candidate list.
    Resolves BRD compound IDs to drug names via sig_info lookup.
    Uses the 'opposite' list — drugs that reverse the disease signature.
    """
    opposite = raw.get("opposite", [])
    if not opposite:
        print("  No opposite signatures found in response")
        return []

    # Resolve drug metadata from sig_ids
    sig_ids = [entry["sig_id"] for entry in opposite]
    print(f"  Resolving drug names for {len(sig_ids)} signatures...")
    drug_info = resolve_drug_names(sig_ids)

    candidates = []
    for entry in opposite:
        sig_id = entry["sig_id"]
        info   = drug_info.get(sig_id, {})

        # Parse cell line and time from sig_id as reliable fallback
        parts     = sig_id.split(":")
        prefix    = parts[0].split("_")
        cell_line = info.get("cell_id",   prefix[1] if len(prefix) > 1 else "")
        timepoint = info.get("pert_time", prefix[2] if len(prefix) > 2 else "")
        brd_id    = parts[1] if len(parts) > 1 else ""
        dose      = parts[2] if len(parts) > 2 else ""
        drug_name = info.get("pert_iname", brd_id)
        score     = entry.get("combined_scores", 0)
        qvalue    = entry.get("qvals", 1.0)

        candidates.append({
            "drug_name":  drug_name,
            "brd_id":     brd_id,
            "sig_id":     sig_id,
            "score":      round(score, 4),
            "pvalue":     entry.get("pvals",   None),
            "qvalue":     entry.get("qvals",   None),
            "zscore":     entry.get("zscores", None),
            "cell_line":  cell_line,
            "dose":       dose,
            "time":       timepoint,
            "moa":        info.get("moa",    ""),
            "target":     info.get("target", ""),
            "phase":      info.get("phase",  ""),
            "confidence": "high" if qvalue < 0.05
                else "medium" if qvalue < 0.20
                else "low"
        })

    # Deduplicate by drug name, keeping highest score per drug
    seen = {}
    for c in candidates:
        drug = c["drug_name"]
        if drug not in seen or c["score"] < seen[drug]["score"]:
            seen[drug] = c

    return sorted(seen.values(), key=lambda x: x["score"], reverse=False)[:top_n]

def query_l1000fwd(signature: dict) -> dict:
    """
    Query L1000FWD (Ma'ayan Lab, Mount Sinai) with an up/down gene signature.

    Two steps:
      1. POST up/down genes → get a result_id
      2. GET result_id      → get ranked drug signatures
    """
    # Step 1: Submit signature
    print("  Submitting signature to L1000FWD...")
    payload = {
        "up_genes":   signature["upregulated"],
        "down_genes": signature["downregulated"]
    }

    resp = requests.post(
        L1000FWD_URL + "sig_search",
        json=payload,
        timeout=60
    )

    if resp.status_code != 200:
        return {
            "error":  f"L1000FWD submission failed: HTTP {resp.status_code}",
            "detail": resp.text[:300]
        }

    result_id = resp.json().get("result_id")
    if not result_id:
        return {
            "error":    "L1000FWD returned no result_id",
            "response": resp.json()
        }

    print(f"  Got result_id: {result_id}")
    print(f"  Fetching top results...")
    time.sleep(2)

    # Step 2: Retrieve ranked results
    resp2 = requests.get(
        L1000FWD_URL + f"result/topn/{result_id}",
        timeout=60
    )

    if resp2.status_code != 200:
        return {
            "error":  f"L1000FWD result fetch failed: HTTP {resp2.status_code}",
            "detail": resp2.text[:300]
        }

    return {
        "source":     "l1000fwd",
        "result_id":  result_id,
        "result_url": f"https://maayanlab.cloud/l1000fwd/#result/{result_id}",
        "data":       resp2.json()
    }

def query_enrichr_fallback(signature: dict) -> list[dict]:
    """
    Fallback using Enrichr LINCS gene sets if L1000FWD is unavailable.
    Queries up and down gene lists separately to find reversal candidates.
    """
    import gseapy as gp
    print("  Falling back to Enrichr/gseapy...")

    candidates = {}

    try:
        up_enr = gp.enrichr(
            gene_list=signature["upregulated"],
            gene_sets=["LINCS_L1000_Chem_Pert_Consensus_Sigs"],
            organism="human",
            outdir=None,
            verbose=False
        ).results

        for _, row in up_enr.head(50).iterrows():
            drug      = row["Term"].split("_")[0]
            direction = "dn" if "_dn_" in row["Term"].lower() else "up"
            if direction == "dn":
                score = 1 - row["Adjusted P-value"]
                if drug not in candidates or score > candidates[drug]["score"]:
                    candidates[drug] = {
                        "drug_name":       drug,
                        "score":           round(score, 4),
                        "pvalue":          row["P-value"],
                        "qvalue":          row["Adjusted P-value"],
                        "overlap_genes":   row["Genes"],
                        "confidence":      "low",
                        "source":          "enrichr_fallback"
                    }
    except Exception as e:
        print(f"  Enrichr up query failed: {e}")

    try:
        dn_enr = gp.enrichr(
            gene_list=signature["downregulated"],
            gene_sets=["LINCS_L1000_Chem_Pert_Consensus_Sigs"],
            organism="human",
            outdir=None,
            verbose=False
        ).results

        for _, row in dn_enr.head(50).iterrows():
            drug      = row["Term"].split("_")[0]
            direction = "up" if "_up_" in row["Term"].lower() else "dn"
            if direction == "up":
                score = 1 - row["Adjusted P-value"]
                if drug not in candidates or score > candidates[drug]["score"]:
                    candidates[drug] = {
                        "drug_name":       drug,
                        "score":           round(score, 4),
                        "pvalue":          row["P-value"],
                        "qvalue":          row["Adjusted P-value"],
                        "overlap_genes":   row["Genes"],
                        "confidence":      "low",
                        "source":          "enrichr_fallback"
                    }
    except Exception as e:
        print(f"  Enrichr down query failed: {e}")

    return sorted(candidates.values(), key=lambda x: x["score"], reverse=True)[:20]

def find_drug_candidates(signature: dict, output_prefix: str = "drug_candidates",
                         data_dir: str = "data", provenance: ProvenanceLog = None) -> dict:
    """
    Main entry point for Tool 6.
    Queries L1000FWD, falls back to Enrichr if needed,
    ranks candidates and saves results to file.
    """
    os.makedirs(data_dir, exist_ok=True)

    raw = query_l1000fwd(signature)

    if "error" in raw:
        print(f"  L1000FWD unavailable: {raw['error']}")
        ranked      = query_enrichr_fallback(signature)
        source_used = "enrichr_fallback"
        result_url  = None
    else:
        ranked      = parse_l1000fwd_results(raw["data"])
        source_used = "l1000fwd"
        result_url  = raw.get("result_url")
        print(f"  Interactive results: {result_url}")

    # Save results
    json_path = os.path.join(data_dir, f"{output_prefix}.json")
    csv_path  = os.path.join(data_dir, f"{output_prefix}.csv")

    with open(json_path, "w") as f:
        json.dump(ranked, f, indent=2)

    if ranked:
        pd.DataFrame(ranked).to_csv(csv_path, index=False)

    print(f"  Results saved to: {json_path}")
    print(f"  Results saved to: {csv_path}")

    result = {
        "source_used":    source_used,
        "result_url":     result_url,
        "n_candidates":   len(ranked),
        "top_candidates": ranked[:10],
        "json_path":      json_path,
        "csv_path":       csv_path
    }


    if provenance:
        n_high = len([d for d in ranked if d.get("confidence") == "high"])
        n_med  = len([d for d in ranked if d.get("confidence") == "medium"])
        provenance.record_step("drug_repurposing", {
            "source":             source_used,
            "n_candidates":       len(ranked),
            "n_high_confidence":  n_high,
            "n_med_confidence":   n_med,
            "result_url":         result.get("result_url", "")
        })
        provenance.record_output_file("drug_candidates_json", json_path)
        provenance.record_output_file("drug_candidates_csv",  csv_path)
        if source_used == "enrichr_fallback":
            provenance.record_warning("drug_repurposing",
                "L1000FWD unavailable — results from Enrichr fallback. "
                "Confidence scores are less reliable than L1000FWD connectivity scores.")
        if n_high == 0:
            provenance.record_warning("drug_repurposing",
                "No high-confidence drug candidates found (q < 0.05). "
                "Results should be interpreted with caution.")

    return result

### TEST ROUTINE

if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()

    SIGNATURE_FILE = "data/GSE226646_signature.json"
    
    # Create the provenance log
    prov = ProvenanceLog.get_or_create("Friedreich ataxia")

    print(f"Loading signature from {SIGNATURE_FILE}...")
    with open(SIGNATURE_FILE) as f:
        signature = json.load(f)

    print(f"Signature: {signature['n_up']} up, {signature['n_down']} down genes")
    print(f"Quality:   {signature['signature_quality']}\n")

    result = find_drug_candidates(signature, output_prefix="GSE226646_drug_candidates", provenance = prov)

    print(f"\nSource used:  {result['source_used']}")
    print(f"Candidates:   {result['n_candidates']}")
    if result.get("result_url"):
        print(f"View online:  {result['result_url']}")

    print(f"\nTop drug candidates (signature reversal):")
    print(f"{'Drug':<30} {'Score':>8}  {'Confidence':<10}  {'Cell':>8}  MoA")
    print("-" * 80)
    for drug in result["top_candidates"]:
        print(
            f"{drug['drug_name']:<30} "
            f"{str(drug.get('score', 'N/A')):>8}  "
            f"{drug.get('confidence',''):.<10}  "
            f"{drug.get('cell_line',''):>8}  "
            f"{drug.get('moa', '')}"
        )

    print(f"\nFull results: {result['csv_path']}")


