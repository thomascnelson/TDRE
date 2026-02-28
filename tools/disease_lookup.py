# tools/disease_lookup.py
import requests
from Bio import Entrez
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.provenance import ProvenanceLog

def lookup_disease(disease_name: str, data_dir: str = "data",
                   provenance: ProvenanceLog = None) -> dict:
    """
    Query NCBI MeSH + Orphanet to assess the disease and return
    structured information about it.
    """
    Entrez.email = os.getenv("EMAIL_ADDRESS")
    Entrez.api_key = os.getenv("NCBI_API_KEY")
    orphanet_api_key = os.getenv("ORPHANET_API_KEY", "")

    result = {
        "disease_name": disease_name,
        "mesh_terms": [],
        "orphanet_info": {},
        "existing_treatments": [],
        "assessment": {}
    }

    # Search MeSH for canonical disease name and synonyms
    handle = Entrez.esearch(db="mesh", term=disease_name)
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        mesh_id = record["IdList"][0]
        handle = Entrez.efetch(db="mesh", id=mesh_id, rettype="full")
        mesh_data = handle.read()
        handle.close()
        result["mesh_terms"] = record["IdList"]
        result["mesh_raw"] = mesh_data[:2000]  # truncate for LLM context

    # Query Orphanet REST API for rare disease info
    print(disease_name)
    orphanet_url = f"https://api.orphacode.org/EN/ClinicalEntity/FindbyName/{requests.utils.quote(disease_name)}"
    try:
        resp = requests.get(orphanet_url, timeout=10,
                       headers={"apiKey": orphanet_api_key})
        if resp.status_code == 200:
            result["orphanet_info"] = resp.json()
        else:
            result["orphanet_info"] = {
                "error": f"HTTP {resp.status_code}",
                "detail": resp.text[:200]
        }
    except Exception as e:
        result["orphanet_info"] = {"error": str(e)}

    # Search PubMed for treatment papers to gauge existing drug landscape
    handle = Entrez.esearch(
        db="pubmed",
        term=f"{disease_name}[Title/Abstract] AND treatment[Title/Abstract]",
        retmax=5,
        sort="relevance"
    )
    pub_record = Entrez.read(handle)
    handle.close()
    result["pubmed_treatment_hits"] = int(pub_record["Count"])

    if provenance:
        provenance.record_step("disease_lookup", {
            "mesh_terms":             result.get("mesh_terms", []),
            "pubmed_treatment_hits":  result.get("pubmed_treatment_hits", 0),
            "orphanet_found":         "error" not in result.get("orphanet_info", {})
        })
        provenance.record_output_file("disease_lookup_raw",
                                      os.path.join(data_dir, f"{disease_name.replace(' ', '_')}_lookup.json"))

    result["feasibility"] = assess_repurposing_feasibility(result, provenance=provenance)

    return result

def assess_repurposing_feasibility(lookup_result: dict, provenance: ProvenanceLog = None) -> dict:
    """
    Rule-based pre-assessment. The LLM will do deeper reasoning,
    but this flags obvious blockers early.
    """
    flags = []
    score = 100  # start optimistic

    if lookup_result["pubmed_treatment_hits"] == 0:
        flags.append("No published treatment literature found — disease may be too poorly characterized")
        score -= 30

    if not lookup_result["mesh_terms"]:
        flags.append("No MeSH term found — disease name may be non-standard, try a synonym")
        score -= 20

    if not lookup_result.get("orphanet_info") or \
       "error" in lookup_result.get("orphanet_info", {}):
        flags.append("Not found in Orphanet — may not be a recognized rare disease")
        score -= 20

    # Store in a variable before recording and returning
    result = {
        "feasibility_score":    max(0, score),
        "flags":                flags,
        "proceed_recommended":  score >= 50
    }

    if provenance:
        provenance.record_step("feasibility_assessment", {
            "feasibility_score":   result["feasibility_score"],
            "flags":               result["flags"],
            "proceed_recommended": result["proceed_recommended"]
        })
        for flag in result["flags"]:
            provenance.record_warning("disease_lookup", flag)

    return result

#### TEST ROUTINE
if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()

    # Create the provenance log
    prov = ProvenanceLog.get_or_create("Friedreich ataxia")
    
    result = lookup_disease("Friedreich ataxia", provenance=prov)
    print(result)
    print(result["feasibility"])

