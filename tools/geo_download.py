# tools/geo_download.py
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import re
import urllib.request

import GEOparse
import pandas as pd

from utils.provenance import ProvenanceLog

DATA_DIR = "data"


# ── Main download function ────────────────────────────────────────────────────

def download_geo_dataset(geo_accession: str,
                         provenance: ProvenanceLog = None) -> dict:
    """
    Download a GEO dataset and extract a count matrix + metadata.
    Returns paths to saved files and a summary of what was found.

    Strategy:
      1. Look for a supplementary count file attached to the GEO record
         (preferred — these are raw integer counts suitable for DESeq2)
      2. Fall back to the GEO matrix file if no supplementary file found
         (often normalized — flagged as a warning)
      3. Return a clear error if neither method works
    """
    os.makedirs(DATA_DIR, exist_ok=True)

    result = {
        "accession":         geo_accession,
        "title":             "",
        "n_samples":         0,
        "platforms":         [],
        "samples":           [],
        "count_matrix_path": None,
        "metadata_path":     None,
        "data_source":       "none",
        "issues":            []
    }

    # Download the GEO soft file
    print(f"  Fetching GEO record {geo_accession}...")
    try:
        gse = GEOparse.get_GEO(
            geo=geo_accession,
            destdir=DATA_DIR,
            silent=True
        )
    except Exception as e:
        result["issues"].append(f"FATAL: Could not fetch GEO record: {e}")
        return result

    result["title"]     = gse.metadata.get("title", [""])[0]
    result["n_samples"] = len(gse.gsms)
    result["platforms"] = list(gse.gpls.keys())

    # ── Extract sample metadata ───────────────────────────────────────────────
    sample_meta = []
    for gsm_name, gsm in gse.gsms.items():
        meta = {
            "sample_id": gsm_name,
            "title":     gsm.metadata.get("title", [""])[0],
            "source":    gsm.metadata.get("source_name_ch1", [""])[0],
            "characteristics": "; ".join([
                item
                for sublist in gsm.metadata.get("characteristics_ch1", [])
                for item in (sublist if isinstance(sublist, list) else [sublist])
            ]),
            "treatment": gsm.metadata.get("treatment_protocol_ch1", [""])[0][:200]
        }
        sample_meta.append(meta)

    meta_df   = pd.DataFrame(sample_meta)
    meta_path = os.path.join(DATA_DIR, f"{geo_accession}_metadata.csv")
    meta_df.to_csv(meta_path, index=False)
    result["metadata_path"] = meta_path
    result["samples"]       = sample_meta

    # ── Try supplementary count file first ───────────────────────────────────
    supp_files = gse.metadata.get("supplementary_file", [])
    count_matrix_path = _try_get_supplementary_counts(geo_accession, supp_files)

    if count_matrix_path:
        result["count_matrix_path"] = count_matrix_path
        result["data_source"]       = "supplementary_counts"
    else:
        # Fall back to GEO matrix file
        matrix_path = _extract_matrix_from_gse(gse, geo_accession)

        if matrix_path is None:
            result["count_matrix_path"] = None
            result["data_source"]       = "none"
            result["issues"].append(
                "FATAL: Could not obtain count matrix from either supplementary "
                "files or GEO matrix file. Cannot proceed with this dataset."
            )
            if provenance:
                provenance.record_warning(
                    "geo_download",
                    f"{geo_accession}: Count matrix could not be obtained — "
                    "pipeline cannot continue with this dataset"
                )
            return result

        result["count_matrix_path"] = matrix_path
        result["data_source"]       = "geo_matrix_file"
        result["issues"].append(
            "No raw count supplementary file found. Using GEO matrix file "
            "which may contain normalized values rather than raw counts. "
            "DESeq2 requires raw integer counts — verify before running DEG."
        )

    # ── Record to provenance ──────────────────────────────────────────────────
    if provenance:
        provenance.record_output_file("count_matrix",
                                      result["count_matrix_path"])
        provenance.record_output_file("sample_metadata",
                                      result["metadata_path"])
        for issue in result["issues"]:
            provenance.record_warning("geo_download", issue)

    return result


# ── Supplementary file download ───────────────────────────────────────────────

def _try_get_supplementary_counts(accession: str,
                                   supp_files: list) -> str | None:
    """
    Look for a raw count matrix in the GEO supplementary files.
    Downloads the file, parses it, and saves a clean CSV.
    Returns the path to the clean CSV, or None if nothing suitable found.
    """
    count_keywords = [
        "count", "raw", "featurecount", "htseq", "rsem",
        "kallisto", "salmon", "star", "expression",
        "read", "matrix", "rnaseq", "rna-seq"
    ]

    print(f"  Supplementary files available:")
    for url in supp_files:
        print(f"    {url}")

    for url in supp_files:
        url_lower = url.lower()
        if not any(k in url_lower for k in count_keywords):
            print(f"  Skipping (no keyword match): {url.split('/')[-1]}")
            continue

        filename = url.split("/")[-1]
        dest     = os.path.join(DATA_DIR, filename)

        # Download if not already cached
        if not os.path.exists(dest):
            print(f"  Downloading: {filename}...")
            try:
                urllib.request.urlretrieve(url, dest)
                print(f"  Downloaded successfully")
            except Exception as e:
                print(f"  Download failed: {e}")
                continue
        else:
            print(f"  Found cached: {filename}")

        # Parse the file — handle quoting issues common in GEO files
        try:
            if dest.endswith(".gz"):
                df = pd.read_csv(dest, sep=",", index_col=0,
                                 compression="gzip", quoting=csv.QUOTE_NONE)
            else:
                df = pd.read_csv(dest, sep=",", index_col=0,
                                 quoting=csv.QUOTE_NONE)

            # Strip stray quote characters from gene ID index
            df.index = df.index.str.strip('"')

            clean_path = os.path.join(DATA_DIR, f"{accession}_counts.csv")
            df.to_csv(clean_path, quoting=csv.QUOTE_NONE)
            print(f"  Count matrix saved: {clean_path} "
                  f"({df.shape[0]} genes x {df.shape[1]} samples)")
            return clean_path

        except Exception as e:
            print(f"  Could not parse {filename}: {e}")
            # Try tab-separated as fallback
            try:
                if dest.endswith(".gz"):
                    df = pd.read_csv(dest, sep="\t", index_col=0,
                                     compression="gzip", quoting=csv.QUOTE_NONE)
                else:
                    df = pd.read_csv(dest, sep="\t", index_col=0,
                                     quoting=csv.QUOTE_NONE)

                df.index = df.index.str.strip('"')
                clean_path = os.path.join(DATA_DIR, f"{accession}_counts.csv")
                df.to_csv(clean_path, quoting=csv.QUOTE_NONE)
                print(f"  Count matrix saved (tab-separated): {clean_path} "
                      f"({df.shape[0]} genes x {df.shape[1]} samples)")
                return clean_path

            except Exception as e2:
                print(f"  Tab-separated parse also failed: {e2}")
                continue

    print("  No supplementary count file found or downloaded successfully")
    return None


# ── GEO matrix fallback ───────────────────────────────────────────────────────

def _extract_matrix_from_gse(gse, accession: str) -> str | None:
    """
    Extract expression matrix from GEO soft file as a last resort.
    Note: RNA-seq data often fails here because GEOparse's pivot_samples
    expects a microarray-style ID_REF column that RNA-seq data doesn't have.
    """
    try:
        pivot = gse.pivot_samples("VALUE")
        path  = os.path.join(DATA_DIR, f"{accession}_matrix.csv")
        pivot.to_csv(path)
        print(f"  GEO matrix extracted: {path}")
        return path

    except KeyError:
        print("  ✗  GEO matrix fallback failed — RNA-seq data does not have")
        print("     the ID_REF column that GEOparse expects.")
        print("     This dataset requires a supplementary count file.")

        # Last resort: try extracting directly from individual sample tables
        try:
            frames = []
            for gsm_name, gsm in gse.gsms.items():
                if gsm.table is not None and not gsm.table.empty:
                    col      = gsm.table.iloc[:, -1]
                    col.name = gsm_name
                    frames.append(col)

            if frames:
                df   = pd.concat(frames, axis=1)
                path = os.path.join(DATA_DIR, f"{accession}_matrix.csv")
                df.to_csv(path)
                print(f"  Sample table matrix extracted: {path}")
                return path

        except Exception as e:
            print(f"  Sample table extraction also failed: {e}")

        return None


# ── Sample group inference ────────────────────────────────────────────────────

def infer_sample_groups(metadata_path: str) -> dict:
    """
    Parse structured key:value pairs from characteristics field.
    Returns best-guess groups for human review.
    """
    meta   = pd.read_csv(metadata_path)
    groups = {}

    for _, row in meta.iterrows():
        characteristics = str(row.get("characteristics", ""))

        parsed = {}
        for item in characteristics.split(";"):
            item = item.strip()
            if ": " in item:
                key, value = item.split(": ", 1)
                parsed[key.strip().lower()] = value.strip().lower()

        disease   = parsed.get("disease",   "")
        treatment = parsed.get("treatment", "")

        if treatment and treatment != "control":
            guess = "exclude"
        elif any(w in disease for w in ("healthy", "normal", "control")):
            guess = "control"
        elif disease:
            guess = "disease"
        else:
            guess = "unknown"

        groups[row["sample_id"]] = {
            "guess":           guess,
            "title":           row.get("title", ""),
            "characteristics": characteristics
        }

    return groups


def assign_sample_groups_interactive(metadata_path: str,
                                     provenance: ProvenanceLog = None) -> dict:
    """
    Present each sample to the user and ask them to confirm or correct
    the automatically inferred group assignment.

    Saves a groups CSV alongside the metadata file for reuse.
    Returns groups keyed by sample title (to match count matrix columns)
    as well as GSM ID, so downstream tools can use either.
    """
    inferred    = infer_sample_groups(metadata_path)
    groups_path = metadata_path.replace("_metadata.csv", "_groups.csv")

    # Check for saved assignment from a previous run
    if os.path.exists(groups_path):
        print(f"\n  Found existing group assignment: {groups_path}")
        response = input("  Use saved assignment? (y/n): ").strip().lower()
        if response == "y":
            saved    = pd.read_csv(groups_path)
            by_title = dict(zip(saved["title"],     saved["group"]))
            by_gsm   = dict(zip(saved["sample_id"], saved["group"]))
            counts   = pd.Series(by_title).value_counts().to_dict()
            print(f"  Loaded {len(by_title)} sample assignments.\n")

            if provenance:
                provenance.record_step("sample_assignment", {
                    "group_counts": counts,
                    "n_excluded":   sum(1 for v in by_title.values()
                                        if v == "exclude"),
                    "ready_for_deg": (counts.get("disease", 0) >= 3 and
                                      counts.get("control", 0) >= 3),
                    "groups_file":  groups_path
                })
                provenance.record_output_file("sample_groups", groups_path)

            return {
                "sample_groups":        by_title,
                "sample_groups_by_gsm": by_gsm,
                "group_counts":         counts,
                "groups_file":          groups_path,
                "ready_for_deg": (counts.get("disease", 0) >= 3 and
                                  counts.get("control", 0) >= 3)
            }

    valid_choices = {"c": "control", "d": "disease", "e": "exclude"}

    print("\n" + "="*65)
    print("  SAMPLE GROUP ASSIGNMENT")
    print("="*65)
    print("  For each sample, confirm or correct the suggested group.")
    print("  Options: [c] control   [d] disease   [e] exclude")
    print("="*65)

    final_groups = {}

    # Sort so similar suggestions appear together
    sorted_samples = sorted(
        inferred.items(),
        key=lambda x: (x[1]["guess"], x[0])
    )

    current_guess = None
    for sample_id, info in sorted_samples:

        if info["guess"] != current_guess:
            current_guess = info["guess"]
            print(f"\n  --- Suggested: {current_guess.upper()} ---")

        print(f"\n  Sample:     {sample_id}  ({info['title']})")
        print(f"  Info:       {info['characteristics']}")
        print(f"  Suggestion: [{info['guess'][0]}] {info['guess']}")

        while True:
            raw = input(
                "  Your choice (c/d/e) or Enter to accept suggestion: "
            ).strip().lower()

            if raw == "":
                chosen = info["guess"]
                if chosen == "unknown":
                    print("  ⚠️  No suggestion — please enter c, d, or e.")
                    continue
                break
            elif raw in valid_choices:
                chosen = valid_choices[raw]
                break
            else:
                print("  Please enter c, d, or e.")

        final_groups[sample_id] = {
            "sample_id": sample_id,
            "title":     info["title"],
            "group":     chosen
        }
        print(f"  → Assigned: {chosen}")

    # Save to CSV
    groups_df = pd.DataFrame(final_groups.values())
    groups_df.to_csv(groups_path, index=False)
    print(f"\n  ✓  Group assignments saved to: {groups_path}")

    by_title = {r["title"]:     r["group"] for r in final_groups.values()}
    by_gsm   = {r["sample_id"]: r["group"] for r in final_groups.values()}
    counts   = pd.Series(by_title).value_counts().to_dict()
    ready    = (counts.get("disease", 0) >= 3 and
                counts.get("control", 0) >= 3)

    # Summary
    print("\n" + "="*65)
    print("  ASSIGNMENT SUMMARY")
    print("="*65)
    for group, count in sorted(counts.items()):
        print(f"  {group:10s}: {count} samples")
    print(f"\n  {'✓  Ready for DESeq2' if ready else '⚠️  Not ready — need at least 3 per group'}")
    print("="*65 + "\n")

    if provenance:
        provenance.record_step("sample_assignment", {
            "group_counts": counts,
            "n_excluded":   len([v for v in final_groups.values()
                                  if v["group"] == "exclude"]),
            "ready_for_deg": ready,
            "groups_file":  groups_path
        })
        provenance.record_output_file("sample_groups", groups_path)
        if not ready:
            provenance.record_warning(
                "sample_assignment",
                f"Insufficient samples for DEG — "
                f"disease: {counts.get('disease', 0)}, "
                f"control: {counts.get('control', 0)}"
            )

    return {
        "sample_groups":        by_title,
        "sample_groups_by_gsm": by_gsm,
        "group_counts":         counts,
        "groups_file":          groups_path,
        "ready_for_deg":        ready
    }


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()

    accession = "GSE226646"
    prov = ProvenanceLog.get_or_create("Friedreich ataxia")

    print(f"Testing GEO download for {accession}...\n")

    result = download_geo_dataset(accession, provenance=prov)

    print(f"\nTitle:        {result.get('title')}")
    print(f"Samples:      {result.get('n_samples')}")
    print(f"Platforms:    {result.get('platforms')}")
    print(f"Data source:  {result.get('data_source')}")
    print(f"Count matrix: {result.get('count_matrix_path')}")
    print(f"Metadata:     {result.get('metadata_path')}")

    if result.get("issues"):
        print("\nIssues flagged:")
        for issue in result["issues"]:
            print(f"  ⚠️  {issue}")

    if result.get("count_matrix_path") is None:
        print("\n✗  Cannot proceed — no count matrix obtained")
        sys.exit(1)

    if result.get("metadata_path"):
        print("\nInferring sample groups...")
        groups = assign_sample_groups_interactive(
            result["metadata_path"], provenance=prov
        )
        print(f"Group counts:  {groups.get('group_counts')}")
        print(f"Ready for DEG: {groups.get('ready_for_deg')}")

