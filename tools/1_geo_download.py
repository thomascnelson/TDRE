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

# ── Main download function ────────────────────────────────────────────────────

def download_geo_dataset(geo_accession: str, data_dir: str = "data",
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
    os.makedirs(data_dir, exist_ok=True)

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
            destdir=data_dir,
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
    meta_path = os.path.join(data_dir, f"{geo_accession}_metadata.csv")
    meta_df.to_csv(meta_path, index=False)
    result["metadata_path"] = meta_path
    result["samples"]       = sample_meta

    # ── Try supplementary count file first ───────────────────────────────────
    supp_files = gse.metadata.get("supplementary_file", [])
    count_matrix_path = _try_get_supplementary_counts(geo_accession, supp_files, data_dir)

    if count_matrix_path:
        result["count_matrix_path"] = count_matrix_path
        result["data_source"]       = "supplementary_counts"
    else:
        # Fall back to GEO matrix file
        matrix_path = _extract_matrix_from_gse(gse, geo_accession, data_dir)

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
        provenance.record_dataset(
            accession   = geo_accession,
            title       = result["title"],
            n_samples   = result["n_samples"],
            n_disease   = 0,   # updated later after sample assignment
            n_control   = 0,   # updated later after sample assignment
            organism    = gse.metadata.get("sample_organism_ch1", [""])[0]
                          or gse.metadata.get("organism_ch1", [""])[0],
            data_source = result["data_source"]
        )
        provenance.record_output_file("count_matrix",
                                      result["count_matrix_path"])
        provenance.record_output_file("sample_metadata",
                                      result["metadata_path"])
        for issue in result["issues"]:
            provenance.record_warning("geo_download", issue)

    return result


# ── Supplementary file download ───────────────────────────────────────────────

def _try_get_supplementary_counts(accession: str,
                                   supp_files: list,
                                   data_dir: str = "data") -> str | None:
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
        dest     = os.path.join(data_dir, filename)

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

        # Parse the file — try comma, tab, and semicolon delimiters.
        # Validate that the result has multiple numeric columns; if not, the
        # wrong delimiter was used (pandas parses silently but produces garbage).
        compression = "gzip" if dest.endswith(".gz") else None
        parsed = False
        for sep_name, sep in [("comma", ","), ("tab", "\t"), ("semicolon", ";")]:
            try:
                df = pd.read_csv(dest, sep=sep, index_col=0,
                                 compression=compression, quoting=csv.QUOTE_NONE)
                df.index = df.index.str.strip('"')
                # Strip trailing commas from column names (artefact of some GEO files)
                df.columns = [c.rstrip(',') for c in df.columns]
                # Drop non-numeric annotation columns (e.g. gene_name, gene_biotype).
                # Count matrix sample columns should be entirely numeric; text
                # annotation columns will fail conversion and are not needed for DESeq2.
                df_num = df.apply(pd.to_numeric, errors='coerce')
                keep   = df_num.columns[df_num.notna().mean() > 0.9]
                df     = df_num[keep]
                if df.shape[1] < 2:
                    print(f"  {sep_name}-separated: only {df.shape[1]} numeric "
                          f"column(s) — trying next delimiter")
                    continue
                clean_path = os.path.join(data_dir, f"{accession}_counts.csv")
                df.to_csv(clean_path)
                print(f"  Count matrix saved ({sep_name}-separated): {clean_path} "
                      f"({df.shape[0]} genes x {df.shape[1]} samples)")
                parsed = True
                break
            except Exception as e:
                print(f"  {sep_name}-separated parse failed: {e}")
                continue

        if parsed:
            return clean_path
        print(f"  Could not parse {filename} with any delimiter")
        continue

    print("  No supplementary count file found or downloaded successfully")
    return None


# ── GEO matrix fallback ───────────────────────────────────────────────────────

def _extract_matrix_from_gse(gse, accession: str, data_dir: str = "data") -> str | None:
    """
    Extract expression matrix from GEO soft file as a last resort.
    Note: RNA-seq data often fails here because GEOparse's pivot_samples
    expects a microarray-style ID_REF column that RNA-seq data doesn't have.
    """
    try:
        pivot = gse.pivot_samples("VALUE")
        path  = os.path.join(data_dir, f"{accession}_matrix.csv")
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
                path = os.path.join(data_dir, f"{accession}_matrix.csv")
                df.to_csv(path)
                print(f"  Sample table matrix extracted: {path}")
                return path

        except Exception as e:
            print(f"  Sample table extraction also failed: {e}")

        return None


# ── Sample group inference ────────────────────────────────────────────────────

def infer_sample_groups(metadata_path: str) -> dict:
    """
    Parse structured key:value pairs from characteristics field and make
    a best-guess group assignment (disease / control / exclude / unknown)
    for each sample.

    GEO metadata is inconsistent: the relevant field could be called
    "disease", "condition", "genotype", "status", "phenotype", "group",
    "diagnosis", "treatment", etc.  We check all of them.
    """
    meta   = pd.read_csv(metadata_path)
    groups = {}

    # Keys that carry disease/condition information (checked in priority order)
    condition_keys = [
        "disease", "condition", "genotype", "genotype/variation",
        "status", "phenotype", "group", "diagnosis", "cell type",
        "cell line", "sample type", "tissue", "subject status"
    ]
    # Keys that indicate an experimental treatment (non-baseline samples)
    treatment_keys = [
        "treatment", "treatment protocol", "drug", "compound",
        "stimulation", "perturbation"
    ]
    # Values that suggest a healthy/control sample
    control_words  = ("healthy", "normal", "control", "wildtype", "wild-type",
                      "wt", "unaffected", "naive", "untreated", "vehicle",
                      "mock", "ctrl")
    # Values that suggest an active treatment (not baseline)
    treated_words  = ("treated", "stimulated", "infected", "transfected",
                      "knockdown", "knockout", "overexpression")

    for _, row in meta.iterrows():
        characteristics = str(row.get("characteristics", ""))

        parsed = {}
        for item in characteristics.split(";"):
            item = item.strip()
            if ": " in item:
                key, value = item.split(": ", 1)
                parsed[key.strip().lower()] = value.strip().lower()

        # Look for a condition value across all candidate keys
        condition_value = ""
        for ck in condition_keys:
            if ck in parsed:
                condition_value = parsed[ck]
                break

        # Look for a treatment value
        treatment_value = ""
        for tk in treatment_keys:
            if tk in parsed:
                treatment_value = parsed[tk]
                break

        # Determine guess
        if treatment_value and not any(cw in treatment_value for cw in control_words):
            # Sample received an active treatment — exclude from disease-vs-control DEG
            guess = "exclude"
        elif any(cw in condition_value for cw in control_words):
            guess = "control"
        elif condition_value:
            # Non-empty condition that doesn't look like a control → disease
            guess = "disease"
        else:
            # Last resort: check title and source for control keywords
            title_source = (str(row.get("title", "")) + " " +
                            str(row.get("source", ""))).lower()
            if any(cw in title_source for cw in control_words):
                guess = "control"
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

    # ── Pre-assignment summary ────────────────────────────────────────────────
    # Parse all unique field names and their values across the dataset so the
    # user can quickly understand what metadata dimensions exist before being
    # asked to classify each sample one at a time.
    all_parsed = {}   # field_name → {value → count}
    meta = pd.read_csv(metadata_path)
    for _, row in meta.iterrows():
        for item in str(row.get("characteristics", "")).split(";"):
            item = item.strip()
            if ": " in item:
                key, value = item.split(": ", 1)
                key   = key.strip().lower()
                value = value.strip().lower()
                all_parsed.setdefault(key, {})
                all_parsed[key][value] = all_parsed[key].get(value, 0) + 1

    print("\n" + "="*65)
    print("  SAMPLE GROUP ASSIGNMENT")
    print("="*65)
    print(f"  Dataset contains {len(inferred)} samples.")
    if all_parsed:
        print("\n  Characteristics fields present in this dataset:")
        for field, values in sorted(all_parsed.items()):
            value_str = ", ".join(
                f"{v} (n={c})" for v, c in sorted(values.items())
            )
            print(f"    {field}: {value_str}")
    print("\n  For each sample, confirm or correct the suggested group.")
    print("  Options: [c] control   [d] disease   [e] exclude")
    print("="*65)

    # ── Batch assignment ──────────────────────────────────────────────────────
    # Offer the user a chance to do field=value → group batch mappings
    # before the per-sample loop.  This handles clean, consistent datasets
    # without requiring 20+ individual confirmations.
    final_groups = {}

    if all_parsed:
        print("\n  BATCH ASSIGNMENT (optional)")
        print("  Assign all samples where a field matches a value to a group.")
        print("  Example: field='condition', value='healthy', group=c")
        print("  Press Enter to skip batch assignment.\n")

        while True:
            raw = input("  Batch rule — field (or Enter to skip): ").strip().lower()
            if not raw:
                break
            if raw not in all_parsed:
                print(f"  Field '{raw}' not found. Available: "
                      f"{', '.join(sorted(all_parsed.keys()))}")
                continue
            val = input(f"  Value for '{raw}': ").strip().lower()
            if val not in all_parsed[raw]:
                print(f"  Value '{val}' not found under '{raw}'. "
                      f"Available: {', '.join(sorted(all_parsed[raw].keys()))}")
                continue
            while True:
                grp_raw = input(
                    f"  Assign '{raw}: {val}' to (c/d/e): "
                ).strip().lower()
                if grp_raw in valid_choices:
                    assigned_group = valid_choices[grp_raw]
                    break
                print("  Please enter c, d, or e.")

            # Apply batch rule to all matching samples
            n_applied = 0
            for sample_id, info in inferred.items():
                char = info["characteristics"].lower()
                # Match "key: value" anywhere in the semicolon-separated string
                if f"{raw}: {val}" in char:
                    final_groups[sample_id] = {
                        "sample_id": sample_id,
                        "title":     info["title"],
                        "group":     assigned_group
                    }
                    n_applied += 1
            print(f"  → Applied to {n_applied} samples as '{assigned_group}'.\n")

        already_assigned = set(final_groups.keys())
        remaining = {k: v for k, v in inferred.items()
                     if k not in already_assigned}
        if remaining:
            print(f"\n  {len(remaining)} sample(s) still need manual assignment.")
        else:
            print("\n  All samples assigned via batch rules.")
    else:
        remaining = dict(inferred)

    # ── Per-sample loop for anything not already batch-assigned ───────────────
    # Sort so similar suggestions appear together
    sorted_remaining = sorted(
        remaining.items(),
        key=lambda x: (x[1]["guess"], x[0])
    )

    current_guess = None
    for sample_id, info in sorted_remaining:

        if info["guess"] != current_guess:
            current_guess = info["guess"]
            print(f"\n  --- Suggested: {current_guess.upper()} ---")

        # Show title/source as fallback if characteristics is empty
        display_info = info["characteristics"].strip()
        if not display_info:
            display_info = (f"title={info['title']}  "
                            f"source={meta.loc[meta['sample_id'] == sample_id, 'source'].values[0] if 'source' in meta.columns else ''}")

        print(f"\n  Sample:     {sample_id}  ({info['title']})")
        print(f"  Info:       {display_info}")
        print(f"  Suggestion: [{info['guess'][0]}] {info['guess']}")

        while True:
            raw = input(
                "  Your choice (c/d/e) or Enter to accept suggestion: "
            ).strip().lower()

            if raw == "":
                chosen = info["guess"]
                if chosen == "unknown":
                    print("  No suggestion — please enter c, d, or e.")
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

    accession = "GSE297335"
    prov = ProvenanceLog.get_or_create("Barth Syndrome")

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

