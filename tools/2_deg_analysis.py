# tools/deg_analysis.py
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv
from utils.provenance import ProvenanceLog

def run_deseq2(count_matrix_path: str, groups_path: str, output_prefix: str = "deg",
               padj_threshold: float = 0.05, lfc_threshold: float = 1.0,
               min_count_filter: int = 10, data_dir: str = "data",
               provenance: ProvenanceLog = None) -> dict:
    """
    Run DESeq2 differential expression analysis.

    groups_path : path to the _groups.csv file produced by Tool 1.
                  Must contain columns: sample_id (GSM ID), title, group.

    Returns path to results file and summary statistics.
    """
    # Load count matrix
    counts = pd.read_csv(count_matrix_path, index_col=0, quoting=csv.QUOTE_NONE)
    counts.index = counts.index.str.strip('"')

    # Load group assignments and build two lookup dicts
    groups_df  = pd.read_csv(groups_path)
    by_gsm     = dict(zip(groups_df["sample_id"], groups_df["group"]))
    by_title   = dict(zip(groups_df["title"],     groups_df["group"]))

    # Determine which key type matches the count matrix columns.
    # Supplementary files often use GSM IDs; GEO matrix pivots use titles.
    gsm_matches   = sum(1 for k in by_gsm   if k in counts.columns)
    title_matches = sum(1 for k in by_title if k in counts.columns)

    if gsm_matches >= title_matches:
        sample_groups = by_gsm
        match_type    = "GSM ID"
    else:
        sample_groups = by_title
        match_type    = "sample title"

    print(f"  Matching samples by {match_type} "
          f"({max(gsm_matches, title_matches)} of "
          f"{len(groups_df)} samples found in count matrix)")

    # Filter to only samples with known groups
    valid_samples = {k: v for k, v in sample_groups.items()
                     if v in ("disease", "control") and k in counts.columns}

    if len([v for v in valid_samples.values() if v == "disease"]) < 3:
        return {"error": "Fewer than 3 disease samples — cannot run DESeq2"}
    if len([v for v in valid_samples.values() if v == "control"]) < 3:
        return {"error": "Fewer than 3 control samples — cannot run DESeq2"}

    # Subset and order
    sample_ids = list(valid_samples.keys())
    counts_subset = counts[sample_ids].T  # PyDESeq2 wants samples as rows

    # Round to integers (required for DESeq2)
    counts_int = counts_subset.round().astype(int)

    # Remove genes with all zeros or very low counts
    counts_int = counts_int.loc[:, counts_int.sum(axis=0) >= min_count_filter]

    # Build metadata dataframe
    meta = pd.DataFrame({
        "condition": [valid_samples[s] for s in sample_ids]
    }, index=sample_ids)

    # Run DESeq2
    dds = DeseqDataSet(
        counts=counts_int,
        metadata=meta,
        design="~condition")
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["condition", "disease", "control"])
    stat_res.summary()

    results = stat_res.results_df.copy()
    results = results.dropna(subset=["padj"])
    results = results.sort_values("padj")

    # Save full results
    os.makedirs(data_dir, exist_ok=True)
    results_path = os.path.join(data_dir, f"{output_prefix}_deseq2_results.csv")
    results.to_csv(results_path)

    # Summary stats
    sig_up   = results[(results["padj"] < padj_threshold) & (results["log2FoldChange"] >  lfc_threshold)]
    sig_down = results[(results["padj"] < padj_threshold) & (results["log2FoldChange"] < -lfc_threshold)]

    result = {
        "results_path": results_path,
        "n_genes_tested": len(results),
        "n_significant_up": len(sig_up),
        "n_significant_down": len(sig_down),
        "top_upregulated": sig_up.head(10).index.tolist(),
        "top_downregulated": sig_down.head(10).index.tolist(),
        "analysis_successful": len(sig_up) + len(sig_down) > 0
    }

    if provenance:
        provenance.record_parameters("deg_analysis", {
            "padj_threshold":   padj_threshold,
            "lfc_threshold":    lfc_threshold,
            "min_count_filter": min_count_filter
        })
        provenance.record_step("deg_analysis", {
            "n_disease_samples":  len([v for v in valid_samples.values()
                                       if v == "disease"]),
            "n_control_samples":  len([v for v in valid_samples.values()
                                       if v == "control"]),
            "n_genes_tested":     result["n_genes_tested"],
            "n_significant_up":   result["n_significant_up"],
            "n_significant_down": result["n_significant_down"],
            "top_upregulated":    result["top_upregulated"],
            "top_downregulated":  result["top_downregulated"]
        })
        provenance.record_output_file("deg_results", result["results_path"])
        if result["n_significant_up"] + result["n_significant_down"] < 20:
            provenance.record_warning("deg_analysis",
                f"Very few significant DEGs found "
                f"({result['n_significant_up']} up, "
                f"{result['n_significant_down']} down) — "
                f"consider relaxing thresholds")

    return result



if __name__ == "__main__":
    import argparse
    from dotenv import load_dotenv
    load_dotenv()

    parser = argparse.ArgumentParser(
        description=(
            "Tool 2 — Run DESeq2 differential gene expression analysis.\n"
            "Requires the count matrix and group assignments produced by Tool 1.\n"
            "Produces: {accession}_deseq2_results.csv"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-a", "--accession", required=True,
                        help="GEO accession number (e.g. GSE297335)")
    parser.add_argument("-d", "--disease", required=True,
                        help="Disease name for the provenance log")
    parser.add_argument("--counts",
                        help="Path to count matrix CSV "
                             "(default: data/{accession}/{accession}_counts.csv)")
    parser.add_argument("--groups",
                        help="Path to sample groups CSV "
                             "(default: data/{accession}/{accession}_groups.csv)")
    parser.add_argument("--data-dir",
                        help="Output directory (default: data/{accession})")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="Adjusted p-value threshold (default: 0.05)")
    parser.add_argument("--lfc", type=float, default=1.0,
                        help="log2 fold-change threshold (default: 1.0)")
    parser.add_argument("--min-count", type=int, default=10,
                        help="Minimum total count filter per gene (default: 10)")
    args = parser.parse_args()

    accession = args.accession.strip().upper()
    data_dir  = args.data_dir or os.path.join("data", accession)
    counts    = args.counts   or os.path.join(data_dir, f"{accession}_counts.csv")
    groups    = args.groups   or os.path.join(data_dir, f"{accession}_groups.csv")

    prov = ProvenanceLog.get_or_create(args.disease)

    print(f"\nRunning DESeq2 for {accession}...")
    print(f"  Count matrix : {counts}")
    print(f"  Groups file  : {groups}")
    print(f"  padj ≤ {args.padj}   |log2FC| ≥ {args.lfc}   min count ≥ {args.min_count}\n")

    result = run_deseq2(
        count_matrix_path = counts,
        groups_path       = groups,
        output_prefix     = accession,
        data_dir          = data_dir,
        padj_threshold    = args.padj,
        lfc_threshold     = args.lfc,
        min_count_filter  = args.min_count,
        provenance        = prov
    )

    if "error" in result:
        print(f"\n  DESeq2 failed: {result['error']}")
        sys.exit(1)

    print(f"  Genes tested       : {result['n_genes_tested']}")
    print(f"  Significant up     : {result['n_significant_up']}")
    print(f"  Significant down   : {result['n_significant_down']}")
    print(f"  Results saved to   : {result['results_path']}")
    print(f"\n  Top upregulated    : {result['top_upregulated'][:5]}")
    print(f"  Top downregulated  : {result['top_downregulated'][:5]}")

