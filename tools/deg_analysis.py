# tools/deg_analysis.py
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv
from utils.provenance import ProvenanceLog

def run_deseq2(count_matrix_path: str, sample_groups: dict, output_prefix: str = "deg", provenance: ProvenanceLog = None) -> dict:
    """
    Run DESeq2 differential expression analysis.
    Returns path to results file and summary statistics.
    """
    # Load count matrix
    counts = pd.read_csv(count_matrix_path, index_col=0, quoting=csv.QUOTE_NONE)
    counts.index = counts.index.str.strip('"')

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
    counts_int = counts_int.loc[:, counts_int.sum(axis=0) >= 10]

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
    os.makedirs("data", exist_ok=True)
    results_path = f"data/{output_prefix}_deseq2_results.csv"
    results.to_csv(results_path)

    # Summary stats
    sig_up = results[(results["padj"] < 0.05) & (results["log2FoldChange"] > 1)]
    sig_down = results[(results["padj"] < 0.05) & (results["log2FoldChange"] < -1)]

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
            "padj_threshold":    0.05,
            "lfc_threshold":     1.0,
            "min_count_filter":  10
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



### TEST ROUTINE
if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()

    # These files were created by Tool 3 — adjust paths if needed
    COUNT_MATRIX = "data/GSE226646_counts.csv"
    GROUPS_FILE  = "data/GSE226646_groups.csv"
    
    # Create the provenance log
    prov = ProvenanceLog.get_or_create("Friedreich ataxia")

    # Load the saved group assignments (title-keyed, to match count matrix)
    groups_df = pd.read_csv(GROUPS_FILE)
    sample_groups = dict(zip(groups_df["title"], groups_df["group"]))

    print("Sample groups loaded:")
    for sample, group in sample_groups.items():
        print(f"  {sample}: {group}")

    print(f"\nRunning DESeq2 on {COUNT_MATRIX}...")
    result = run_deseq2(COUNT_MATRIX, sample_groups, output_prefix="GSE226646", provenance = prov)

    if "error" in result:
        print(f"\n✗  DESeq2 failed: {result['error']}")
    else:
        print("Analysis complete")
        print(f"   Genes tested:       {result['n_genes_tested']}")
        print(f"   Significant up:     {result['n_significant_up']}")
        print(f"   Significant down:   {result['n_significant_down']}")
        print(f"   Results saved to:   {result['results_path']}")
        print(f"\n   Top upregulated:   {result['top_upregulated']}")
        print(f"   Top downregulated: {result['top_downregulated']}")

