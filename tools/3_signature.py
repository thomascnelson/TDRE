# tools/signature.py
import json
import re
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mygene
import pandas as pd

from utils.provenance import ProvenanceLog

def build_signature(deg_results_path: str,
                    padj_threshold: float = 0.05,
                    lfc_threshold: float = 1.0,
                    max_genes: int = 150,
                    provenance: ProvenanceLog = None) -> dict:

    results = pd.read_csv(deg_results_path, index_col=0)

    sig = results[
        (results["padj"] < padj_threshold) &
        (results["log2FoldChange"].abs() > lfc_threshold)
    ].copy()

    sig_up   = sig[sig["log2FoldChange"] > 0].nsmallest(max_genes // 2, "padj")
    sig_down = sig[sig["log2FoldChange"] < 0].nsmallest(max_genes // 2, "padj")

    mg = mygene.MyGeneInfo()

    def map_genes(gene_list):
        if not gene_list:
            return []
        result = mg.querymany(
            gene_list,
            scopes="symbol,ensembl.gene,entrezgene",
            fields="symbol,entrezgene",
            species="human",
            returnall=True
        )
        mapped = []
        for hit in result.get("out", []):
            if "symbol" in hit:
                mapped.append(hit["symbol"])
        return list(set(mapped))

    up_genes   = map_genes(sig_up.index.tolist())
    down_genes = map_genes(sig_down.index.tolist())

    # Filter pseudogenes, keep named lncRNAs
    pseudogene_patterns = r'(P\d+$|^\d+P\d+$)'
    up_genes   = [g for g in up_genes   if not re.search(pseudogene_patterns, g)]
    down_genes = [g for g in down_genes if not re.search(pseudogene_patterns, g)]

    signature = {
        "upregulated": up_genes,
        "downregulated": down_genes,
        "n_up": len(up_genes),
        "n_down": len(down_genes),
        "thresholds_used": {
            "padj": padj_threshold,
            "log2fc": lfc_threshold
        },
        "signature_quality": "good" if len(up_genes) >= 10 and len(down_genes) >= 10
                             else "weak — consider relaxing thresholds"
    }

    # Derive output paths from the input path
    base = deg_results_path.replace("_deseq2_results.csv", "")

    # JSON — for pipeline reuse
    json_path = f"{base}_signature.json"
    with open(json_path, "w") as f:
        json.dump(signature, f, indent=2)

    # CSV — for the investigator to inspect
    csv_path = f"{base}_signature.csv"
    max_len  = max(len(up_genes), len(down_genes))
    sig_df   = pd.DataFrame({
        "upregulated":   up_genes   + [""] * (max_len - len(up_genes)),
        "downregulated": down_genes + [""] * (max_len - len(down_genes))
    })
    sig_df.to_csv(csv_path, index=False)

    signature["json_path"] = json_path
    signature["csv_path"]  = csv_path

    print(f"  Signature saved to: {json_path}")
    print(f"  Signature saved to: {csv_path}")

    if provenance:
        provenance.record_parameters("signature_building", {
            "padj_threshold": padj_threshold,
            "lfc_threshold":  lfc_threshold,
            "max_genes":      max_genes
        })
        provenance.record_step("signature_building", {
            "n_up":              signature["n_up"],
            "n_down":            signature["n_down"],
            "signature_quality": signature["signature_quality"],
            "top_up":            signature["upregulated"][:15],
            "top_down":          signature["downregulated"][:15]
        })
        provenance.record_output_file("signature_json", signature["json_path"])
        provenance.record_output_file("signature_csv",  signature["csv_path"])
        if signature["signature_quality"].startswith("weak"):
            provenance.record_warning("signature_building",
                f"Weak signature: only {signature['n_up']} up and "
                f"{signature['n_down']} down genes — "
                f"consider relaxing padj or lfc thresholds")

    return signature


if __name__ == "__main__":
    import argparse
    from dotenv import load_dotenv
    load_dotenv()

    parser = argparse.ArgumentParser(
        description=(
            "Tool 3 — Build a disease transcriptomic signature from DESeq2 results.\n"
            "Maps DEG gene IDs to canonical symbols via MyGeneInfo.\n"
            "Produces: {accession}_signature.json and _signature.csv"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-a", "--accession", required=True,
                        help="GEO accession number (e.g. GSE297335)")
    parser.add_argument("-d", "--disease", required=True,
                        help="Disease name for the provenance log")
    parser.add_argument("--deg-results",
                        help="Path to DESeq2 results CSV "
                             "(default: data/{accession}/{accession}_deseq2_results.csv)")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="Adjusted p-value threshold (default: 0.05)")
    parser.add_argument("--lfc", type=float, default=1.0,
                        help="log2 fold-change threshold (default: 1.0)")
    parser.add_argument("--max-genes", type=int, default=150,
                        help="Maximum signature size — split evenly up/down (default: 150)")
    args = parser.parse_args()

    accession   = args.accession.strip().upper()
    data_dir    = os.path.join("data", accession)
    deg_results = args.deg_results or os.path.join(
        data_dir, f"{accession}_deseq2_results.csv"
    )

    prov = ProvenanceLog.get_or_create(args.disease)

    print(f"\nBuilding signature from {deg_results}...")
    print(f"  padj ≤ {args.padj}   |log2FC| ≥ {args.lfc}   max genes: {args.max_genes}\n")

    signature = build_signature(
        deg_results_path = deg_results,
        padj_threshold   = args.padj,
        lfc_threshold    = args.lfc,
        max_genes        = args.max_genes,
        provenance       = prov
    )

    print(f"\n  Quality            : {signature['signature_quality']}")
    print(f"  Thresholds used    : padj < {signature['thresholds_used']['padj']},"
          f"  |log2FC| > {signature['thresholds_used']['log2fc']}")
    print(f"  Upregulated genes  : {signature['n_up']}")
    print(f"  Downregulated genes: {signature['n_down']}")
    print(f"\n  Top upregulated    : {signature['upregulated'][:10]}")
    print(f"  Top downregulated  : {signature['downregulated'][:10]}")
    print(f"\n  JSON : {signature['json_path']}")
    print(f"  CSV  : {signature['csv_path']}")



