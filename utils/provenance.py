# utils/provenance.py
import json
import os
from datetime import datetime


class ProvenanceLog:
    """
    Incrementally builds a machine-readable record of a drug repurposing analysis.
    Each tool appends its own section when it completes.

    The log is written to disk after every update so that if the pipeline
    crashes partway through, nothing is lost. The orchestrator reads this
    log at the end to generate the human-readable Markdown report.

    Usage:
        # At the start of the pipeline (in the orchestrator)
        prov = ProvenanceLog("Friedreich ataxia")

        # In each tool, passed as an optional argument
        prov.record_step("disease_lookup", {...})
        prov.record_warning("geo_download", "No raw count file found...")
        prov.record_parameters("deg_analysis", {"padj": 0.05, "lfc": 1.0})
        prov.record_output_file("deg_results", "data/GSE226646_deseq2_results.csv")

        # At the end of the pipeline
        report_path = prov.generate_report(drug_candidates, llm_summary)
    """

    def __init__(self, disease_name: str, output_dir: str = "reports"):
        os.makedirs(output_dir, exist_ok=True)

        self.disease_name = disease_name
        self.output_dir   = output_dir

        # Sanitise disease name for use in filenames
        safe_name = disease_name.replace(" ", "_").replace("/", "-")

        self.log_path    = os.path.join(output_dir, f"{safe_name}_provenance.json")
        self.report_path = os.path.join(output_dir, f"{safe_name}_report.md")

        self.data = {
            "disease_name":   disease_name,
            "created":        datetime.now().isoformat(),
            "last_updated":   datetime.now().isoformat(),
            "pipeline_steps": {},   # step_name → {completed, ...details}
            "warnings":       [],   # list of {step, time, message}
            "parameters":     {},   # step_name → {param: value}
            "output_files":   {},   # label → filepath
            "dataset":        {},   # GEO dataset metadata
            "summary_stats":  {}    # key numbers for the report header
        }

        self._save()
        print(f"  Provenance log: {self.log_path}")


    # ── Recording methods ─────────────────────────────────────────────────────

    def record_step(self, step_name: str, details: dict):
        """
        Record the completion of a pipeline step with its key details.

        step_name should be one of:
            disease_lookup, geo_search, dataset_selection, geo_download,
            sample_assignment, deg_analysis, signature_building, drug_repurposing
        """
        self.data["pipeline_steps"][step_name] = {
            "completed": datetime.now().isoformat(),
            **details
        }
        self._update_summary_stats(step_name, details)
        self._save()

    def record_warning(self, step: str, message: str):
        """Record a QC warning raised during the analysis."""
        self.data["warnings"].append({
            "step":    step,
            "time":    datetime.now().isoformat(),
            "message": message
        })
        print(f"  ⚠️  Warning recorded ({step}): {message[:80]}")
        self._save()

    def record_parameters(self, step: str, params: dict):
        """Record the parameters used at a given step."""
        self.data["parameters"][step] = params
        self._save()

    def record_output_file(self, label: str, path: str):
        """Record a file produced by the pipeline."""
        self.data["output_files"][label] = path
        self._save()

    def record_dataset(self, accession: str, title: str, n_samples: int,
                       n_disease: int, n_control: int, organism: str,
                       data_source: str, tissue: str = ""):
        """Record the GEO dataset chosen for analysis."""
        self.data["dataset"] = {
            "accession":   accession,
            "title":       title,
            "n_samples":   n_samples,
            "n_disease":   n_disease,
            "n_control":   n_control,
            "organism":    organism,
            "tissue":      tissue,
            "data_source": data_source,
            "geo_url":     f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
        }
        self._save()

    def record_rejected_dataset(self, accession: str, title: str, reason: str):
        """Record a dataset that was considered but rejected."""
        rejected = self.data["pipeline_steps"].get("geo_search", {})
        rejected_list = rejected.get("rejected_datasets", [])
        rejected_list.append({
            "accession": accession,
            "title":     title,
            "reason":    reason
        })
        if "geo_search" not in self.data["pipeline_steps"]:
            self.data["pipeline_steps"]["geo_search"] = {}
        self.data["pipeline_steps"]["geo_search"]["rejected_datasets"] = rejected_list
        self._save()


    # ── Convenience loaders ───────────────────────────────────────────────────

    @classmethod
    def load(cls, log_path: str) -> "ProvenanceLog":
        """
        Load an existing provenance log from disk.
        Useful for resuming a crashed pipeline or generating a report
        from a completed run.
        """
        with open(log_path) as f:
            data = json.load(f)

        instance = cls.__new__(cls)
        instance.disease_name = data["disease_name"]
        instance.output_dir   = os.path.dirname(log_path)
        instance.log_path     = log_path

        safe_name = data["disease_name"].replace(" ", "_").replace("/", "-")
        instance.report_path  = os.path.join(
            instance.output_dir, f"{safe_name}_report.md"
        )
        instance.data = data
        return instance

    def get_step(self, step_name: str) -> dict:
        """Return the recorded details for a given step, or empty dict."""
        return self.data["pipeline_steps"].get(step_name, {})

    def has_warnings(self) -> bool:
        return len(self.data["warnings"]) > 0

    def steps_completed(self) -> list[str]:
        return list(self.data["pipeline_steps"].keys())

    @classmethod
    def get_or_create(cls, disease_name: str, output_dir: str = "reports") -> "ProvenanceLog":
        """
        Load an existing provenance log for this disease if one exists,
        otherwise create a new one. This is the right call to use in
        standalone tool testing so runs accumulate rather than overwrite.
        """
        safe_name = disease_name.replace(" ", "_").replace("/", "-")
        log_path  = os.path.join(output_dir, f"{safe_name}_provenance.json")
    
        if os.path.exists(log_path):
            print(f"  Resuming provenance log: {log_path}")
            return cls.load(log_path)
        else:
            return cls(disease_name, output_dir)

    # ── Report generation ─────────────────────────────────────────────────────

    def generate_report(
        self,
        drug_candidates: list[dict],
        llm_narrative:   str = ""
    ) -> str:
        """
        Generate a human-readable Markdown report from the provenance log.

        drug_candidates : ranked list of dicts from Tool 6
        llm_narrative   : optional free-text summary written by the LLM
                          orchestrator to be inserted after the header

        Returns the path to the saved report file.
        """
        lines = []
        d     = self.data
        ds    = d.get("dataset", {})
        stats = d.get("summary_stats", {})
        steps = d.get("pipeline_steps", {})
        params = d.get("parameters", {})

        # ── Header ────────────────────────────────────────────────────────────
        lines += [
            f"# Drug Repurposing Report: {d['disease_name']}",
            "",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}  ",
            f"**Pipeline version:** 0.1.0  ",
            "",
            "---",
            ""
        ]

        # ── LLM narrative summary (if provided) ───────────────────────────────
        if llm_narrative:
            lines += [
                "## Executive Summary",
                "",
                llm_narrative,
                "",
                "---",
                ""
            ]

        # ── Disease overview ──────────────────────────────────────────────────
        disease_step = steps.get("disease_lookup", {})
        lines += [
            "## Disease Overview",
            "",
            f"**Disease:** {d['disease_name']}  ",
        ]
        if disease_step.get("mesh_terms"):
            lines.append(f"**MeSH ID(s):** {', '.join(disease_step['mesh_terms'])}  ")
        if disease_step.get("pubmed_treatment_hits") is not None:
            lines.append(
                f"**PubMed treatment literature:** "
                f"{disease_step['pubmed_treatment_hits']} publications  "
            )
        feasibility = disease_step.get("feasibility", {})
        if feasibility:
            lines += [
                f"**Repurposing feasibility score:** "
                f"{feasibility.get('feasibility_score', 'N/A')}/100  ",
            ]
            for flag in feasibility.get("flags", []):
                lines.append(f"- ⚠️ {flag}")
        lines.append("")

        # ── Dataset ───────────────────────────────────────────────────────────
        if ds:
            lines += [
                "## Dataset",
                "",
                f"**GEO Accession:** [{ds.get('accession')}]({ds.get('geo_url')})  ",
                f"**Title:** {ds.get('title')}  ",
                f"**Organism:** {ds.get('organism')}  ",
            ]
            if ds.get("tissue"):
                lines.append(f"**Tissue:** {ds.get('tissue')}  ")
            lines += [
                f"**Total samples:** {ds.get('n_samples')}  ",
                f"**Disease samples:** {ds.get('n_disease')}  ",
                f"**Control samples:** {ds.get('n_control')}  ",
                f"**Count data source:** {ds.get('data_source')}  ",
                ""
            ]

        # Rejected datasets
        rejected = steps.get("geo_search", {}).get("rejected_datasets", [])
        if rejected:
            lines += ["**Datasets considered but rejected:**", ""]
            for r in rejected:
                lines.append(f"- `{r['accession']}` — {r['title']}: {r['reason']}")
            lines.append("")

        # ── Sample groups ─────────────────────────────────────────────────────
        assign_step = steps.get("sample_assignment", {})
        if assign_step:
            lines += [
                "## Sample Group Assignment",
                "",
                "| Group | N |",
                "|-------|---|",
            ]
            for group, count in assign_step.get("group_counts", {}).items():
                lines.append(f"| {group} | {count} |")
            lines.append("")

        # ── DEG analysis ──────────────────────────────────────────────────────
        deg_step   = steps.get("deg_analysis", {})
        deg_params = params.get("deg_analysis", {})
        if deg_step:
            lines += [
                "## Differential Gene Expression Analysis",
                "",
                "**Method:** DESeq2 (PyDESeq2)  ",
                f"**Genes tested:** {deg_step.get('n_genes_tested', 'N/A')}  ",
                f"**Significantly upregulated:** {deg_step.get('n_significant_up', 'N/A')}  ",
                f"**Significantly downregulated:** "
                f"{deg_step.get('n_significant_down', 'N/A')}  ",
                ""
            ]
            if deg_params:
                lines += [
                    "**Parameters:**",
                    "",
                    f"- Adjusted p-value threshold: {deg_params.get('padj_threshold', 0.05)}",
                    f"- log₂ fold-change threshold: {deg_params.get('lfc_threshold', 1.0)}",
                    f"- Minimum count filter: {deg_params.get('min_count_filter', 10)}",
                    ""
                ]

        # ── Signature ─────────────────────────────────────────────────────────
        sig_step   = steps.get("signature_building", {})
        sig_params = params.get("signature_building", {})
        if sig_step:
            lines += [
                "## Disease Transcriptomic Signature",
                "",
                f"**Quality:** {sig_step.get('signature_quality', 'N/A')}  ",
                f"**Upregulated genes:** {sig_step.get('n_up', 'N/A')}  ",
                f"**Downregulated genes:** {sig_step.get('n_down', 'N/A')}  ",
                ""
            ]
            if sig_step.get("top_up"):
                lines += [
                    "**Top upregulated:**  ",
                    "`" + "`, `".join(sig_step["top_up"][:15]) + "`",
                    ""
                ]
            if sig_step.get("top_down"):
                lines += [
                    "**Top downregulated:**  ",
                    "`" + "`, `".join(sig_step["top_down"][:15]) + "`",
                    ""
                ]

        # ── Drug candidates ───────────────────────────────────────────────────
        drug_step = steps.get("drug_repurposing", {})
        if drug_candidates:
            lines += [
                "## Drug Repurposing Candidates",
                "",
                "Drugs are ranked by L1000FWD combined score "
                "(more negative = stronger reversal of disease signature).  ",
                ""
            ]

            if drug_step.get("result_url"):
                lines += [
                    f"**Interactive results:** [{drug_step['result_url']}]"
                    f"({drug_step['result_url']})  ",
                    ""
                ]

            lines += [
                "| Rank | Drug | Score | Q-value | Confidence | "
                "Cell Line | MoA | Target |",
                "|------|------|-------|---------|------------|"
                "-----------|-----|--------|"
            ]

            for i, drug in enumerate(drug_candidates, 1):
                name  = drug.get("drug_name", "")
                # Truncate long IUPAC names in the table
                if len(name) > 35:
                    name = name[:32] + "..."
                lines.append(
                    f"| {i} "
                    f"| {name} "
                    f"| {drug.get('score', '')} "
                    f"| {drug.get('qvalue', '')} "
                    f"| {drug.get('confidence', '')} "
                    f"| {drug.get('cell_line', '')} "
                    f"| {drug.get('moa', '')} "
                    f"| {drug.get('target', '')} |"
                )

            lines.append("")
            lines += [
                "> **Interpretation note:** Candidates should be treated as "
                "hypotheses for further investigation, not validated drug targets. "
                "Experimental validation in disease-relevant cell or animal models "
                "is required before drawing conclusions.",
                ""
            ]

        # ── QC warnings ───────────────────────────────────────────────────────
        if d["warnings"]:
            lines += [
                "## Quality Control Warnings",
                "",
                "The following issues were flagged during the analysis:",
                ""
            ]
            for w in d["warnings"]:
                lines.append(f"- **{w['step']}** ({w['time'][:10]}): {w['message']}")
            lines.append("")

        # ── Output files ──────────────────────────────────────────────────────
        if d["output_files"]:
            lines += [
                "## Output Files",
                "",
                "| File | Path |",
                "|------|------|"
            ]
            for label, path in d["output_files"].items():
                lines.append(f"| {label} | `{path}` |")
            lines.append("")

        # ── Methods ───────────────────────────────────────────────────────────
        lines += [
            "## Methods",
            "",
            "**Transcriptomic data:** Retrieved from NCBI Gene Expression Omnibus (GEO). "
            "Raw count matrices were downloaded and quality-assessed prior to analysis.  ",
            "",
            "**Differential expression:** Performed using PyDESeq2, a Python "
            "implementation of the DESeq2 negative binomial model. Genes were "
            "filtered to a minimum count sum of 10 across all samples. "
            "Significance was assessed using Wald tests with Benjamini-Hochberg "
            "FDR correction.  ",
            "",
            "**Drug repurposing:** Disease signature genes were submitted to "
            "L1000FWD (Ma'ayan Laboratory, Mount Sinai) for connectivity scoring "
            "against the LINCS L1000 chemical perturbagen database. Drugs with "
            "negatively concordant signatures (i.e. those that reverse the disease "
            "gene expression pattern) were ranked by combined score.  ",
            ""
        ]

        # Write the report
        report_text = "\n".join(lines)
        with open(self.report_path, "w") as f:
            f.write(report_text)

        print(f"  Report saved to: {self.report_path}")
        return self.report_path


    # ── Internal helpers ──────────────────────────────────────────────────────

    def _save(self):
        self.data["last_updated"] = datetime.now().isoformat()
        with open(self.log_path, "w") as f:
            json.dump(self.data, f, indent=2, default=str)

    def _update_summary_stats(self, step_name: str, details: dict):
        """Keep a flat summary of key numbers for easy access."""
        stats = self.data["summary_stats"]

        if step_name == "deg_analysis":
            stats["n_genes_tested"]       = details.get("n_genes_tested")
            stats["n_significant_up"]     = details.get("n_significant_up")
            stats["n_significant_down"]   = details.get("n_significant_down")

        elif step_name == "signature_building":
            stats["n_signature_up"]   = details.get("n_up")
            stats["n_signature_down"] = details.get("n_down")

        elif step_name == "drug_repurposing":
            stats["n_drug_candidates"]      = details.get("n_candidates")
            stats["n_high_confidence_drugs"] = details.get("n_high_confidence")


# ── Standalone test ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    import json

    print("Testing ProvenanceLog with simulated pipeline run...\n")

    # Simulate a complete pipeline run
    prov = ProvenanceLog("Friedreich ataxia", output_dir="reports")

    prov.record_step("disease_lookup", {
        "mesh_terms": ["D005621"],
        "pubmed_treatment_hits": 847,
        "feasibility": {
            "feasibility_score": 80,
            "flags": [],
            "proceed_recommended": True
        }
    })

    prov.record_step("geo_search", {
        "n_datasets_found": 3,
        "dataset_selected": "GSE226646"
    })
    prov.record_rejected_dataset(
        "GSE83496", "FRDA iPSC study",
        "Only 2 control samples — insufficient for DESeq2"
    )

    prov.record_dataset(
        accession   = "GSE226646",
        title       = "Skeletal muscle transcriptome in Friedreich ataxia",
        n_samples   = 20,
        n_disease   = 7,
        n_control   = 6,
        organism    = "Homo sapiens",
        tissue      = "Skeletal muscle",
        data_source = "supplementary_counts"
    )

    prov.record_warning(
        "geo_download",
        "7 samples excluded from analysis (treatment: rhuEPO) — "
        "retained disease vs healthy control comparison only"
    )

    prov.record_parameters("deg_analysis", {
        "padj_threshold":  0.05,
        "lfc_threshold":   1.0,
        "min_count_filter": 10
    })
    prov.record_step("deg_analysis", {
        "n_disease_samples": 7,
        "n_control_samples": 6,
        "n_genes_tested":    14592,
        "n_significant_up":  560,
        "n_significant_down": 619
    })
    prov.record_output_file("deg_results", "data/GSE226646_deseq2_results.csv")

    prov.record_parameters("signature_building", {
        "padj_threshold": 0.05,
        "lfc_threshold":  1.0,
        "max_genes":      150
    })
    prov.record_step("signature_building", {
        "n_up":              68,
        "n_down":            39,
        "signature_quality": "good",
        "top_up":   ["MALAT1", "S100A4", "SYNE2", "TGFB2", "WNK1"],
        "top_down": ["NDUFS7", "ATP5F1A", "IMPA2", "EFHD1", "DECR2"]
    })
    prov.record_output_file("signature_json", "data/GSE226646_signature.json")
    prov.record_output_file("signature_csv",  "data/GSE226646_signature.csv")

    prov.record_step("drug_repurposing", {
        "source":            "l1000fwd",
        "n_candidates":      48,
        "n_high_confidence": 4,
        "result_url": "https://maayanlab.cloud/l1000fwd/#result/699c933ed33550002d5858cd"
    })
    prov.record_output_file(
        "drug_candidates_csv", "data/GSE226646_drug_candidates.csv"
    )

    # Load the real drug candidates if available, otherwise use a stub
    try:
        with open("data/GSE226646_drug_candidates.json") as f:
            drug_candidates = json.load(f)
    except FileNotFoundError:
        drug_candidates = [
            {"drug_name": "PP1 analogue", "score": -11.99, "qvalue": 0.004,
             "confidence": "high", "cell_line": "PC3", "moa": "", "target": ""},
            {"drug_name": "Lovastatin",   "score": -5.55,  "qvalue": 0.171,
             "confidence": "low",  "cell_line": "HA1E", "moa": "", "target": ""}
        ]

    llm_narrative = (
        "This analysis identified 48 drug candidates that reverse the "
        "transcriptomic signature of Friedreich's ataxia in skeletal muscle. "
        "The top hit, a pyrazolopyrimidine compound structurally related to the "
        "Src kinase inhibitor PP1, showed the strongest signature reversal "
        "(combined score -11.99, q=0.004). Notably, Lovastatin also appeared "
        "in the candidate list, consistent with published reports of statin "
        "effects on frataxin expression. Mitochondrial complex I subunits "
        "(NDUFS7, NDUFB3) featured prominently among downregulated genes, "
        "confirming the expected mitochondrial dysfunction signature."
    )

    report_path = prov.generate_report(drug_candidates, llm_narrative)

    print(f"\nProvenance log: {prov.log_path}")
    print(f"Report:         {report_path}")
    print(f"Steps recorded: {prov.steps_completed()}")
    print(f"Warnings:       {len(prov.data['warnings'])}")


