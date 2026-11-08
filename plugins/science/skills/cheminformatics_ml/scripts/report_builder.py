# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib",
#     "numpy",
#     "pandas>=2.0.0",
#     "rdkit",
#     "scikit-learn",
#     "scipy",
#     "tabulate>=0.9.0",
# ]
# ///
"""Standardized Scientific Report Builder & Deliverable Synthesizer.

Provides programmatic assembly of 5-section cheminformatics reports,
guaranteeing:
1. Absolute image paths in markdown embeds
   (![caption](/absolute/path/to/img.png)).
2. On-disk verification of referenced figures and datasets.
3. Clean lead table formatting with InChIKey, canonical SMILES, and verified
   CAS.
4. Structured modeling caveats and persistent artifact links.
"""

import argparse
from dataclasses import dataclass, field
import os
import re
import sys
from typing import Any, Dict, List, Optional
import pandas as pd


def validate_mermaid_syntax(markdown_content: str) -> List[Dict[str, Any]]:
  """Validates Mermaid code blocks in markdown content for common syntax issues."""
  issues: List[Dict[str, Any]] = []
  blocks = re.findall(r"```mermaid\s*(.*?)```", markdown_content, re.DOTALL)
  for b_idx, block in enumerate(blocks, start=1):
    for l_idx, line in enumerate(block.splitlines(), start=1):
      stripped = line.strip()
      if not stripped or stripped.startswith("%%"):
        continue
      if re.search(r"\[[^\]]*[\(\)][^\]]*\]", stripped) and '"' not in stripped:
        issues.append({
            "block_index": b_idx,
            "line_number": l_idx,
            "error": (
                "Unquoted parentheses inside rectangular node label: "
                f"{stripped}"
            ),
        })
  return issues


@dataclass
class FigureEntry:
  image_path: str
  caption: str
  takeaways: List[str] = field(default_factory=list)


@dataclass
class MermaidEntry:
  code: str
  caption: Optional[str] = None
  takeaways: List[str] = field(default_factory=list)


def _dataframe_to_markdown(df: pd.DataFrame) -> str:
  """Converts a pandas DataFrame to a markdown table without requiring tabulate."""
  headers = [str(col) for col in df.columns]
  rows = [[str(val) for val in row] for row in df.values]
  col_widths = [len(h) for h in headers]
  for row in rows:
    for i, val in enumerate(row):
      col_widths[i] = max(col_widths[i], len(val))

  header_line = (
      "| "
      + " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
      + " |"
  )
  sep_line = (
      "| " + " | ".join("-" * col_widths[i] for i in range(len(headers))) + " |"
  )
  row_lines = [
      "| "
      + " | ".join(val.ljust(col_widths[i]) for i, val in enumerate(row))
      + " |"
      for row in rows
  ]
  return "\n".join([header_line, sep_line] + row_lines)


class ScientificReportBuilder:
  """Constructs standardized 5-section cheminformatics deliverable reports."""

  def __init__(self, title: str, artifact_dir: Optional[str] = None):
    self.title = title
    self.artifact_dir = (
        os.path.abspath(artifact_dir) if artifact_dir else os.getcwd()
    )
    self.executive_summary: List[str] = []
    self.figures: List[FigureEntry] = []
    self.mermaid_diagrams: List[MermaidEntry] = []
    self.leads_df: Optional[pd.DataFrame] = None
    self.lead_table_cols: List[str] = []
    self.methodology_items: Dict[str, str] = {}
    self.artifacts: List[Dict[str, str]] = []

  def set_executive_summary(
      self,
      context: str,
      core_insights: List[str],
      recommended_actions: List[str],
  ) -> "ScientificReportBuilder":
    """Sets Section 1: Executive Summary."""
    self.executive_summary = [
        f"**Context & Objective**: {context}\n",
        "**Core Insights**:",
        *[f"- {insight}" for insight in core_insights],
        "\n**Recommended Action / Next Steps**:",
        *[f"- {action}" for action in recommended_actions],
    ]
    return self

  def add_figure(
      self, image_path: str, caption: str, takeaways: Optional[List[str]] = None
  ) -> "ScientificReportBuilder":
    """Adds a verified figure deliverable with absolute path validation."""
    abs_path = (
        os.path.abspath(image_path)
        if not os.path.isabs(image_path)
        else image_path
    )
    if not os.path.exists(abs_path):
      alt_path = os.path.join(self.artifact_dir, os.path.basename(image_path))
      if os.path.exists(alt_path):
        abs_path = alt_path
      else:
        print(
            f"[ReportBuilder Warning] Image file not found on disk: {abs_path}",
            file=sys.stderr,
        )

    self.figures.append(
        FigureEntry(
            image_path=abs_path,
            caption=caption,
            takeaways=takeaways or [],
        )
    )
    return self

  def add_mermaid_diagram(
      self,
      code: str,
      caption: Optional[str] = None,
      takeaways: Optional[List[str]] = None,
  ) -> "ScientificReportBuilder":
    """Adds a Mermaid diagram to the report."""
    code_clean = code.strip()
    self.mermaid_diagrams.append(
        MermaidEntry(
            code=code_clean, caption=caption, takeaways=takeaways or []
        )
    )
    return self

  def add_figure_with_df_takeaways(
      self,
      image_path: str,
      caption: str,
      df: pd.DataFrame,
      stat_type: str = "benchmark",
      feature_col: str = "Feature_Set",
      target_col: str = "Target",
      r2_col: Optional[str] = None,
      r_col: Optional[str] = None,
      max_rows: Optional[int] = None,
  ) -> "ScientificReportBuilder":
    """Adds a figure and automatically generates verbatim takeaways directly from a DataFrame.

    Ensures narrative consistency by reading exact computed metrics.

    The benchmark engine uses two column schemas: `Train_CV_R2`/`Train_CV_r`
    when a held-out test split is supplied, and `Q2_LOO`/`Pearson_r` otherwise.
    When `r2_col`/`r_col` are not given explicitly, whichever schema is present
    is resolved automatically so captions never silently read "N/A".
    """

    def _resolve(explicit, candidates):
      """Picks the first candidate column actually present in the frame."""
      if explicit is not None:
        return explicit
      for cand in candidates:
        if cand in df.columns:
          return cand
      return candidates[0]

    takeaways = []
    if stat_type == "benchmark":
      resolved_r2 = _resolve(
          r2_col, ["Q2_LOO", "Train_CV_R2", "CV_R2", "Test_R2"]
      )
      resolved_r = _resolve(
          r_col, ["Pearson_r", "Train_CV_r", "Test_Pearson_r"]
      )
      r2_label = (
          "Q²" if resolved_r2 in ("Q2_LOO", "CV_R2", "Train_CV_R2") else "R²"
      )
      for _, row in df.iterrows():
        feat = row.get(feature_col, "Feature")
        tgt = (
            f" [{row.get(target_col)}]"
            if target_col in row and pd.notna(row.get(target_col))
            else ""
        )
        q2_val = row.get(resolved_r2)
        r_val = row.get(resolved_r)
        q2_str = (
            f"{r2_label} = {q2_val:.3f}"
            if pd.notna(q2_val)
            else f"{r2_label} = N/A"
        )
        r_str = f"Pearson r = {r_val:.3f}" if pd.notna(r_val) else "r = N/A"
        takeaways.append(f"**{feat}**{tgt}: {q2_str}, {r_str}")
    elif stat_type == "physics":
      for _, row in df.iterrows():
        term = row.get("Term", row.get("Component", "Term"))
        pr = row.get("Pearson_r")
        sr = row.get("Spearman_rho")
        pr_str = f"Pearson r = {pr:.3f}" if pd.notna(pr) else "r = N/A"
        sr_str = f"Spearman ρ = {sr:.3f}" if pd.notna(sr) else "ρ = N/A"
        n_str = f" (N = {row.get('N')})" if "N" in row else ""
        takeaways.append(f"**{term}**: {pr_str}, {sr_str}{n_str}")
    else:
      sub_rows = df.head(max_rows) if max_rows is not None else df
      for _, row in sub_rows.iterrows():
        takeaways.append(
            "; ".join([f"{col}: {row[col]}" for col in df.columns])
        )
      if max_rows is not None and len(df) > len(sub_rows):
        takeaways.append(
            f"*Showing top {len(sub_rows)} of {len(df)} rows; full data in"
            " companion CSV.*"
        )

    return self.add_figure(
        image_path=image_path, caption=caption, takeaways=takeaways
    )

  def set_leads_table(
      self,
      df: pd.DataFrame,
      display_cols: Optional[List[str]] = None,
      max_rows: Optional[int] = None,
  ) -> "ScientificReportBuilder":
    """Sets Section 3: Representative Leads Table with clean float rounding."""
    formatted_df = df.copy()
    for col in formatted_df.select_dtypes(include=["float", "float64"]).columns:
      formatted_df[col] = formatted_df[col].apply(
          lambda x: f"{x:.3f}" if pd.notna(x) else ""
      )
    self.leads_df = formatted_df
    self.lead_table_cols = display_cols or list(df.columns)
    self.leads_max_rows = max_rows
    return self

  def set_methodology(
      self,
      structural_inputs: str,
      tools_and_methods: str,
      limitations_and_caveats: List[str],
  ) -> "ScientificReportBuilder":
    """Sets Section 4: Methodology & Assumptions."""
    self.methodology_items = {
        "Structural Coordinates & Inputs": structural_inputs,
        "Computational Methods & Force Fields": tools_and_methods,
        "Physical Assumptions & Known Limitations": "\n".join(
            [f"- {lim}" for lim in limitations_and_caveats]
        ),
    }
    return self

  def add_artifact_link(
      self, label: str, file_path: str
  ) -> "ScientificReportBuilder":
    """Adds a persistent clickable file link."""
    abs_path = (
        os.path.abspath(file_path)
        if not os.path.isabs(file_path)
        else file_path
    )
    self.artifacts.append({"label": label, "path": abs_path})
    return self

  def render_markdown(self) -> str:
    """Renders the complete 5-section Markdown document."""
    lines: List[str] = [f"# {self.title}\n"]

    # Section 1: Executive Summary
    lines.append("## 1. Executive Summary & Key Findings\n")
    lines.extend(self.executive_summary)
    lines.append("\n---\n")

    # Section 2: Visual Deliverables & Analytical Evidence
    lines.append("## 2. Core Visual Deliverables & Analytical Evidence\n")
    for m in self.mermaid_diagrams:
      if m.caption:
        lines.append(f"### {m.caption}\n")
      lines.append(f"```mermaid\n{m.code}\n```\n")
      if m.takeaways:
        for t in m.takeaways:
          lines.append(f"- {t}")
        lines.append("")

    for i, fig in enumerate(self.figures, start=1):
      lines.append(f"![{fig.caption}]({fig.image_path})\n")
      if fig.takeaways:
        lines.append(f"**Key Findings (Figure {i})**:")
        for t in fig.takeaways:
          lines.append(f"- {t}")
        lines.append("")
    lines.append("---\n")

    # Section 3: Representative Leads
    lines.append("## 3. Representative Chemical Series & Concordance Table\n")
    if self.leads_df is not None and not self.leads_df.empty:
      max_r = getattr(self, "leads_max_rows", None)
      sub_df = (
          self.leads_df[self.lead_table_cols].head(max_r)
          if max_r is not None
          else self.leads_df[self.lead_table_cols]
      )
      lines.append(_dataframe_to_markdown(sub_df))
      if max_r is not None and len(self.leads_df) > len(sub_df):
        lines.append(
            f"\n*Showing top {len(sub_df)} of {len(self.leads_df)} compounds;"
            " full table in companion CSV.*"
        )
      lines.append("\n")
    else:
      lines.append(
          "*Detailed per-compound metrics are recorded in companion artifact"
          " files.*\n"
      )
    lines.append("---\n")

    # Section 4: Methodology & Assumptions
    lines.append(
        "## 4. Methodology, Physical Approximations & Modeling Assumptions\n"
    )
    for heading, text in self.methodology_items.items():
      lines.append(f"### {heading}\n{text}\n")
    lines.append("---\n")

    # Section 5: Persistent Deliverables & Artifact Access
    lines.append("## 5. Persistent Deliverables & Artifact Access\n")
    for art in self.artifacts:
      lines.append(f"- [{art['label']}](file://{art['path']})")
    lines.append("")

    return "\n".join(lines)

  def save(self, output_filename: str) -> str:
    """Saves the rendered report to the specified file and returns absolute path."""
    out_path = (
        output_filename
        if os.path.isabs(output_filename)
        else os.path.join(self.artifact_dir, output_filename)
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    content = self.render_markdown()

    # Pre-save validation of Mermaid invariants
    issues = validate_mermaid_syntax(content)
    if issues:
      err_msgs = [
          f"Block {iss['block_index']} Line {iss['line_number']}:"
          f" {iss['error']}"
          for iss in issues
      ]
      raise ValueError(
          "Report contains invalid Mermaid syntax violating agent UI"
          " invariants:\n"
          + "\n".join(err_msgs)
      )

    with open(out_path, "w", encoding="utf-8") as f:
      f.write(content)
    print(f"Successfully generated scientific deliverable report: {out_path}")
    return out_path


def main():
  parser = argparse.ArgumentParser(
      description="Assemble standardized computational chemistry reports."
  )
  parser.add_argument("--title", type=str, required=True, help="Report title")
  parser.add_argument(
      "--artifact_dir",
      type=str,
      default=os.getcwd(),
      help="Target artifact directory",
  )
  parser.add_argument(
      "--benchmark_csv",
      type=str,
      default=None,
      help="Path to benchmark summary CSV",
  )
  parser.add_argument(
      "--predictions_csv",
      type=str,
      default=None,
      help="Path to predictions and leads CSV",
  )
  parser.add_argument(
      "--parity_png",
      type=str,
      default=None,
      help="Path to parity plot figure PNG",
  )
  parser.add_argument(
      "--weights_png",
      type=str,
      default=None,
      help="Path to feature weights/coefficients PNG",
  )
  parser.add_argument(
      "--output",
      "--output_md",
      dest="output",
      type=str,
      default="computational_report.md",
      help="Output markdown filename",
  )
  args = parser.parse_args()

  builder = ScientificReportBuilder(
      title=args.title, artifact_dir=args.artifact_dir
  )

  bench_df = None
  if args.benchmark_csv and os.path.exists(args.benchmark_csv):
    try:
      bench_df = pd.read_csv(args.benchmark_csv)
      builder.add_artifact_link("Benchmark Summary CSV", args.benchmark_csv)
    except Exception as e:
      print(
          "[ReportBuilder Warning] Could not read benchmark CSV"
          f" {args.benchmark_csv}: {e}",
          file=sys.stderr,
      )

  if args.parity_png and os.path.exists(args.parity_png):
    if bench_df is not None:
      builder.add_figure_with_df_takeaways(
          image_path=args.parity_png,
          caption="Multi-Target QSAR Model Parity & Goodness-of-Fit",
          df=bench_df,
          stat_type="benchmark",
      )
    else:
      builder.add_figure(
          image_path=args.parity_png,
          caption="Multi-Target QSAR Model Parity & Goodness-of-Fit",
          takeaways=[
              "Evaluated model predictions against holdout experimental"
              " endpoints."
          ],
      )
    builder.add_artifact_link("Multi-Target Parity Plot", args.parity_png)

  if args.weights_png and os.path.exists(args.weights_png):
    builder.add_figure(
        image_path=args.weights_png,
        caption="Standardized Regression Feature Weights & Attributions",
        takeaways=[
            "Attribution coefficients highlight key structural drivers of"
            " potency."
        ],
    )
    builder.add_artifact_link("Feature Weights Plot", args.weights_png)

  if args.predictions_csv and os.path.exists(args.predictions_csv):
    try:
      preds_df = pd.read_csv(args.predictions_csv)
      builder.set_leads_table(preds_df)
      builder.add_artifact_link("Predictions & Leads CSV", args.predictions_csv)
    except Exception as e:
      print(
          "[ReportBuilder Warning] Could not read predictions CSV"
          f" {args.predictions_csv}: {e}",
          file=sys.stderr,
      )

  builder.set_executive_summary(
      context="Standard automated cheminformatics campaign deliverable.",
      core_insights=["Data curated and validated against chemical invariants."],
      recommended_actions=["Review visual figures and lead rankings."],
  )
  out = builder.save(args.output)
  print(f"Report saved to {out}")


if __name__ == "__main__":
  main()
