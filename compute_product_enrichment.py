#!/usr/bin/env python3
"""
compute_product_enrichment.py

Compute subclass-specific relative enrichment of product annotations from the TSV
produced by summarize_mibig_bgc_genes.py and write results to an Excel workbook.

Example:
    python compute_product_enrichment.py \
        --input mibig_bgc_gene_summary.tsv \
        --output top_products_by_subclass.xlsx \
        --top-n 20 \
        --min-count 10
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import pandas as pd


LOGGER = logging.getLogger("compute_product_enrichment")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Identify product annotations overrepresented within each BGC subclass "
            "relative to all other subclasses."
        )
    )
    parser.add_argument(
        "--input",
        default="mibig_bgc_gene_summary.tsv",
        help="Input TSV from summarize_mibig_bgc_genes.py (default: mibig_bgc_gene_summary.tsv).",
    )
    parser.add_argument(
        "--output",
        default="top_products_by_subclass.xlsx",
        help="Output Excel workbook (default: top_products_by_subclass.xlsx).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of enriched products to report per subclass (default: 20).",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help="Minimum subclass count required for a product to be considered (default: 10).",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-6,
        help="Small constant added to denominator proportions to avoid division by zero (default: 1e-6).",
    )
    parser.add_argument(
        "--unknown-label",
        default="Unknown",
        help="Label to use for missing subclass values (default: Unknown).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s: %(message)s",
    )


def sanitize_sheet_name(name: str) -> str:
    """
    Remove or replace characters invalid in Excel sheet names.
    """
    name = str(name).strip()
    name = re.sub(r'[:\\/?*\[\]]', "_", name)
    return name[:31] or "Sheet"


def compute_enrichment(
    df: pd.DataFrame,
    top_n: int,
    min_count: int,
    epsilon: float,
    unknown_label: str,
) -> dict[str, pd.DataFrame]:
    if "BGC_subclass" not in df.columns:
        raise ValueError("Input TSV must contain a 'BGC_subclass' column.")

    df = df.copy()
    df["BGC_subclass"] = df["BGC_subclass"].fillna(unknown_label).replace("", unknown_label)

    product_cols = [col for col in df.columns if col.startswith("product:")]
    if not product_cols:
        raise ValueError("No product columns found. Expected columns beginning with 'product:'.")

    results: dict[str, pd.DataFrame] = {}

    for subclass in sorted(df["BGC_subclass"].unique()):
        subclass_df = df[df["BGC_subclass"] == subclass]
        other_df = df[df["BGC_subclass"] != subclass]

        sub_counts = subclass_df[product_cols].sum(numeric_only=True)
        other_counts = other_df[product_cols].sum(numeric_only=True)

        sub_counts = sub_counts[sub_counts >= min_count]
        if sub_counts.empty:
            LOGGER.warning(
                "Skipping subclass '%s' because no products met min_count=%d",
                subclass,
                min_count,
            )
            continue

        other_counts = other_counts.reindex(sub_counts.index).fillna(0)

        total_sub = float(sub_counts.sum()) or 1.0
        total_other = float(other_counts.sum()) or 1.0

        enrichment = (sub_counts / total_sub) / ((other_counts / total_other) + epsilon)

        result_df = pd.DataFrame(
            {
                "Product": [idx.replace("product:", "") for idx in enrichment.index],
                "Subclass_Count": sub_counts.values,
                "Other_Count": other_counts.values,
                "Relative_Enrichment": enrichment.values,
            }
        ).sort_values("Relative_Enrichment", ascending=False)

        results[str(subclass)] = result_df.head(top_n).reset_index(drop=True)

    return results


def write_excel(results: dict[str, pd.DataFrame], output_path: Path) -> None:
    if not results:
        raise ValueError("No enrichment results were generated.")

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary_rows = []
        for subclass, result_df in results.items():
            sheet_name = sanitize_sheet_name(subclass)
            result_df.to_excel(writer, sheet_name=sheet_name, index=False)
            summary_rows.append({"BGC_subclass": subclass, "Reported_Products": len(result_df)})

        pd.DataFrame(summary_rows).sort_values("BGC_subclass").to_excel(
            writer,
            sheet_name="Summary",
            index=False,
        )


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input TSV not found: {input_path}")

    df = pd.read_csv(input_path, sep="\t")
    results = compute_enrichment(
        df=df,
        top_n=args.top_n,
        min_count=args.min_count,
        epsilon=args.epsilon,
        unknown_label=args.unknown_label,
    )

    output_path = Path(args.output)
    write_excel(results, output_path)
    LOGGER.info("Wrote enrichment workbook to %s", output_path.resolve())


if __name__ == "__main__":
    main()
