#!/usr/bin/env python3
"""
summarize_mibig_bgc_genes.py

Parse MIBiG JSON and GenBank files to summarize gene_kind and product counts
across BGC classes/subclasses.

Example:
    python summarize_mibig_bgc_genes.py \
        --json-folder mibig_json_4.0 \
        --gbk-folder mibig_gbk_4.0 \
        --output mibig_bgc_gene_summary.tsv
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from Bio import SeqIO


LOGGER = logging.getLogger("summarize_mibig_bgc_genes")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize gene_kind and product annotations from MIBiG GenBank files, "
            "grouped by BGC class and subclass from companion JSON files."
        )
    )
    parser.add_argument(
        "--json-folder",
        default="mibig_json_4.0",
        help="Folder containing MIBiG JSON files (default: mibig_json_4.0).",
    )
    parser.add_argument(
        "--gbk-folder",
        default="mibig_gbk_4.0",
        help="Folder containing MIBiG GenBank files (default: mibig_gbk_4.0).",
    )
    parser.add_argument(
        "--output",
        default="mibig_bgc_gene_summary.tsv",
        help="Output TSV filename (default: mibig_bgc_gene_summary.tsv).",
    )
    parser.add_argument(
        "--keep-product-case",
        action="store_true",
        help="Keep original product capitalization instead of lowercasing/stripping.",
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


def validate_input_folder(path_str: str, label: str) -> Path:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(f"{label} is not a directory: {path}")
    return path


def normalize_text(value: str) -> str:
    return " ".join(str(value).strip().split()).lower()


def infer_kingdom(organism_name: str) -> str:
    """
    Heuristic kingdom inference from organism name.

    This is intentionally conservative and should be described as heuristic in any
    downstream README/manuscript text.
    """
    org = organism_name.lower()

    bacteria_markers = (
        "streptomyces",
        "bacillus",
        "escherichia",
        "pseudomonas",
        "salinispora",
        "micromonospora",
        "kitasatospora",
        "actinoplanes",
        "actinosynnema",
        "myxococcus",
        "cyanobacter",
    )
    fungi_markers = (
        "aspergillus",
        "penicillium",
        "fusarium",
        "candida",
        "saccharomyces",
        "trichoderma",
        "claviceps",
        "monascus",
    )

    if any(marker in org for marker in bacteria_markers):
        return "Bacteria"
    if any(marker in org for marker in fungi_markers):
        return "Fungi"
    return "Other"


def load_bgc_info(json_folder: Path) -> Dict[str, Dict[str, str]]:
    bgc_info: Dict[str, Dict[str, str]] = {}
    json_files = sorted(json_folder.glob("*.json"))
    if not json_files:
        raise FileNotFoundError(f"No JSON files found in {json_folder}")

    for json_file in json_files:
        with json_file.open("r", encoding="utf-8") as handle:
            data = json.load(handle)

        accession = data.get("accession")
        if not accession:
            LOGGER.warning("Skipping JSON with no accession: %s", json_file.name)
            continue

        classes = data.get("biosynthesis", {}).get("classes", [])
        if classes:
            main_class = classes[0].get("class") or "Unknown"
            subclass = classes[0].get("subclass") or "None"
        else:
            main_class = "Unknown"
            subclass = "None"

        bgc_info[accession] = {"class": main_class, "subclass": subclass}

    LOGGER.info("Loaded class/subclass annotations for %d accessions", len(bgc_info))
    return bgc_info


def parse_genbank_features(
    gbk_folder: Path, normalize_products: bool = True
) -> Tuple[Dict[str, List[Dict[str, str]]], Dict[str, str]]:
    bgc_gene_data: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    bgc_kingdom: Dict[str, str] = {}

    gbk_files = sorted(gbk_folder.glob("*.gbk"))
    if not gbk_files:
        raise FileNotFoundError(f"No GenBank files found in {gbk_folder}")

    for gbk_file in gbk_files:
        accession = gbk_file.stem

        for record in SeqIO.parse(str(gbk_file), "genbank"):
            if accession not in bgc_kingdom:
                organism = record.annotations.get("organism", "")
                bgc_kingdom[accession] = infer_kingdom(organism)

            for feature in record.features:
                if feature.type != "CDS":
                    continue

                gene_kind = feature.qualifiers.get("gene_kind", ["unknown"])[0]
                product = feature.qualifiers.get("product", ["unknown"])[0]

                gene_kind = " ".join(str(gene_kind).strip().split()) or "unknown"
                product = " ".join(str(product).strip().split()) or "unknown"
                if normalize_products:
                    product = normalize_text(product)

                bgc_gene_data[accession].append(
                    {"gene_kind": gene_kind, "product": product}
                )

    LOGGER.info("Parsed CDS annotations from %d GenBank files", len(gbk_files))
    return bgc_gene_data, bgc_kingdom


def build_summary(
    bgc_info: Dict[str, Dict[str, str]],
    bgc_gene_data: Dict[str, List[Dict[str, str]]],
    bgc_kingdom: Dict[str, str],
) -> defaultdict:
    summary = defaultdict(Counter)

    for accession, genes in bgc_gene_data.items():
        cls = bgc_info.get(accession, {}).get("class", "Unknown")
        subclass = bgc_info.get(accession, {}).get("subclass", "None")
        kingdom = bgc_kingdom.get(accession, "Unknown")

        # Count each accession once per class/subclass row for kingdom prevalence.
        summary[(cls, subclass)][("BGC_count", kingdom)] += 1

        for gene in genes:
            summary[(cls, subclass)][("gene_kind", gene["gene_kind"])] += 1
            summary[(cls, subclass)][("product", gene["product"])] += 1

    return summary


def collect_feature_levels(summary: defaultdict) -> Tuple[List[str], List[str]]:
    gene_kinds = set()
    products = set()

    for counts in summary.values():
        for feature_type, value in counts:
            if feature_type == "gene_kind":
                gene_kinds.add(value)
            elif feature_type == "product":
                products.add(value)

    return sorted(gene_kinds), sorted(products)


def write_summary_tsv(
    output_path: Path,
    summary: defaultdict,
    all_gene_kinds: List[str],
    all_products: List[str],
) -> None:
    with output_path.open("w", encoding="utf-8") as out:
        header = (
            ["BGC_class", "BGC_subclass"]
            + [f"gene_kind:{g}" for g in all_gene_kinds]
            + [f"product:{p}" for p in all_products]
            + ["Bacteria_count", "Fungi_count", "Other_count"]
        )
        out.write("\t".join(header) + "\n")

        for (cls, subclass), counts in sorted(summary.items()):
            row = [cls, subclass]
            row += [str(counts.get(("gene_kind", g), 0)) for g in all_gene_kinds]
            row += [str(counts.get(("product", p), 0)) for p in all_products]
            row += [
                str(counts.get(("BGC_count", k), 0))
                for k in ["Bacteria", "Fungi", "Other"]
            ]
            out.write("\t".join(row) + "\n")


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    json_folder = validate_input_folder(args.json_folder, "JSON folder")
    gbk_folder = validate_input_folder(args.gbk_folder, "GenBank folder")
    output_path = Path(args.output)

    bgc_info = load_bgc_info(json_folder)
    bgc_gene_data, bgc_kingdom = parse_genbank_features(
        gbk_folder, normalize_products=not args.keep_product_case
    )
    summary = build_summary(bgc_info, bgc_gene_data, bgc_kingdom)
    all_gene_kinds, all_products = collect_feature_levels(summary)
    write_summary_tsv(output_path, summary, all_gene_kinds, all_products)

    LOGGER.info("Wrote summary TSV to %s", output_path.resolve())


if __name__ == "__main__":
    main()
