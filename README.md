# Analysis of Gene Content and Product Enrichment in MIBiG 4.0 Biosynthetic Gene Clusters

## Overview

This repository contains Python scripts for analyzing biosynthetic gene clusters (BGCs) from the MIBiG 4.0 database. The workflow parses GenBank and JSON files to summarize gene annotations and product distributions across BGC classes and subclasses, and identifies products that are overrepresented within specific subclasses.

This analysis is intended to support the identification of patterns in biosynthetic gene content and highlight potentially informative or unusual enzymatic functions in natural product biosynthesis.

## Workflow

The analysis consists of two main steps:

### 1. BGC Gene and Product Summary

Script:
summarize_mibig_bgc_genes.py

This script:
- Parses MIBiG JSON files to extract BGC class and subclass annotations
- Parses GenBank files to extract CDS-level annotations (gene_kind, product)
- Aggregates counts across BGC classes and subclasses

Output:
mibig_bgc_gene_summary.tsv

### 2. Product Enrichment Analysis

Script:
compute_product_enrichment.py

This script:
- Takes the summary TSV as input
- Computes relative enrichment of products within each subclass compared to all other subclasses
- Reports the top enriched products per subclass

Output:
top_products_by_subclass.xlsx

## Requirements

Python ≥ 3.8

Packages:
pandas
biopython
openpyxl

Install:
pip install pandas biopython openpyxl

## Input Data

Requires MIBiG 4.0 dataset:

mibig_json_4.0/
mibig_gbk_4.0/

## Usage

python summarize_mibig_bgc_genes.py --json-folder mibig_json_4.0 --gbk-folder mibig_gbk_4.0 --output mibig_bgc_gene_summary.tsv

python compute_product_enrichment.py --input mibig_bgc_gene_summary.tsv --output top_products_by_subclass.xlsx

## Notes

- BGC classification uses first JSON class entry
- Kingdom assignment is heuristic
- Product annotations come from GenBank and may vary
- Products normalized by default
- Default min_count = 10

## Reproducibility

Tested on Python 3.x (Linux/WSL)
