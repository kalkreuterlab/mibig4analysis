#!/usr/bin/env python3
import pandas as pd
import re

# --- Settings ---
tsv_file = "mibig_bgc_gene_summary_organized.tsv"
output_file = "top_products_by_subclass.xlsx"
epsilon = 1e-6      # small value to avoid division by zero
min_count = 10        # optional: filter products with very few counts

# --- Helper function to sanitize sheet names ---
def sanitize_sheet_name(name):
    """
    Remove or replace characters invalid in Excel sheet names.
    """
    name = str(name).strip()
    name = re.sub(r'[:\\/?*\[\]]', '_', name)
    return name[:30]  # Excel sheet name limit

# --- Load TSV ---
df = pd.read_csv(tsv_file, sep="\t")

# Fill missing subclass values
df['BGC_subclass'] = df['BGC_subclass'].fillna("Unknown")

# Identify product columns
product_cols = [c for c in df.columns if c.startswith("product:")]

# --- Compute top 20 overrepresented products per subclass ---
results = {}

for subclass in df['BGC_subclass'].unique():
    sub_df = df[df['BGC_subclass'] == subclass]
    other_df = df[df['BGC_subclass'] != subclass]

    # Sum product counts
    sub_counts = sub_df[product_cols].sum()
    other_counts = other_df[product_cols].sum()

    # Optional filter: ignore very low counts
    sub_counts = sub_counts[sub_counts >= min_count]

    # Avoid division by zero
    total_sub = sub_counts.sum() if sub_counts.sum() > 0 else 1
    total_other = other_counts.sum() if other_counts.sum() > 0 else 1

    # Compute relative enrichment
    enrichment = (sub_counts / total_sub) / ((other_counts / total_other) + epsilon)

    # Take top 20
    top20 = enrichment.sort_values(ascending=False).head(20)
    results[subclass] = top20

# --- Write to Excel ---
try:
    import openpyxl
except ImportError:
    raise ImportError("Please install openpyxl: pip install openpyxl")

with pd.ExcelWriter(output_file) as writer:
    for subclass, top20 in results.items():
        # Remove 'product:' prefix for readability
        top20.index = [p.replace("product:", "") for p in top20.index]
        # Sanitize sheet name
        sheet_name = sanitize_sheet_name(subclass)
        top20.to_frame(name="Relative_Enrichment").to_excel(writer, sheet_name=sheet_name)

print(f"✅ Excel file written to {output_file}")
