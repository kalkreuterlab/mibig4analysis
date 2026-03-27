#!/usr/bin/env python3
import json
import glob
import os
from collections import defaultdict, Counter
from Bio import SeqIO

# --- Settings ---
json_folder = "mibig_json_4.0"
gbk_folder = "mibig_gbk_4.0"
output_file = "mibig_bgc_gene_summary.tsv"

# --- Step 1: load BGC type from JSONs ---
bgc_info = {}  # BGC accession -> {class, subclass}

for json_file in glob.glob(os.path.join(json_folder, "*.json")):
    with open(json_file) as f:
        data = json.load(f)

    accession = data.get("accession")
    if not accession:
        continue

    # --- BGC class ---
    classes = data.get("biosynthesis", {}).get("classes", [])
    if classes:
        main_class = classes[0].get("class", "Unknown")
        subclass = classes[0].get("subclass", "None")
    else:
        main_class = "Unknown"
        subclass = "None"

    bgc_info[accession] = {
        "class": main_class,
        "subclass": subclass
    }

# --- Step 2: Parse GBK files ---
bgc_gene_data = defaultdict(list)  # BGC accession -> list of gene dicts
bgc_kingdom = {}  # BGC accession -> bacteria/fungi/other

for gbk_file in glob.glob(os.path.join(gbk_folder, "*.gbk")):
    accession = os.path.basename(gbk_file).split(".")[0]

    for record in SeqIO.parse(gbk_file, "genbank"):
        # --- Kingdom inference ---
        if accession not in bgc_kingdom:
            org_name = record.annotations.get("organism", "").lower()
            if "streptomyces" in org_name or "bacillus" in org_name or "escherichia" in org_name:
                bgc_kingdom[accession] = "Bacteria"
            elif "aspergillus" in org_name or "penicillium" in org_name or "fusarium" in org_name:
                bgc_kingdom[accession] = "Fungi"
            else:
                bgc_kingdom[accession] = "Other"

        # --- Extract CDS features ---
        for feature in record.features:
            if feature.type == "CDS":
                gene_kind = feature.qualifiers.get("gene_kind", ["unknown"])[0]
                product = feature.qualifiers.get("product", ["unknown"])[0]

                bgc_gene_data[accession].append({
                    "gene_kind": gene_kind,
                    "product": product
                })

# --- Step 3: Summarize counts by BGC class ---
summary = defaultdict(Counter)
for accession, genes in bgc_gene_data.items():
    cls = bgc_info.get(accession, {}).get("class", "Unknown")
    subclass = bgc_info.get(accession, {}).get("subclass", "None")
    kingdom = bgc_kingdom.get(accession, "Unknown")

    for g in genes:
        summary[(cls, subclass)][("gene_kind", g["gene_kind"])] += 1
        summary[(cls, subclass)][("product", g["product"])] += 1
        summary[(cls, subclass)][("BGC_count", kingdom)] += 1  # keep kingdom prevalence

# --- Step 4: Write TSV ---
all_gene_kinds = set()
all_products = set()
for counts in summary.values():
    for (ftype, val) in counts:
        if ftype == "gene_kind":
            all_gene_kinds.add(val)
        elif ftype == "product":
            all_products.add(val)
all_gene_kinds = sorted(all_gene_kinds)
all_products = sorted(all_products)

with open(output_file, "w") as out:
    header = ["BGC_class", "BGC_subclass"] + \
             [f"gene_kind:{g}" for g in all_gene_kinds] + \
             [f"product:{p}" for p in all_products] + \
             ["Bacteria_count", "Fungi_count", "Other_count"]
    out.write("\t".join(header) + "\n")

    for (cls, subclass), counts in sorted(summary.items()):
        row = [cls, subclass]
        row += [str(counts.get(("gene_kind", g), 0)) for g in all_gene_kinds]
        row += [str(counts.get(("product", p), 0)) for p in all_products]
        row += [str(counts.get(("BGC_count", k), 0)) for k in ["Bacteria", "Fungi", "Other"]]
        out.write("\t".join(row) + "\n")

print(f"✅ TSV summary written to {output_file}")
