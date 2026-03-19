#!/usr/bin/env python3
"""Download miRNA-lncRNA and miRNA-mRNA interactions from ENCORI API."""

import urllib.request
import urllib.parse
import time
import os
import sys

PROJECT_DIR = "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
REF_DIR = os.path.join(PROJECT_DIR, "data/reference")

# Our 16 DE-lncRNAs
DE_LNCRNAS = [
    "H19", "LINC01445", "FP671120.6", "LINC01936", "DNM3OS",
    "LINC00707", "ZFPM2-AS1", "LYPLAL1-AS1", "CYP1B1-AS1",
    "AC046143.1", "AC083864.5", "LINC01605", "LINC01711",
    "GRASLND", "AD000090.1", "AL606500.1"
]

BASE_URL = "https://rnasysu.com/encori/api/miRNATarget/"


def query_encori_mirna_target(target, gene_type="lncRNA", clip_exp=0, program_num=0):
    """Query ENCORI API for miRNA-target interactions."""
    params = {
        "assembly": "hg38",
        "geneType": gene_type,
        "miRNA": "all",
        "clipExpNum": str(clip_exp),
        "degraExpNum": "0",
        "pancancerNum": "0",
        "programNum": str(program_num),
        "program": "None",
        "target": target,
        "cellType": "all"
    }
    url = BASE_URL + "?" + urllib.parse.urlencode(params)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read().decode("utf-8")
            lines = [l for l in data.strip().split("\n") if l.strip()]
            if len(lines) > 1:
                return lines  # header + data rows
            return []
    except Exception as e:
        print(f"  Error querying {target}: {e}", file=sys.stderr)
        return []


def query_encori_cerna(lncrna):
    """Query ENCORI ceRNA API for lncRNA-mRNA ceRNA pairs."""
    cerna_url = "https://rnasysu.com/encori/api/ceRNA/"
    params = {
        "assembly": "hg38",
        "geneType": "lncRNA",
        "ceRNA": lncrna,
        "miRNAnum": "1",
        "family": "all",
        "pval": "0.05",
        "fdr": "0.05",
        "pancancerNum": "0"
    }
    url = cerna_url + "?" + urllib.parse.urlencode(params)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read().decode("utf-8")
            lines = [l for l in data.strip().split("\n") if l.strip()]
            if len(lines) > 1:
                return lines
            return []
    except Exception as e:
        print(f"  Error querying ceRNA for {lncrna}: {e}", file=sys.stderr)
        return []


# ── 1. Download miRNA-lncRNA interactions ──
print("=" * 60)
print(" Downloading miRNA-lncRNA interactions from ENCORI")
print("=" * 60)

all_mirna_lncrna = []
header_mirna_lncrna = None

for lncrna in DE_LNCRNAS:
    print(f"  Querying: {lncrna}...", end=" ")

    # Try with relaxed parameters first (clipExpNum=0, programNum=0)
    lines = query_encori_mirna_target(lncrna, gene_type="lncRNA",
                                       clip_exp=0, program_num=0)
    if not lines:
        # Try alternate name format
        alt_name = lncrna.replace(".", "-")
        if alt_name != lncrna:
            lines = query_encori_mirna_target(alt_name, gene_type="lncRNA",
                                               clip_exp=0, program_num=0)
    if lines:
        if header_mirna_lncrna is None:
            header_mirna_lncrna = lines[0]
        data_lines = [l for l in lines[1:] if not l.startswith("miRNAid")]
        all_mirna_lncrna.extend(data_lines)
        print(f"found {len(data_lines)} interactions")
    else:
        print("no data")

    time.sleep(0.5)  # Rate limiting

# Save miRNA-lncRNA interactions
output_file = os.path.join(REF_DIR, "mirna_lncrna_interactions.tsv")
if header_mirna_lncrna and all_mirna_lncrna:
    with open(output_file, "w") as f:
        f.write(header_mirna_lncrna + "\n")
        for line in all_mirna_lncrna:
            f.write(line + "\n")
    print(f"\n  Total miRNA-lncRNA interactions: {len(all_mirna_lncrna)}")
    print(f"  Saved to: {output_file}")
else:
    print("\n  No miRNA-lncRNA interactions found.")

# ── 2. Get ceRNA pairs from ENCORI ──
print("\n" + "=" * 60)
print(" Downloading ceRNA pairs from ENCORI")
print("=" * 60)

all_cerna = []
header_cerna = None

for lncrna in DE_LNCRNAS:
    print(f"  Querying ceRNA: {lncrna}...", end=" ")
    lines = query_encori_cerna(lncrna)
    if lines:
        if header_cerna is None:
            header_cerna = lines[0]
        data_lines = [l for l in lines[1:] if not l.startswith("geneID")]
        all_cerna.extend(data_lines)
        print(f"found {len(data_lines)} ceRNA pairs")
    else:
        print("no data")
    time.sleep(0.5)

# Save ceRNA pairs
cerna_file = os.path.join(REF_DIR, "encori_cerna_pairs.tsv")
if header_cerna and all_cerna:
    with open(cerna_file, "w") as f:
        f.write(header_cerna + "\n")
        for line in all_cerna:
            f.write(line + "\n")
    print(f"\n  Total ceRNA pairs: {len(all_cerna)}")
    print(f"  Saved to: {cerna_file}")
else:
    print("\n  No ceRNA pairs found.")

# ── 3. Download miRNA-mRNA interactions for co-expressed mRNAs ──
print("\n" + "=" * 60)
print(" Downloading miRNA-mRNA interactions (for ceRNA mRNA targets)")
print("=" * 60)

# Get unique mRNA targets from ceRNA results
mrna_targets = set()
if all_cerna:
    for line in all_cerna:
        fields = line.split("\t")
        if len(fields) >= 5:
            mrna_targets.add(fields[1])  # geneName column

# Also get unique miRNAs from lncRNA interactions
mirna_set = set()
if all_mirna_lncrna:
    for line in all_mirna_lncrna:
        fields = line.split("\t")
        if len(fields) >= 2:
            mirna_set.add(fields[1])  # miRNAname

print(f"  Unique miRNAs targeting our lncRNAs: {len(mirna_set)}")
print(f"  Unique mRNA ceRNA partners: {len(mrna_targets)}")

# Query miRNA-mRNA interactions for our specific miRNAs
# Instead of querying each mRNA individually, get bulk miRNA-mRNA data
# for each miRNA that targets our lncRNAs
all_mirna_mrna = []
header_mirna_mrna = None

# Query miRNA-mRNA in bulk (all targets for specific miRNAs)
mirna_list = sorted(mirna_set)
print(f"  Querying miRNA-mRNA interactions for {len(mirna_list)} miRNAs...")

batch_size = 20
for i in range(0, min(len(mirna_list), 200), 1):
    mirna = mirna_list[i]
    params = {
        "assembly": "hg38",
        "geneType": "mRNA",
        "miRNA": mirna,
        "clipExpNum": "1",
        "degraExpNum": "0",
        "pancancerNum": "0",
        "programNum": "1",
        "program": "None",
        "target": "all",
        "cellType": "all"
    }
    url = BASE_URL + "?" + urllib.parse.urlencode(params)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read().decode("utf-8")
            lines = [l for l in data.strip().split("\n") if l.strip()]
            if len(lines) > 1:
                if header_mirna_mrna is None:
                    header_mirna_mrna = lines[0]
                all_mirna_mrna.extend(lines[1:])
    except Exception as e:
        pass

    if (i + 1) % 20 == 0:
        print(f"    Processed {i + 1} / {min(len(mirna_list), 200)} miRNAs...")
    time.sleep(0.3)

# Save miRNA-mRNA interactions
mrna_file = os.path.join(REF_DIR, "mirna_mrna_interactions.tsv")
if header_mirna_mrna and all_mirna_mrna:
    with open(mrna_file, "w") as f:
        f.write(header_mirna_mrna + "\n")
        for line in all_mirna_mrna:
            f.write(line + "\n")
    print(f"\n  Total miRNA-mRNA interactions: {len(all_mirna_mrna)}")
    print(f"  Saved to: {mrna_file}")
else:
    print("\n  No miRNA-mRNA interactions found.")

print("\n" + "=" * 60)
print(" ENCORI download complete.")
print("=" * 60)
