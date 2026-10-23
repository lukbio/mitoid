"""Build the invertebrate COX1 reference set from RefSeq mitogenomes.

Usage:
    python scripts/build_invertebrate_db.py --email you@example.org

Same output format as the vertebrate set (public/db/cox1.json.gz), but each
record carries its own genetic code: invertebrate mitochondria don't share one.
Most use table 5, echinoderms and flatworms table 9, cnidarians and sponges
table 4, and a few others something else. The code comes from NCBI Taxonomy
(each taxon's mitochondrial genetic code), so the app translates every
reference the way the animal itself does.

Fetches only coding sequences (fasta_cds_na), like build_plant_db.py, and
caches downloads in data/cache_inv/.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path

from build_plant_db import cached, eutil, search_ids, summaries, usable

DEFAULT_QUERY = ("Metazoa[Organism] NOT Vertebrata[Organism] AND mitochondrion[Filter] "
                 "AND refseq[Filter] AND complete genome[Title]")
BATCH = 100
MIN_LEN = 1200
COX1_NAMES = {"COX1", "CO1", "COI", "COXI", "MTCO1", "COXA"}


def parse_cox1(fasta: str) -> dict[str, str]:
    """{accession: COX1 sequence} from a fasta_cds_na dump. Many invertebrate
    records name the gene only through [protein=...], so check both."""
    out: dict[str, str] = {}
    for block in fasta.split("\n>"):
        header, _, seq = block.lstrip(">").partition("\n")
        m = re.match(r"lcl\|(.+?)_cds_", header)
        if not m:
            continue
        g = re.search(r"\[gene=([^\]]+)\]", header)
        p = re.search(r"\[protein=([^\]]+)\]", header)
        gene = re.sub(r"[^A-Z0-9]", "", g.group(1).upper()) if g else ""
        protein = p.group(1).lower() if p else ""
        if gene not in COX1_NAMES and not re.search(r"cytochrome c oxidase subunit (i|1)\b", protein):
            continue
        seq = re.sub(r"\s+", "", seq).upper()
        if len(seq) > len(out.get(m.group(1), "")):
            out[m.group(1)] = seq
    return out


def lineages(taxids: list[str], email: str, cache: Path) -> dict[str, dict]:
    """{taxid: {"phylum", "class", "order", "family", "code"}}; code is the
    taxon's mitochondrial genetic code."""
    out = {}
    for i in range(0, len(taxids), 300):
        chunk = taxids[i:i + 300]
        xml = cached(cache / f"tax_{chunk[0]}_{len(chunk)}.xml",
                     lambda: eutil("efetch", email, db="taxonomy", id=",".join(chunk)))
        # Split per top-level taxon; the lineage inside also has <Taxon> tags.
        for body in re.split(r"\n<Taxon>\n", xml)[1:]:
            tid = re.search(r"<TaxId>(\d+)</TaxId>", body)
            code = re.search(r"<MGCId>(\d+)</MGCId>", body)
            if not tid:
                continue
            ranks = dict((r, n) for n, r in re.findall(
                r"<ScientificName>([^<]+)</ScientificName>\s*<Rank>([^<]+)</Rank>", body))
            out[tid.group(1)] = {
                "phylum": ranks.get("phylum", ""), "class": ranks.get("class", ""),
                "order": ranks.get("order", ""), "family": ranks.get("family", ""),
                "code": int(code.group(1)) if code else 5,
            }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True)
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--out", default="public/db/cox1_inv.json.gz")
    ap.add_argument("--cache", default="data/cache_inv")
    args = ap.parse_args()
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)

    ids = search_ids(args.query, args.email)
    print(f"{len(ids)} invertebrate mitogenomes match", file=sys.stderr)
    cox1: dict[str, str] = {}
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        f = cache / f"cds_{chunk[0]}_{chunk[-1]}_{len(chunk)}.fa"
        if not f.exists():
            print(f"  fetching CDS {i + 1}-{i + len(chunk)} of {len(ids)}", file=sys.stderr)
        cox1.update(parse_cox1(cached(f, lambda: eutil("efetch", args.email, db="nuccore", id=",".join(chunk),
                                                       rettype="fasta_cds_na", retmode="text"))))
    info = summaries(ids, args.email, cache)
    tax = lineages(sorted({v["taxid"] for v in info.values() if v["taxid"]}), args.email, cache)

    entries, seen = [], set()
    for acc, seq in cox1.items():
        meta = info.get(acc)
        if not meta or not usable(seq, MIN_LEN) or meta["species"] in seen:
            continue
        seen.add(meta["species"])
        lin = tax.get(meta["taxid"], {})
        entries.append({
            "acc": acc,
            "species": meta["species"],
            # Some invertebrate groups have no class in NCBI; the phylum still says what it is.
            "class": lin.get("class") or lin.get("phylum") or "Other",
            "order": lin.get("order", ""),
            "family": lin.get("family", ""),
            "code": lin.get("code", 5),
            "seq": seq,
        })
    entries.sort(key=lambda e: (e["class"], e["order"], e["family"], e["species"]))
    db = {
        "gene": "COX1",
        "code": 5,
        "built": date.today().isoformat(),
        "source": "NCBI RefSeq complete invertebrate mitochondrial genomes",
        "query": args.query,
        "count": len(entries),
        "records": entries,
    }
    out = Path(args.out)
    with gzip.open(out, "wt", encoding="utf-8", compresslevel=9) as fh:
        json.dump(db, fh, separators=(",", ":"))
    print(f"wrote {out}: {len(entries)} species, {out.stat().st_size / 1e6:.1f} MB", file=sys.stderr)
    print("codes:", dict(Counter(e["code"] for e in entries)), file=sys.stderr)
    print("classes:", Counter(e["class"] for e in entries).most_common(15), file=sys.stderr)


if __name__ == "__main__":
    main()
