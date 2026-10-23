"""Build the fungal reference set: ITS from NCBI's RefSeq Targeted Loci project.

Usage:
    python scripts/build_fungi_db.py --email you@example.org

Fungi are barcoded with ITS (the nuclear ribosomal internal transcribed
spacers, ITS1-5.8S-ITS2), not COX1. NCBI curates one ITS reference per type
specimen in BioProject PRJNA177353 ("Fungal ITS RefSeq Targeted Loci"). ITS is
not protein-coding, so the set has "code": 0 and the app skips the protein
screen for it. Downloads are cached in data/cache_fungi/.
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

from build_plant_db import cached, eutil, lineages, search_ids, summaries

DEFAULT_QUERY = "PRJNA177353[BioProject]"
BATCH = 500
MIN_LEN, MAX_LEN = 300, 1500


def parse_fasta(text: str) -> dict[str, str]:
    out = {}
    for block in text.split("\n>"):
        header, _, seq = block.lstrip(">").partition("\n")
        if header:
            out[header.split()[0]] = re.sub(r"\s+", "", seq).upper()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True)
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--out", default="public/db/its.json.gz")
    ap.add_argument("--cache", default="data/cache_fungi")
    args = ap.parse_args()
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)

    ids = search_ids(args.query, args.email)
    print(f"{len(ids)} ITS records match", file=sys.stderr)
    seqs: dict[str, str] = {}
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        f = cache / f"its_{chunk[0]}_{chunk[-1]}_{len(chunk)}.fa"
        if not f.exists():
            print(f"  fetching {i + 1}-{i + len(chunk)} of {len(ids)}", file=sys.stderr)
        seqs.update(parse_fasta(cached(f, lambda: eutil("efetch", args.email, db="nuccore", id=",".join(chunk),
                                                         rettype="fasta", retmode="text"))))
    info = summaries(ids, args.email, cache)
    tax = lineages(sorted({v["taxid"] for v in info.values() if v["taxid"]}), args.email, cache)

    entries, seen = [], set()
    for acc, seq in seqs.items():
        meta = info.get(acc)
        if not meta or meta["species"] in seen or " sp." in meta["species"]:
            continue
        if not MIN_LEN <= len(seq) <= MAX_LEN or set(seq) - set("ACGTN") or seq.count("N") > 10:
            continue
        seen.add(meta["species"])
        lin = tax.get(meta["taxid"], {})
        entries.append({
            "acc": acc,
            "species": meta["species"],
            "class": lin.get("class") or "Other",
            "order": lin.get("order", ""),
            "family": lin.get("family", ""),
            "seq": seq,
        })
    entries.sort(key=lambda e: (e["class"], e["order"], e["family"], e["species"]))
    db = {
        "gene": "ITS",
        "code": 0,
        "built": date.today().isoformat(),
        "source": "NCBI RefSeq Targeted Loci, fungal ITS (PRJNA177353)",
        "query": args.query,
        "count": len(entries),
        "records": entries,
    }
    out = Path(args.out)
    with gzip.open(out, "wt", encoding="utf-8", compresslevel=9) as fh:
        json.dump(db, fh, separators=(",", ":"))
    print(f"wrote {out}: {len(entries)} species, {out.stat().st_size / 1e6:.1f} MB", file=sys.stderr)
    print("classes:", Counter(e["class"] for e in entries).most_common(12), file=sys.stderr)


if __name__ == "__main__":
    main()
