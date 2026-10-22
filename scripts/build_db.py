"""Build MitoID's reference database: COX1 from every RefSeq vertebrate mitogenome.

Usage:
    python scripts/build_db.py --email you@example.org [--query "..."] [--out public/db/cox1.json.gz]

Downloads GenBank records from NCBI in batches (cached in data/cache/, so reruns
only fetch what's new), pulls out the COX1 coding sequence and the taxonomy of
each species, and writes one gzipped JSON file the web app loads.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
import time
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_QUERY = ("Vertebrata[Organism] AND mitochondrion[Filter] AND refseq[Filter] "
                 "AND complete genome[Title]")
COX1_NAMES = {"COX1", "CO1", "COI", "COXI", "MTCO1"}
BATCH = 200


def eutil(endpoint: str, email: str, **params) -> bytes:
    url = f"{EUTILS}/{endpoint}.fcgi?" + urlencode({**params, "email": email, "tool": "mitoid"})
    for attempt in range(4):
        try:
            with urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:  # NCBI drops connections now and then
            if attempt == 3:
                raise
            print(f"  retrying after {e}", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    raise AssertionError


def search_ids(query: str, email: str) -> list[str]:
    xml = eutil("esearch", email, db="nuccore", term=query, retmax=100000).decode()
    return re.findall(r"<Id>(\d+)</Id>", xml)


def fetch_batches(ids: list[str], email: str, cache: Path) -> list[Path]:
    cache.mkdir(parents=True, exist_ok=True)
    files = []
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        f = cache / f"batch_{chunk[0]}_{chunk[-1]}_{len(chunk)}.gb"
        if not f.exists() or f.stat().st_size == 0:
            print(f"  fetching {i + 1}-{i + len(chunk)} of {len(ids)}", file=sys.stderr)
            f.write_bytes(eutil("efetch", email, db="nuccore", id=",".join(chunk),
                                rettype="gb", retmode="text"))
            time.sleep(0.4)  # stay under 3 requests/second
        files.append(f)
    return files


def is_cox1(feature) -> bool:
    q = feature.qualifiers
    gene = re.sub(r"[^A-Z0-9]", "", q.get("gene", [""])[0].upper())
    if gene in COX1_NAMES:
        return True
    product = q.get("product", [""])[0].lower()
    return bool(re.search(r"cytochrome c oxidase subunit (i|1)\b", product))


def lineage_rank(taxonomy: list[str], suffixes: tuple[str, ...]) -> str:
    """Pick a rank out of NCBI's unranked lineage by its Latin ending
    (-formes / -idae); good enough for vertebrates."""
    for name in taxonomy:
        if name.endswith(suffixes):
            return name
    return ""


CLASSES = ["Mammalia", "Aves", "Amphibia", "Actinopterygii", "Actinistia", "Chondrichthyes", "Lepidosauria",
           "Testudines", "Crocodylia", "Cladistia", "Dipnoi", "Coelacanthimorpha",
           "Hyperoartia", "Myxini"]


def record_entry(rec) -> dict | None:
    for feat in rec.features:
        if feat.type == "CDS" and is_cox1(feat):
            try:
                seq = str(feat.extract(rec.seq)).upper()
            except UndefinedSequenceError:
                return None  # record came without sequence (e.g. a CON entry)
            if len(seq) < 1200 or set(seq) - set("ACGTN"):
                return None  # partial or ambiguous; skip
            tax = rec.annotations.get("taxonomy", [])
            return {
                "acc": rec.id,
                "species": rec.annotations.get("organism", ""),
                "class": next((c for c in CLASSES if c in tax), "Other"),
                "order": lineage_rank(tax, ("formes",)),
                "family": lineage_rank(tax, ("idae",)),
                "seq": seq,
            }
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True, help="contact address NCBI asks for")
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--out", default="public/db/cox1.json.gz")
    ap.add_argument("--cache", default="data/cache")
    args = ap.parse_args()

    ids = search_ids(args.query, args.email)
    print(f"{len(ids)} records match", file=sys.stderr)
    entries, seen = [], set()
    for f in fetch_batches(ids, args.email, Path(args.cache)):
        for rec in SeqIO.parse(io.StringIO(f.read_text()), "genbank"):
            e = record_entry(rec)
            if e and e["species"] not in seen:  # one record per species
                seen.add(e["species"])
                entries.append(e)
    entries.sort(key=lambda e: (e["class"], e["order"], e["family"], e["species"]))

    db = {
        "gene": "COX1",
        "built": date.today().isoformat(),
        "source": "NCBI RefSeq complete mitochondrial genomes",
        "query": args.query,
        "count": len(entries),
        "records": entries,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8", compresslevel=9) as fh:
        json.dump(db, fh, separators=(",", ":"))
    by_class: dict[str, int] = {}
    for e in entries:
        by_class[e["class"]] = by_class.get(e["class"], 0) + 1
    print(f"wrote {out}: {len(entries)} species, {out.stat().st_size / 1e6:.1f} MB", file=sys.stderr)
    print(", ".join(f"{k} {v}" for k, v in sorted(by_class.items(), key=lambda kv: -kv[1])), file=sys.stderr)


if __name__ == "__main__":
    main()
