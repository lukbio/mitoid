"""Build the plant reference sets: rbcL and matK from RefSeq chloroplast genomes.

Usage:
    python scripts/build_plant_db.py --email you@example.org

Plants are identified with chloroplast barcodes (rbcL + matK, the CBOL standard)
rather than COX1, which evolves too slowly in plant mitochondria to tell species
apart. To keep the download small this fetches only each genome's coding
sequences (FASTA, not full GenBank records), then takes names and lineages from
NCBI Taxonomy. Downloads are cached in data/cache_plants/.

Writes public/db/rbcl.json.gz and public/db/matk.json.gz, in the same format as
the COX1 set, with "code": 11 (the plastid/bacterial genetic code).
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_QUERY = ("Viridiplantae[Organism] AND chloroplast[Filter] AND refseq[Filter] "
                 "AND complete genome[Title]")
BATCH = 100
MARKERS = {
    # gene name -> (output file, minimum length of a usable full-length CDS)
    "rbcL": ("rbcl", 1200),
    "matK": ("matk", 1000),
}


def eutil(endpoint: str, email: str, **params) -> bytes:
    data = urlencode({**params, "email": email, "tool": "mitoid"}).encode()
    for attempt in range(5):
        try:
            # POST, so long id lists don't hit URL length limits.
            with urlopen(Request(f"{EUTILS}/{endpoint}.fcgi", data=data), timeout=300) as r:
                return r.read()
        except Exception as e:
            if attempt == 4:
                raise
            print(f"  retrying after {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    raise AssertionError


def search_ids(query: str, email: str) -> list[str]:
    xml = eutil("esearch", email, db="nuccore", term=query, retmax=100000).decode()
    return re.findall(r"<Id>(\d+)</Id>", xml)


def cached(path: Path, fetch) -> str:
    if not path.exists() or path.stat().st_size == 0:
        path.write_bytes(fetch())
        time.sleep(0.4)  # stay under NCBI's 3 requests/second
    return path.read_text()


def parse_cds(fasta: str) -> dict[str, dict[str, str]]:
    """{accession: {gene: sequence}} for the marker genes in a fasta_cds_na dump.
    Headers look like >lcl|NC_000932.1_cds_NP_051067.1_51 [gene=rbcL] [protein=...] ..."""
    out: dict[str, dict[str, str]] = {}
    for block in fasta.split("\n>"):
        header, _, seq = block.lstrip(">").partition("\n")
        m = re.match(r"lcl\|(.+?)_cds_", header)
        g = re.search(r"\[gene=([^\]]+)\]", header)
        if not m or not g:
            continue
        gene = next((k for k in MARKERS if k.lower() == g.group(1).strip().lower()), None)
        if gene is None:
            continue
        seq = re.sub(r"\s+", "", seq).upper()
        acc = m.group(1)
        # Keep the longest copy (a few plastomes list a gene twice).
        if len(seq) > len(out.setdefault(acc, {}).get(gene, "")):
            out[acc][gene] = seq
    return out


def summaries(ids: list[str], email: str, cache: Path) -> dict[str, dict]:
    """{accession: {"species", "taxid"}} from esummary."""
    out = {}
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        js = json.loads(cached(cache / f"sum_{chunk[0]}_{len(chunk)}.json",
                               lambda: eutil("esummary", email, db="nuccore", id=",".join(chunk), retmode="json")))
        for uid in js["result"]["uids"]:
            d = js["result"][uid]
            out[d["accessionversion"]] = {"species": d.get("organism", ""), "taxid": str(d.get("taxid", ""))}
    return out


def lineages(taxids: list[str], email: str, cache: Path) -> dict[str, dict[str, str]]:
    """{taxid: {"class", "order", "family"}} from NCBI Taxonomy's ranked lineage."""
    out = {}
    for i in range(0, len(taxids), 300):
        chunk = taxids[i:i + 300]
        xml = cached(cache / f"tax_{chunk[0]}_{len(chunk)}.xml",
                     lambda: eutil("efetch", email, db="taxonomy", id=",".join(chunk)))
        for taxon in re.findall(r"<Taxon>\s*<TaxId>(\d+)</TaxId>(.*?)</LineageEx>", xml, re.S):
            tid, body = taxon
            ranks = dict((r, n) for n, r in re.findall(r"<ScientificName>([^<]+)</ScientificName>\s*<Rank>([^<]+)</Rank>", body))
            out[tid] = {"class": ranks.get("class", ""), "order": ranks.get("order", ""), "family": ranks.get("family", "")}
    return out


def usable(seq: str, min_len: int) -> bool:
    return len(seq) >= min_len and not (set(seq) - set("ACGTN")) and seq.count("N") < 20


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True)
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--outdir", default="public/db")
    ap.add_argument("--cache", default="data/cache_plants")
    args = ap.parse_args()
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)

    ids = search_ids(args.query, args.email)
    print(f"{len(ids)} chloroplast genomes match", file=sys.stderr)

    genes: dict[str, dict[str, str]] = {}
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        f = cache / f"cds_{chunk[0]}_{chunk[-1]}_{len(chunk)}.fa"
        if not f.exists():
            print(f"  fetching CDS {i + 1}-{i + len(chunk)} of {len(ids)}", file=sys.stderr)
        genes.update(parse_cds(cached(f, lambda: eutil("efetch", email=args.email, db="nuccore", id=",".join(chunk),
                                                       rettype="fasta_cds_na", retmode="text"))))
    info = summaries(ids, args.email, cache)
    tax = lineages(sorted({v["taxid"] for v in info.values() if v["taxid"]}), args.email, cache)

    for gene, (name, min_len) in MARKERS.items():
        entries, seen = [], set()
        for acc, by_gene in genes.items():
            seq = by_gene.get(gene, "")
            meta = info.get(acc)
            if not meta or not usable(seq, min_len) or meta["species"] in seen:
                continue
            seen.add(meta["species"])  # one reference per species
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
            "gene": gene,
            "code": 11,
            "built": date.today().isoformat(),
            "source": "NCBI RefSeq complete chloroplast genomes",
            "query": args.query,
            "count": len(entries),
            "records": entries,
        }
        out = Path(args.outdir) / f"{name}.json.gz"
        out.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(out, "wt", encoding="utf-8", compresslevel=9) as fh:
            json.dump(db, fh, separators=(",", ":"))
        print(f"wrote {out}: {len(entries)} species, {out.stat().st_size / 1e6:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
