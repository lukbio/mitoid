# MitoID

**Identify an animal or plant DNA sequence by its barcode genes, and see where it sits on a tree, entirely in the browser.**

**Live:** https://yeezy1218.github.io/mitoid/

Two modes:

- **Animals:** COX1, the standard animal barcode, from **7,687 vertebrate species** (every RefSeq complete mitochondrial genome). Paste a barcode, any stretch of COX1, or a whole mitogenome.
- **Plants:** the chloroplast genes **rbcL (12,495 species) and matK (12,168 species)**, the plant barcode chosen by the Consortium for the Barcode of Life, from every RefSeq complete chloroplast genome. Paste either gene, both, or a whole chloroplast genome. Plant mitochondrial genes change too slowly to separate species, which is why plants use chloroplast markers instead of COX1.

MitoID screens the sequence against the reference set, aligns it to the closest candidates, and draws a neighbor-joining tree of the query among its relatives. In plant mode it searches both genes and says whether they agree. Nothing is uploaded. A barcode takes about a second, and a whole 154 kb chloroplast genome takes about 3.5 s.

## How it works

1. **Two-stage screen.** The query and all references are reduced to k-mer sets.
   - *DNA:* FracMinHash sketches of canonical 11-mers (both strands, 1 in 4 kept) rank near-identical sequences.
   - *Protein:* the query is translated in all six frames (vertebrate mitochondrial code, NCBI table 2, for COX1; the plastid code, table 11, for rbcL and matK), and amino-acid 5-mers are compared with each reference's protein. Most differences between related species are silent third-codon changes, so protein k-mers still find relatives at 80–90% DNA identity, where shared DNA k-mers are mostly noise. On an American toad barcode, the DNA screen alone kept 2 of the true top 10 relatives; with the protein screen, 9 of the top 10 hits are toads.
2. **Locate.** Exact k-mer hits against the best reference show which strand carries the gene, and where it sits in a longer sequence. In a whole genome, short k-mers also match by chance far from the gene, so the densest gene-length cluster of hits is used.
3. **Align.** Every candidate is aligned to that region with Needleman–Wunsch with free end gaps, so a 658 bp barcode isn't penalised for the parts of the 1.5 kb gene it doesn't cover. Hits are ranked by identity over the overlap. Overlaps under half the shorter sequence are discounted, because a few bases lining up end to end can be "100% identical" by chance.
4. **Place.** Jukes–Cantor distances between the query and its top 10 hits feed a neighbor-joining tree, rooted at its midpoint.

The verdict uses rough rules of thumb that differ by gene. For COX1, ≥98% identity suggests the same species and ≥92% the same genus. rbcL is so conserved that species in one genus are often identical, so its thresholds are stricter and it warns you to check matK. When several species tie for the best match (seven oak species share an identical matK, because oaks hybridize and share chloroplasts), MitoID reports the genus rather than picking one.

Performance note: comparing k-mer sets walks two sorted arrays, which is fine for a barcode. For a whole chloroplast genome (~300,000 protein k-mers) against 12,000 references it took 12 s, so when one set is much larger, the smaller one is binary-searched in it instead. The screen now takes 0.5 s.

All algorithms (k-mer sketching, translation, alignment, neighbor-joining, midpoint rooting, tree layout) are written from scratch in TypeScript in `src/lib/`. They run in a Web Worker so the page stays responsive, and they're covered by unit tests (`npm test`).

## Reference databases

- `scripts/build_db.py` queries NCBI for RefSeq vertebrate mitogenomes, downloads them in batches (cached), extracts the COX1 coding sequence and the taxonomy (class, order, family) of each species, and writes `public/db/cox1.json.gz` (1.9 MB).
- `scripts/build_plant_db.py` does the same for RefSeq chloroplast genomes. To keep the download small it fetches only coding sequences (`fasta_cds_na`, not full GenBank records), then takes names from esummary and ranked lineages from NCBI Taxonomy. It writes `public/db/rbcl.json.gz` (0.6 MB) and `public/db/matk.json.gz` (1.0 MB). The plant script has no dependencies beyond Python.

```bash
pip install -r scripts/requirements.txt
python scripts/build_db.py --email you@example.org
python scripts/build_plant_db.py --email you@example.org
```

## Examples

- **Toad COX1 barcode:** the standard 658 bp barcode region of *Anaxyrus americanus*.
- **English oak matK:** an 889 bp matK barcode region. Seven *Quercus* species match it equally, so MitoID reports *Quercus* sp.
- **Thale cress chloroplast genome:** the full 154 kb *Arabidopsis thaliana* plastome (NC_000932.1). MitoID finds rbcL and matK inside it (matK on the reverse strand), and both identify *A. thaliana*.
- **Gibbon mitogenome:** a 16.5 kb white-cheeked gibbon mitochondrial genome assembled from raw Illumina reads (SRA run SRR1425539) with my [mitogenome assembly pipeline](https://github.com/yeezy1218/mitogenome-pipeline). MitoID finds COX1 inside it and identifies *Nomascus leucogenys* at 99.9%.

## Development

```bash
npm install
npm run dev      # http://localhost:5173
npm test         # unit tests (Vitest)
npm run build    # production build in dist/
```

Pushing to `main` runs lint, tests and the build, then deploys to GitHub Pages.

## Limits

- Animals: vertebrates only, one gene (COX1). Invertebrates use a different mitochondrial genetic code and aren't in the reference set.
- Plants: land plants and green algae with a complete chloroplast genome in RefSeq. ITS2, the nuclear plant barcode, isn't included. Fungi aren't covered.
- One reference per species, so variation within a species isn't represented.
- The tree is a quick neighbor-joining guide built from one gene, not a published phylogeny.

## License

MIT
