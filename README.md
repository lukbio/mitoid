# MitoID

**Identify a vertebrate mitochondrial sequence by its COX1 barcode, and see where it sits on a tree, entirely in the browser.**

**Live:** https://yeezy1218.github.io/mitoid/

Paste a COX1 barcode, any stretch of COX1, or a whole mitogenome. MitoID screens it against COX1 from **7,687 vertebrate species** (every RefSeq complete mitochondrial genome), aligns it to the closest candidates, and draws a neighbor-joining tree of the query among its relatives. Nothing is uploaded; a search takes about a second.

## How it works

1. **Two-stage screen.** The query and all references are reduced to k-mer sets.
   - *DNA:* FracMinHash sketches of canonical 11-mers (both strands, 1 in 4 kept) rank near-identical sequences.
   - *Protein:* the query is translated in all six frames with the vertebrate mitochondrial code (NCBI table 2), and amino-acid 5-mers are compared with each reference's COX1 protein. Most differences between related species are silent third-codon changes, so protein k-mers still find relatives at 80–90% DNA identity, where shared DNA k-mers are mostly noise. On an American toad barcode, the DNA screen alone kept 2 of the true top 10 relatives; with the protein screen, 9 of the top 10 hits are toads.
2. **Locate.** Exact k-mer hits against the best reference show which strand carries COX1, and where the gene sits in a longer sequence such as a whole mitogenome.
3. **Align.** Every candidate is aligned to that region with Needleman–Wunsch with free end gaps, so a 658 bp barcode isn't penalised for the parts of the 1.5 kb gene it doesn't cover. Hits are ranked by identity over the overlap. Overlaps under half the shorter sequence are discounted, because a few bases lining up end to end can be "100% identical" by chance.
4. **Place.** Jukes–Cantor distances between the query and its top 10 hits feed a neighbor-joining tree, rooted at its midpoint.

The verdict uses common rules of thumb for vertebrate COX1: ≥98% identity suggests the same species, ≥92% the same genus.

All algorithms (k-mer sketching, translation, alignment, neighbor-joining, midpoint rooting, tree layout) are written from scratch in TypeScript in `src/lib/`. They run in a Web Worker so the page stays responsive, and they're covered by unit tests (`npm test`).

## Reference database

`scripts/build_db.py` queries NCBI for RefSeq vertebrate mitogenomes, downloads them in batches (cached), extracts the COX1 coding sequence and the taxonomy (class, order, family) of each species, and writes `public/db/cox1.json.gz` (1.9 MB).

```bash
pip install -r scripts/requirements.txt
python scripts/build_db.py --email you@example.org
```

## Examples

- **Toad COX1 barcode:** the standard 658 bp barcode region of *Anaxyrus americanus*.
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

- Vertebrates only, and one gene (COX1). Invertebrates use a different genetic code and aren't in the reference set.
- One reference per species, so variation within a species isn't represented.
- The tree is a quick neighbor-joining guide built from one gene, not a published phylogeny.

## License

MIT
