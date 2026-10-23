# MitoID

**Identify an animal, plant or fungus from a DNA sequence by its barcode genes, and see where it sits on a tree, entirely in the browser.**

**Live:** https://yeezy1218.github.io/mitoid/

Three modes:

- **Animals:** COX1, the standard animal barcode, from **15,612 species**: 7,687 vertebrates and 7,925 invertebrates (insects, crustaceans, mollusks, spiders, corals, worms...), from every RefSeq complete mitochondrial genome. Paste a barcode, any stretch of COX1, or a whole mitogenome.
- **Plants:** the chloroplast genes **rbcL (12,495 species) and matK (12,168 species)**, the plant barcode chosen by the Consortium for the Barcode of Life, from every RefSeq complete chloroplast genome. Paste either gene, both, or a whole chloroplast genome. Plant mitochondrial genes change too slowly to separate species, which is why plants use chloroplast markers instead of COX1.
- **Fungi:** ITS, the fungal barcode (nuclear ribosomal ITS1-5.8S-ITS2), from **19,216 species** in NCBI's curated RefSeq ITS collection (one sequence per type specimen).

MitoID screens the sequence against the reference set, aligns it to the closest candidates, and draws a neighbor-joining tree of the query among its relatives. In plant mode it searches both genes and says whether they agree. Nothing is uploaded. A barcode takes about a second, and a whole 154 kb chloroplast genome takes about 3.5 s.

## How it works

1. **Two-stage screen.** The query and all references are reduced to k-mer sets.
   - *DNA:* FracMinHash sketches of canonical 11-mers (both strands, 1 in 4 kept) rank near-identical sequences.
   - *Protein:* the query is translated in all six frames and amino-acid 5-mers are compared with each reference's protein. Most differences between related species are silent third-codon changes, so protein k-mers still find relatives at 80–90% DNA identity, where shared DNA k-mers are mostly noise. On an American toad barcode, the DNA screen alone kept 2 of the true top 10 relatives; with the protein screen, 9 of the top 10 hits are toads.
   - *Genetic codes:* animal mitochondria don't all read DNA the same way. Vertebrates use NCBI table 2, most invertebrates table 5 (AGA/AGG = Ser), echinoderms and flatworms table 9 (AAA = Asn), cnidarians and sponges table 4, tunicates table 13; plastids use table 11. Each reference stores its own code (taken from NCBI Taxonomy when the set is built), and the query is translated once per code and compared with each reference under that reference's code.
   - ITS doesn't code for a protein, so fungi are screened with DNA k-mers only.
2. **Locate.** Exact k-mer hits against the best reference show which strand carries the gene, and where it sits in a longer sequence. In a whole genome, short k-mers also match by chance far from the gene, so the densest gene-length cluster of hits is used.
3. **Align.** Every candidate is aligned to that region with Needleman–Wunsch with free end gaps, so a 658 bp barcode isn't penalised for the parts of the 1.5 kb gene it doesn't cover. Hits are ranked by identity over the overlap. Overlaps under half the shorter sequence are discounted, because a few bases lining up end to end can be "100% identical" by chance.
4. **Place.** Jukes–Cantor distances between the query and its top 10 hits feed a neighbor-joining tree, rooted at its midpoint.

The verdict uses rough rules of thumb that differ by gene. For COX1, ≥98% identity suggests the same species and ≥92% the same genus. rbcL is so conserved that species in one genus are often identical, so its thresholds are stricter and it warns you to check matK. Fungal ITS uses 98.5% for species and 94.3% for genus (Vu et al. 2019). When several species tie for the best match (seven oak species share an identical matK, because oaks hybridize and share chloroplasts), MitoID reports the genus rather than picking one.

Performance note: comparing k-mer sets walks two sorted arrays, which is fine for a barcode. For a whole chloroplast genome (~300,000 protein k-mers) against 12,000 references it took 12 s, so when one set is much larger, the smaller one is binary-searched in it instead. The screen now takes 0.5 s.

All algorithms (k-mer sketching, translation, alignment, neighbor-joining, midpoint rooting, tree layout) are written from scratch in TypeScript in `src/lib/`. They run in a Web Worker so the page stays responsive, and they're covered by unit tests (`npm test`).

## Accuracy

`scripts/benchmark.ts` measures it (`npx vite-node scripts/benchmark.ts`; full table in [docs/benchmark.md](docs/benchmark.md)). It takes 200 random species per group, cuts a barcode-length window from each one's reference gene, adds 1% random substitutions, and searches it. In a second run the species is removed from the reference set first, to see what MitoID says about a species it has never seen.

| Gene | Group | Species in the set: right | Species missing: right genus/family | names a sister species | wrong |
|---|---|---:|---:|---:|---:|
| COX1 | vertebrates | 98.5% | 69.5% | 17.0% | 7.0% |
| COX1 | invertebrates | 98.5% | 70.0% | 10.0% | 6.0% |
| matK | plants | 96.5% | 78.5% | 12.0% | 9.0% |
| rbcL | plants | 93.0% | 83.0% | 2.0% | 14.0% |
| ITS | fungi | 100.0% | 74.5% | 8.5% | 2.5% |

When the species is in the set, MitoID never named the wrong species; the rest were less specific answers (mostly rbcL ties). When it's missing, most mistakes are sister species that really are 98-100% identical in the barcode gene (Humboldt and Galápagos penguins, for example, have identical COX1), so the genus is still right. The remaining "wrong" answers are mostly ties across several genera and genera NCBI has renamed or split. Searches take about 1 s for COX1 and the plant genes and 0.2 s for ITS.

## Reference databases

- `scripts/build_db.py` queries NCBI for RefSeq vertebrate mitogenomes, downloads them in batches (cached), extracts the COX1 coding sequence and the taxonomy (class, order, family) of each species, and writes `public/db/cox1.json.gz` (1.9 MB).
- `scripts/build_plant_db.py` does the same for RefSeq chloroplast genomes. To keep the download small it fetches only coding sequences (`fasta_cds_na`, not full GenBank records), then takes names from esummary and ranked lineages from NCBI Taxonomy. It writes `public/db/rbcl.json.gz` (0.6 MB) and `public/db/matk.json.gz` (1.0 MB).
- `scripts/build_invertebrate_db.py` builds `public/db/cox1_inv.json.gz` (2.2 MB) from RefSeq invertebrate mitogenomes the same way, and stores each species' mitochondrial genetic code from NCBI Taxonomy. The app searches it together with the vertebrate set.
- `scripts/build_fungi_db.py` builds `public/db/its.json.gz` (1.6 MB) from BioProject PRJNA177353.

Only the vertebrate script needs Biopython; the others use plain Python.

```bash
pip install -r scripts/requirements.txt
python scripts/build_db.py --email you@example.org
python scripts/build_plant_db.py --email you@example.org
python scripts/build_invertebrate_db.py --email you@example.org
python scripts/build_fungi_db.py --email you@example.org
```

## Examples

- **Toad COX1 barcode:** the standard 658 bp barcode region of *Anaxyrus americanus*.
- **Fruit fly COX1 barcode:** a 658 bp *Drosophila melanogaster* barcode from GenBank (MW122888.1), identified at 100%.
- **Oyster mushroom ITS:** a 675 bp *Pleurotus ostreatus* ITS from GenBank (QB181443.1), not the RefSeq reference; 99.7% to *P. ostreatus*.
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

- Only species with a complete mitochondrial or chloroplast genome (or a RefSeq ITS) are included. That's most vertebrate families but a small fraction of insects, so many invertebrates will only match to genus or family.
- Sister species can be 98-100% identical in a barcode gene (see the benchmark); no barcode-only method can separate them.
- Plants: ITS2, the nuclear plant barcode, isn't included. Bacteria (16S) and protists aren't covered.
- One sequence per sample: mixed samples (eDNA metabarcoding) and raw reads aren't supported.
- One reference per species, so variation within a species isn't represented.
- The tree is a quick neighbor-joining guide built from one gene, not a published phylogeny.

## License

MIT
