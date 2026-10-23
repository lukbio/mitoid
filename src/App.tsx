import { useEffect, useRef, useState } from "react";
import { ResultView } from "./components/ResultView";
import type { Result } from "./lib/identify";
import { gcPercent, parseInput } from "./lib/seq";
import { useIdentifier } from "./useIdentifier";
import type { Mode } from "./worker";

const EXAMPLES: Record<Mode, { file: string; label: string }[]> = {
  animals: [
    { file: "toad_barcode.fasta", label: "American toad, COX1 barcode (658 bp)" },
    { file: "gibbon_mitogenome.fasta", label: "gibbon, whole mitogenome (16.5 kb)" },
    { file: "fruitfly_barcode.fasta", label: "fruit fly, COX1 barcode (658 bp)" },
  ],
  plants: [
    { file: "oak_matk.fasta", label: "English oak, matK (889 bp)" },
    { file: "arabidopsis_chloroplast.fasta", label: "thale cress, whole chloroplast genome (154 kb)" },
  ],
  fungi: [{ file: "fungi_its.fasta", label: "oyster mushroom, ITS (675 bp)" }],
};
const MODES: { mode: Mode; label: string; markers: string }[] = [
  { mode: "animals", label: "Animals (COX1)", markers: "COX1" },
  { mode: "plants", label: "Plants (rbcL + matK)", markers: "rbcL and matK" },
  { mode: "fungi", label: "Fungi (ITS)", markers: "ITS" },
];
const GITHUB = "https://github.com/yeezy1218/mitoid";
const PIPELINE = "https://github.com/yeezy1218/mitogenome-pipeline";

function savedMode(): Mode {
  try {
    const m = localStorage.getItem("mitoid.mode");
    return m === "plants" || m === "fungi" ? m : "animals";
  } catch {
    return "animals";
  }
}

/** A marker counts as found when its best hit lines up over at least half the gene. */
const found = (r: Result) => r.verdict.level !== "none" || (r.hits[0]?.alignment.coverage ?? 0) >= 0.5;

/** One sentence on whether two markers agree (plants: rbcL and matK). */
function agreement(results: Result[]): string | null {
  if (results.length < 2) return null;
  const [a, b] = results.map((r) => r.hits[0].ref.species);
  if (a === b) return `Both ${results[0].gene} and ${results[1].gene} point to ${a}.`;
  const [ga, gb] = [a, b].map((s) => s.split(" ")[0]);
  if (ga === gb) return `${results[0].gene} and ${results[1].gene} agree on the genus (${ga}) but not the species.`;
  return `${results[0].gene} and ${results[1].gene} disagree (${a} vs. ${b}); the sequence may be mixed or misidentified.`;
}

export default function App() {
  const [mode, setMode] = useState<Mode>(savedMode);
  const { sets, loadError, busy, results, error, elapsed, identify, reset } = useIdentifier(mode);
  const [text, setText] = useState("");
  const parsed = parseInput(text);
  const resultsRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    try {
      localStorage.setItem("mitoid.mode", mode);
    } catch {
      /* private mode etc.; not important */
    }
  }, [mode]);

  useEffect(() => {
    if (results) resultsRef.current?.scrollIntoView({ block: "start" });
  }, [results]);

  const switchMode = (m: Mode) => {
    setMode(m);
    reset();
  };

  const loadExample = async (file: string) => {
    const res = await fetch(`${import.meta.env.BASE_URL}examples/${file}`);
    setText(await res.text());
    reset();
  };

  const submit = (e: React.FormEvent) => {
    e.preventDefault();
    if (parsed.seq) identify(parsed.seq);
  };

  const markers = MODES.find((m) => m.mode === mode)!.markers;
  const speciesCount = sets ? sets.map((s) => `${s.count.toLocaleString()} (${s.gene})`).join(" and ") : "...";

  let inputInfo = "";
  if (parsed.seq) {
    inputInfo = `${parsed.seq.length.toLocaleString()} bp, ${gcPercent(parsed.seq).toFixed(1)}% GC`;
    if (parsed.ambiguous) inputInfo += `, ${parsed.ambiguous} ambiguous`;
    if (parsed.dropped) inputInfo += `, ${parsed.dropped} non-DNA characters ignored`;
  }

  let status = "";
  if (loadError) status = `Couldn't load the reference set (${loadError}).`;
  else if (!sets) status = "Loading reference set...";
  else if (parsed.seq.length > 0 && parsed.seq.length < 100) status = "Needs at least 100 bp.";

  const shown = results ? (results.filter(found).length ? results.filter(found) : results.slice(0, 1)) : [];
  const agree = agreement(results ? results.filter(found) : []);

  return (
    <>
      <header className="topbar">
        <div className="wrap">
          <span className="brand">MitoID</span>
          <nav>
            <a href="#about">about</a>
            <a href={GITHUB}>source</a>
          </nav>
        </div>
      </header>

      <main className="wrap">
        <fieldset className="mode">
          <legend className="visually-hidden">Organisms</legend>
          {MODES.map(({ mode: m, label }) => (
            <label key={m} className={m === mode ? "on" : undefined}>
              <input type="radio" name="mode" value={m} checked={m === mode} onChange={() => switchMode(m)} />
              {label}
            </label>
          ))}
        </fieldset>

        <p className="intro">
          {mode === "animals" && (
            <>
              Identify an animal mitochondrial sequence, vertebrate or invertebrate, by its COX1 gene. Paste a barcode, part of COX1, or a whole
              mitogenome, in FASTA or plain text.
            </>
          )}
          {mode === "plants" && (
            <>
              Identify a land plant or green alga by its chloroplast barcodes, rbcL and matK. Paste either gene, both, or
              a whole chloroplast genome, in FASTA or plain text.
            </>
          )}
          {mode === "fungi" && (
            <>
              Identify a fungus (mushroom, mold or yeast) by its ITS region, the fungal barcode. Paste ITS1, ITS2 or the
              whole ITS1-5.8S-ITS2 stretch, in FASTA or plain text.
            </>
          )}{" "}
          It's compared against {markers} from NCBI RefSeq, locally in your browser.
        </p>

        <form onSubmit={submit}>
          <label htmlFor="seq">Sequence</label>
          <textarea id="seq" value={text} onChange={(e) => setText(e.target.value)} spellCheck={false} rows={9} />
          <div className="under-input">
            <span aria-live="polite">{inputInfo}</span>
            <span>
              Example:{" "}
              {EXAMPLES[mode].map((ex, i) => (
                <span key={ex.file}>
                  {i > 0 && " / "}
                  <button type="button" className="linklike" onClick={() => loadExample(ex.file)}>
                    {ex.label}
                  </button>
                </span>
              ))}
            </span>
          </div>
          <button type="submit" className="submit" disabled={!sets || busy || parsed.seq.length < 100}>
            {busy ? "Searching..." : "Search"}
          </button>
          <span className="status" aria-live="polite">
            {status}
          </span>
          {error && <p className="error">{error}</p>}
        </form>

        {results && results.length > 0 && (
          <div ref={resultsRef} className="results-all">
            <p className="small muted searched">
              Searched {results[0].queryLength.toLocaleString()} bp against {markers}
              {elapsed !== null && <> in {(elapsed / 1000).toFixed(1)} s</>}.
              {mode === "plants" && shown.length === 1 && results.length > 1 && found(shown[0]) && (
                <> Only {shown[0].gene} was found in your sequence.</>
              )}
            </p>
            {agree && <p className="agree">{agree}</p>}
            {shown.map((r) => (
              <ResultView key={`${r.gene}-${results.indexOf(r)}`} result={r} heading={results.length > 1 ? `${r.gene} result` : "Result"} />
            ))}
          </div>
        )}

        <section id="about" className="about">
          <h2>About</h2>
          <p>
            MitoID compares your sequence to standard DNA barcode genes. For animals that's COX1, from every
            mitogenome in NCBI RefSeq, vertebrate and invertebrate. Plant mitochondrial genes change too slowly to tell species apart, so for plants it
            uses the chloroplast genes rbcL and matK (the plant barcode chosen by the Consortium for the Barcode of Life),
            from every chloroplast genome in RefSeq. Fungi use ITS, a stretch of nuclear ribosomal DNA, from NCBI's
            curated RefSeq ITS collection (one sequence per type specimen). Current reference sets: {speciesCount} species.
          </p>
          <p>
            The search first narrows the reference set down with k-mers. It uses both DNA 11-mers and amino-acid 5-mers
            from a six-frame translation, since related species often differ only at silent third codon positions. Animal
            mitochondria don't all read DNA the same way, so each reference is translated with its own genetic code
            (vertebrates use NCBI table 2, most invertebrates table 5, echinoderms and flatworms table 9, corals and
            sponges table 4), and plastid genes with the bacterial code. ITS doesn't code for a protein, so fungi are
            searched with DNA k-mers only. The top candidates are then aligned to your
            sequence (Needleman-Wunsch with free end gaps), and a neighbor-joining tree is built from the pairwise
            distances.
          </p>
          <p>
            As a rough rule for COX1, 98% identity or more usually means the same species and 92% or more the
            same genus. rbcL is much more conserved: different species in one genus are often identical, so it's best at
            genus and family level, and matK is the better marker for species. For fungal ITS the usual cutoffs are 98.5%
            for species and 94.3% for genus. Only species with a complete genome in
            RefSeq are included, so many insects and other invertebrates will only match to genus or family. Nothing you paste leaves your computer.
          </p>
          <p>
            The gibbon example was assembled from raw Illumina reads with my{" "}
            <a href={PIPELINE}>mitogenome assembly pipeline</a>. The code and the scripts that build the reference sets
            are on <a href={GITHUB}>GitHub</a>.
          </p>
        </section>
      </main>
    </>
  );
}
