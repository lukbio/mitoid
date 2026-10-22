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
  ],
  plants: [
    { file: "oak_matk.fasta", label: "English oak, matK (889 bp)" },
    { file: "arabidopsis_chloroplast.fasta", label: "thale cress, whole chloroplast genome (154 kb)" },
  ],
};
const GITHUB = "https://github.com/yeezy1218/mitoid";
const PIPELINE = "https://github.com/yeezy1218/mitogenome-pipeline";

function savedMode(): Mode {
  try {
    return localStorage.getItem("mitoid.mode") === "plants" ? "plants" : "animals";
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

  const animal = mode === "animals";
  const markers = animal ? "COX1" : "rbcL and matK";
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
          {(["animals", "plants"] as const).map((m) => (
            <label key={m} className={m === mode ? "on" : undefined}>
              <input type="radio" name="mode" value={m} checked={m === mode} onChange={() => switchMode(m)} />
              {m === "animals" ? "Animals (COX1)" : "Plants (rbcL + matK)"}
            </label>
          ))}
        </fieldset>

        <p className="intro">
          {animal ? (
            <>
              Identify a vertebrate mitochondrial sequence by its COX1 gene. Paste a barcode, part of COX1, or a whole
              mitogenome, in FASTA or plain text.
            </>
          ) : (
            <>
              Identify a land plant or green alga by its chloroplast barcodes, rbcL and matK. Paste either gene, both, or
              a whole chloroplast genome, in FASTA or plain text.
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
              {!animal && shown.length === 1 && results.length > 1 && found(shown[0]) && (
                <> Only {shown[0].gene} was found in your sequence.</>
              )}
            </p>
            {agree && <p className="agree">{agree}</p>}
            {shown.map((r) => (
              <ResultView key={`${r.gene}-${results.indexOf(r)}`} result={r} heading={animal ? "Result" : `${r.gene} result`} />
            ))}
          </div>
        )}

        <section id="about" className="about">
          <h2>About</h2>
          <p>
            MitoID compares your sequence to standard DNA barcode genes. For vertebrates that's COX1, from every
            mitogenome in NCBI RefSeq. Plant mitochondrial genes change too slowly to tell species apart, so for plants it
            uses the chloroplast genes rbcL and matK (the plant barcode chosen by the Consortium for the Barcode of Life),
            from every chloroplast genome in RefSeq. Current reference sets: {speciesCount} species.
          </p>
          <p>
            The search first narrows the reference set down with k-mers. It uses both DNA 11-mers and amino-acid 5-mers
            from a six-frame translation (vertebrate mitochondrial code for COX1, the plastid code for rbcL and matK), since
            related species often differ only at silent third codon positions. The top candidates are then aligned to your
            sequence (Needleman-Wunsch with free end gaps), and a neighbor-joining tree is built from the pairwise
            distances.
          </p>
          <p>
            As a rough rule for vertebrate COX1, 98% identity or more usually means the same species and 92% or more the
            same genus. rbcL is much more conserved: different species in one genus are often identical, so it's best at
            genus and family level, and matK is the better marker for species. Invertebrates and fungi aren't in the
            reference sets. Nothing you paste leaves your computer.
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
