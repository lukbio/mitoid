import { useCallback, useEffect, useRef, useState } from "react";
import { AlignmentView } from "./components/AlignmentView";
import { TreeView } from "./components/TreeView";
import type { Hit } from "./lib/identify";
import { gcPercent, parseInput } from "./lib/seq";
import { useIdentifier } from "./useIdentifier";

const EXAMPLES = [
  { file: "toad_barcode.fasta", label: "American toad, COX1 barcode (658 bp)" },
  { file: "gibbon_mitogenome.fasta", label: "gibbon, whole mitogenome (16.5 kb)" },
];
const GITHUB = "https://github.com/yeezy1218/mitoid";
const PIPELINE = "https://github.com/yeezy1218/mitogenome-pipeline";

const pct = (x: number) => `${(x * 100).toFixed(1)}%`;

/** Tracks an element's content width; the callback ref works for elements
 * that mount later (the tree only exists once there's a result). */
function useWidth<T extends HTMLElement>() {
  const [width, setWidth] = useState(720);
  const observer = useRef<ResizeObserver | null>(null);
  const ref = useCallback((el: T | null) => {
    observer.current?.disconnect();
    if (!el) return;
    observer.current = new ResizeObserver(([e]) => setWidth(Math.round(e.contentRect.width)));
    observer.current.observe(el);
  }, []);
  return [ref, width] as const;
}

export default function App() {
  const { db, loadError, busy, result, error, elapsed, identify, reset } = useIdentifier();
  const [text, setText] = useState("");
  const [selected, setSelected] = useState(0);
  const [treeRef, treeWidth] = useWidth<HTMLDivElement>();
  const parsed = parseInput(text);
  const resultsRef = useRef<HTMLElement>(null);

  useEffect(() => {
    if (result) {
      setSelected(0);
      resultsRef.current?.scrollIntoView({ block: "start" });
    }
  }, [result]);

  const loadExample = async (file: string) => {
    const res = await fetch(`${import.meta.env.BASE_URL}examples/${file}`);
    setText(await res.text());
    reset();
  };

  const submit = (e: React.FormEvent) => {
    e.preventDefault();
    if (parsed.seq) identify(parsed.seq);
  };

  const top: Hit | undefined = result?.hits[0];
  const hit: Hit | undefined = result?.hits[selected];
  const speciesCount = db ? db.count.toLocaleString() : "about 7,700";

  let inputInfo = "";
  if (parsed.seq) {
    inputInfo = `${parsed.seq.length.toLocaleString()} bp, ${gcPercent(parsed.seq).toFixed(1)}% GC`;
    if (parsed.ambiguous) inputInfo += `, ${parsed.ambiguous} ambiguous`;
    if (parsed.dropped) inputInfo += `, ${parsed.dropped} non-DNA characters ignored`;
  }

  let status = "";
  if (loadError) status = `Couldn't load the reference set (${loadError}).`;
  else if (!db) status = "Loading reference set...";
  else if (parsed.seq.length > 0 && parsed.seq.length < 100) status = "Needs at least 100 bp.";

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
        <p className="intro">
          Identify a vertebrate mitochondrial sequence by its COX1 gene. Paste a barcode, part of COX1, or a whole
          mitogenome, in FASTA or plain text. It's compared against COX1 from {speciesCount} species in NCBI RefSeq,
          locally in your browser.
        </p>

        <form onSubmit={submit}>
          <label htmlFor="seq">Sequence</label>
          <textarea id="seq" value={text} onChange={(e) => setText(e.target.value)} spellCheck={false} rows={9} />
          <div className="under-input">
            <span aria-live="polite">{inputInfo}</span>
            <span>
              Example:{" "}
              {EXAMPLES.map((ex, i) => (
                <span key={ex.file}>
                  {i > 0 && " / "}
                  <button type="button" className="linklike" onClick={() => loadExample(ex.file)}>
                    {ex.label}
                  </button>
                </span>
              ))}
            </span>
          </div>
          <button type="submit" className="submit" disabled={!db || busy || parsed.seq.length < 100}>
            {busy ? "Searching..." : "Search"}
          </button>
          <span className="status" aria-live="polite">
            {status}
          </span>
          {error && <p className="error">{error}</p>}
        </form>

        {result && top && hit && (
          <section ref={resultsRef} className="results" aria-label="Results">
            <h2>Result</h2>
            <p className="verdict">
              {result.verdict.level === "none" ? (
                <b>{result.verdict.label}.</b>
              ) : (
                <>
                  Best match: <i className="sp">{top.ref.species}</i>, {pct(top.alignment.identity)} identical over{" "}
                  {top.alignment.length.toLocaleString()} bp.
                </>
              )}{" "}
              {result.verdict.text}
            </p>
            <p className="small muted">
              Searched {result.queryLength.toLocaleString()} bp
              {result.region && (
                <>
                  ; COX1 is at {(result.region.start + 1).toLocaleString()}-{result.region.end.toLocaleString()}
                </>
              )}
              {result.reversed && <> on the reverse strand</>}
              {elapsed !== null && <>; took {(elapsed / 1000).toFixed(1)} s</>}.
            </p>

            <h3>Top {result.hits.length} hits</h3>
            <div className="table-wrap">
              <table className="hits">
                <thead>
                  <tr>
                    <th scope="col"></th>
                    <th scope="col">Species</th>
                    <th scope="col" className="hide-narrow">Family</th>
                    <th scope="col" className="hide-narrow">Accession</th>
                    <th scope="col" className="num">Identity</th>
                    <th scope="col" className="num">Aligned</th>
                  </tr>
                </thead>
                <tbody>
                  {result.hits.map((h, i) => (
                    <tr
                      key={h.ref.acc}
                      className={i === selected ? "selected" : undefined}
                      onClick={() => setSelected(i)}
                      tabIndex={0}
                      onKeyDown={(e) => (e.key === "Enter" || e.key === " ") && setSelected(i)}
                      aria-selected={i === selected}
                    >
                      <td className="num muted">{i + 1}</td>
                      <td>
                        <i>{h.ref.species}</i>
                      </td>
                      <td className="muted hide-narrow">{h.ref.family || h.ref.order || h.ref.class}</td>
                      <td className="hide-narrow">
                        <a href={`https://www.ncbi.nlm.nih.gov/nuccore/${h.ref.acc}`} onClick={(e) => e.stopPropagation()}>
                          {h.ref.acc}
                        </a>
                      </td>
                      <td className="num">{pct(h.alignment.identity)}</td>
                      <td className="num">{h.alignment.length.toLocaleString()} bp</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className="small muted">Click a row to show its alignment below.</p>

            {result.tree && (
              <>
                <h3>Neighbor-joining tree</h3>
                <div ref={treeRef}>
                  <TreeView tree={result.tree} width={treeWidth} />
                </div>
                <p className="small muted">
                  Jukes-Cantor distances, midpoint rooted. It's built from one gene, so treat it as a rough guide.
                </p>
              </>
            )}

            <h3>
              Alignment with <i>{hit.ref.species}</i>
            </h3>
            <p className="small muted">
              {hit.alignment.matches.toLocaleString()} identical, {hit.alignment.mismatches} substitutions,{" "}
              {hit.alignment.gaps} gap positions.
            </p>
            <AlignmentView alignment={hit.alignment} refName={hit.ref.species} offset={result.region?.start ?? 0} />
          </section>
        )}

        <section id="about" className="about">
          <h2>About</h2>
          <p>
            I made this as a side project while working on mitochondrial genome assembly in an undergraduate research lab.
            It compares your sequence to COX1, the gene used for animal DNA barcoding, from every vertebrate mitogenome in
            NCBI RefSeq ({speciesCount} species).
          </p>
          <p>
            The search first narrows the reference set down with k-mers. It uses both DNA 11-mers and amino-acid 5-mers
            from a six-frame translation with the vertebrate mitochondrial code, since related species often differ only
            at silent third codon positions. The top candidates are then aligned to your sequence (Needleman-Wunsch with
            free end gaps), and a neighbor-joining tree is built from the pairwise distances.
          </p>
          <p>
            As a rough rule for vertebrate COX1, 98% identity or more usually means the same species and 92% or more the
            same genus, but it varies between groups. Invertebrates aren't in the reference set. Nothing you paste leaves
            your computer.
          </p>
          <p>
            The gibbon example was assembled from raw Illumina reads with my{" "}
            <a href={PIPELINE}>mitogenome assembly pipeline</a>. The code and the script that builds the reference set
            are on <a href={GITHUB}>GitHub</a>.
          </p>
          <p className="small muted">
            Garrett Dungan
            {db && <>, reference set built {db.built}</>}
          </p>
        </section>
      </main>
    </>
  );
}
