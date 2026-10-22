import { useCallback, useEffect, useRef, useState } from "react";
import { AlignmentView } from "./components/AlignmentView";
import { TreeView } from "./components/TreeView";
import { matchScore, type Hit } from "./lib/identify";
import { gcPercent, parseInput } from "./lib/seq";
import { useIdentifier } from "./useIdentifier";

const EXAMPLES = [
  { file: "toad_barcode.fasta", label: "Toad COX1 barcode", note: "658 bp, the standard barcode region" },
  { file: "gibbon_mitogenome.fasta", label: "Gibbon mitogenome", note: "16.5 kb, assembled from raw reads" },
];
const GITHUB = "https://github.com/yeezy1218/mitoid";
const PIPELINE = "https://github.com/yeezy1218/mitogenome-pipeline";

const pct = (x: number) => `${(x * 100).toFixed(1)}%`;

/** Tracks an element's content width; the callback ref works for elements
 * that mount later (the tree panel only exists once there's a result). */
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
      resultsRef.current?.scrollIntoView({ behavior: "smooth", block: "start" });
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

  const hit: Hit | undefined = result?.hits[selected];

  return (
    <div className="page">
      <header className="masthead">
        <h1>
          Mito<span>ID</span>
        </h1>
        <p className="tagline">
          Paste a vertebrate mitochondrial sequence. MitoID finds its closest species by COX1, the standard animal DNA barcode,
          and places it on a tree, right in your browser.
        </p>
      </header>

      <form className="panel input-panel" onSubmit={submit}>
        <label htmlFor="seq" className="label">
          Sequence <span className="hint">raw or FASTA: a COX1 barcode, a single gene, or a whole mitogenome</span>
        </label>
        <textarea
          id="seq"
          value={text}
          onChange={(e) => setText(e.target.value)}
          placeholder={">my_sample\nACTCTATACCTAATCTTCGGCGCATGAGCTGGAATAGTGGGCACCGCCCTAAGTCTA..."}
          spellCheck={false}
          rows={8}
        />
        <div className="input-meta" aria-live="polite">
          {parsed.seq ? (
            <>
              {parsed.seq.length.toLocaleString()} bases · GC {gcPercent(parsed.seq).toFixed(1)}%
              {parsed.ambiguous > 0 && <> · {parsed.ambiguous} ambiguous</>}
              {parsed.dropped > 0 && <> · {parsed.dropped} non-DNA characters ignored</>}
            </>
          ) : (
            "Nothing pasted yet."
          )}
        </div>
        <div className="actions">
          <button type="submit" className="primary" disabled={!db || busy || parsed.seq.length < 100}>
            {busy ? "Searching…" : "Identify"}
          </button>
          <span className="or">or try</span>
          {EXAMPLES.map((ex) => (
            <button type="button" key={ex.file} className="chip" onClick={() => loadExample(ex.file)}>
              {ex.label} <span className="chip-note">{ex.note}</span>
            </button>
          ))}
        </div>
        <p className="db-status" aria-live="polite">
          {loadError
            ? `Couldn't load the reference set: ${loadError}`
            : db
              ? `Reference set: COX1 from ${db.count.toLocaleString()} vertebrate species (${db.source}, built ${db.built}).`
              : "Loading the reference set…"}
        </p>
        {parsed.seq.length > 0 && parsed.seq.length < 100 && <p className="warn">Paste at least 100 bases.</p>}
        {error && <p className="error">{error}</p>}
      </form>

      {result && hit && (
        <section ref={resultsRef} className="results" aria-label="Results">
          <div className={`panel verdict level-${result.verdict.level}`}>
            <div className="verdict-label">{result.verdict.level === "none" ? "Result" : `Best guess · ${result.verdict.level}`}</div>
            <div className="verdict-name">{result.verdict.label}</div>
            <p>{result.verdict.text}</p>
            <ul className="facts">
              <li>{result.queryLength.toLocaleString()} bases searched</li>
              {result.region && (
                <li>
                  COX1 found at {(result.region.start + 1).toLocaleString()}–{result.region.end.toLocaleString()}
                </li>
              )}
              {result.reversed && <li>read on the reverse strand</li>}
              {elapsed !== null && <li>{(elapsed / 1000).toFixed(2)} s</li>}
            </ul>
          </div>

          <div className="panel">
            <h2>Closest species</h2>
            <div className="table-wrap">
              <table className="hits">
                <thead>
                  <tr>
                    <th scope="col">#</th>
                    <th scope="col">Species</th>
                    <th scope="col">Group</th>
                    <th scope="col" className="num">Identity</th>
                    <th scope="col" className="num">Overlap</th>
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
                      <td>{i + 1}</td>
                      <td>
                        <i>{h.ref.species}</i>{" "}
                        <a href={`https://www.ncbi.nlm.nih.gov/nuccore/${h.ref.acc}`} target="_blank" rel="noreferrer" onClick={(e) => e.stopPropagation()}>
                          {h.ref.acc}
                        </a>
                      </td>
                      <td className="group">{[h.ref.family, h.ref.order, h.ref.class].filter(Boolean)[0]}</td>
                      <td className="num">
                        <span className="bar" style={{ "--w": pct(matchScore(h.alignment)) } as React.CSSProperties} />
                        {pct(h.alignment.identity)}
                      </td>
                      <td className="num">{h.alignment.length.toLocaleString()} bp</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className="note">Click a row to see its alignment below. Identity is measured over the overlapping region.</p>
          </div>

          {result.tree && (
            <div className="panel" ref={treeRef}>
              <h2>Where it sits</h2>
              <TreeView tree={result.tree} width={treeWidth} />
              <p className="note">
                Neighbor-joining tree from Jukes–Cantor distances, midpoint-rooted. A quick guide to relationships, not a
                published phylogeny.
              </p>
            </div>
          )}

          <div className="panel">
            <h2>
              Alignment with <i>{hit.ref.species}</i>
            </h2>
            <p className="note">
              {hit.alignment.matches.toLocaleString()} identical, {hit.alignment.mismatches} substitutions,{" "}
              {hit.alignment.gaps} gap positions. Differences are highlighted.
            </p>
            <AlignmentView alignment={hit.alignment} refName={hit.ref.species} offset={result.region?.start ?? 0} />
          </div>
        </section>
      )}

      <footer className="panel about">
        <h2>How it works</h2>
        <ol>
          <li>
            <b>Screen.</b> Your sequence and all {db ? db.count.toLocaleString() : "7,000+"} references are reduced to
            FracMinHash sketches of 11-base k-mers (both strands), and the 25 references sharing the most are kept.
          </li>
          <li>
            <b>Locate.</b> Exact k-mer matches with the best reference show which strand carries COX1 and, for a long
            sequence like a whole mitogenome, where the gene is.
          </li>
          <li>
            <b>Align.</b> Each candidate is aligned to that region (Needleman–Wunsch with free end gaps) to get exact
            identity.
          </li>
          <li>
            <b>Place.</b> Pairwise distances (Jukes–Cantor) feed a neighbor-joining tree, rooted at its midpoint.
          </li>
        </ol>
        <p>
          Everything runs in your browser; nothing you paste is uploaded. The reference set is built from NCBI RefSeq with
          a Python script, and the same COX1 genes come out of my{" "}
          <a href={PIPELINE} target="_blank" rel="noreferrer">
            mitogenome assembly pipeline
          </a>
          . Identity thresholds (≥98% same species, ≥92% same genus) are rules of thumb for vertebrate COX1 and vary by
          group. <a href={GITHUB} target="_blank" rel="noreferrer">Source on GitHub</a>.
        </p>
      </footer>
    </div>
  );
}
