import { useCallback, useRef, useState } from "react";
import type { Result } from "../lib/identify";
import { AlignmentView } from "./AlignmentView";
import { TreeView } from "./TreeView";

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

/** Everything found for one marker gene: verdict, hit table, tree, alignment. */
export function ResultView({ result, heading }: { result: Result; heading: string }) {
  const [selected, setSelected] = useState(0);
  const [treeRef, treeWidth] = useWidth<HTMLDivElement>();
  const top = result.hits[0];
  const hit = result.hits[selected] ?? top;
  if (!top || !hit) return null;
  const gene = result.gene;

  return (
    <section className="results" aria-label={heading}>
      <h2>{heading}</h2>
      <p className="verdict">
        {result.verdict.level === "none" ? (
          <b>{result.verdict.label}.</b>
        ) : (
          <>
            Best match: <i className="sp">{result.verdict.label}</i>, {pct(top.alignment.identity)} identical over{" "}
            {top.alignment.length.toLocaleString()} bp.
          </>
        )}{" "}
        {result.verdict.text}
      </p>
      <p className="small muted">
        {result.region && (
          <>
            {/* Report positions on the sequence as pasted; the search worked on the
                reverse complement when the gene was on the other strand. */}
            {gene} is at{" "}
            {result.reversed
              ? `${(result.queryLength - result.region.end + 1).toLocaleString()}-${(result.queryLength - result.region.start).toLocaleString()}`
              : `${(result.region.start + 1).toLocaleString()}-${result.region.end.toLocaleString()}`}
            {result.reversed ? " on the reverse strand" : ""}.
          </>
        )}
        {!result.region && result.reversed && <>Read on the reverse strand.</>}
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
          <h3>Neighbor-joining tree ({gene})</h3>
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
        {result.reversed && " Your positions below count along the reverse complement of what you pasted."}
      </p>
      <AlignmentView alignment={hit.alignment} refName={hit.ref.species} offset={result.region?.start ?? 0} />
    </section>
  );
}
