import type { Alignment } from "../lib/align";

const WIDTH = 60;

/** Pairwise alignment in blocks of 60 columns; differences are highlighted.
 * `offset` shifts your positions when only a window of your sequence was aligned. */
export function AlignmentView({ alignment, refName, offset = 0 }: { alignment: Alignment; refName: string; offset?: number }) {
  const { aRow, bRow, aStart, bStart } = alignment;
  const blocks = [];
  let aPos = aStart + offset;
  let bPos = bStart;
  for (let i = 0; i < aRow.length; i += WIDTH) {
    const a = aRow.slice(i, i + WIDTH);
    const b = bRow.slice(i, i + WIDTH);
    const aFrom = aPos + 1;
    const bFrom = bPos + 1;
    aPos += a.replace(/-/g, "").length;
    bPos += b.replace(/-/g, "").length;
    blocks.push(
      <div className="aln-block" key={i}>
        <div className="aln-row">
          <span className="aln-label">You</span>
          <span className="aln-pos">{aFrom}</span>
          <code>
            {[...a].map((c, k) => (
              <span key={k} className={c === b[k] ? undefined : c === "-" || b[k] === "-" ? "gap" : "diff"}>
                {c}
              </span>
            ))}
          </code>
        </div>
        <div className="aln-row">
          <span className="aln-label" title={refName}>Ref</span>
          <span className="aln-pos">{bFrom}</span>
          <code>
            {[...b].map((c, k) => (
              <span key={k} className={c === a[k] ? undefined : c === "-" || a[k] === "-" ? "gap" : "diff"}>
                {c}
              </span>
            ))}
          </code>
        </div>
      </div>,
    );
  }
  return (
    <div className="aln" aria-label={`Alignment of your sequence with ${refName}`}>
      {blocks}
    </div>
  );
}
