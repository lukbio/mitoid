import { useMemo } from "react";
import { QUERY_NAME } from "../lib/identify";
import { layout, type TreeNode } from "../lib/tree";

const ROW = 26;
const LABEL_W = 230;
const PAD = 12;

/** Rectangular phylogram of the query among its closest references. */
export function TreeView({ tree, width = 720 }: { tree: TreeNode; width?: number }) {
  const { leaves, lines, depth } = useMemo(() => layout(tree), [tree]);
  const plotW = Math.max(120, width - LABEL_W - 2 * PAD);
  const sx = (x: number) => PAD + (depth ? (x / depth) * plotW : 0);
  const sy = (y: number) => PAD + y * ROW + ROW / 2;
  const height = leaves.length * ROW + 2 * PAD + 28;

  // Scale bar: a round number near a fifth of the tree depth.
  const raw = depth / 5 || 0.01;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const bar = [1, 2, 5, 10].map((m) => m * mag).find((v) => v >= raw) ?? raw;

  return (
    <svg
      className="tree"
      viewBox={`0 0 ${width} ${height}`}
      role="img"
      aria-label={`Neighbor-joining tree of your sequence and ${leaves.length - 1} closest references`}
    >
      <g className="tree-lines">
        {lines.map((l, i) => (
          <line key={i} x1={sx(l.x1)} y1={sy(l.y1)} x2={sx(l.x2)} y2={sy(l.y2)} />
        ))}
      </g>
      {leaves.map((leaf) => {
        const me = leaf.name === QUERY_NAME;
        return (
          <g key={leaf.name} className={me ? "tree-leaf me" : "tree-leaf"}>
            <circle cx={sx(leaf.x)} cy={sy(leaf.y)} r={me ? 5 : 3} />
            <text x={sx(leaf.x) + 9} y={sy(leaf.y)} dominantBaseline="middle" fontStyle={me ? "normal" : "italic"}>
              {leaf.name}
            </text>
          </g>
        );
      })}
      <g className="tree-scale" transform={`translate(${PAD}, ${height - 18})`}>
        <line x1={0} x2={(bar / (depth || 1)) * plotW} y1={0} y2={0} />
        <text x={(bar / (depth || 1)) * plotW + 8} y={0} dominantBaseline="middle">
          {bar} substitutions per site
        </text>
      </g>
    </svg>
  );
}
