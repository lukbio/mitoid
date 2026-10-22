/** Neighbor-joining trees, midpoint rooting and a simple rectangular layout. */

export interface TreeNode {
  name?: string;
  /** Branch length to the parent. */
  length: number;
  children: TreeNode[];
}

/** Saitou & Nei neighbor-joining on a symmetric distance matrix. Unrooted
 * result returned with the last three nodes joined at a central node. */
export function neighborJoining(names: string[], dist: number[][]): TreeNode {
  if (names.length === 0) throw new Error("no taxa");
  let nodes: TreeNode[] = names.map((name) => ({ name, length: 0, children: [] }));
  if (nodes.length === 1) return nodes[0];
  let d = dist.map((row) => row.slice());

  while (nodes.length > 2) {
    const n = nodes.length;
    const r = d.map((row) => row.reduce((s, x) => s + x, 0));
    let bi = 0;
    let bj = 1;
    let bq = Infinity;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const q = (n - 2) * d[i][j] - r[i] - r[j];
        if (q < bq) {
          bq = q;
          bi = i;
          bj = j;
        }
      }
    }
    const li = Math.max(0, 0.5 * d[bi][bj] + (r[bi] - r[bj]) / (2 * (n - 2)));
    const lj = Math.max(0, d[bi][bj] - li);
    const joined: TreeNode = {
      length: 0,
      children: [
        { ...nodes[bi], length: li },
        { ...nodes[bj], length: lj },
      ],
    };
    const keep = nodes.map((_, k) => k).filter((k) => k !== bi && k !== bj);
    const newRow = keep.map((k) => Math.max(0, 0.5 * (d[bi][k] + d[bj][k] - d[bi][bj])));
    d = keep.map((k, x) => [...keep.map((k2) => d[k][k2]), newRow[x]]);
    d.push([...newRow, 0]);
    nodes = [...keep.map((k) => nodes[k]), joined];
  }
  const [a, b] = nodes;
  const len = d[0][1];
  return { length: 0, children: [{ ...a, length: len / 2 }, { ...b, length: len / 2 }] };
}

interface Edge {
  to: number;
  length: number;
}

/** Re-root at the midpoint of the longest leaf-to-leaf path. */
export function midpointRoot(tree: TreeNode): TreeNode {
  // Flatten to an undirected graph.
  const names: (string | undefined)[] = [];
  const adj: Edge[][] = [];
  const add = (node: TreeNode, parent: number, length: number): void => {
    const id = names.length;
    names.push(node.name);
    adj.push([]);
    if (parent >= 0) {
      adj[id].push({ to: parent, length });
      adj[parent].push({ to: id, length });
    }
    for (const c of node.children) add(c, id, c.length);
  };
  add(tree, -1, 0);
  const leaves = adj.map((_, i) => i).filter((i) => names[i] !== undefined);
  if (leaves.length < 3) return tree;

  const farthest = (start: number) => {
    const dist = new Map<number, number>([[start, 0]]);
    const prev = new Map<number, number>();
    const stack = [start];
    while (stack.length) {
      const u = stack.pop()!;
      for (const e of adj[u]) {
        if (!dist.has(e.to)) {
          dist.set(e.to, dist.get(u)! + e.length);
          prev.set(e.to, u);
          stack.push(e.to);
        }
      }
    }
    let best = start;
    for (const l of leaves) if (dist.get(l)! > dist.get(best)!) best = l;
    return { best, dist, prev };
  };
  const a = farthest(leaves[0]).best;
  const { best: b, dist, prev } = farthest(a);
  const half = dist.get(b)! / 2;

  // Walk from b back towards a until we pass the halfway point.
  let v = b;
  while (prev.has(v) && dist.get(prev.get(v)!)! >= half) v = prev.get(v)!;
  const u = prev.get(v)!; // the edge u-v straddles the midpoint
  const fromU = half - dist.get(u)!;
  const edgeLen = dist.get(v)! - dist.get(u)!;

  const build = (node: number, from: number, length: number): TreeNode => ({
    name: names[node],
    length,
    children: adj[node].filter((e) => e.to !== from).map((e) => build(e.to, node, e.length)),
  });
  return {
    length: 0,
    children: [build(u, v, Math.max(0, fromU)), build(v, u, Math.max(0, edgeLen - fromU))],
  };
}

export interface LaidOutLeaf {
  name: string;
  x: number;
  y: number;
}
export interface LaidOutLine {
  x1: number;
  y1: number;
  x2: number;
  y2: number;
}

/** Rectangular phylogram: x = distance from root, y = leaf order (0..n-1). */
export function layout(tree: TreeNode): { leaves: LaidOutLeaf[]; lines: LaidOutLine[]; depth: number } {
  const leaves: LaidOutLeaf[] = [];
  const lines: LaidOutLine[] = [];
  let depth = 0;
  const walk = (node: TreeNode, x: number): number => {
    depth = Math.max(depth, x);
    if (!node.children.length) {
      const y = leaves.length;
      leaves.push({ name: node.name ?? "", x, y });
      return y;
    }
    const ys = node.children.map((c) => {
      const cy = walk(c, x + c.length);
      lines.push({ x1: x, y1: cy, x2: x + c.length, y2: cy });
      return cy;
    });
    lines.push({ x1: x, y1: Math.min(...ys), x2: x, y2: Math.max(...ys) });
    return (Math.min(...ys) + Math.max(...ys)) / 2;
  };
  walk(tree, 0);
  return { leaves, lines, depth };
}
