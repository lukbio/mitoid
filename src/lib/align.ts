/**
 * Semi-global pairwise alignment (Needleman–Wunsch with free end gaps).
 *
 * Leading and trailing gaps cost nothing, so a short barcode fragment aligns
 * inside the full-length reference (and vice versa) without being penalised
 * for the parts the other sequence doesn't cover. Identity is measured over
 * the overlapping region only, which is how barcode identity is usually quoted.
 */

const MATCH = 2;
const MISMATCH = -3;
const GAP = -5;

export interface Alignment {
  /** Identical columns / aligned columns in the overlap (gaps inside it count). */
  identity: number;
  /** Columns in the overlap, including internal gaps. */
  length: number;
  /** Overlap as a fraction of the shorter sequence. A tiny end-to-end overlap
   * can be 100% identical by chance, so low coverage means "not a real match". */
  coverage: number;
  matches: number;
  mismatches: number;
  gaps: number;
  /** Where the overlap lies on each sequence (0-based, end exclusive). */
  aStart: number;
  aEnd: number;
  bStart: number;
  bEnd: number;
  /** Gapped strings of the overlap, for display. */
  aRow: string;
  bRow: string;
}

const DIAG = 1;
const UP = 2;
const LEFT = 3;

export function align(a: string, b: string): Alignment {
  const n = a.length;
  const m = b.length;
  const cols = m + 1;
  const trace = new Uint8Array((n + 1) * cols);
  let prev = new Int32Array(cols); // first row 0: free leading gaps in a
  let cur = new Int32Array(cols);
  let best = -Infinity;
  let bestI = 0;
  let bestJ = 0;

  for (let i = 1; i <= n; i++) {
    cur[0] = 0; // free leading gaps in b
    const ai = a.charCodeAt(i - 1);
    for (let j = 1; j <= m; j++) {
      const bj = b.charCodeAt(j - 1);
      const s = ai === bj && ai !== 78 /* N */ ? MATCH : MISMATCH;
      let v = prev[j - 1] + s;
      let t = DIAG;
      const up = prev[j] + GAP;
      if (up > v) {
        v = up;
        t = UP;
      }
      const left = cur[j - 1] + GAP;
      if (left > v) {
        v = left;
        t = LEFT;
      }
      cur[j] = v;
      trace[i * cols + j] = t;
    }
    // free trailing gaps: the alignment may end on the last column...
    if (cur[m] > best) {
      best = cur[m];
      bestI = i;
      bestJ = m;
    }
    [prev, cur] = [cur, prev];
  }
  // ...or on the last row (prev holds row n after the final swap).
  for (let j = 1; j <= m; j++) {
    if (prev[j] > best) {
      best = prev[j];
      bestI = n;
      bestJ = j;
    }
  }

  let i = bestI;
  let j = bestJ;
  let matches = 0;
  let mismatches = 0;
  let gaps = 0;
  const aRow: string[] = [];
  const bRow: string[] = [];
  while (i > 0 && j > 0) {
    const t = trace[i * cols + j];
    if (t === DIAG) {
      const x = a[i - 1];
      const y = b[j - 1];
      if (x === y && x !== "N") matches++;
      else mismatches++;
      aRow.push(x);
      bRow.push(y);
      i--;
      j--;
    } else if (t === UP) {
      gaps++;
      aRow.push(a[i - 1]);
      bRow.push("-");
      i--;
    } else {
      gaps++;
      aRow.push("-");
      bRow.push(b[j - 1]);
      j--;
    }
  }
  const length = matches + mismatches + gaps;
  return {
    identity: length ? matches / length : 0,
    length,
    coverage: Math.min(n, m) ? Math.min(bestI - i, bestJ - j) / Math.min(n, m) : 0,
    matches,
    mismatches,
    gaps,
    aStart: i,
    aEnd: bestI,
    bStart: j,
    bEnd: bestJ,
    aRow: aRow.reverse().join(""),
    bRow: bRow.reverse().join(""),
  };
}

/** Jukes–Cantor corrected distance from a proportion of differing sites. */
export function jukesCantor(pDistance: number): number {
  if (pDistance <= 0) return 0;
  const x = 1 - (4 / 3) * pDistance;
  return x > 0 ? -0.75 * Math.log(x) : 3; // saturated: cap instead of Infinity
}
