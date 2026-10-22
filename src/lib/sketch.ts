/**
 * FracMinHash sketches of canonical k-mers.
 *
 * Every k-mer (and its reverse complement, whichever encodes smaller) is
 * hashed; only k-mers whose hash falls under 1/SCALE of the hash space are
 * kept. Two sketches then share roughly 1/SCALE of the k-mers their full sets
 * share, so overlap between sketches estimates overlap between sequences,
 * without storing every k-mer of ~7,700 references.
 */

export const K = 11;
export const SCALE = 4;
const MASK = (1 << (2 * K)) - 1;
const THRESHOLD = Math.floor(0x100000000 / SCALE);
const CODE: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

/** 32-bit integer mix (murmur3 finalizer). */
export function mix32(x: number): number {
  x ^= x >>> 16;
  x = Math.imul(x, 0x85ebca6b);
  x ^= x >>> 13;
  x = Math.imul(x, 0xc2b2ae35);
  x ^= x >>> 16;
  return x >>> 0;
}

/** Calls fn(canonicalCode, forwardCode, position) for every k-mer without N. */
export function forEachKmer(seq: string, fn: (canon: number, fwd: number, pos: number) => void) {
  let fwd = 0;
  let rev = 0;
  let valid = 0;
  const shift = 2 * (K - 1);
  for (let i = 0; i < seq.length; i++) {
    const c = CODE[seq[i]];
    if (c === undefined) {
      valid = 0;
      continue;
    }
    fwd = ((fwd << 2) | c) & MASK;
    rev = (rev >>> 2) | ((3 - c) << shift);
    if (++valid >= K) fn(Math.min(fwd, rev), fwd, i - K + 1);
  }
}

/** Sorted, de-duplicated hashes of the kept k-mers. */
export function sketch(seq: string): Uint32Array {
  const kept = new Set<number>();
  forEachKmer(seq, (canon) => {
    const h = mix32(canon);
    if (h < THRESHOLD) kept.add(h);
  });
  return Uint32Array.from(kept).sort();
}

export const PROTEIN_K = 5;
const AMINO = "ACDEFGHIKLMNPQRSTVWY";

/**
 * Every amino-acid k-mer of the given proteins (no subsampling: proteins are
 * short), as sorted unique integers. Protein k-mers survive the silent
 * third-position changes that make distant relatives share few DNA k-mers.
 */
export function proteinKmers(proteins: string[]): Uint32Array {
  const kept = new Set<number>();
  for (const p of proteins) {
    let code = 0;
    let valid = 0;
    for (let i = 0; i < p.length; i++) {
      const c = AMINO.indexOf(p[i]);
      if (c < 0) {
        valid = 0;
        code = 0;
        continue;
      }
      code = (code * 20 + c) % 20 ** PROTEIN_K;
      if (++valid >= PROTEIN_K) kept.add(code);
    }
  }
  return Uint32Array.from(kept).sort();
}

/** Number of shared values between two sorted arrays. */
export function intersectSize(a: Uint32Array, b: Uint32Array): number {
  let i = 0;
  let j = 0;
  let n = 0;
  while (i < a.length && j < b.length) {
    if (a[i] === b[j]) {
      n++;
      i++;
      j++;
    } else if (a[i] < b[j]) i++;
    else j++;
  }
  return n;
}

/**
 * Overlap relative to the smaller sketch, so a COX1 reference is scored fully
 * against a whole mitogenome that contains it, and a short barcode fragment
 * against the full-length gene.
 */
export function containment(query: Uint32Array, ref: Uint32Array): number {
  const smaller = Math.min(query.length, ref.length);
  return smaller ? intersectSize(query, ref) / smaller : 0;
}
