/** The search: sketch screen, strand and region finding, alignment, tree. */

import { align, jukesCantor, type Alignment } from "./align";
import { reverseComplement, sixFrames, translate } from "./seq";
import { containment, forEachKmer, intersectSize, proteinKmers, sketch } from "./sketch";
import { midpointRoot, neighborJoining, type TreeNode } from "./tree";

export interface RefRecord {
  acc: string;
  species: string;
  class: string;
  order: string;
  family: string;
  seq: string;
}

export interface Database {
  gene: string;
  built: string;
  source: string;
  count: number;
  records: RefRecord[];
}

export interface IndexedDatabase extends Database {
  sketches: Uint32Array[];
  proteins: Uint32Array[];
}

export interface Hit {
  ref: RefRecord;
  /** Sketch overlap used for the first screen (0-1). */
  screen: number;
  alignment: Alignment;
}

export interface Result {
  queryLength: number;
  /** Whether the query had to be reverse-complemented to match COX1. */
  reversed: boolean;
  /** Where COX1 sits in the (oriented) query, when the query is longer than the gene. */
  region: { start: number; end: number } | null;
  hits: Hit[];
  tree: TreeNode | null;
  verdict: Verdict;
  timings: Record<string, number>;
}

export interface Verdict {
  level: "species" | "genus" | "family" | "order" | "none";
  label: string;
  text: string;
}

export const QUERY_NAME = "Your sequence";
const SCREEN_DNA = 25;
const SCREEN_PROTEIN = 50;
const SHOW_HITS = 10;
const TREE_TAXA = 10;
const MIN_QUERY = 100;

export function indexDatabase(db: Database): IndexedDatabase {
  return {
    ...db,
    sketches: db.records.map((r) => sketch(r.seq)),
    proteins: db.records.map((r) => proteinKmers([translate(r.seq)])),
  };
}

/** Rough identity thresholds used in vertebrate COX1 barcoding. */
export function verdictFor(hits: Hit[]): Verdict {
  const top = hits[0];
  if (!top) return { level: "none", label: "No match", text: "Nothing in the reference set looks like this sequence." };
  const id = top.alignment.identity * 100;
  if (top.alignment.coverage < 0.5) {
    return { level: "none", label: "No match", text: "Only a short stretch lines up with anything in the reference set. This may not be vertebrate COX1." };
  }
  const sp = top.ref.species;
  const genus = sp.split(" ")[0];
  if (id >= 98) {
    return { level: "species", label: sp, text: "That's in the range usually seen within one species." };
  }
  if (id >= 92) {
    return { level: "genus", label: `${genus} sp.`, text: `Probably the same genus (${genus}), possibly a species that isn't in the reference set.` };
  }
  if (id >= 85 && top.ref.family) {
    return { level: "family", label: top.ref.family, text: `Probably the same family (${top.ref.family}), but not a close match.` };
  }
  if (id >= 78 && top.ref.order) {
    return { level: "order", label: top.ref.order, text: `Only a distant match; maybe the same order (${top.ref.order}).` };
  }
  return { level: "none", label: "No close match", text: `The best hit is only ${id.toFixed(1)}% identical. This may not be vertebrate COX1.` };
}

/** Which strand of the query carries the gene, and where it lies. Counts exact
 * k-mer matches against the best reference on each strand. */
export function orient(query: string, ref: string): { seq: string; reversed: boolean; region: { start: number; end: number } | null } {
  const refFwd = new Set<number>();
  forEachKmer(ref, (_c, fwd) => refFwd.add(fwd));
  const hitsOn = (s: string) => {
    const pos: number[] = [];
    forEachKmer(s, (_c, fwd, p) => {
      if (refFwd.has(fwd)) pos.push(p);
    });
    return pos;
  };
  const fwdPos = hitsOn(query);
  const rc = reverseComplement(query);
  const revPos = hitsOn(rc);
  const reversed = revPos.length > fwdPos.length;
  const seq = reversed ? rc : query;
  const pos = reversed ? revPos : fwdPos;
  if (query.length <= ref.length + 300 || pos.length < 5) return { seq, reversed, region: null };
  pos.sort((a, b) => a - b);
  // Ignore stray matches outside the gene: trim 2% at each end of the hit cloud.
  const lo = pos[Math.floor(pos.length * 0.02)];
  const hi = pos[Math.ceil(pos.length * 0.98) - 1];
  const pad = 150;
  return { seq, reversed, region: { start: Math.max(0, lo - pad), end: Math.min(seq.length, hi + pad + 11) } };
}

/** Identity, discounted when the overlap covers under half the shorter sequence. */
export function matchScore(a: Alignment): number {
  return a.coverage >= 0.5 ? a.identity : a.identity * a.coverage;
}

function pDistance(a: Alignment): number {
  const sites = a.matches + a.mismatches;
  return sites ? a.mismatches / sites : 0.75;
}

export function identify(query: string, db: IndexedDatabase, now: () => number = () => performance.now()): Result {
  if (query.length < MIN_QUERY) throw new Error(`Sequence too short: need at least ${MIN_QUERY} bases, got ${query.length}.`);
  const timings: Record<string, number> = {};
  let t = now();

  // Two screens, then align the union: DNA k-mers rank near-identical
  // sequences well; protein k-mers (all six frames of the query) still find
  // relatives at 80-90% DNA identity, where shared DNA k-mers are mostly noise.
  const qs = sketch(query);
  const qp = proteinKmers(sixFrames(query));
  const byDna = db.sketches
    .map((s, i) => ({ i, score: containment(qs, s) }))
    .sort((a, b) => b.score - a.score)
    .slice(0, SCREEN_DNA);
  const byProtein = db.proteins
    .map((p, i) => ({ i, score: p.length ? intersectSize(qp, p) / p.length : 0 }))
    .sort((a, b) => b.score - a.score)
    .slice(0, SCREEN_PROTEIN);
  const screenScore = new Map<number, number>();
  for (const { i, score } of [...byDna, ...byProtein]) screenScore.set(i, Math.max(score, screenScore.get(i) ?? 0));
  const screened = [...screenScore].map(([i, score]) => ({ i, score })).sort((a, b) => b.score - a.score);
  timings.screen = now() - t;
  t = now();

  const best = db.records[byDna[0].score >= byProtein[0].score ? byDna[0].i : byProtein[0].i];
  const { seq, reversed, region } = orient(query, best.seq);
  const target = region ? seq.slice(region.start, region.end) : seq;

  const hits: Hit[] = screened
    .map(({ i, score }) => ({ ref: db.records[i], screen: score, alignment: align(target, db.records[i].seq) }))
    .sort((a, b) => matchScore(b.alignment) - matchScore(a.alignment) || b.screen - a.screen)
    .slice(0, SHOW_HITS);
  timings.align = now() - t;
  t = now();

  let tree: TreeNode | null = null;
  const treeHits = hits.slice(0, TREE_TAXA);
  if (treeHits.length >= 2) {
    const names = [QUERY_NAME, ...treeHits.map((h) => h.ref.species)];
    const n = names.length;
    const d = Array.from({ length: n }, () => new Array<number>(n).fill(0));
    treeHits.forEach((h, k) => {
      d[0][k + 1] = d[k + 1][0] = jukesCantor(pDistance(h.alignment));
    });
    for (let a = 0; a < treeHits.length; a++) {
      for (let b = a + 1; b < treeHits.length; b++) {
        const x = jukesCantor(pDistance(align(treeHits[a].ref.seq, treeHits[b].ref.seq)));
        d[a + 1][b + 1] = d[b + 1][a + 1] = x;
      }
    }
    tree = midpointRoot(neighborJoining(names, d));
  }
  timings.tree = now() - t;

  return { queryLength: query.length, reversed, region, hits, tree, verdict: verdictFor(hits), timings };
}
