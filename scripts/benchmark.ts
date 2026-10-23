/**
 * Accuracy benchmark: how often MitoID names the right species, and what it
 * says when the species isn't in the reference set.
 *
 * Usage:  npx vite-node scripts/benchmark.ts [-- --n 200 --seed 1 --only ITS]
 *
 * --only reruns one gene and keeps the other genes' rows from the last run.
 *
 * For each reference set it samples species and makes a test query from each:
 * a barcode-length window (658 bp, the standard COX1 barcode; 800 bp for the
 * plant genes; all of ITS, ~600 bp) with 1% random substitutions, standing in for sequencing error
 * and variation within a species. Two experiments:
 *
 *  1. Species in the set: the query's own species stays in the reference set.
 *     Right = the verdict names it, or names its genus when several species
 *     tie (the gene genuinely can't separate them). Wrong = another species.
 *  2. Species missing: the species is removed first (only species whose genus
 *     has another species in the set are sampled). Right = the verdict stops
 *     at the correct genus or family. Sister = it names another species of
 *     the right genus (sister species can be 98-100% identical in a barcode
 *     gene, which no barcode-only method can fix). Wrong = wrong genus or
 *     family.
 *
 * Writes docs/benchmark.md and docs/benchmark.json.
 */

import { existsSync, readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { identify, indexDatabase, mergeDatabases, type Database, type IndexedDatabase, type Result } from "../src/lib/identify";

const args = process.argv.slice(2);
const opt = (name: string, def: number) => {
  const i = args.indexOf(`--${name}`);
  return i >= 0 ? Number(args[i + 1]) : def;
};
const N = opt("n", 200);
const SEED = opt("seed", 1);
const SUB_RATE = 0.01;
const onlyAt = args.indexOf("--only");
const ONLY = onlyAt >= 0 ? args[onlyAt + 1] : null;

/** Small seeded PRNG (mulberry32), so a run can be repeated exactly. */
function rng(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function loadSet(names: string[]): IndexedDatabase {
  const dbs = names.map((n) => JSON.parse(gunzipSync(readFileSync(`public/db/${n}.json.gz`)).toString("utf8")) as Database);
  return indexDatabase(mergeDatabases(dbs));
}

function makeQuery(seq: string, len: number, rand: () => number): string {
  const start = seq.length > len ? Math.floor(rand() * (seq.length - len)) : 0;
  const window = seq.slice(start, start + len);
  let out = "";
  for (const ch of window) {
    if (rand() < SUB_RATE) {
      const others = "ACGT".replace(ch, "");
      out += others[Math.floor(rand() * others.length)];
    } else out += ch;
  }
  return out;
}

/** The set without record i (arrays are shallow-copied, nothing re-indexed). */
function without(db: IndexedDatabase, i: number): IndexedDatabase {
  const drop = <T,>(xs: T[]) => [...xs.slice(0, i), ...xs.slice(i + 1)];
  return { ...db, records: drop(db.records), sketches: drop(db.sketches), proteins: drop(db.proteins), count: db.count - 1 };
}

const genusOf = (sp: string) => sp.split(" ")[0];

type Outcome = "right" | "vague" | "sister" | "wrong";
type Tally = Record<Outcome, number>;
const tally = (): Tally => ({ right: 0, vague: 0, sister: 0, wrong: 0 });

/** Experiment 1: species present. */
function judgePresent(r: Result, species: string): Outcome {
  const v = r.verdict;
  if (v.level === "species") return v.label === species ? "right" : "wrong";
  if (v.level === "genus" && v.label === `${genusOf(species)} sp.`) return "right"; // a tie inside the genus
  return "vague";
}

/** Experiment 2: species missing. */
function judgeMissing(r: Result, species: string, family: string): Outcome {
  const v = r.verdict;
  if (v.level === "species") return genusOf(v.label) === genusOf(species) ? "sister" : "wrong";
  if (v.level === "genus") return v.label === `${genusOf(species)} sp.` ? "right" : "wrong";
  if (v.level === "family") return family && v.label === family ? "right" : "wrong";
  return "vague";
}

function sample<T>(xs: T[], n: number, rand: () => number): T[] {
  const a = [...xs];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(rand() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a.slice(0, n);
}

interface Row {
  set: string;
  group: string;
  experiment: "present" | "missing";
  n: number;
  right: number;
  vague: number;
  sister: number;
  wrong: number;
  ms: number;
}

const SETS: { name: string; files: string[]; len: number; groups: (code: number | undefined) => string }[] = [
  { name: "COX1", files: ["cox1", "cox1_inv"], len: 658, groups: (c) => (c === 2 ? "vertebrates" : "invertebrates") },
  { name: "matK", files: ["matk"], len: 800, groups: () => "plants" },
  { name: "rbcL", files: ["rbcl"], len: 800, groups: () => "plants" },
  { name: "ITS", files: ["its"], len: 2000, groups: () => "fungi" },
];

const previous: Row[] = ONLY && existsSync("docs/benchmark.json") ? JSON.parse(readFileSync("docs/benchmark.json", "utf8")).rows : [];
const rows: Row[] = [];
for (const set of SETS) {
  if (ONLY && set.name !== ONLY) {
    rows.push(...previous.filter((r) => r.set === set.name));
    continue;
  }
  const t0 = performance.now();
  const db = loadSet(set.files);
  console.error(`${set.name}: ${db.count} species, indexed in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
  const byGroup = new Map<string, number[]>();
  db.records.forEach((r, i) => {
    const g = set.groups(r.code ?? db.code);
    byGroup.set(g, [...(byGroup.get(g) ?? []), i]);
  });
  const genusCount = new Map<string, number>();
  for (const r of db.records) genusCount.set(genusOf(r.species), (genusCount.get(genusOf(r.species)) ?? 0) + 1);

  for (const [group, idx] of byGroup) {
    for (const experiment of ["present", "missing"] as const) {
      const rand = rng(SEED * 1000 + rows.length);
      const pool = experiment === "present" ? idx : idx.filter((i) => (genusCount.get(genusOf(db.records[i].species)) ?? 0) > 1);
      const picks = sample(pool, N, rand);
      const t = tally();
      const start = performance.now();
      for (const i of picks) {
        const rec = db.records[i];
        const query = makeQuery(rec.seq, set.len, rand);
        const target = experiment === "present" ? db : without(db, i);
        const r = identify(query, target);
        t[experiment === "present" ? judgePresent(r, rec.species) : judgeMissing(r, rec.species, rec.family)]++;
      }
      const ms = (performance.now() - start) / picks.length;
      rows.push({ set: set.name, group, experiment, n: picks.length, ...t, ms });
      console.error(`  ${group} / ${experiment}: ${t.right}/${picks.length} right, ${t.sister} sister, ${t.wrong} wrong (${ms.toFixed(0)} ms per query)`);
    }
  }
}

const pct = (x: number, n: number) => (n ? `${((100 * x) / n).toFixed(1)}%` : "-");
const lines = [
  "# MitoID accuracy benchmark",
  "",
  `Generated by \`scripts/benchmark.ts\` (seed ${SEED}, up to ${N} species per row). Each query is a barcode-length`,
  `window of a reference gene (COX1 658 bp, matK/rbcL 800 bp, all of ITS) with ${SUB_RATE * 100}% random substitutions.`,
  "",
  "**Species in the set:** right = the right species, or its genus when several species tie.",
  "**Species missing:** the species is removed first; right = stops at the correct genus or family,",
  "sister = names another species of the right genus, wrong = wrong genus or family. Vague = a less specific answer",
  "than possible (order level, or no call).",
  "",
  "| Gene | Group | Experiment | n | Right | Vague | Sister species | Wrong | ms/query |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|",
  ...rows.map(
    (r) =>
      `| ${r.set} | ${r.group} | ${r.experiment === "present" ? "species in the set" : "species missing"} | ${r.n} | ${pct(r.right, r.n)} | ${pct(r.vague, r.n)} | ${r.experiment === "present" ? "-" : pct(r.sister, r.n)} | ${pct(r.wrong, r.n)} | ${r.ms.toFixed(0)} |`,
  ),
  "",
];
mkdirSync("docs", { recursive: true });
writeFileSync("docs/benchmark.md", lines.join("\n"));
writeFileSync("docs/benchmark.json", JSON.stringify({ seed: SEED, n: N, subRate: SUB_RATE, rows }, null, 2));
console.log(lines.join("\n"));
