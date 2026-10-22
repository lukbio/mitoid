import { describe, expect, it } from "vitest";
import { align, jukesCantor } from "./align";
import { identify, indexDatabase, orient, verdictFor, type Database, type Hit } from "./identify";
import { gcPercent, parseInput, reverseComplement, sixFrames, translate } from "./seq";
import { containment, intersectSize, proteinKmers, sketch } from "./sketch";
import { layout, midpointRoot, neighborJoining, type TreeNode } from "./tree";

/** Deterministic pseudo-random DNA. */
function randomDna(n: number, seed = 1): string {
  let x = seed;
  let s = "";
  for (let i = 0; i < n; i++) {
    x = (Math.imul(x, 1103515245) + 12345) >>> 0;
    s += "ACGT"[(x >>> 16) & 3];
  }
  return s;
}

/** Copy with every k-th base changed. */
function mutate(seq: string, every: number): string {
  const swap: Record<string, string> = { A: "G", G: "A", C: "T", T: "C" };
  return [...seq].map((c, i) => (i % every === every - 1 ? swap[c] : c)).join("");
}

const leafNames = (t: TreeNode): string[] => (t.children.length ? t.children.flatMap(leafNames) : [t.name!]);

describe("parseInput", () => {
  it("reads FASTA, RNA, ambiguity codes and junk", () => {
    const p = parseInput(">my sample\nacgu ny\n12acgt\n>second\nGGGG");
    expect(p.name).toBe("my sample");
    expect(p.seq).toBe("ACGTNNACGT");
    expect(p.ambiguous).toBe(2);
    expect(p.dropped).toBe(2);
  });
  it("accepts a bare sequence", () => {
    expect(parseInput("AC GT\n-AC").seq).toBe("ACGTAC");
  });
});

describe("seq helpers", () => {
  it("reverse-complements", () => expect(reverseComplement("AACGTN")).toBe("NACGTT"));
  it("computes GC ignoring N", () => expect(gcPercent("GGCCAATTNN")).toBe(50));
  it("translates with the vertebrate mitochondrial code", () => {
    // ATA = Met, TGA = Trp, AGA = stop in table 2 (Ile, stop, Arg in the standard code)
    expect(translate("ATGATATGAAGANNN")).toBe("MMW*X");
  });
  it("gives six frames, three per strand", () => {
    const f = sixFrames("ATGGCC");
    expect(f).toHaveLength(6);
    expect(f[0]).toBe("MA");
    expect(f[1]).toBe(translate(reverseComplement("ATGGCC")));
  });
});

describe("protein k-mers", () => {
  it("are unchanged by synonymous (third-position) substitutions", () => {
    // Leu codons CTT/CTC/CTA/CTG and Gly codons GGT/GGC/GGA/GGG: same protein.
    const a = "CTTGGTCTTGGTCTTGGTCTTGGT";
    const b = "CTCGGCCTAGGACTGGGGCTTGGC";
    expect(a).not.toBe(b);
    expect(Array.from(proteinKmers([translate(a)]))).toEqual(Array.from(proteinKmers([translate(b)])));
  });
});

describe("sketch", () => {
  const a = randomDna(1500, 7);
  it("is strand-independent", () => {
    expect(Array.from(sketch(a))).toEqual(Array.from(sketch(reverseComplement(a))));
  });
  it("keeps about 1/4 of k-mers", () => {
    const n = sketch(a).length;
    expect(n).toBeGreaterThan(250);
    expect(n).toBeLessThan(500);
  });
  it("scores a gene inside a longer genome as fully contained", () => {
    const genome = randomDna(6000, 3) + a + randomDna(6000, 4);
    expect(containment(sketch(a), sketch(genome))).toBe(1);
    expect(containment(sketch(a), sketch(randomDna(1500, 99)))).toBeLessThan(0.05);
  });
  it("counts shared values of sorted arrays", () => {
    expect(intersectSize(Uint32Array.from([1, 3, 5, 9]), Uint32Array.from([2, 3, 9, 10]))).toBe(2);
  });
  it("gives the same count when one array is much larger (binary-search path)", () => {
    const big = Uint32Array.from({ length: 5000 }, (_, i) => i * 3);
    const small = Uint32Array.from([0, 4, 9, 300, 14997, 20000]);
    expect(intersectSize(small, big)).toBe(4);
    expect(intersectSize(big, small)).toBe(4);
  });
});

describe("align", () => {
  it("finds 100% identity for identical sequences", () => {
    const s = randomDna(400, 2);
    const r = align(s, s);
    expect(r.identity).toBe(1);
    expect(r.length).toBe(400);
  });
  it("measures substitutions", () => {
    const s = randomDna(1000, 5);
    const r = align(mutate(s, 10), s);
    expect(r.mismatches).toBe(100);
    expect(r.identity).toBeCloseTo(0.9, 5);
  });
  it("doesn't penalise a fragment for the reference's overhangs", () => {
    const ref = randomDna(1500, 6);
    const frag = ref.slice(300, 950);
    const r = align(frag, ref);
    expect(r.identity).toBe(1);
    expect([r.bStart, r.bEnd]).toEqual([300, 950]);
  });
  it("handles an internal deletion as one gap run", () => {
    const s = randomDna(600, 8);
    const r = align(s.slice(0, 300) + s.slice(303), s);
    expect(r.gaps).toBe(3);
    expect(r.matches).toBe(597);
  });
  it("Jukes–Cantor grows faster than p and caps when saturated", () => {
    expect(jukesCantor(0)).toBe(0);
    expect(jukesCantor(0.1)).toBeGreaterThan(0.1);
    expect(jukesCantor(0.8)).toBe(3);
  });
});

describe("trees", () => {
  // Additive distances from ((A:1,B:2):1,(C:3,D:1)) — NJ recovers them exactly.
  const names = ["A", "B", "C", "D"];
  const d = [
    [0, 3, 5, 3],
    [3, 0, 6, 4],
    [5, 6, 0, 4],
    [3, 4, 4, 0],
  ];
  it("neighbor-joining groups the right pairs", () => {
    const t = midpointRoot(neighborJoining(names, d));
    const sides = t.children.map((c) => leafNames(c).sort().join(""));
    expect(new Set(leafNames(t))).toEqual(new Set(names));
    // Midpoint of the longest path (C to B, length 6) is on C's side here, so
    // check the pairings through the cherry structure instead of exact sides.
    const cherries = (n: TreeNode): string[] =>
      n.children.length === 2 && n.children.every((c) => !c.children.length)
        ? [leafNames(n).sort().join("")]
        : n.children.flatMap(cherries);
    expect(cherries(t).concat(sides)).toEqual(expect.arrayContaining(["AB"]));
  });
  it("lays out every leaf on its own row", () => {
    const { leaves, lines } = layout(midpointRoot(neighborJoining(names, d)));
    expect(leaves.map((l) => l.y).sort()).toEqual([0, 1, 2, 3]);
    expect(lines.length).toBeGreaterThan(4);
  });
  it("copes with two taxa", () => {
    const t = neighborJoining(["X", "Y"], [[0, 0.2], [0.2, 0]]);
    expect(leafNames(t).sort()).toEqual(["X", "Y"]);
  });
});

describe("identify", () => {
  // A toy database: one "true" gene, near relatives and unrelated decoys.
  const gene = randomDna(1500, 42);
  const db: Database = {
    gene: "COX1",
    built: "test",
    source: "test",
    count: 6,
    records: [
      { acc: "T1", species: "Bufo truus", class: "Amphibia", order: "Anura", family: "Bufonidae", seq: gene },
      { acc: "T2", species: "Bufo cousinus", class: "Amphibia", order: "Anura", family: "Bufonidae", seq: mutate(gene, 20) },
      { acc: "T3", species: "Rana farus", class: "Amphibia", order: "Anura", family: "Ranidae", seq: mutate(gene, 7) },
      { acc: "D1", species: "Decoy one", class: "Aves", order: "", family: "", seq: randomDna(1500, 100) },
      { acc: "D2", species: "Decoy two", class: "Aves", order: "", family: "", seq: randomDna(1500, 101) },
      { acc: "D3", species: "Decoy three", class: "Aves", order: "", family: "", seq: randomDna(1500, 102) },
    ],
  };
  const idx = indexDatabase(db);

  it("identifies an exact barcode fragment", () => {
    const r = identify(gene.slice(100, 760), idx);
    expect(r.hits[0].ref.species).toBe("Bufo truus");
    expect(r.hits[0].alignment.identity).toBe(1);
    expect(r.verdict.level).toBe("species");
    expect(r.reversed).toBe(false);
  });
  it("finds the gene on the reverse strand inside a whole genome", () => {
    const genome = randomDna(5000, 9) + mutate(gene, 50) + randomDna(8000, 10);
    const r = identify(reverseComplement(genome), idx);
    expect(r.reversed).toBe(true);
    expect(r.region).not.toBeNull();
    expect(r.region!.start).toBeLessThanOrEqual(5000);
    expect(r.region!.end).toBeGreaterThanOrEqual(6500);
    expect(r.hits[0].ref.species).toBe("Bufo truus");
    expect(r.hits[0].alignment.identity).toBeCloseTo(0.98, 2);
    expect(r.tree && leafNames(r.tree)).toContain("Your sequence");
  });
  it("ranks relatives by identity", () => {
    const r = identify(mutate(gene, 25), idx);
    expect(r.hits.slice(0, 3).map((h) => h.ref.acc)).toEqual(["T1", "T2", "T3"]);
  });
  it("rejects very short input", () => {
    expect(() => identify("ACGT", idx)).toThrow(/too short/);
  });
  it("finds the gene in a long sequence despite chance matches elsewhere", () => {
    // Sprinkle short pieces of the gene far away from the real copy.
    const decoys = [0, 1, 2, 3, 4, 5].map((k) => randomDna(9000, 200 + k) + gene.slice(k * 200, k * 200 + 40)).join("");
    const genome = decoys + randomDna(2000, 300) + gene + randomDna(20000, 301);
    const { region } = orient(genome, gene);
    const start = decoys.length + 2000;
    expect(region).not.toBeNull();
    expect(region!.start).toBeGreaterThan(start - 400);
    expect(region!.end).toBeLessThan(start + gene.length + 400);
  });
  it("orient leaves short queries unwindowed", () => {
    expect(orient(gene.slice(0, 500), gene).region).toBeNull();
  });
  it("verdict thresholds", () => {
    const hit = (id: number): Hit[] => [
      { ref: db.records[0], screen: 1, alignment: { ...align("A", "A"), identity: id } },
    ];
    expect(verdictFor(hit(0.99)).level).toBe("species");
    expect(verdictFor(hit(0.95)).label).toBe("Bufo sp.");
    expect(verdictFor(hit(0.88)).level).toBe("family");
    expect(verdictFor(hit(0.8)).level).toBe("order");
    expect(verdictFor(hit(0.6)).level).toBe("none");
    expect(verdictFor([]).level).toBe("none");
  });
});
