import { describe, expect, it } from "vitest";
import { align } from "./align";
import { geneLocation, MARKER_RULES, mergeDatabases, rankLabels, verdictFor, type Hit, type RefRecord, type Result } from "./identify";
import { translate } from "./seq";

describe("genetic codes", () => {
  it("table 11 (plastid) differs from table 2 (vertebrate mitochondrial)", () => {
    // ATA: Met in table 2, Ile in 11. TGA: Trp vs stop. AGA: stop vs Arg.
    expect(translate("ATATGAAGA", 2)).toBe("MW*");
    expect(translate("ATATGAAGA", 11)).toBe("I*R");
  });
  it("rejects codes it doesn't know", () => {
    expect(() => translate("ATG", 3)).toThrow(/genetic code 3/);
  });
});

describe("marker thresholds", () => {
  const ref: RefRecord = { acc: "X", species: "Quercus alba", class: "Magnoliopsida", order: "Fagales", family: "Fagaceae", seq: "" };
  const hit = (id: number): Hit[] => [{ ref, screen: 1, alignment: { ...align("A", "A"), identity: id } }];

  it("99% is species level for COX1 and matK but only genus level for rbcL", () => {
    expect(verdictFor(hit(0.99), "COX1").level).toBe("species");
    expect(verdictFor(hit(0.99), "matK").level).toBe("species");
    expect(verdictFor(hit(0.99), "rbcL").level).toBe("genus");
  });
  it("warns that rbcL can't always separate species", () => {
    expect(verdictFor(hit(1), "rbcL").text).toMatch(/matK/);
  });
  it("names the expected marker when nothing matches", () => {
    expect(verdictFor(hit(0.5), "rbcL").text).toMatch(/plant rbcL/);
  });
  it("falls back to COX1 rules for an unknown gene", () => {
    expect(verdictFor(hit(0.985), "ND5").level).toBe("species");
  });
  it("reports the genus when several species tie for the best match", () => {
    const tie = (species: string, id: number): Hit => ({
      ref: { ...ref, species, acc: species },
      screen: 1,
      alignment: { ...align("A", "A"), identity: id },
    });
    const hits = [tie("Quercus robur", 1), tie("Quercus aliena", 1), tie("Castanea sativa", 0.998)];
    const v = verdictFor(hits, "matK");
    expect(v.level).toBe("genus");
    expect(v.label).toBe("Quercus sp.");
    expect(v.text).toMatch(/2 Quercus species/);
    const mixed = verdictFor([tie("Quercus robur", 1), tie("Castanea sativa", 1)], "matK");
    expect(mixed.level).toBe("family");
    expect(verdictFor([tie("Quercus robur", 1), tie("Quercus aliena", 0.99)], "matK").level).toBe("species");
  });
  it("orders each marker's thresholds sensibly", () => {
    for (const r of Object.values(MARKER_RULES)) {
      expect(r.species).toBeGreaterThan(r.genus);
      expect(r.genus).toBeGreaterThan(r.family);
      expect(r.family).toBeGreaterThan(r.order);
    }
  });
});

describe("result helpers", () => {
  const ref: RefRecord = { acc: "X", species: "A b", class: "", order: "", family: "", seq: "" };
  const hitAt = (aStart: number, aEnd: number, identity = 1): Hit => ({
    ref,
    screen: 1,
    alignment: { ...align("A", "A"), aStart, aEnd, identity },
  });
  const result = (reversed: boolean, region: Result["region"], hits: Hit[]): Result => ({
    gene: "matK", queryLength: 10000, reversed, region, hits, tree: null,
    verdict: { level: "species", label: "", text: "" }, timings: {},
  });

  it("gene location undoes the window and the strand", () => {
    // Window starts at 1,000; the gene is bases 150-1,665 of the window.
    expect(geneLocation(result(false, { start: 1000, end: 3000 }, [hitAt(150, 1665)]))).toEqual({ start: 1151, end: 2665, reverse: false });
    // Same on the reverse complement of a 10 kb sequence: positions flip.
    expect(geneLocation(result(true, { start: 1000, end: 3000 }, [hitAt(150, 1665)]))).toEqual({ start: 7336, end: 8850, reverse: true });
    expect(geneLocation(result(false, null, [hitAt(0, 500)]))).toBeNull();
  });

  it("marks tied ranks with =", () => {
    const hits = [1, 1, 1, 0.998, 0.99, 0.99, 0.95].map((id) => hitAt(0, 1, id));
    expect(rankLabels(hits)).toEqual(["=1", "=1", "=1", "4", "=5", "=5", "7"]);
  });
});

describe("animal genetic codes", () => {
  it("translates the codons that differ between mitochondrial codes", () => {
    // AGA, ATA, TGA, AAA under vertebrate, invertebrate, echinoderm and cnidarian codes.
    expect(translate("AGAATATGAAAA", 2)).toBe("*MWK");
    expect(translate("AGAATATGAAAA", 5)).toBe("SMWK");
    expect(translate("AGAATATGAAAA", 9)).toBe("SIWN");
    expect(translate("AGAATATGAAAA", 4)).toBe("RIWK");
  });

  it("merged sets keep each record's code", () => {
    const rec = (acc: string, code?: number): RefRecord => ({ acc, species: acc, class: "", order: "", family: "", code, seq: "" });
    const base = { gene: "COX1", built: "2026-01-01", source: "x" };
    const merged = mergeDatabases([
      { ...base, code: 2, count: 1, records: [rec("vert")] },
      { ...base, code: 5, built: "2026-02-01", count: 2, records: [rec("fly", 5), rec("urchin", 9)] },
    ]);
    expect(merged.records.map((r) => r.code)).toEqual([2, 5, 9]);
    expect(merged.count).toBe(3);
    expect(merged.built).toBe("2026-02-01");
  });
});
