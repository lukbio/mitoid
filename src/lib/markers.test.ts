import { describe, expect, it } from "vitest";
import { align } from "./align";
import { MARKER_RULES, verdictFor, type Hit, type RefRecord } from "./identify";
import { translate } from "./seq";

describe("genetic codes", () => {
  it("table 11 (plastid) differs from table 2 (vertebrate mitochondrial)", () => {
    // ATA: Met in table 2, Ile in 11. TGA: Trp vs stop. AGA: stop vs Arg.
    expect(translate("ATATGAAGA", 2)).toBe("MW*");
    expect(translate("ATATGAAGA", 11)).toBe("I*R");
  });
  it("rejects codes it doesn't know", () => {
    expect(() => translate("ATG", 5)).toThrow(/genetic code 5/);
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
