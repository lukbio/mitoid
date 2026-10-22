/** Parsing and basic DNA helpers. */

export interface ParsedInput {
  name: string;
  seq: string;
  /** Characters dropped because they aren't nucleotides (digits, spaces, etc.). */
  dropped: number;
  ambiguous: number;
}

const COMPLEMENT: Record<string, string> = { A: "T", C: "G", G: "C", T: "A", N: "N" };

/** Accepts raw sequence or FASTA (first record only). RNA (U) is read as DNA,
 * IUPAC ambiguity codes become N, everything else is ignored. */
export function parseInput(text: string): ParsedInput {
  const lines = text.split(/\r?\n/);
  let name = "Your sequence";
  const body: string[] = [];
  let seenHeader = false;
  for (const line of lines) {
    if (line.startsWith(">")) {
      if (seenHeader) break;
      seenHeader = true;
      name = line.slice(1).trim() || name;
      continue;
    }
    body.push(line);
  }
  let seq = "";
  let dropped = 0;
  let ambiguous = 0;
  for (const ch of body.join("").toUpperCase()) {
    if ("ACGT".includes(ch)) seq += ch;
    else if (ch === "U") seq += "T";
    else if ("RYSWKMBDHVN".includes(ch)) {
      seq += "N";
      ambiguous++;
    } else if (!/\s|-|\*|\./.test(ch)) dropped++;
  }
  return { name, seq, dropped, ambiguous };
}

export function reverseComplement(seq: string): string {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]] ?? "N";
  return out;
}

// Vertebrate mitochondrial code (NCBI table 2), codons in TCAG order:
// AGA/AGG are stops, ATA is Met, TGA is Trp.
const BASES = "TCAG";
const TABLE2 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";

/** Translate in frame 0; codons with N become X. */
export function translate(seq: string): string {
  let out = "";
  for (let i = 0; i + 3 <= seq.length; i += 3) {
    const a = BASES.indexOf(seq[i]);
    const b = BASES.indexOf(seq[i + 1]);
    const c = BASES.indexOf(seq[i + 2]);
    out += a < 0 || b < 0 || c < 0 ? "X" : TABLE2[a * 16 + b * 4 + c];
  }
  return out;
}

/** All six reading frames. */
export function sixFrames(seq: string): string[] {
  const rc = reverseComplement(seq);
  return [0, 1, 2].flatMap((f) => [translate(seq.slice(f)), translate(rc.slice(f))]);
}

export function gcPercent(seq: string): number {
  let gc = 0;
  let acgt = 0;
  for (const ch of seq) {
    if (ch === "G" || ch === "C") gc++;
    if (ch !== "N") acgt++;
  }
  return acgt ? (100 * gc) / acgt : 0;
}
