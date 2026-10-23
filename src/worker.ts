/// <reference lib="webworker" />
/** Loads and indexes reference sets on first use, then answers queries off the main thread. */

import { identify, indexDatabase, mergeDatabases, type Database, type IndexedDatabase, type Result } from "./lib/identify";

export type Mode = "animals" | "plants" | "fungi";

/** Reference sets per mode, searched in this order; files in one group are the
 * same gene and are searched as one set. */
export const MODE_FILES: Record<Mode, string[][]> = {
  animals: [["cox1", "cox1_inv"]],
  plants: [["matk"], ["rbcl"]],
  fungi: [["its"]],
};

export interface SetInfo {
  gene: string;
  count: number;
  built: string;
  source: string;
}

export type WorkerRequest = { type: "load"; mode: Mode } | { type: "identify"; id: number; seq: string; mode: Mode };
export type WorkerResponse =
  | { type: "ready"; mode: Mode; sets: SetInfo[]; ms: number }
  | { type: "result"; id: number; results: Result[] }
  | { type: "error"; id?: number; mode?: Mode; message: string };

const ctx = self as unknown as DedicatedWorkerGlobalScope;
const loading = new Map<Mode, Promise<IndexedDatabase[]>>();

async function fetchSet(name: string): Promise<Database> {
  const res = await fetch(`${import.meta.env.BASE_URL}db/${name}.json.gz`);
  if (!res.ok) throw new Error(`could not load ${name} (${res.status})`);
  // Servers may or may not decompress for us; check the gzip magic bytes.
  const bytes = new Uint8Array(await res.arrayBuffer());
  const gz = bytes[0] === 0x1f && bytes[1] === 0x8b;
  const text = gz
    ? await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text()
    : new TextDecoder().decode(bytes);
  return JSON.parse(text) as Database;
}

async function fetchGroup(names: string[]): Promise<IndexedDatabase> {
  return indexDatabase(mergeDatabases(await Promise.all(names.map(fetchSet))));
}

function load(mode: Mode): Promise<IndexedDatabase[]> {
  let p = loading.get(mode);
  if (!p) {
    const t = performance.now();
    p = Promise.all(MODE_FILES[mode].map(fetchGroup));
    loading.set(mode, p);
    p.then(
      (dbs) =>
        ctx.postMessage({
          type: "ready",
          mode,
          sets: dbs.map((d) => ({ gene: d.gene, count: d.count, built: d.built, source: d.source })),
          ms: performance.now() - t,
        } satisfies WorkerResponse),
      (err) => {
        loading.delete(mode); // allow a retry
        ctx.postMessage({ type: "error", mode, message: (err as Error).message } satisfies WorkerResponse);
      },
    );
  }
  return p;
}

ctx.onmessage = async (e: MessageEvent<WorkerRequest>) => {
  const msg = e.data;
  if (msg.type === "load") {
    load(msg.mode).catch(() => {});
    return;
  }
  try {
    const dbs = await load(msg.mode);
    ctx.postMessage({ type: "result", id: msg.id, results: dbs.map((db) => identify(msg.seq, db)) } satisfies WorkerResponse);
  } catch (err) {
    ctx.postMessage({ type: "error", id: msg.id, message: (err as Error).message } satisfies WorkerResponse);
  }
};
