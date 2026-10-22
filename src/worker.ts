/// <reference lib="webworker" />
/** Loads and indexes the reference database once, then answers queries off the main thread. */

import { identify, indexDatabase, type Database, type IndexedDatabase } from "./lib/identify";

export type WorkerRequest = { type: "identify"; id: number; seq: string };
export type WorkerResponse =
  | { type: "ready"; count: number; built: string; source: string; ms: number }
  | { type: "result"; id: number; result: ReturnType<typeof identify> }
  | { type: "error"; id?: number; message: string };

const ctx = self as unknown as DedicatedWorkerGlobalScope;
let db: IndexedDatabase | null = null;

async function load() {
  const t = performance.now();
  const res = await fetch(`${import.meta.env.BASE_URL}db/cox1.json.gz`);
  if (!res.ok || !res.body) throw new Error(`Could not load the reference database (${res.status}).`);
  // Servers may or may not decompress for us; check the gzip magic bytes.
  const bytes = new Uint8Array(await res.arrayBuffer());
  const gz = bytes[0] === 0x1f && bytes[1] === 0x8b;
  const text = gz
    ? await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text()
    : new TextDecoder().decode(bytes);
  db = indexDatabase(JSON.parse(text) as Database);
  ctx.postMessage({ type: "ready", count: db.count, built: db.built, source: db.source, ms: performance.now() - t } satisfies WorkerResponse);
}

ctx.onmessage = (e: MessageEvent<WorkerRequest>) => {
  const msg = e.data;
  if (msg.type !== "identify") return;
  if (!db) {
    ctx.postMessage({ type: "error", id: msg.id, message: "Still loading the reference database." } satisfies WorkerResponse);
    return;
  }
  try {
    ctx.postMessage({ type: "result", id: msg.id, result: identify(msg.seq, db) } satisfies WorkerResponse);
  } catch (err) {
    ctx.postMessage({ type: "error", id: msg.id, message: (err as Error).message } satisfies WorkerResponse);
  }
};

load().catch((err) => ctx.postMessage({ type: "error", message: (err as Error).message } satisfies WorkerResponse));
