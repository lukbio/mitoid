import { useCallback, useEffect, useRef, useState } from "react";
import type { Result } from "./lib/identify";
import type { WorkerRequest, WorkerResponse } from "./worker";

export interface DbInfo {
  count: number;
  built: string;
  source: string;
  ms: number;
}

/** Owns the search worker: loads the database once, runs one query at a time. */
export function useIdentifier() {
  const worker = useRef<Worker | null>(null);
  const nextId = useRef(0);
  const [db, setDb] = useState<DbInfo | null>(null);
  const [loadError, setLoadError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [result, setResult] = useState<Result | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [elapsed, setElapsed] = useState<number | null>(null);
  const started = useRef(0);

  useEffect(() => {
    const w = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
    worker.current = w;
    w.onmessage = (e: MessageEvent<WorkerResponse>) => {
      const msg = e.data;
      if (msg.type === "ready") setDb({ count: msg.count, built: msg.built, source: msg.source, ms: msg.ms });
      else if (msg.type === "error" && msg.id === undefined) setLoadError(msg.message);
      else if (msg.id === nextId.current) {
        setBusy(false);
        setElapsed(performance.now() - started.current);
        if (msg.type === "result") setResult(msg.result);
        else if (msg.type === "error") setError(msg.message);
      }
    };
    return () => w.terminate();
  }, []);

  const identify = useCallback((seq: string) => {
    if (!worker.current) return;
    nextId.current += 1;
    started.current = performance.now();
    setBusy(true);
    setError(null);
    setResult(null);
    worker.current.postMessage({ type: "identify", id: nextId.current, seq } satisfies WorkerRequest);
  }, []);

  const reset = useCallback(() => {
    nextId.current += 1;
    setBusy(false);
    setResult(null);
    setError(null);
  }, []);

  return { db, loadError, busy, result, error, elapsed, identify, reset };
}
