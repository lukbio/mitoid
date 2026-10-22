import { useCallback, useEffect, useRef, useState } from "react";
import type { Result } from "./lib/identify";
import type { Mode, SetInfo, WorkerRequest, WorkerResponse } from "./worker";

/** Owns the search worker: loads a mode's reference sets once, runs one query at a time. */
export function useIdentifier(mode: Mode) {
  const worker = useRef<Worker | null>(null);
  const nextId = useRef(0);
  const [sets, setSets] = useState<Partial<Record<Mode, SetInfo[]>>>({});
  const [loadError, setLoadError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [results, setResults] = useState<Result[] | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [elapsed, setElapsed] = useState<number | null>(null);
  const started = useRef(0);

  useEffect(() => {
    const w = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
    worker.current = w;
    w.onmessage = (e: MessageEvent<WorkerResponse>) => {
      const msg = e.data;
      if (msg.type === "ready") setSets((s) => ({ ...s, [msg.mode]: msg.sets }));
      else if (msg.type === "error" && msg.id === undefined) setLoadError(msg.message);
      else if (msg.id === nextId.current) {
        setBusy(false);
        setElapsed(performance.now() - started.current);
        if (msg.type === "result") setResults(msg.results);
        else if (msg.type === "error") setError(msg.message);
      }
    };
    return () => w.terminate();
  }, []);

  // Load the current mode's sets as soon as it's chosen, not on first search.
  useEffect(() => {
    setLoadError(null);
    worker.current?.postMessage({ type: "load", mode } satisfies WorkerRequest);
  }, [mode]);

  const identify = useCallback(
    (seq: string) => {
      if (!worker.current) return;
      nextId.current += 1;
      started.current = performance.now();
      setBusy(true);
      setError(null);
      setResults(null);
      worker.current.postMessage({ type: "identify", id: nextId.current, seq, mode } satisfies WorkerRequest);
    },
    [mode],
  );

  const reset = useCallback(() => {
    nextId.current += 1;
    setBusy(false);
    setResults(null);
    setError(null);
  }, []);

  return { sets: sets[mode] ?? null, loadError, busy, results, error, elapsed, identify, reset };
}
