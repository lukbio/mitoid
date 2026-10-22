import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// GitHub Pages serves the site from /mitoid/; locally it's at the root.
export default defineConfig(({ command }) => ({
  base: command === "build" ? "/mitoid/" : "/",
  plugins: [react()],
  worker: { format: "es" },
}));
