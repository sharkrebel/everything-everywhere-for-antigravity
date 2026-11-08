// Launcher of the sidecars of the Android CLI agent plugin, run by the agent host with its bundled Node:
//
//   "command": "node", "args": ["../launch.mjs", "<sidecar name>"]
//
// The implementation of the sidecars is in the Android CLI (`android sidecar run <name>`, which downloads the "sidecars" pack), so the
// agent host could run it directly, but it does not look up commands outside of its PATH and shows nothing if the command cannot run.
// This launcher finds the Android CLI and runs the sidecar with it; if it is not installed, or too old to run the sidecar, it serves a
// page that explains how to install or update it, on the port of the sidecar, until the user asks to retry.
//
// It only runs the `android` launcher installed by the Android CLI installers (which downloads the CLI itself on first use, and runs
// it), never the CLI (`android-cli`) directly. The launcher is, in order:
// - `android` from the PATH,
// - `android` where the installers put it, for an agent host started before the installation (the installers add its directory to the
//   PATH of new processes only).
//
// Only uses the Node standard library: the plugin ships this file as-is.

import { spawn } from "node:child_process";
import { randomUUID } from "node:crypto";
import { existsSync } from "node:fs";
import { createServer } from "node:http";
import { homedir } from "node:os";
import { join } from "node:path";

const WEB_PORT_ENV_VAR = "ANTIGRAVITY_SIDECAR_WEB_PORT";
const IS_WINDOWS = process.platform === "win32";

/** Exit code of the Android CLI for an invalid command line, e.g. a CLI without the `sidecar` command. */
const USAGE_ERROR_EXIT_CODE = 2;

/** How much of the end of the output of a failed CLI run is shown on the page. */
const MAX_OUTPUT_CHARS = 4000;

/** The official page of the Android CLI, with its installation instructions. */
const INSTALL_URL = "https://goo.gle/android-cli";

const SIDECAR_NAMES = {
  "device-hub": "Android Device Hub",
  "compose-preview": "Android Compose Preview",
};

const name = process.argv[2];
if (!name) {
  console.error("Usage: node launch.mjs <sidecar name>");
  process.exit(USAGE_ERROR_EXIT_CODE);
}
const title = SIDECAR_NAMES[name] ?? name;

/** The process of the running Android CLI, killed when this process is terminated. */
let child = null;
/** The server of the page shown when the Android CLI cannot run the sidecar. */
let pageServer = null;
/** Whether the agent host asked this process to stop: exit as soon as the Android CLI has stopped the sidecar. */
let stopping = false;

// How the agent host stops the sidecar:
// - On macOS and Linux, it sends SIGINT to the process group of this process, which includes the Android CLI and the sidecar server:
//   forwarding the signal to the CLI is redundant there, but makes this process stoppable on its own too. The host kills this process
//   (only) 5s later if it is still running; the CLI then finishes stopping the sidecar on its own.
// - On Windows, it runs `taskkill /T /F` on this process, which kills the whole process tree (no signal handler runs), like it did when
//   the host ran the CLI directly.
for (const signal of ["SIGINT", "SIGTERM", "SIGHUP"]) {
  process.on(signal, () => {
    stopping = true;
    if (child) {
      // The CLI stops the sidecar server gracefully (e.g. the mirroring sessions of the Device Hub), then exits; so does this process.
      child.kill(signal);
    } else {
      process.exit(0);
    }
  });
}

/** The `android` launcher commands to try, in order (see the top of this file). */
function cliCandidates() {
  const candidates = [IS_WINDOWS ? "android.exe" : "android"];
  const installed = IS_WINDOWS
    ? join(homedir(), "AppData", "AndroidCLI", "android.exe") // install.cmd
    : join(homedir(), ".local", "bin", "android"); // install.sh
  if (existsSync(installed)) {
    candidates.push(installed);
  }
  return candidates;
}

/**
 * Runs `<command> sidecar run <name>`, forwarding its output to the log of the sidecar. Resolves with the outcome:
 * - { notFound: true } if the command does not exist,
 * - { exitCode, output } when it exits, with the end of its output.
 */
function runSidecar(command) {
  return new Promise((resolve) => {
    let output = "";
    const keep = (chunk) => {
      output = (output + chunk.toString()).slice(-MAX_OUTPUT_CHARS);
    };
    console.error(`Starting the ${title} with: ${command} sidecar run ${name}`);
    const proc = spawn(command, ["sidecar", "run", name], { stdio: ["ignore", "pipe", "pipe"], windowsHide: true });
    child = proc;
    proc.stdout.on("data", (chunk) => {
      process.stdout.write(chunk);
      keep(chunk);
    });
    proc.stderr.on("data", (chunk) => {
      process.stderr.write(chunk);
      keep(chunk);
    });
    proc.on("error", (error) => {
      child = null;
      resolve(error.code === "ENOENT" ? { notFound: true } : { exitCode: -1, output: `${output}\n${error.message}` });
    });
    proc.on("close", (code, signal) => {
      child = null;
      resolve({ exitCode: code ?? -1, signal, output });
    });
  });
}

/** Whether a CLI run failed because the CLI does not have the `sidecar` command (or its `run` subcommand). */
function isMissingSidecarCommand(result) {
  return result.exitCode === USAGE_ERROR_EXIT_CODE && /is not a valid command/.test(result.output);
}

function escapeHtml(text) {
  return text.replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
}

/** The page shown for [problem], one of "not-installed", "outdated", "failed". */
function pageHtml(problem, details) {
  const terminal = IS_WINDOWS ? "a Command Prompt" : "a terminal";
  let body;
  if (problem === "not-installed") {
    body =
      `<h1>The Android CLI is not installed</h1>` +
      `<p>The ${escapeHtml(title)} runs with the Android CLI (<code>android</code>), which was not found.</p>` +
      `<p>Install it by following the instructions at ` +
      `<a href="${INSTALL_URL}" target="_blank" rel="noopener noreferrer">${INSTALL_URL}</a>, then click Retry.</p>` +
      `<p>If the Android CLI is installed elsewhere, add it to the PATH and restart Antigravity.</p>`;
  } else if (problem === "outdated") {
    const update = details.command.includes("/") || details.command.includes("\\") ? `"${details.command}" update` : "android update";
    body =
      `<h1>The Android CLI needs to be updated</h1>` +
      `<p>This version of the Android CLI (<code>${escapeHtml(details.command)}</code>) cannot run the ${escapeHtml(title)}.</p>` +
      `<p>Update it by running this command in ${terminal}:</p><pre>${escapeHtml(update)}</pre>` +
      `<p>Then click Retry.</p>`;
  } else {
    body =
      `<h1>The ${escapeHtml(title)} stopped</h1>` +
      `<p><code>${escapeHtml(details.command)} sidecar run ${escapeHtml(name)}</code> exited with code ${details.exitCode}.</p>` +
      (details.output.trim() ? `<pre class="output">${escapeHtml(details.output.trim())}</pre>` : "") +
      `<p>Click Retry to start it again.</p>`;
  }
  // Identifies this page, see the polling after Retry below.
  const pageMarker = `<meta name="android-sidecar-launcher" content="${randomUUID()}">`;
  return `<!doctype html>
<html>
<head>
<meta charset="utf-8">
${pageMarker}
<title>${escapeHtml(title)}</title>
<style>
  :root { color-scheme: light dark; font-family: system-ui, sans-serif; font-size: 14px; }
  body { margin: 0; padding: 24px; max-width: 720px; line-height: 1.5; }
  h1 { font-size: 18px; margin: 0 0 12px; }
  pre { padding: 8px 12px; border-radius: 6px; background: rgba(127, 127, 127, 0.15); white-space: pre-wrap; word-break: break-all; }
  pre.output { max-height: 240px; overflow: auto; font-size: 12px; }
  button { font: inherit; padding: 6px 16px; border-radius: 6px; border: 1px solid rgba(127, 127, 127, 0.5); cursor: pointer; }
  button:disabled { cursor: default; opacity: 0.6; }
</style>
</head>
<body>
${body}
<button id="retry">Retry</button>
<script>
  const button = document.getElementById("retry");
  button.addEventListener("click", async () => {
    button.disabled = true;
    button.textContent = "Starting\\u2026";
    try { await fetch("/__launcher/retry", { method: "POST" }); } catch (e) {}
    // The launcher frees the port and starts the sidecar: reload once something else serves it (the setup page of the Android CLI,
    // the sidecar, or a new page of the launcher if the retry failed). Like the setup page of the Android CLI:
    // - Fetches location.href: the pane may be loaded through a proxy of the agent host (e.g. a plugin:// URL), which maps it to the port.
    // - Does not reload on errors: while the port is not served, the agent host answers with its own error page (e.g. "plugin
    //   unavailable"), which would replace this page for good.
    // - Does not reload while this same page is still served (identified by its page ID), before the launcher has freed the port.
    const poll = async () => {
      try {
        const response = await fetch(location.href, { cache: "no-store" });
        if (response.ok && !(await response.text()).includes(${JSON.stringify(pageMarker)})) {
          location.reload();
          return;
        }
      } catch (e) {}
      setTimeout(poll, 500);
    };
    setTimeout(poll, 500);
  });
</script>
</body>
</html>
`;
}

/** Serves the page for [problem] on the port of the sidecar until the user clicks Retry. */
function showPageUntilRetry(problem, details) {
  const port = Number(process.env[WEB_PORT_ENV_VAR] || 0);
  return new Promise((resolve) => {
    const html = pageHtml(problem, details);
    pageServer = createServer((req, res) => {
      if (req.method === "POST" && req.url === "/__launcher/retry") {
        res.writeHead(204).end();
        // Free the port for the sidecar before starting it.
        pageServer.close(() => resolve());
        pageServer.closeAllConnections?.();
        pageServer = null;
        return;
      }
      res.writeHead(200, { "Content-Type": "text/html; charset=utf-8", "Cache-Control": "no-store" }).end(html);
    });
    pageServer.on("error", (error) => {
      console.error(`Cannot show the ${problem} page on port ${port}: ${error.message}`);
      // Exit, so that the agent host restarts the launcher (with a backoff).
      process.exit(1);
    });
    pageServer.listen(port, "127.0.0.1", () => {
      console.error(`Showing the ${problem} page on http://localhost:${pageServer.address().port}`);
    });
  });
}

async function main() {
  for (;;) {
    let ran = null;
    for (const command of cliCandidates()) {
      const result = await runSidecar(command);
      if (result.notFound) {
        continue;
      }
      ran = { command, ...result };
      break;
    }
    if (ran === null) {
      await showPageUntilRetry("not-installed", {});
    } else if (stopping || ran.signal) {
      // Stopped by the agent host (the signal was forwarded to the CLI), or killed.
      process.exit(0);
    } else if (isMissingSidecarCommand(ran)) {
      await showPageUntilRetry("outdated", ran);
    } else if (ran.exitCode === 0) {
      process.exit(0);
    } else {
      await showPageUntilRetry("failed", ran);
    }
  }
}

main();
