---
name: compose-preview
description: Instructions for the Android Compose Preview sidecar, which renders Jetpack Compose @Preview functions with `android preview render` and shows them in the Android Compose Preview pane. Use it when the user wants to see a Compose preview, or to review the result of UI changes — re-render after changing the code so that the pane shows the result.
---

# Android Compose Preview Skill

The **Android Compose Preview** sidecar (plugin `android-cli`) is a local web server run by `android sidecar run compose-preview`.
It renders Compose `@Preview` functions with the `android preview render` CLI, and shows the result in the **Android Compose Preview**
pane, one preview below the other. A function can produce several images (one per `@Preview` annotation, or per value of a
`PreviewParameterProvider`): the pane shows them side by side. Configuring a new session replaces the previews shown.

## Rendering previews

1. **Agree on what to render.** Find the candidate `@Preview` functions (e.g. `grep -rn "@Preview" --include=*.kt <module>`, or by
   running `<android_bin> preview list --reference <reference> [--all] [--format=text|json] [-- <build_args>...]` in `project_dir`
   using `android_bin` from `GET /api/status` to also discover custom MultiPreview functions) and **ask the user** which previews to
   render, and confirm the arguments of `android preview render` (see [Arguments](#arguments)). Do not guess silently when several
   modules or previews match.
2. **Configure and render** with `POST /api/session` (see [HTTP API](#http-api)), with `"wait": true` so you get the result. The
   sidecar runs `android preview render`; do not run it yourself.
3. **Check the result.** If the job `status` is `failed`, read its `output` (build errors, render errors, missing classes), fix
   the problem, and render again with `POST /api/render`. If it is `completed`, tell the user to open the
   **Android Compose Preview** pane. You can look at the rendered images yourself: their paths are in the job `images`, and per
   preview in `GET /api/state`.
4. **After changing the code** of a rendered preview, **re-render** with `POST /api/render` and `"wait": true`, **only after
   all your code edits are saved**: never start the render in the same batch of parallel tool calls as an edit, or it may build
   the old code. Fix any render failure. Then tell the user that the pane shows the result, summarizing what you changed.

## Arguments

**Do not run `android preview render` (or `android preview render --help`) yourself.** The sidecar runs it, and it may not even be available
in your shell. To check the arguments, configure the session: a wrong reference fails within seconds, and the error is in the job
`output` (e.g. `No supported build system found for path: ...`). Fix the arguments and call `POST /api/session` again.

The sidecar runs, in `project_dir`, a single command that builds the project once and renders all the previews of the session (or,
if the command line would be too long, several commands one after the other, each with some of the previews):

```bash
android preview render [--daemon] [--component-tree] <reference> --fqn <fqn1> --fqn <fqn2>... [-- <build_args>...]
```

`--daemon` and `--component-tree` are settings of the pane (`--component-tree` on by default, `--daemon` off by default, see `/api/settings`). With `--daemon`, the renderer
keeps running in the background until it is idle for 10 minutes, which makes the next renders faster. With `--component-tree`, it
writes next to each image a JSON file with its composables (see `bounds_url` below), which the pane uses to highlight the composable
under the mouse.

- `reference`: the previewable unit, `[<system>://]<path>?<target>`:
  - `<path>` is a directory relative to `project_dir` (`.` for `project_dir` itself) from which the build system is detected;
  - `<target>` identifies the module or target in that build system: for Gradle, the project path and the build variant,
    `:<project>[@<variant>]` (e.g. `:app@debug`, `:feature:home@demoDebug`). The variant defaults to `debug`: a module with product
    flavors has no `debug` variant, so give the full variant name, flavors first (e.g. `demoDebug`). For Bazel, a label.
  - `<system>` is the build system name. When it is omitted, the build system is detected from `<path>`; if that fails, or
    several build systems match, give it explicitly. Built-in systems: `gradle`, `bazel`, `soong`, `lightbuild`; the user may
    have configured custom ones (e.g. `blaze`).
  - Examples: `gradle://.?:app@debug`, `gradle://.?:feature:home@demoDebug`, `bazel://.?//app/src/main:ui`,
    `blaze://.?//java/com/example/ui:ui` (google3, with `project_dir` the `google3` directory).
  - Ask the user if you are not sure which build system or target to use.
- `fqn`: the fully qualified name of a `@Preview` function to render (`fqns` in the API): the package, the class containing the
  function (for a top-level function in `HomeScreen.kt`, the file class `HomeScreenKt`), and the function name, e.g.
  `com.example.ui.HomeScreenKt.HomeScreenPreview`. The target must contain the files that declare the functions. You can list the
  available `@Preview` (and MultiPreview) FQNs in a target without rendering by running
  `<android_bin> preview list --reference <reference> [--all] [--format=text|json] [-- <build_args>...]` in `project_dir` (`--all`
  includes local workspace module dependencies; `--format=json` returns `@Preview` annotation parameters and
  `PreviewParameterProvider` metadata).
- `build_args`: arguments passed verbatim to the build system (optional).

## Finding the server

The port is chosen by the sidecar host and **changes every time the sidecar restarts**. The server listens on the loopback interface
only. Find its URL, in this order:

1. **The server info file**, written when the server starts and deleted when it stops:

   ```bash
   cat ~/.android/cli/sidecars/compose-preview.json
   # {"name":"compose-preview","url":"http://127.0.0.1:<PORT>","port":<PORT>,"pid":...,"started_at":...}
   ```

   (`$ANDROID_USER_HOME/cli/sidecars/compose-preview.json` if `ANDROID_USER_HOME` is set.)
2. As a last resort, the last "listening" line of the sidecar log, with a shell command (the file is not readable with file-viewing
   tools, and the command may need to run outside the sandbox):
   `tail -n 5 ~/.gemini/jetski/sidecar_data/android-cli/compose-preview/logs/sidecar.log`.

Then check it responds: `curl -s <URL>/api/status`. If it does not, the URL is stale or the sidecar is not running: ask the user
to open the **Android Compose Preview** pane (from the plugins menu of the conversation), which starts the sidecar, then look
for its URL again.

Call the server directly on `127.0.0.1:<PORT>`. Do not use the pane's proxied URL, which requires the user's browser credentials.

## HTTP API

All responses are JSON. Errors are `{"error": "..."}` with a 4xx/5xx status. Requests with an `Origin` header from another site are
rejected with 403; plain `curl` requests are accepted.

**Always send the header `X-Sidecar-Conversation: <your conversation ID>` with `POST /api/session`, `POST /api/render` and `GET
/api/state`.** The sidecar keeps the session, the previews and the last render of each conversation apart, and the pane of a
conversation only shows the previews of that conversation: without the header, the user does not see your previews. Use the header,
not a query parameter, and keep the URLs without a query string.

- `POST /api/session` `{"project_dir": "<absolute path>", "reference": "...", "fqns": ["..."], "build_args": [],
  "output_dir": "<absolute path>", "render": true, "wait": true}`: sets what to render, and renders it unless `"render": false`.
  `build_args` and `output_dir` (default: `<project_dir>/build/compose-preview`) are optional. The pane stops showing the previous
  previews. Returns `{"session": {...}, "job": {...}}`.
- `POST /api/render` `{"wait": true}`: re-renders all the previews of the current session (e.g. after changing code). If a
  render is already running, it returns that render instead of starting another one. With `"wait": true`, the request returns
  when the render is done; otherwise poll `GET /api/jobs/<id>`. The render builds the project first, so it can take minutes.
- `GET /api/jobs/<id>`: a render job: `{"id", "commands": [[...], ...], "status": "running"|"completed"|"failed", "error",
  "output": [lines], "images": ["<path>", ...]}`, with the command that renders all the previews (or the commands run one after the
  other, each with some of the previews). The job fails if any
  preview fails: the `error` lists the failed previews with their render error (e.g. a missing class), and the other previews
  are still rendered. If the build fails, all the previews fail with `android preview failed with exit code <n>`: the reason is in
  the `output`.
- `GET /api/state`: everything the pane shows: `{"session", "render": <last job>, "previews": [{"fqn", "name",
  "status": "pending"|"rendering"|"completed"|"failed", "error", "images": [{"name", "path", "image_url", "bounds_url", "configuration"}]}]}`. A
  preview whose render failed keeps the images of its previous render. `bounds_url` (`GET /api/bounds/<preview>/<image>`, or `null`
  when the component tree setting is off) is the component tree of the image (also copied next to the image as
  `<image name>.json`): `{"previewId", "methodFQN", "width", "height", "components": [...]}`, a tree of `{"name", "bounds":
  {"left", "top", "right", "bottom"}, "sourceLocation": {"fileName", "lineNumber"}, "semantics", "parameters", "modifiers",
  "children"}` composables, with bounds in pixels of the image (see [Reading the component tree](#reading-the-component-tree)).
  `configuration` (or `null` without a component tree) is the configuration the image was rendered with, from its component tree:
  `{"display_name", "width", "height", "density_dpi", "font_scale", "ui_mode", "device_id"}` (each field only if known);
  `ui_mode & 0x30 == 0x20` is night mode, and `ui_mode & 0x0f` the type of device (1 normal, 4 TV, 6 watch...).
- `GET /api/status`: `{"status": "ok", "android_bin", "rendering": bool}`.
- `GET /api/settings`: `{"use_daemon": false, "use_daemon_default": false, "component_tree": true, "component_tree_default": true}`:
  whether `android preview render` runs with `--daemon` and `--component-tree`. The user changes them in the pane; `POST
  /api/settings` with the settings to change (e.g. `{"use_daemon": true}`) changes them too, from the next render on, and keeps the
  others. Do not change them unless the user asks.

Example:

```bash
curl -s -X POST http://127.0.0.1:$PORT/api/session -H "X-Sidecar-Conversation: $CONVERSATION_ID" \
  -H 'Content-Type: application/json' -d '{
  "project_dir": "/home/me/MyApp", "reference": "<reference>",
  "fqns": ["com.example.ui.HomeScreenKt.HomeScreenPreview"], "wait": true}'
```

## Files

The images of each preview are copied to `<output_dir>/<fqn>/`, and replaced by the next render of that preview. They are named
`<function>_<n>.png` (`n` from 0), or `<function>_<display name>.png` for the values of a `PreviewParameterProvider`. The
`@Preview(name = ...)` is not used, and the order of the numbered images does not tell their configuration: read it from the
`configuration` of each image in `GET /api/state` (e.g. night mode, font scale, device), or open the images, and do not assume an
order.

## Reading the component tree

The component tree describes what each image shows, without looking at its pixels. For each composable: `name`, `bounds` (in
pixels of the image), `sourceLocation` (`fileName` and `lineNumber` of the call), and, when the renderer could extract them:

- `semantics`: what it means to the user and to accessibility services, e.g. `text`, `contentDescription`, `role`, `testTag`,
  `stateDescription`, `toggleableState`, `selected`, `disabled`, `heading`, `error`, `progress`, `actions`.
- `modifiers`: e.g. `padding`, `size`, `width`, `height`, `fill`, `background` (color or brush, and shape), `shadow`, `alpha`,
  `clip`, `clickable`, `scroll`, `weight`, `align`.
- `parameters`: the arguments of the call that are not their default value, e.g. the `color` and `fontSize` of a `Text`.

Use it to find the code that draws a part of an image (`sourceLocation`), to check a change (e.g. the `padding` or `text` you
expect), and to review accessibility (clickable composables without `contentDescription` or text, missing `role`). To get only the
trees, without rendering images, run `android preview render --component-tree --no-image` yourself (the pane does not show trees
without images).

## Guidelines

1. Render only the previews the user asked for: rendering more previews makes each render longer.
2. After each render, report failures with the relevant `output` lines instead of silently retrying.
