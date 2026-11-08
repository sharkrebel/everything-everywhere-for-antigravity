---
name: android-device-hub
description: Instructions for working with the Android Device Hub sidecar, which lists, starts, stops, creates and removes Android Virtual Devices (AVDs), lists physical devices connected to adb, and lets the user mirror and control any of them (video, touch, keyboard, mouse wheel) in the Device Hub pane. Also explains how to inspect and interact with the screen the user is mirroring, with `android layout`, `android screen` and `adb`.
---

# Android Device Hub Skill

The **Android Device Hub** sidecar (plugin `android-cli`) is a local web server run by `android sidecar run device-hub`. Its pane
shows the user's devices and lets them:

- start, stop, create and remove **emulators** (AVDs), through the `android emulator` CLI;
- **mirror** any running emulator or **physical device** connected to adb: live video of the display, with touch, keyboard and
  mouse-wheel input forwarded to the device.

Physical devices include USB devices and devices attached with `adb connect` (e.g. a phone tunneled from another machine). They can
only be mirrored: start, stop, create and remove only apply to emulators.

## When to use it

- The user asks to **see or interact with a device screen**: tell them to open the **Android Device Hub** pane and click **Mirror**
  on the device (it must be online). You cannot see the mirrored video yourself; inspect the screen with the commands in
  [Inspecting the mirrored screen](#inspecting-the-mirrored-screen).
- You need to know **which devices the user sees** in the pane, or to start/stop an emulator so that the pane reflects it.

For scripted work, the `android emulator` CLI (`android emulator list|start|stop|create|remove`) and `adb` do the same operations
directly, and are usually simpler than the HTTP API. Operations started from the CLI also show up in the pane, since the pane reads
the same AVD and emulator files.

## Inspecting the mirrored screen

The mirror shows the user the live screen, but you can only inspect it with the `android` CLI and `adb`. They run on the same adb
server as the sidecar, so they reach every device listed by the pane, including physical devices attached with `adb connect` (e.g.
`127.0.0.1:6555`). Any input you send shows up live in the user's mirror view. Taking a screenshot or a layout does not disturb an
open mirror.

1. **Pick the device** the user is looking at: its `serial` in `GET /api/devices` (only online devices can be inspected), or
   `adb devices`.
2. **Read the UI structure** (start here; it is fast and text-based):
   ```bash
   android layout --device <serial> --pretty        # add --full for hidden and non-interactive elements, --flat for a list
   ```
   It returns the visible text, resource IDs, content descriptions, states and bounds of the UI elements.
3. **Look at the pixels** when the layout is not enough (custom-drawn UI, games, WebViews, images, or to check how the screen
   actually looks):
   ```bash
   android screen capture --device <serial> --annotate --output /tmp/screen.png
   ```
   `--annotate` draws numbered boxes around the UI elements detected in the screenshot. Open the PNG to view it. Omit
   `--annotate` for a plain screenshot.
4. **Act on an element**: tap the center of the element bounds reported by `android layout`, or substitute the center of an
   annotated box `#N` with `android screen resolve`:
   ```bash
   adb -s <serial> shell $(android screen resolve --screenshot /tmp/screen.png --string "input tap #7")
   ```
   Other input: `adb -s <serial> shell input text '<text>'`, `input keyevent KEYCODE_BACK`, `input swipe x1 y1 x2 y2`.

A screen capture takes a few seconds (more for a device attached over a tunnel), so re-capture only after actions that change the
screen. Secure windows (e.g. banking apps, DRM video) capture as black, as they appear in the mirror.

## Finding the server

The port is chosen by the sidecar host and **changes every time the sidecar restarts**. The server listens on the loopback interface
only. Find the current port from the last "listening" line of the sidecar log:

```bash
grep -o 'Device Hub listening on http://localhost:[0-9]*' \
  ~/.gemini/jetski/sidecar_data/android-cli/device-hub/logs/sidecar.log | tail -1
```

Then check it responds: `curl -s http://127.0.0.1:$PORT/api/status`.

Call the server directly on `127.0.0.1:$PORT` from the machine it runs on. Do not use the pane's proxied URL, which requires the
user's browser credentials.

## HTTP API

All responses are JSON. Errors are `{"error": "..."}` with a 4xx status.

### Read

- `GET /api/status`: `{"status": "ok", "port": <loopback port>, "android_bin": ..., "active_jobs": N, "server_time": ms}`
- `GET /api/devices`: `{"devices": [...]}`, emulators first, then physical devices:
  - Emulator: `{"kind": "emulator", "id": "<avd id>", "name", "target": "<API level>", "path", "status": "Online"|"Offline",
    "serial": "emulator-5554"|null, "active_job": {"id", "action", "status"}|null, "error": "<error of the last job on the
    device, if it failed>"|null}`
  - Physical device: `{"kind": "physical", "id": "<serial>", "name": "<model>", "target": "Android 15 (API 35)", "status":
    "Online"|"Unauthorized"|"Offline", "state": "<adb state>", "serial": "<serial>", "active_job": null}`
  - Both have `"api_level": "36"|null`, `"abi": "x86_64"|null`, `"android_version": "Android 16.0 (\"Baklava\")"|null` and
    `"form_factor": "phone"|"wear"|"tv"|"automotive"|"desktop"|"xr"`.
  - "Online" for an emulator means its process is running; it may still be booting. Wait for boot with
    `adb -s <serial> wait-for-device shell 'while [ -z "$(getprop sys.boot_completed)" ]; do sleep 1; done'`.
  - "Unauthorized" means the user must accept the "Allow USB debugging?" prompt on the device.
- `GET /api/profiles`: `{"profiles": [...]}`, the device profiles accepted by `/api/create`.
- `GET /api/jobs`: `{"jobs": [...]}`, the 50 most recent jobs, newest first: `{"id", "cmd", "action", "target", "status":
  "running"|"completed"|"failed", "output": [lines], "error", "created_at", "completed_at"}`.
- `GET /api/logs`: `{"logs": [...]}`, the server log shown in the pane.

### Actions

Actions run as background jobs. They return immediately with `{"job_id", "status", ...}`; poll `GET /api/jobs` for the result. The
body is JSON:

- `POST /api/start` `{"device": "<avd id>", "cold": false}`: start an emulator. The job completes once the device has booted.
- `POST /api/stop` `{"device": "<avd id>"}`: stop an emulator.
- `POST /api/create` `{"profile": "<profile>"}`: create an AVD from a profile listed by `/api/profiles`.
- `POST /api/remove` `{"device": "<avd id>", "force": false}`: delete an AVD. Stop it first; this is not checked. `force` deletes
  the AVD's `.ini` file even if its `.avd` folder is missing or invalid.

Example:

```bash
curl -s -X POST http://127.0.0.1:$PORT/api/start -H 'Content-Type: application/json' -d '{"device": "medium_phone"}'
curl -s http://127.0.0.1:$PORT/api/jobs
```

Requests with an `Origin` header from another site are rejected with 403; plain `curl` requests are accepted.

### Mirroring

The pane mirrors a device through the WebSocket `/mirror/ws`, which is for the web UI only. To send input to a device yourself, use
`adb shell input` instead. Opening a mirror of a device that is already mirrored in another view takes over the session: the other
view is disconnected.

## Guidelines

1. Check `GET /api/devices` (or `android emulator list` and `adb devices`) before starting a device, to avoid starting one that is
   already running.
2. Do not start several jobs on the same device at once (e.g. start and remove): the server does not serialize them.
3. After starting an emulator for the user to look at, tell them they can mirror it from the Device Hub pane.
