"""Record a looping GIF of the 3-D view of a ``geomech.viz`` report.

The report is opened with ``?view`` in headless Chrome, driven over the DevTools protocol:
the playback is paused, the scrub bar is set to each GIF frame's time in turn, and the
viewer is captured. ffmpeg then builds the GIF with a palette made for it.

    python scripts/report_gif.py report.html out.gif --start 0 --end 5 --fps 15 --width 560

Needs Google Chrome (or Chromium), ffmpeg and ``websocket-client``
(``python -m pip install websocket-client``); none of them is a dependency of geomech.
The report loads three.js from its CDN, so the machine needs network access.
"""

from __future__ import annotations

import argparse
import base64
import json
import shutil
import subprocess
import tempfile
import time
import urllib.request
from pathlib import Path

CHROME = ["google-chrome", "chromium", "chromium-browser"]


class Page:
    """A minimal DevTools-protocol client over one page's websocket."""

    def __init__(self, ws_url: str):
        import websocket

        self.ws = websocket.create_connection(ws_url, timeout=60, suppress_origin=True)
        self.n = 0

    def call(self, method: str, **params):
        self.n += 1
        self.ws.send(json.dumps({"id": self.n, "method": method, "params": params}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == self.n:
                if "error" in msg:
                    raise RuntimeError(f"{method}: {msg['error']}")
                return msg.get("result", {})

    def js(self, expr: str, await_promise: bool = False):
        r = self.call(
            "Runtime.evaluate", expression=expr, returnByValue=True, awaitPromise=await_promise
        )
        if "exceptionDetails" in r:
            raise RuntimeError(f"JS error: {r['exceptionDetails']}")
        return r["result"].get("value")


def _chrome() -> str:
    for name in CHROME:
        path = shutil.which(name)
        if path:
            return path
    raise SystemExit("no Chrome or Chromium on PATH")


def _wait_json(port: int, path: str, timeout: float = 20.0):
    t0 = time.time()
    while True:
        try:
            with urllib.request.urlopen(f"http://127.0.0.1:{port}{path}", timeout=2) as r:
                return json.loads(r.read())
        except OSError:
            if time.time() - t0 > timeout:
                raise
            time.sleep(0.2)


# The camera and its controls are local to the report's script; a copy of the page made for
# the recording hands them out (the published report is not changed).
_HOOK = "var controls = new THREE.OrbitControls(camera, renderer.domElement);"
_EXPOSE = " window.__gmCam = { camera: camera, controls: controls, toThree: toThree };"

# Pause, put the scrub bar on the last stored state at or before time t, and place the
# camera: its default offset from the target scaled by `dist`, the target on the body named
# `follow` (default: the framed centre).
_SEEK = """
(function (t, follow, dist) {
  var S = JSON.parse(document.getElementById("gm-data").textContent).scene;
  var k = 0;
  while (k + 1 < S.t.length && S.t[k + 1] <= t + 1e-9) k++;
  var G = window.__gmCam;
  if (!G.off) {
    G.c0 = G.controls.target.clone();
    G.off = G.camera.position.clone().sub(G.c0);
  }
  var tgt = G.c0.clone();
  if (follow) {
    var b = S.bodies.filter(function (b) { return b.name === follow; })[0];
    if (!b) throw new Error("no body named " + follow);
    tgt = G.toThree([b.pos[3 * k], b.pos[3 * k + 1], b.pos[3 * k + 2]]);
  }
  G.controls.target.copy(tgt);
  G.camera.position.copy(tgt).add(G.off.clone().multiplyScalar(dist));
  G.controls.update();
  var scrub = document.getElementById("gm-scrub");
  scrub.value = String(k);
  scrub.dispatchEvent(new Event("input"));
  return new Promise(function (ok) {
    requestAnimationFrame(function () { requestAnimationFrame(function () { ok(k); }); });
  });
})(%s, %s, %s)
"""


def record(
    report: Path,
    out: Path,
    start: float = 0.0,
    end: float | None = None,
    fps: float = 15.0,
    width: int = 560,
    size: int = 720,
    theme: str = "light",
    colors: int = 128,
    port: int = 9333,
    hold: float = 0.6,
    distance: float = 1.0,
    follow: str | None = None,
    speed: float = 1.0,
) -> Path:
    """Write the GIF of *report* between simulated times *start* and *end* (s).

    *size* is the viewport width in CSS pixels used to render (the viewer is square), *width*
    the GIF width; *hold* seconds of the last frame end the loop. *distance* scales the
    report's default camera distance; *follow* names a body the camera tracks; *speed* is
    simulated seconds per second of GIF.
    """
    report, out = Path(report).resolve(), Path(out).resolve()
    tmp = Path(tempfile.mkdtemp(prefix="gm-gif-"))
    text = report.read_text(encoding="utf-8")
    if _HOOK not in text:
        raise SystemExit(f"{report}: not a geomech.viz report (camera setup not found)")
    copy = tmp / "report.html"
    copy.write_text(text.replace(_HOOK, _HOOK + _EXPOSE, 1), encoding="utf-8")
    proc = subprocess.Popen(
        [
            _chrome(),
            "--headless=new",
            f"--remote-debugging-port={port}",
            f"--user-data-dir={tmp / 'profile'}",
            "--no-first-run",
            "--hide-scrollbars",
            "--enable-unsafe-swiftshader",
            "--use-angle=swiftshader",
            f"--window-size={size},{size + 200}",
            "about:blank",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        targets = _wait_json(port, "/json/list")
        ws = next(t["webSocketDebuggerUrl"] for t in targets if t.get("type") == "page")
        page = Page(ws)
        page.call(
            "Emulation.setDeviceMetricsOverride",
            width=size,
            height=size + 200,
            deviceScaleFactor=1,
            mobile=False,
        )
        page.call("Page.enable")
        page.call("Page.navigate", url=copy.as_uri() + "?view")
        t0 = time.time()
        ready = (
            "document.documentElement && document.documentElement.getAttribute('data-gm-viewer')"
        )
        while page.js(ready) != "ready":
            if time.time() - t0 > 60:
                raise RuntimeError("the viewer did not get ready (no network for three.js?)")
            time.sleep(0.3)
        if theme in ("light", "dark"):
            page.js(f"document.documentElement.setAttribute('data-theme', {theme!r})")
        t_last = page.js(
            "(function(){var S=JSON.parse(document.getElementById('gm-data').textContent)"
            ".scene; return S.t[S.t.length-1];})()"
        )
        end = t_last if end is None else min(end, t_last)
        rect = page.js(
            "(function(){var r=document.getElementById('gm-viewer').getBoundingClientRect();"
            "return [r.left, r.top, r.width, r.height];})()"
        )
        clip = dict(x=rect[0], y=rect[1], width=rect[2], height=rect[3], scale=1)
        n = max(2, int(round((end - start) / speed * fps)) + 1)
        frames = tmp / "frames"
        frames.mkdir()
        for i in range(n):
            t = start + (end - start) * i / (n - 1)
            page.js(
                _SEEK % (json.dumps(t), json.dumps(follow), json.dumps(distance)),
                await_promise=True,
            )
            shot = page.call("Page.captureScreenshot", format="png", clip=clip)
            (frames / f"f{i:04d}.png").write_bytes(base64.b64decode(shot["data"]))
        # hold the last frame a little before the loop restarts
        last = frames / f"f{n - 1:04d}.png"
        for j in range(int(round(hold * fps))):
            shutil.copy(last, frames / f"f{n + j:04d}.png")
        vf = f"fps={fps},scale={width}:-1:flags=lanczos"
        subprocess.run(
            [
                "ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
                "-i", str(frames / "f%04d.png"),
                "-vf", f"{vf},split[a][b];[a]palettegen=max_colors={colors}:stats_mode=diff[p];"
                "[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle",
                "-loop", "0", str(out),
            ],
            check=True,
        )  # fmt: skip
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("report", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--start", type=float, default=0.0, help="first simulated time, s")
    ap.add_argument("--end", type=float, default=None, help="last simulated time, s")
    ap.add_argument("--fps", type=float, default=15.0)
    ap.add_argument("--width", type=int, default=560, help="GIF width, px")
    ap.add_argument("--size", type=int, default=720, help="render width, CSS px")
    ap.add_argument("--theme", default="light", choices=["light", "dark", "auto"])
    ap.add_argument("--colors", type=int, default=128)
    ap.add_argument("--distance", type=float, default=1.0, help="camera distance factor")
    ap.add_argument("--speed", type=float, default=1.0, help="simulated s per GIF s")
    ap.add_argument("--follow", default=None, help="name of a body the camera tracks")
    args = ap.parse_args(argv)
    path = record(
        args.report, args.out, args.start, args.end, args.fps, args.width, args.size,
        args.theme, args.colors, distance=args.distance, follow=args.follow, speed=args.speed,
    )  # fmt: skip
    print(f"wrote {path} ({path.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
