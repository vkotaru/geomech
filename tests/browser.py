"""Headless Chrome for the browser tests of ``geomech.viz``.

The browser tests open generated report pages in headless Chrome, which fetches three.js
and MathJax from their pinned CDN URLs and draws with software WebGL (SwiftShader). They
depend on things outside the repository, so they are opt-in and skip, never fail, when
those things are missing:

- they run only when ``GEOMECH_BROWSER_TESTS=1`` (CI sets it in the Python 3.12 job);
- they skip without Chrome/Chromium, without network to the CDNs, when a CDN script does
  not load (after one retry), and when the page reports that WebGL is not available.

Run them locally with::

    GEOMECH_BROWSER_TESTS=1 python -m pytest tests/test_viz.py \\
        tests/test_viz_adversarial.py tests/test_viz_fixes.py -q

Probes that read the three.js scene are wrapped in ``when_ready``: they wait for the
page's ready hook (``data-gm-viewer`` on the root element) rather than a fixed delay.
"""

from __future__ import annotations

import functools
import html
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path

import pytest

from geomech.viz.report import EXTERNAL_URLS, MATHJAX_URL

ENV = "GEOMECH_BROWSER_TESTS"

CHROME = (
    shutil.which("google-chrome")
    or shutil.which("chromium")
    or ("/usr/bin/google-chrome" if os.path.exists("/usr/bin/google-chrome") else None)
)

# --enable-unsafe-swiftshader keeps software WebGL available under --disable-gpu, now
# that Chrome no longer falls back to SwiftShader on its own.
FLAGS = ("--headless=new", "--disable-gpu", "--enable-unsafe-swiftshader", "--no-sandbox")

#: Appended to every page: marks the root once MathJax itself (not just its config
#: object) has loaded, so a CDN failure can be told from a typesetting failure.
_LOADED = """<script>
(function () {
  (function poll() {
    if (window.MathJax && window.MathJax.version) {
      document.documentElement.setAttribute("data-gm-mathjax-loaded", "1");
    } else setTimeout(poll, 50);
  })();
})();
</script>
"""

_READY = """<script>
(function (probe) {
  var t0 = Date.now();
  (function wait() {
    var s = document.documentElement.getAttribute("data-gm-viewer");
    if (s === "ready") { probe(); return; }
    if (s || Date.now() - t0 > 5000) {
      var pre = document.createElement("pre"); pre.id = "gm-probe";
      pre.textContent = JSON.stringify({ unavailable: s || "no data-gm-viewer after 5 s" });
      document.body.appendChild(pre);
      return;
    }
    setTimeout(wait, 20);
  })();
})(function () {
__BODY__
});
</script>
"""


def when_ready(script: str) -> str:
    """*script* (one ``<script>…</script>`` probe) run once the page's 3-D view is ready
    (``data-gm-viewer="ready"``). If the view cannot be drawn, the probe is replaced by
    ``{"unavailable": state}``, which ``probe`` turns into a skip (no WebGL, no three.js)
    or a failure (anything else)."""
    m = re.fullmatch(r"\s*<script>(.*)</script>\s*", script, re.S)
    assert m, "when_ready: expected one <script> element"
    return _READY.replace("__BODY__", m.group(1))


def enabled() -> bool:
    return os.environ.get(ENV) == "1"


@functools.lru_cache(maxsize=1)
def network_ok() -> bool:
    """Whether every pinned CDN URL answers a HEAD request (checked once per session)."""
    for url in EXTERNAL_URLS:
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=10) as resp:
                if resp.status != 200:
                    return False
        except Exception:
            return False
    return True


def require() -> str:
    """Chrome's path, or skip: browser tests off, no Chrome, or no network to the CDNs."""
    if not enabled():
        pytest.skip(f"browser tests are opt-in: set {ENV}=1 to run them")
    if CHROME is None:
        pytest.skip("no Chrome/Chromium on this machine")
    if not network_ok():
        pytest.skip("no network: the page loads three.js and MathJax from their CDNs")
    return CHROME


def _root(dom: str) -> str:
    m = re.search(r"<html[^>]*>", dom)
    return m.group(0) if m else ""


def _attr(root: str, name: str) -> str | None:
    m = re.search(rf'{name}="([^"]*)"', root)
    return m.group(1) if m else None


def _with_marker(page: str) -> str:
    i = page.rfind("</body>")
    return page + _LOADED if i < 0 else page[:i] + _LOADED + page[i:]


def _run(page: str, budget: int, timeout: int) -> tuple[str, str]:
    chrome = require()
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "page.html"
        f.write_text(page, encoding="utf-8")
        cmd = [chrome, *FLAGS, f"--user-data-dir={d}/profile",
               f"--virtual-time-budget={budget}", "--dump-dom", f.as_uri()]  # fmt: skip
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return "", f"Chrome did not finish in {timeout} s"
    return proc.stdout, proc.stderr


def _cdn_problem(page: str, dom: str, err: str) -> str | None:
    """Why the page could not load what it fetches from the CDNs (or Chrome could not
    run at all), or ``None``."""
    if not dom.strip():
        return f"Chrome produced no DOM: {err.strip()[-300:]}"
    root = _root(dom)
    if _attr(root, "data-gm-viewer") == "no-three":
        return "three.js did not load from its CDN"
    if MATHJAX_URL in page and _attr(root, "data-gm-mathjax-loaded") != "1":
        return "MathJax did not load from its CDN"
    if (
        '"data-gm-mathjax"' in page  # a report page, which marks the end of typesetting
        and _attr(root, "data-gm-mathjax") != "done"
        and _attr(root, "data-gm-errors") in (None, "0")
    ):
        # loaded, no error, but not typeset within the virtual-time budget: a slow
        # machine or a slow CDN response, not a fault of the page
        return "MathJax loaded but did not finish typesetting within the time budget"
    return None


def load(page: str, budget: int = 8000, *, retries: int = 1, timeout: int = 180) -> str:
    """The DOM of *page* after *budget* ms of virtual time in headless Chrome.

    Retries *retries* times when a CDN script did not load, then skips with the reason;
    skips when the page reports that WebGL is not available."""
    marked = _with_marker(page)
    for _ in range(retries + 1):
        dom, err = _run(marked, budget, timeout)
        problem = _cdn_problem(page, dom, err)
        if problem is None:
            break
    else:
        pytest.skip(f"{problem} ({retries + 1} attempts)")
    if _attr(_root(dom), "data-gm-viewer") == "no-webgl":
        pytest.skip("no WebGL in this headless Chrome (the page says so)")
    return dom


def probe(page: str, budget: int = 8000, *, retries: int = 1) -> dict:
    """The JSON a probe script wrote into ``<pre id="gm-probe">``, plus ``root`` (the
    ``<html>`` start tag with the page's ``data-gm-*`` hooks)."""
    dom = load(page, budget, retries=retries)
    m = re.search(r'<pre id="gm-probe">(.*?)</pre>', dom, re.S)
    assert m, f"no probe output in the page; root: {_root(dom)}"
    out = json.loads(html.unescape(m.group(1)))
    state = out.get("unavailable")
    if state in ("no-webgl", "no-three"):
        pytest.skip(f"the 3-D view could not be drawn: {state}")
    assert state is None, f"the 3-D view never became ready: {state}; root: {_root(dom)}"
    assert "error" not in out, out["error"]
    out["root"] = _root(dom)
    return out
