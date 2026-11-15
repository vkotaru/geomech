"""The agent skill (``skills/geomech``): front matter and every fenced Python block run in a
subprocess."""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "skills" / "geomech"
FILES = ["SKILL.md", "reference.md"]

_BLOCK = re.compile(r"^```python\n(.*?)^```", re.S | re.M)


def _blocks():
    out = []
    for name in FILES:
        text = (SKILL / name).read_text(encoding="utf-8")
        for i, m in enumerate(_BLOCK.finditer(text)):
            out.append(pytest.param(m.group(1), id=f"{name}-{i}"))
    return out


def test_front_matter():
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    assert m, "SKILL.md must start with YAML front matter"
    fields = {}
    for line in m.group(1).splitlines():
        key, sep, value = line.partition(":")
        assert sep and re.fullmatch(r"[a-z_-]+", key), f"not a 'key: value' line: {line!r}"
        fields[key] = value.strip()
    assert set(fields) == {"name", "description"}
    assert fields["name"] == "geomech"
    assert re.fullmatch(r"[a-z0-9-]{1,64}", fields["name"])
    desc = fields["description"]
    assert 50 < len(desc) <= 1024 and "Use when" in desc
    assert "<" not in desc and ">" not in desc


def test_no_scala_reference_and_book_name():
    for name in FILES:
        text = (SKILL / name).read_text(encoding="utf-8")
        assert not re.search(r"\bscala\b", text, re.I)
        assert "LLM" not in text


def test_every_worked_example_is_present():
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    for title in ("spherical pendulum", "double rigid pendulum", "quadrotor with a cable"):
        assert f"Worked example: {title}" in text


@pytest.mark.parametrize("code", _blocks())
def test_snippet_runs(code, tmp_path):
    env = dict(os.environ, PYTHONPATH=str(ROOT), MPLBACKEND="Agg")
    proc = subprocess.run(
        [sys.executable, "-W", "error", "-c", code],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr[-4000:]
