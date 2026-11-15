"""The agent skill is usable from Claude Code, Codex and Gemini CLI.

``skills/geomech`` is the canonical copy. Claude Code reads ``.claude/skills`` (checked in
``test_viz_skill.py``); Codex and Gemini CLI read ``.agents/skills``, which this file keeps
identical. The front matter also has to satisfy the Agent Skills specification
(https://agentskills.io/specification), which all three tools follow.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "skills" / "geomech"
AGENTS_COPY = ROOT / ".agents" / "skills" / "geomech"
FILES = ["SKILL.md", "reference.md"]


def _front_matter(path: Path) -> dict[str, str]:
    m = re.match(r"^---\n(.*?)\n---\n", path.read_text(encoding="utf-8"), re.S)
    assert m, "SKILL.md must start with YAML front matter"
    out = {}
    for line in m.group(1).splitlines():
        key, _, value = line.partition(":")
        out[key.strip()] = value.strip()
    return out


def test_agents_copy_is_identical():
    for name in FILES:
        assert (AGENTS_COPY / name).read_bytes() == (SKILL / name).read_bytes(), (
            f".agents/skills/geomech/{name} differs from skills/geomech/{name}; copy it"
        )
    assert sorted(p.name for p in AGENTS_COPY.iterdir()) == sorted(p.name for p in SKILL.iterdir())


def test_front_matter_follows_the_agent_skills_spec():
    fields = _front_matter(AGENTS_COPY / "SKILL.md")
    name, desc = fields["name"], fields["description"]
    # 1-64 characters, lowercase letters, digits and single hyphens, equal to the directory
    assert re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", name) and len(name) <= 64
    assert name == AGENTS_COPY.name == SKILL.name
    assert 0 < len(desc) <= 1024 and "Use when" in desc


def test_skill_stays_short_enough_to_load_whole():
    """The specification recommends under 500 lines; the rest goes in reference.md."""
    assert len((SKILL / "SKILL.md").read_text(encoding="utf-8").splitlines()) < 500


def test_skill_is_agent_neutral():
    """Nothing in the skill names one agent, its tools or its directories."""
    for name in FILES:
        text = (SKILL / name).read_text(encoding="utf-8")
        for word in ("Claude", "Anthropic", "Codex", "Gemini", "Artifact", "SendUserFile"):
            assert word not in text, f"{name} mentions {word}"


def test_setup_section_tells_the_agent_how_to_install():
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    assert "## Setup" in text and "pip install geomech" in text
