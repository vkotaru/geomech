"""The agent skill is usable from Claude Code, Codex and Gemini CLI.

``skills/geomech`` is the only copy of the skill. Agents working in this repository reach it
through ``AGENTS.md`` (read by Codex and, when no ``CLAUDE.md`` exists, by Claude Code),
``CLAUDE.md`` and ``GEMINI.md``, which only import ``AGENTS.md``. The front matter has to
satisfy the Agent Skills specification (https://agentskills.io/specification), which all three
tools follow, so that the folder can be copied into any of their skill directories.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "skills" / "geomech"
FILES = ["SKILL.md", "reference.md"]


def _front_matter(path: Path) -> dict[str, str]:
    m = re.match(r"^---\n(.*?)\n---\n", path.read_text(encoding="utf-8"), re.S)
    assert m, "SKILL.md must start with YAML front matter"
    out = {}
    for line in m.group(1).splitlines():
        key, _, value = line.partition(":")
        out[key.strip()] = value.strip()
    return out


def test_front_matter_follows_the_agent_skills_spec():
    fields = _front_matter(SKILL / "SKILL.md")
    name, desc = fields["name"], fields["description"]
    # 1-64 characters, lowercase letters, digits and single hyphens, equal to the directory
    assert re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", name) and len(name) <= 64
    assert name == SKILL.name
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


def test_agents_md_points_at_the_skill():
    text = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    for target in ("skills/geomech/SKILL.md", "skills/geomech/reference.md", "CONTRIBUTING.md"):
        assert target in text and (ROOT / target).is_file(), target


def test_claude_and_gemini_files_only_import_agents_md():
    """One source of instructions: the per-agent files import it and add nothing."""
    for name in ("CLAUDE.md", "GEMINI.md"):
        assert (ROOT / name).read_text(encoding="utf-8").strip() == "@AGENTS.md", name
