"""``geomech.__version__`` is the installed package version, the one in ``pyproject.toml``,
and ``CITATION.cff`` names the same version."""

import pathlib
import re

import geomech

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _pyproject_version():
    # tomllib needs Python 3.11; CI also runs 3.10, so read the one line directly.
    text = (ROOT / "pyproject.toml").read_text()
    (match,) = re.findall(r'^version\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    return match


def test_version_matches_pyproject():
    assert geomech.__version__ == _pyproject_version()


def test_citation_version_matches_pyproject():
    lines = (ROOT / "CITATION.cff").read_text().splitlines()
    (line,) = [ln for ln in lines if ln.startswith("version:")]
    assert line.split(":", 1)[1].strip() == _pyproject_version()
