"""Sphinx configuration for the geomech documentation.

Build locally with ``python -m sphinx -W --keep-going -b html docs _build/html``.
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

from __future__ import annotations

import geomech

# ── Project ───────────────────────────────────────────────────────────
project = "geomech"
author = "Prasanth Kotaru"
copyright = f"2020–2026, {author}"

release = geomech.__version__
version = ".".join(release.split(".")[:2])

# ── Extensions ────────────────────────────────────────────────────────
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # Google-style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",  # the guides are Markdown
    "sphinx_copybutton",
]

# MyST: $…$ / $$…$$ math, \begin{align} environments, ::: fences, definition lists.
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
    "deflist",
]
# Anchors for [text](page.md#heading) links, as GitHub generates them.
myst_heading_anchors = 3

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# ── autodoc ───────────────────────────────────────────────────────────
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autodoc_class_signature = "separated"
# Optional extras, imported only by some functions; never needed to read the docs.
autodoc_mock_imports = ["casadi", "scipy"]

napoleon_google_docstring = True
napoleon_numpy_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

# ── HTML (furo) ───────────────────────────────────────────────────────
html_theme = "furo"
html_title = "geomech"  # the version is shown at the bottom of the sidebar instead
templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["sidebar-version.css"]
html_sidebars = {
    "**": [
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/scroll-start.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
        "sidebar/version.html",
        "sidebar/variant-selector.html",
    ]
}
html_theme_options = {
    "navigation_with_keys": True,
    "source_repository": "https://github.com/vkotaru/geomech/",
    "source_branch": "main",
    "source_directory": "docs/",
}
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True


# ── Docstrings that are not valid reStructuredText ────────────────────
# Many module docstrings are plain-text tables (rule names, dispatch grids) laid out with
# indentation. Rather than let docutils reflow them into broken block quotes, a docstring
# that does not parse cleanly is shown as preformatted text, after its first paragraph
# when that paragraph parses on its own. A dataclass's generated ``Name(field=…)``
# docstring repeats the signature and is dropped.
_IGNORED = (
    "Unknown interpreted text role",  # Sphinx roles (:class:, :func:, …)
    "Unknown directive type",  # Sphinx directives (seealso, versionadded, …)
)


def _rst_is_clean(text: str, env) -> bool:
    import io

    from docutils.core import publish_doctree

    stream = io.StringIO()
    publish_doctree(
        text,
        settings_overrides={
            "report_level": 2,
            "halt_level": 5,
            "warning_stream": stream,
            "file_insertion_enabled": False,
            "raw_enabled": False,
            "env": env,  # Sphinx roles look it up
        },
    )
    report = stream.getvalue()
    problems = [
        line
        for line in report.splitlines()
        if any(level in line for level in ("(WARNING/2)", "(ERROR/3)", "(SEVERE/4)"))
        and not any(ignored in line for ignored in _IGNORED)
    ]
    return not problems


def _preformat_unparsable(app, what, name, obj, options, lines):
    import dataclasses

    if (
        what == "class"
        and dataclasses.is_dataclass(obj)
        and lines
        and lines[0].startswith(f"{obj.__name__}(")
    ):
        lines[:] = []
        return
    if not lines or _rst_is_clean("\n".join(lines), app.env):
        return
    head: list[str] = []
    if "" in lines:
        cut = lines.index("")
        if _rst_is_clean("\n".join(lines[:cut]), app.env):
            head, lines_rest = lines[:cut] + [""], lines[cut + 1 :]
        else:
            lines_rest = list(lines)
    else:
        lines_rest = list(lines)
    while lines_rest and not lines_rest[-1].strip():
        lines_rest.pop()
    lines[:] = head + ["::", ""] + ["   " + line for line in lines_rest] + [""]


def setup(app):
    # before napoleon (priority 500): check the docstring as written, not its rewrite
    app.connect("autodoc-process-docstring", _preformat_unparsable, priority=400)
