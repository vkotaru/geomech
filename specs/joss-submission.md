# JOSS submission plan

Target: Journal of Open Source Software (paper.md in this repository), with an arXiv
preprint of the longer technical paper at submission. Requirements checked against the
JOSS submission and paper guides on 2026-11-14.

## Checklist

- [ ] 1. Merge the visualizer, report and skill PR into `main`
  - [ ] final-review fixes (browser tests skip cleanly, controller called once, skill gaps, polish)
  - [ ] verify, commit, PII scan, push `feature/viz-skill`, open PR, CI green, merge
- [ ] 2. Repository readiness
  - [ ] `CONTRIBUTING.md`, code of conduct, support statement
  - [ ] `CITATION.cff`, issue templates, license declared in `pyproject.toml`
  - [ ] installation checked in a clean environment
  - [ ] docs moved to Sphinx + Furo on Read the Docs, matching udaan: left sidebar grouped into
        Getting started / Guides / Theory / API reference, short page titles, MyST Markdown,
        autodoc API pages; `.readthedocs.yaml`; links updated (author: import the repo on readthedocs.org)
  - [ ] PyPI publishing like udaan: `publish.yml` on `v*` tags, `release` environment, trusted
        publishing, GitHub release with notes (author: add the pending publisher on pypi.org)
- [ ] 3. Verified systems library: one function per paper system, citation, stored reference output, `check_eom`
- [ ] 4. `paper/paper.md` + `paper/paper.bib` (750-1750 words; Summary, Statement of need,
      State of the field, Software design, Research impact statement, AI usage disclosure,
      Acknowledgements, References) + one figure; JOSS draft-PDF GitHub Action
- [ ] 5. arXiv preprint: rework the IEEE draft (verified systems only, equations generated
      from the systems library, S² formula fix, new URLs, check_eom evidence)
- [ ] 6. Release v0.2.0: tag, PyPI upload; Zenodo archive DOI at acceptance
- [ ] 7. Final adversarial review of the paper against the code

## Open questions for the author

- Evidence of research use (papers, theses, projects, other users).
- Development gap May-October 2026 versus JOSS's "active development across six months".
