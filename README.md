# Suppression of electronic pair order for slow lattice phonons

Davis Rens, SuperMatics Inc. (ORCID 0000-0001-5839-5197).

This repository accompanies the paper [10.5281/zenodo.24366508](https://doi.org/10.5281/zenodo.24366508) and its replay supplement [10.5281/zenodo.24366511](https://doi.org/10.5281/zenodo.24366511). The paper studies the spinful Holstein and bond SSH models without direct electron interactions. The proven slow-phonon bound concerns the decay of fixed local electronic correlations under the stated conditions. It is not a general ceiling on superconductivity, a finite-mass helicity theorem, or a prediction of a nonzero transition.

## Reproduce the checks

With Python 3.12 or newer, run `python -B replay.py --all-mutants`. This standard-library replay checks the stated finite arithmetic and confirms rejection of 25 sign and coefficient mutants. The optional proxy modes of `replay/base.py` and `replay/linear.py` require an explicitly supplied external CSV; these digitized proxies are not needed for the exact replay and are not thermodynamic results.

For the finite Lean lemmas, enter `lean/`, install the pinned toolchain and run `lake exe cache get` followed by `lake build`. The pinned Mathlib revision and transitive dependencies are recorded in `lakefile.toml` and `lake-manifest.json`. The Lean lemmas and rational replay do not formalize the paper's analytic, infinite-volume or ensemble arguments. The expanded proof is in `RIGOR.md`.

To build the paper, run `tectonic main.tex` from `note/`, or run `pdflatex main.tex`, `bibtex main`, and `pdflatex main.tex` twice more. The source includes citations only, not third-party papers.

Code is licensed under Apache-2.0; the paper is available under CC BY 4.0. DOI links resolve once the records are published. The GitHub Zenodo integration is not enabled for this repository.
