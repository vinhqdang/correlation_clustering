# Submission package: Algorithmica, Springer Nature template (sn-jnl)

## Build
    pdflatex main && bibtex main && pdflatex main && pdflatex main
Class: `sn-jnl.cls` (Springer Nature template v3.1), options `pdflatex,sn-mathphys-num` (numbered references, `sn-mathphys-num.bst`).
Current state: 113 pages, no undefined references, no overfull boxes, no hyperref warnings.
Underfull-vbox reports are suppressed with `\vbadness=10000` (they come from the class output routine and also occur on the unmodified template).

## Files to upload
- `main.tex`
- `sections/*.tex` (all 28 files, incl. appendix files)
- `refs.bib` and the generated `main.bbl` (portals usually need the .bbl)
- `sn-jnl.cls`, `sn-mathphys-num.bst`
- `main.pdf`
- Not needed: `main.aux/.log/.out/.blg`, `README_submission.md`.
(If the portal wants one flat file, inline the `\input` files; no figures are used.)

## Structure notes
- Theorem-like environments: theorem, lemma, proposition, corollary, conjecture in `thmstyleone`; definition, remark, example, observation in `thmstylethree`; all share one counter numbered within sections (as before). Open Problem has its own counter within sections. Appendix sections are lettered (A, B, ...).
- Order follows the template: Statements and Declarations, appendices, bibliography.
- Headings: at most three levels (section, subsection; no subsubsection). The former `\paragraph` run-in headings are now `\runinhead` (italic run-in text, not a sectioning level).
- No footnotes, no figures. Acknowledgements: none (a title-page section would be added if needed).
- The appendices live in the `appendices` environment of main.tex (`\input` of `app_*.tex` and the appendix sections `optimal`, `components`, `nl`, `locality`, `repack`). They can be split off as an "Online Resource" PDF by moving that block to a second document that `\input`s the same files and uses `xr` for cross-references.
- Declarations contain a commented placeholder `% LLM-use statement (author decision pending)`; no statement about AI tools was added.

## Checklist for the Algorithmica portal (guidelines pages were not readable; verify)
- [ ] Springer Nature LaTeX template is accepted/recommended; PDF accepted for first submission.
- [ ] Abstract 150-250 words (now about 240 incl. formulas), no undefined abbreviations, no citations. Keywords 4-6 (now 6). MSC codes accepted/required.
- [ ] Numbered citations in square brackets, DOIs as https://doi.org links (present for all but the arXiv-only item, marked as arXiv).
- [ ] Declarations section before the reference list: Funding, Competing interests, Data/Code availability, Author contributions; edit Funding if a grant applies.
- [ ] Decimal headings, max three levels.
- [ ] Whether the portal wants a separate title page / cover letter (cover_letter.md and highlights.md exist in paper_algorithmica/), and whether suggested reviewers are requested.
- [ ] Supplementary material limits (appendices are 70+ pages; check whether an Online Resource is preferred).
- [ ] No simultaneous submission elsewhere; no submission or publication fee; open access optional after acceptance.
- [ ] Author data: ORCID, affiliation, e-mail, corresponding author flag; confirm whether the portal needs a postcode for the affiliation.
- [ ] Author decision on the LLM-use statement.
