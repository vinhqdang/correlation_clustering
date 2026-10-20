# Algorithmica (Springer, journal 453): submission requirements

Checked 2026-10-20. Legend: V = read on an official Springer page; U = unverified.

## Access status

The journal pages are not readable from this environment. `link.springer.com/journal/453`,
`.../submission-guidelines` and `www.springer.com/journal/453/...` redirect to `idp.springer.com`
and then return a JavaScript "Client Challenge" page (about 300 characters of text, no content), both
through WebFetch and through curl with a cookie jar and a browser user agent. No journal-specific item
below could be read, so every journal-specific item is U. Two generic Springer Nature pages were
readable (marked V (generic)); they apply to Springer Nature journals in general, not necessarily to
Algorithmica's own instructions.

## Items

| Item | Finding | Status |
|---|---|---|
| LaTeX template | Springer Nature recommends its "Springer Nature LaTeX authoring template" (zip or Overleaf) for any Springer Nature journal (generic page `springernature.com/gp/authors/campaigns/latex-author-support`; it does not name sn-jnl or svjour3 and does not say whether svjour3 is still supported). Algorithmica's own template instruction not seen. Used: `svjour3` (options `smallextended,envcountsame,envcountsect`), copied from `paper_mpc/`, as for paper 1. | V (generic) / U (journal) |
| Reference style | Not read. Used numeric `plainnat` with natbib `[n]`, as in paper 1. The Springer `spmpsci.bst` is not installed. | U |
| Abstract length | Not read. Abstract was shortened to about 230 words (the full-version abstract had about 350). | U |
| Keywords | Not read. 6 keywords given (svjour3 `\keywords`). | U |
| MSC codes | Not read whether required. Given: 68Q27 (parameterized complexity), 68W40 (analysis of algorithms), 05C85 (graph algorithms), MSC 2020. I set the label in the local `svjour3.cls` copy to "(2020)" (the stock class prints "(2000)"). | U |
| Length policy / page cap | Not read. No cap known. The PDF is 100 pages: main text and declarations on pp. 1-33, references p. 34, appendices pp. 35-100. | U |
| Supplementary / appendix policy | Not read. Appendices are placed after the references. | U |
| Open access / APC | Not read. Prior note (`venues.md`, secondary source): hybrid journal, CC BY option, subscription route without fee. APC amount not found. | U |
| Review model | Not read (single-blind assumed, not verified). | U |
| arXiv / preprint policy | Generic Springer Nature policy: "a unified policy that encourages posting of preprints on preprint servers" and supports citation and open licensing of preprints (page `springernature.com/gp/policies/editorial-policies`). Journal-specific wording not read. | V (generic) |
| Declarations | Not read for this journal. Included, following paper 1: competing interests, funding, data and code availability (plus the author contribution is trivial, single author). | U |
| AI-use disclosure | Generic page says AI use should be transparent and that human accountability cannot be transferred to AI; the page names separate guidance on "AI use in manuscript preparation" and "AI use in research practice" that was not readable. The dedicated AI policy URL tried returned 404. Manuscript contains no statement about AI tools; the author must check the current policy and add any statement it requires before submitting. | V (generic) / U |
| Submission portal | Not read. Springer journals normally use the Editorial Manager system; the journal's own link was not seen. | U |
| Typical review time | Not read (journal metrics page blocked). | U |

## To do by the author before submitting

1. Open the journal's "Submission guidelines" in a normal browser and check: template (sn-jnl vs svjour3), reference style, abstract limit, page or appendix policy, required declarations, AI-disclosure wording.
2. If sn-jnl is required, the sources in `paper_algorithmica/sections/` can be reused unchanged; only `main.tex` needs the new class and the theorem environments.
3. Funding: the declaration says "Not applicable"; a commented placeholder for a grant statement is in `main.tex`.
