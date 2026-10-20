# Venue research: "Cluster Editing above a P3-packing"

Research date: 2026-10-20 (all "verified on" entries use this date). Tools: web search and page fetch only.
Legend: **V** = read on the official or primary page on the verified date. **S** = seen only in a search snippet or a third-party listing, treat as unconfirmed. **U** = unverified, inference or from memory, do not rely on it.
Paper profile: about 27 pages main text plus 69 pages appendix, parameterized dichotomy (W[1]-hardness, XP algorithm, ETH lower bound), partial Lean verification.

## 0. Headline findings

1. Three of the candidate deadlines have already passed as of 2026-10-20:
   - IPEC 2026: 2026-06-29.
   - SODA 2027: 2026-07-09 (S).
   - STACS 2027: 2026-10-11, 23:59 AoE. It closed 9 days ago.
2. ICALP 2027 Track A (the only Track A round) and ESA 2027 have no published dates yet.
3. WG 2027 is the only open conference with a verified date: paper due 2027-02-05 (abstract 2027-01-31).
4. Journals (TALG, Algorithmica, JCSS, TCS, SIDMA, JGAA, DAM, I&C, TOCS) have no deadline. Page-level policies for most were not readable (HTTP 403 or login redirects). Those rows are marked U or S.
5. Related work: a WG 2026 LIPIcs paper (vol. 376, paper 17) appears to study Correlation Clustering with kernels by structural parameters of a "fuzzy" edge graph. Check it for overlap before submitting (see section 5).

## 1. Conference table

| Venue | Deadline / notification | Page limit and appendix | Template | Anonymity | arXiv | OA / cost | Source | Verified on |
|---|---|---|---|---|---|---|---|---|
| **IPEC 2026** (L'Aquila, 2-4 Sep 2026, part of ALGO) | Deadline 2026-06-29 23:59 AoE, **passed**. Notification 2026-08-02. | No hard limit. First 500 lines (about 12 pages) should suffice for review, the rest is read at PC discretion. Final version: 15 pages plus bibliography, **no appendix**. | LIPIcs via socg-lipics-v2021.cls, `anonymous` option | Lightweight double-blind | Allowed (drafts on arXiv or web) | LIPIcs, open access (post-proceedings, camera-ready end of Sep). APC not stated (V: none mentioned; LIPIcs is normally free for authors, U). | https://algo-conference.org/2026/ipec/ , https://easychair.org/cfp/IPEC2026 | 2026-10-20 (V) |
| **IPEC 2027** | Not announced. Pattern from 2026 suggests about late June 2027 (U). | Expect same as 2026 (U) | LIPIcs (U) | Expect double-blind (U) | Expect allowed (U) | LIPIcs (U) | https://algo-conference.org/2026 (no ALGO 2027 page yet) | 2026-10-20 (not announced) |
| **ESA 2026, Track A** (reference for ESA 2027) | Abstract 2026-04-21, paper 2026-04-23 AoE, notification 2026-06-26 | 500 lines excluding front matter, references and a clearly marked appendix. Appendix may be the full version, read at PC discretion, **not published**. | socg-lipics-v2021.cls | **Not anonymous** (single-blind) | Not mentioned on page (others OK, dual submission to proceedings or journals barred) | LIPIcs, open access | https://algo-conference.org/2026/esa/ | 2026-10-20 (V) |
| **ESA 2027** | **No ALGO 2027 or ESA 2027 page found.** If ESA 2026 timing is repeated, expect about April 2027 (U). | Expect as 2026 (U) | LIPIcs (U) | Expect single-blind (U) | Expect fine (U) | LIPIcs (U) | https://algo-conference.org/2026 (nav has no 2027 entry) | 2026-10-20 |
| **ICALP 2027** (Krakow, 5-9 Jul 2027) | Track A: single round, **dates all TBA**. Track B round 1: 2026-09-15, notification 2026-12-08. Deadlines firm. | **No page limit.** First 10 pages after the title page carry merits, contributions and key ideas. Everything else (appendix included) is read at committee discretion. | LIPIcs class recommended, pdflatex mandatory, 11pt, single column, 1-inch margins. Large deviations risk desk rejection. | **Lightweight double-blind** (third-person self-reference) | **Allowed** (arXiv, web, talks) | LIPIcs (S) | https://icalp.uj.edu.pl/track-b/ (page shared by both tracks) | 2026-10-20 (V for rules, dates TBA) |
| **STACS 2027** (Goettingen, 8-12 Mar 2027) | Deadline 2026-10-11 23:59 AoE, **passed**. Rebuttal 30 Nov-2 Dec. Notification 2026-12-22. Track A includes parameterized algorithms. | Rules not on the page read. Earlier editions: 15 pages (2024), 12 pages (2021) excluding refs and appendix, appendix at PC discretion (S, other years). | LIPIcs (earlier editions, S) | Not stated (U) | Not stated (U) | LIPIcs (U) | https://events.gwdg.de/event/1460/ | 2026-10-20 (V for dates) |
| **SODA 2027** (Philadelphia, 24-27 Jan 2027) | Deadline 2026-07-09 (UTC-12), **passed**. Aggregators conflict (one says 2026-06-10). Notification not found. | Not found (official SIAM page not reached). | SIAM format (U) | Not found (U) | Not found (U) | Not found | https://www.getpaperpilot.com/deadlines/soda-2027.html , https://www.myhuiban.com/conference/383.md | 2026-10-20 (S, third-party only) |
| **MFCS 2027** | **No 2027 call found.** MFCS 2026 deadline was 2026-04-24 (S), so expect about April 2027 (U). | Not found | LIPIcs in recent years (U) | Not found | Not found | LIPIcs (U) | https://lics.siglog.org/newsletters/272.html (2026 listing) | 2026-10-20 (S) |
| **WG 2027** (Schloss Schney, 15-18 Jun 2027) | Abstract 2027-01-31, **paper 2027-02-05**, notification 2027-04-09, final 2027-05-10 (all 23:59 AoE) | **500 lines** (cls count), about 12 pages, excluding references. Violating or evading the limit may mean desk rejection. Clearly marked appendix allowed, **no length limit**, read at reviewers' discretion. | LIPIcs via socg-lipics-v2021.cls (`anonymous` option) | Lightweight double-blind | Allowed (arXiv, web, talks) | LIPIcs (S: the page mentions LIPIcs guidelines only) | https://www.uni-bamberg.de/en/wg2027/submissions/ , https://www.uni-bamberg.de/en/wg2027/call-for-papers/ | 2026-10-20 (V) |
| **SWAT** (2026 edition as reference) | SWAT 2026: submission 2026-02-20, notification 2026-04-10, Copenhagen 17-19 Jun 2026. 2028 edition not announced. | Not found for 2026 | LIPIcs (S) | Not found | Not found | LIPIcs | https://swat2026.compute.dtu.dk/ | 2026-10-20 (S) |
| **FCT 2027** | Nothing found | - | - | - | - | - | - | **Unverified** |

## 2. Journals

Most journal pages were not readable by the tools (403 or login redirect). Treat everything marked S or U as needing a check on the publisher's page before submission.

| Journal | Template / length | Review model | OA / APC | arXiv preprint | Fit and precedent | Source | Verified on |
|---|---|---|---|---|---|---|---|
| **ACM Trans. on Algorithms (TALG)** | ACM format. Page cap not found (U). The Li-Pilipczuk-Sorge paper is 43 pages. | Not found (U) | ACM went fully open access from 2026-01-01. Options are ACM Open (institutional agreement) or an APC. Waivers are rare. Check whether BUV or the corresponding author's institution is an ACM Open member. APC amount not found. | Sherpa summary says the submitted version may go in repositories under conditions (arXiv not named specifically). | **Precedent for exactly this problem line** (see section 5). Long-form papers accepted. | https://v2.sherpa.ac.uk/id/publisher_policy/21 , https://www.acm.org/articles/pubs-newsletter/2025/blue-diamond-january-2024 | 2026-10-20 (S) |
| **Algorithmica** (Springer) | Springer template. No fixed page cap found (U). | Not found (U) | Hybrid OA, CC BY option. APC amount not found. Subscription route has no fee. | Springer generally permits preprints (U) | Regular IPEC special issues exist (IPEC 2016, 2021, 2022). Cao-Chen "Cluster editing: kernelization based on edge cuts" appeared here (2012). No 2026 special issue call found. | https://oapublishing.mpdl.mpg.de/journals/Record/Springer-453 , https://www.springer.com/journal/453/volumes-and-issues/83-8 , https://link.springer.com/doi/10.1007/s00453-011-9595-1 | 2026-10-20 (S) |
| **JCSS** (Elsevier) | Elsevier. Page cap not found. | Not found (U) | Hybrid. One aggregator lists APC USD 2,520, another lists N/A, so conflicting. | Not found (Elsevier normally allows preprints, U) | Fomin et al. "Tight bounds ... Cluster Editing with a small number of clusters", JCSS 80(7), 2014. Good for long dichotomy papers. | https://www.journalmetrics.org/journal/journal-of-computer-and-system-sciences , https://core.ac.uk/display/29192941 | 2026-10-20 (S) |
| **Theory of Computing Systems** (Springer) | Springer | Not found | Hybrid (U), APC not found | Not found (U) | Van Bevern-Froese-Komusiewicz journal version is here (see section 5). Runs conference special issues (CSR 2019 and 2020, STACS 2021). No IPEC or MFCS 2026-27 call found. | https://dblp1.uni-trier.de/db/journals/mst/mst67.html (search listing) | 2026-10-20 (S) |
| **SIAM J. Discrete Math (SIDMA)** | SIAM format. Length cap not found for SIDMA (sister journals: 20 recommended, 40 hard cap, S). | Not found | Optional one-time OA fee of USD 3,750 on SIAM author page (not SIDMA-specific, S). | Sherpa: submitted version may go to a preprint repository. | Plausible but graph-theoretic and discrete-math orientation, less typical for FPT dichotomies. | https://epubs.siam.org/journal-authors , https://v2.sherpa.ac.uk/id/publisher_policy/71 | 2026-10-20 (S) |
| **J. Graph Algorithms and Applications (JGAA)** | Not found | Peer-reviewed, model not stated | **Diamond OA, no author fees** (V) | Not found | Good fit for graph algorithms, accepts Regular, Concise and Survey articles. No cluster-editing precedent checked. | https://www.jgaa.info/ | 2026-10-20 (V for OA and fees) |
| **Discrete Applied Mathematics** (Elsevier) | - | - | - | - | **Unverified, not researched** | - | - |
| **Information and Computation** (Elsevier) | - | - | - | - | **Unverified, not researched** | - | - |

## 3. Fit analysis for this paper

- **Appendix tolerance.** Best fits are ICALP (no limit, first 10 pages must carry the story), WG (500 lines plus unlimited appendix), ESA (500 lines plus appendix that may be the full version, but not published). IPEC reviews "the full version" but the final version has no appendix. A 69-page appendix with the Lean material will be read at most partially at any conference, so the main text must be self-contained.
- **Compression need.** ICALP needs a 10-page core, ESA and WG need about 12 pages (500 lines). The current 27-page main text must be cut to about 12 for either.
- **Anonymity.** ICALP, IPEC and WG require a lightweight double-blind version (no names, affiliations or emails, third-person self-citation). ESA is not anonymous. Author name: Quang-Vinh Dang, British University Vietnam (per researcher information). Note the account email domain (rmit.edu.vn) differs from the BUV affiliation; confirm the intended affiliation line before filling submission forms.
- **arXiv.** Allowed at ICALP, IPEC and WG (V). ESA page is silent (V), dual-submission rule only bars other proceedings or journals.
- **Topical fit.** A parameterized dichotomy paper with W[1]-hardness, XP, ETH and CSP representation is a natural IPEC paper, and plausible for ICALP Track A, ESA Track A and STACS Track A. The soft-triangle-free polynomial case (majority closure) adds a CSP flavor useful for ICALP.

## 4. Ranking and recommendation

1. **ICALP 2027, Track A** (recommended). Reasons: strongest appendix policy of the major venues (no limit, 10-page core, discretionary rest), arXiv explicitly allowed, July 2027 timing, and good prestige for a dichotomy result. Risk: the Track A deadline is still TBA. By analogy with ICALP 2026 (deadline in early February, from a snippet, S) expect about early February 2027 (U). Single round only, so no second chance. Act: watch https://icalp.uj.edu.pl and the EATCS announcements.
2. **ESA 2027, Track A.** Reasons: parameterized algorithms are welcome, non-anonymous review, appendix may hold the full version. Risks: appendix unpublished, 500-line limit forces heavy cuts, no 2027 call yet (expected about April 2027, U).
3. **ACM Transactions on Algorithms** (journal, no deadline). Reasons: precedent for this exact problem line (43-page TALG paper), tolerates long text, does not force a 12-page cut. Risks: 96 pages in total is far beyond normal even for journals, so the appendix would need trimming or splitting. APC and ACM Open status of the author's institution must be checked.

Fallbacks: **WG 2027** (verified deadline 2027-02-05, unlimited appendix, but lower prestige and a safe place if the ICALP deadline slips). **IPEC 2027** (best topical fit and a likely Algorithmica special issue, but the final version excludes the appendix, and the date is unknown).
Do not plan on STACS 2027, SODA 2027 or IPEC 2026: all closed.

Suggested sequence: post to arXiv now (allowed at ICALP, IPEC, WG), prepare an anonymous 10-12 page core, submit to ICALP 2027 Track A, and keep WG or ESA as the next target. Do not submit the same results to two conferences at once (all check dual submission). Then send the full version to TALG or a similar journal after the conference decision.

## 5. Prior work: where it appeared

| Work | Venues | Source | Verified on |
|---|---|---|---|
| Li, Pilipczuk, Sorge, "Cluster Editing Parameterized above Modification-disjoint P3-packings" (arXiv 1910.08517, v1 Oct 2019, v3 Oct 2023) | STACS 2021, LIPIcs vol. 187, 49:1-49:16. Journal: ACM Trans. Algorithms 20(1), Art. 3, 43 pages (published Dec 2023, volume year 2024), DOI 10.1145/3626526. | https://arxiv.org/abs/1910.08517 , https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.STACS.2021.49 , https://repositum.tuwien.at/handle/20.500.12708/208560 | 2026-10-20 (S, from search snippets, DOI not opened) |
| van Bevern, Froese, Komusiewicz, "Parameterizing edge modification problems above lower bounds" (arXiv 1512.04047) | CSR 2016, LNCS 9691, pp. 57-72. Theory of Computing Systems 62(3), 2018, pp. 739-770, DOI 10.1007/s00224-016-9746-5. | https://arxiv.org/abs/1512.04047 , https://pure.nsu.ru/portal/en/publications/parameterizing-edge-modification-problems-above-lower-bounds(d9fe170b-52aa-4b46-9915-fdb99cc1a3d0)/export.html | 2026-10-20 (S) |

Other cluster editing and correlation clustering papers at candidate venues:

| Paper | Venue | Source | Verified on |
|---|---|---|---|
| Fomin, Kratsch, Pilipczuk, Pilipczuk, Villanger, Tight bounds for parameterized complexity of Cluster Editing (with a small number of clusters) | STACS 2013 (LIPIcs vol. 20, pp. 32-43), JCSS 80(7), 2014, 1430-1447 | https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.STACS.2013.32 | 2026-10-20 (S) |
| Cao, Chen, Cluster Editing: Kernelization based on Edge Cuts | Algorithmica 64(1), 2012, 152-169 (arXiv 1008.4250; conference version in IPEC 2010 is likely but unconfirmed, U) | https://link.springer.com/doi/10.1007/s00453-011-9595-1 | 2026-10-20 (S) |
| PACE 2021 Challenge: Cluster Editing report | IPEC 2021, LIPIcs 214, paper 26 | https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.IPEC.2021.26 | 2026-10-20 (S) |
| Correlation Clustering with fuzzy-edge-graph structural parameters (kernels: fuzzy degeneracy, fuzzy closure) | WG 2026, LIPIcs vol. 376, paper 17. Title and authors not extracted (PDF text compressed). **Read before submitting, possible overlap.** | https://drops.dagstuhl.de/storage/00lipics/lipics-vol376-wg2026/LIPIcs.WG.2026.17/LIPIcs.WG.2026.17.pdf (local copy: /root/.claude/projects/-home-user-correlation-clustering/404aaf8c-18f1-5f19-bc11-1c2160450d1f/tool-results/webfetch-1792487737294-ixqeoe.pdf) | 2026-10-20 (partly read) |
| Zhou, Zhang, Guo, FPT constant-factor approximation for Correlation Clustering (arXiv 2503.00281); venue not found | - | https://arxiv.org/pdf/2503.00281 | 2026-10-20 (S, venue **unverified**) |
| Fair Correlation Clustering and graph parameters (arXiv 2602.15683) says FPT by treewidth for Correlation Clustering was previously open (relevant to related-work claims about treewidth); venue not found | - | https://arxiv.org/pdf/2602.15683 | 2026-10-20 (S, venue **unverified**) |
| IPEC 2024: no cluster editing paper in the accepted list (a Modularity Clustering paper only). IPEC 2025 list not found. | IPEC 2024 | https://www.fmi.uni-jena.de/en/12914/ipec-adt-2024 | 2026-10-20 (S) |

## 6. Not verified / open items

- ICALP 2027 Track A deadline and notification (TBA on the official page).
- ESA 2027 and ALGO 2027 (location, dates, call) and IPEC 2027 (no page).
- MFCS 2027, FCT 2027, SWAT 2028: no calls found.
- SODA 2027: official SIAM page not reached, notification date and page rules unknown, deadline conflicts across aggregators.
- STACS 2027 page limit, anonymity and arXiv policy (only dates were readable).
- All journal review models, page limits, current APCs and arXiv terms (publisher pages blocked). DAM and I&C not researched.
- APC for LIPIcs proceedings (assumed free for authors, U).
- Whether the ICALP "first 10 pages" rule is identical for Track A (page is shared and shows Track B dates; the 2026 Track A CFP at https://icalppodcspaa2026.cs.rhul.ac.uk/icalp-cfp has the same wording, V).
