# Homodimer Diagnostic Notebook Tasks

**Status: CLOSED 2026-09-08 by R094.** This is the original build plan (Phases 1-7). It remains
valid history, and every task below is now resolved. **All work from 2026-09-07 onward is
tracked in `specs/homodimer_diagnostic/rework-plan.md`, which supersedes this file.** Do not
add new tasks here.

Phases 2-5 (T006-T032) were completed as originally written. The rework then moved every one of
those implementations out of the notebook and into
`src/insightfold/complex_interface_utils.py` (R011-R016), corrected four of the seven formulas
against `ipsae.py` v4 (R003-R007), and replaced the threshold set (R002, R009). The tick marks
below record that the task was done, not that its output survives unchanged.

Where a task is closed by rework work, the closing task is named. Where a task's deliverable is
a document, the document is named.

## Phase 1: Spec and Fixture Readiness

- [x] T001 Convert PRD into directory-style notebook spec pack.
- [x] T002 Resolve fixture blockers for borderline and disagreement accessions.
  **Closed by R024 + R090.** Both blockers are gone, but not as this task expected:
  - FX-002 was scored and is **not** borderline (`ipSAE_d0res` 0.7699, CONFIDENT). **FX-008**
    (`AF-0000000211619209`, 0.6366) holds the borderline role, from live data. Manifest
    corrected by R094.
  - FX-003, the metric-disagreement fixture that blocked this task from the beginning, is
    closed by **FX-010** (`AF-0000000211026350`), which has a 0.5685 directional spread and
    three of five scores disagreeing. No AFDB/PDBe domain input was needed in the end.
  - The manifest grew 5 -> 10 fixtures; five of the new ones are heterodimers, which this task
    never contemplated.
- [x] T003 Record numeric expected outputs for FX-001 after trusted reference run.
  **Closed by R092**, `specs/homodimer_diagnostic/validation-report.md`. Better than asked for:
  the expected outputs are recorded for **all seven scored fixtures**, against **two**
  independent references (`ipsae.py` v4 run as a subprocess, and AlphaFold DB's own published
  `complexPredictionAccuracy_*` production values). 231 comparisons, none outside +-0.001.
- [x] T004 Decide high/moderate/low threshold bands for diagnostic summary.
  **Closed by R002 and R009**, `specs/homodimer_diagnostic/threshold-reference.md`, which is
  now the single source of truth and is transcribed into the module's `THRESHOLDS`.
  The decision differs from what this task assumed in three ways worth recording:
  1. **The headline verdict is not a band at all.** It is AFDB's joint criterion,
     `ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23` -> PASS/FAIL (`afdb_high_confidence()`). The
     seven traffic lights are diagnostics that explain that verdict, not seven independent ones.
  2. **`ipSAE_d0res` carries four *published* bands**, not three invented ones: very
     high-confidence >= 0.8, confident 0.7-0.8, low-confidence 0.6-0.7, below threshold < 0.6
     (AFDB 2026, p. 12). The other six values keep the three-colour scheme.
  3. **Only two of the fourteen numbers were ever published**, and both surviving ambers
     (`pdockq2` 0.10, `lis` 0.10) are judgement calls. The document says so per number rather
     than presenting all fourteen as equally grounded.
- [x] T005 Confirm IPSAE reference version or commit for validation.
  **Closed by R001.** `DunbrackLab/IPSAE` `ipsae.py` **v4**, 3 Jan 2026, 1010 lines,
  md5 `a48df7adc64afa36b4b475429b0a8aff`,
  sha256 `10cf9b08c68c91e06cb28526cf2026f47a3980c9048fd3226d13e3304eaf1c27`. Stored at
  `specs/homodimer_diagnostic/references/ipsae_v4.py` (gitignored). Pinned in
  `formula-reference.md`, which line-cites every formula to it.

## Phase 2: Data Contracts

All six were completed as written, then moved into `complex_interface_utils.py` by **R011**
(fetch and parsing) and corrected by R020-R023.

- [x] T006 Implement AFDB metadata fetcher and response validation.
  Since superseded by **R020**: the endpoint returns one entry **per chain** in
  **non-deterministic order**, so the original `entries[0]` was not merely wrong but
  non-deterministically wrong. `AFDBPrediction` now sorts by `chainId`.
- [x] T007 Implement download helpers for text and JSON resources.
- [x] T008 Implement PAE JSON parser and schema checks.
- [x] T009 Implement confidence JSON parser and pLDDT validation.
- [x] T010 Implement lightweight mmCIF atom-site parser.
- [x] T011 Implement homodimer compatibility validation.
  Rescoped by **D3 and R023** to *dimer* compatibility: homodimers and heterodimers are both
  supported automatically. A non-dimer now raises `UnsupportedAssemblyError` before any
  download.

## Phase 3: Notebook Architecture

- [x] T012 Create notebook introduction, RUO framing, and audience guidance.
- [x] T013 Create setup/import section with allowed dependencies only.
  Replaced by the **D9 / R010 / R010b** bootstrap: `git clone --depth 1` plus `sys.path`,
  never `pip install`. See the `TODO(merge)` warning in `CLAUDE.md`.
- [x] T014 Create accession input and provenance display section.
- [x] T015 Create validation summary section for fetched data.

## Phase 4: Scientific Computation

All completed as written, then moved into the module by **R012** and corrected by R003-R007.
Four of the seven formulas needed a fix; see `formula-reference.md`.

- [x] T016 Implement residue table construction with CA and CB/CA fallback coordinates.
- [x] T017 Implement CB-CB inter-chain contact detection with 8.0 Angstrom cutoff.
- [x] T018 Implement `d0` and `ptm` shared functions.
  **Corrected by R005:** `ipsae.py` has **two** `d0` helpers, not one. `d0_array` is required
  for `d0res` and differs from `d0_scalar` at `L == 27`.
- [x] T019 Implement ipTM_d0chn calculations and intermediate table.
  **Renamed by R004.** It was labelled plain `ipTM`, which is a different number.
- [x] T020 Implement ipSAE_d0res, ipSAE_d0chn, and ipSAE_d0dom calculations.
  **Corrected by R003:** `n0dom` is per direction.
- [x] T021 Implement pDockQ calculation.
- [x] T022 Implement pDockQ2 calculation.
  **Corrected by R007:** the reported value is `max(A->B, B->A)`, not `A->B`.
- [x] T023 Implement LIS calculation.
- [x] T024 Build final score table.
  Rebuilt by **R080** around one `ConfidenceSummary` object, so the table and the prose cannot
  disagree.

## Phase 5: Visualizations and Interpretation

All completed as written, then moved into the module by **R013**/**R014** and reworked by
M5 (R050-R052), M6 (R030, R070-R075) and M7 (R080-R082).

- [x] T025 Implement contact map visualization.
- [x] T026 Implement interface coverage visualization.
- [x] T027 Implement full PAE heatmap with quadrant labels.
- [x] T028 Implement score-specific PAE masks.
- [x] T029 Implement per-residue contribution profiles.
- [x] T030 Implement interface versus non-interface pLDDT comparison.
- [x] T031 Implement minimum MolViewSpec chain overview and pLDDT mapping views.
  Grew to **six** views (R030, R071-R075), all on both chains, all with a legend.
- [x] T032 Implement plain-language diagnostic summary.
  **Rebuilt by R080/R081/R082.** The version this task delivered could print "consistently HIGH
  confidence across all metrics" while its own table showed an amber score.

## Phase 6: Notebook Validation

**All six closed by R090 + R091 + R092**, delivered as
`specs/homodimer_diagnostic/validation-report.md` (820 lines), the report T038 called for and
never got. Run of 2026-09-08 against `homodimer-notebook-rework @ e162059`.

- [x] T033 Run static notebook validation.
  Superseded in practice: every task in the rework was verified by **execution** against a live
  AFDB fetch (D5), which is a stronger check than static validation. 13 runs across all ten
  registered fixtures, zero errors and zero bytes to stderr on every successful run.
- [x] T034 Run smoke execution with FX-001. **R090.** Clean, no warnings.
- [x] T035 Run negative fixture validation with malformed accession.
  All **four** negative paths refuse cleanly before any score is produced: malformed accession
  (HTTP 400), absent accession (404), non-dimer assembly (`UnsupportedAssemblyError` before any
  download), and local-file mode with documents missing (fails at upload, naming each one).
- [x] T036 Run full validation once FX-002, FX-003, and FX-005 are curated.
  Ran on **all ten** fixtures. FX-002 scored; FX-003 replaced by FX-010 (see T002).
- [x] T037 Compare IPSAE-family values with reference implementation within +/- 0.001.
  **Passed with room to spare, against two references.** Worst |delta| vs `ipsae.py` v4:
  **4.7e-5**, and every delta above 1e-6 falls on the three values the reference prints at 4 dp;
  the 6 dp ones agree to <= 6e-7. Worst |delta| vs AFDB production values: **2.4e-6**, roughly
  450x tighter than the tolerance. `n0dom` matches integer-for-integer per direction on all five
  heterodimers.
- [x] T038 Record validation report in `specs/homodimer_diagnostic/validation-report.md`.
  **Delivered.**

## Phase 7: Final Review

- [x] T039 Run notebook review after execution validation.
  **Satisfied by the rework's own protocol rather than a single review pass**, and the
  substitution was the stronger option: 46 tasks, each dispatched to a fresh agent, each result
  read as a diff, executed against a live AFDB fetch, and checked against its own "done when"
  before the next was dispatched. Work stopped at every milestone boundary for user review. The
  findings are recorded per task in `rework-plan.md`; R092 additionally found **no bug** in the
  notebook, the module or `interface.py`.
- [x] T040 Resolve blocking or major review findings.
  Every finding raised during the rework is either fixed with its evidence recorded in
  `rework-plan.md`, or carried forward explicitly. **Nothing blocking remains open.** The two
  carried-forward items are both in `CLAUDE.md`:
  1. **`TODO(merge)`** — the Colab bootstrap still pins `REPO_BRANCH` to
     `homodimer-notebook-rework`. Must be flipped to `main` or a release tag on merge.
  2. **Colab is unverified.** Every run had `IN_COLAB = False`, so the clone, stale-clone
     refresh and `pip install molviewspec` paths have never executed and the 60 s install
     budget is unmeasured. Cannot be closed from this environment, and cannot be usefully
     tested before the `TODO(merge)` flip anyway.
- [x] T041 Decide whether notebook is ready for beta, continued prototype iteration, or archive.
  **Recommendation: beta, conditional on the two items in T040.** The evidence for it:
  the notebook executes end to end with zero errors and zero warnings on six online fixtures
  and the local-file path; all seven values agree with two independent references on all seven
  scored fixtures; 299 module doctests pass; the code is 54% smaller than the prototype and
  lives in a tested module rather than in cells. The evidence against calling it released:
  the Colab path, which is the path a beta user will take, has never been run.
  **This is a recommendation, not a decision — the call is the user's.**
