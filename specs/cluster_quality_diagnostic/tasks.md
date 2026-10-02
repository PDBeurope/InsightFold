# Tasks

`[P]` marks tasks that can run in parallel without touching the same files or sharing assumptions.

## Phase A: Spec Preparation

- [ ] T001 Run `$notebook-spec-review` over this spec pack and resolve its findings
- [x] T002 Run `$fixture-selection` to confirm and freeze the fixture manifest. **Done 2026-09-21**: 12 cluster responses, 15 prediction responses and 5 synthetic files pinned under `specs/cluster_quality_diagnostic/fixtures/`, 3.15 MB gzipped, with `_expected.json` and `_alignments.json` holding measured expected values
- [x] T003 Re-test `P12345`. **Done**: reproduced HTTP 500 at 2.6 s spacing, so OQ-4 resolves to a genuine per-accession backend fault, not rate limiting

## Phase B: Data And Fixtures

**Phase B is complete.** Remaining open item: no cluster below 30 members has been located for the low-n caveat; construct one synthetically during implementation.

- [x] T010 Pin cluster API JSON. **Done**: 12 responses gzipped under `fixtures/cluster/`
- [x] T011 Pin prediction API responses. **Done**: 15 responses under `fixtures/prediction/`
- [x] T012 Pin refusal-gate responses. **Done**: FX-08, FX-09, FX-10
- [x] T013 Record expected output snapshots. **Done**: `_expected.json` and `_alignments.json`, all values computed from the pins
- [ ] T014 Implement a **built-in default residue colour scheme** in the module that needs no external file. `cleancolors.json` is an optional enhancement loaded only when present and readable, pending a licence check
- [x] T015 Build FX-11 `SYN-COLLAPSE`. **Done**
- [x] T016 Build FX-13 `SYN-A3M` and its wrong-sequence variant. **Done**
- [x] T017 Investigate the deleted-member path. **Done, and the premise was wrong**: all 20 of `Q9I1F6`'s top members were checked; 14 are deleted from UniProtKB yet AFDB serves a sequence for every one. FX-12 is therefore synthetic, and the manifest records that no real trigger is known

## Phase C: Notebook Architecture

- [x] T020 Create `src/insightfold/cluster_quality_utils.py` with the section structure from `notebook-design.md` **Done.** module written
- [x] T021 Define the dataclasses in the variable handoff table: `ClusterMembers`, `QueryMeta`, `FamilySummary`, `Position`, `ViewTarget` **Done.** QueryMeta, ClusterMembers, FamilySummary, Position, ViewTarget, AlignmentResult, StructureSource, ColourRun
- [x] T022 Define API base URLs as single constants (NFR-008) with per-endpoint rate policy attached to each fetch helper, never to shared calling code (NFR-003) **Done.** single constants; spacing inside each fetch helper
- [x] T023 Create the notebook skeleton `notebooks/cluster_quality_diagnostic.ipynb` with the cell IDs from `cell-blueprint.md` as markdown placeholders **Done.** notebooks/cluster_quality_diagnostic.ipynb, 37 cells

## Phase D: Implementation

- [x] T030 Implement accession normalisation and validation: **normalise** the `AF-{ACC}-F1` form and print the interpretation, reject malformed input, before any network call (REQ-001) **Done.** normalises the AF form and prints the interpretation
- [x] T031 Implement `fetch_prediction()` against the prediction endpoint, no rate spacing (REQ-001, NFR-003) **Done.** no spacing
- [x] T032 Implement `fetch_cluster()` against the cluster endpoint with enforced 2.5 s spacing and 5xx backoff (REQ-002, REQ-014, NFR-003) **Done.** 2.5 s spacing enforced internally, 5xx backoff
- [x] T033 Implement the refusal gates: fragmented, too small, malformed, absent, 5xx. **Degrade per clustering**: refuse the run only when the sequence cluster is unusable (REQ-014, REQ-014a) **Done.** degrades per clustering
- [ ] T033a Implement the prediction-endpoint failure paths from `data-contracts.md` §2: 404, 200-without-sequence, 5xx (REQ-011)
- [ ] T033b Implement the low-n caveat for clusters below 30 members (REQ-004)
- [ ] T033c Implement isoform substitution reporting: print which sequence is aligned and which structure rendered (REQ-001a)
- [x] T034 Implement `family_summary()`: mean, median, sd, band composition, p25-p75 width (REQ-003) **Done.** verified against pins
- [x] T035 Implement `locate_query()` handling both presence and absence, returning a category rather than only a percentile (REQ-004, REQ-006) **Done.** handles absence; p5/p95 category
- [ ] T036 Implement `resolve_representative()` for the Foldseek hop (REQ-006)
- [x] T037 Implement `identify_extremes()`, asserting rather than assuming the descending sort, and sorting defensively with a logged warning if the assumption fails (REQ-009) **Done.** sort asserted, defensive re-sort
- [x] T038 Implement the collapse rule once, in `identify_extremes()`, returning 2 or 3 `ViewTarget`s consumed by both the 3D and alignment sections (REQ-010, REQ-011) **Done.** computed once in identify_view_targets
- [x] T039 [P] Implement `plot_distribution()` as a seaborn raincloud, replacing the source notebook's ptitprince version (REQ-V1) **Done.**
- [x] T040 [P] Implement `plot_structure_cluster_distribution()`. Axis labelled "structurally similar proteins; your protein is not among them". **No query reference line and no representative marker**: the member sets are disjoint from the sequence cluster, so the query is normally in neither (REQ-V2, REQ-006) **Done.** no query marker
- [x] T041 [P] Implement `plot_length_histogram()` and `plot_length_vs_plddt()` (REQ-V3, REQ-V4) **Done.**
- [x] T044 Implement `align_pair()` using Biopython `PairwiseAligner(mode='local')`, BLOSUM62, gaps −11/−1, returning identity, aligned length and query coverage, with the convention recorded (REQ-011) **Done.** local, reports identity, columns, coverage and aligned range
- [x] T044a Implement the mandatory text output of identity, aligned length and coverage, outside any try/except that guards the figure (REQ-011a) **Done.** outside the figure try/except
- [x] T044b Implement `plot_plddt_profile()` from the downloaded mmCIF, one per view target (REQ-010a, REQ-V8) **Done.**
- [x] T045 Implement `plot_alignment()` via pyMSAviz, falling back to the built-in colour scheme when `cleancolors.json` is absent, and degrading to a named message when pyMSAviz is unavailable (REQ-V7) **Done.** pyMSAviz with fallback
- [x] T046 Implement the Mol\* views through `$molviewspec-rendering`, pLDDT-coloured, lazy import, at most 3 viewers (REQ-V6, REQ-010) **Done.** copies the complex_interface_utils.py pattern
- [x] T046a Implement `superpose()` over `tmtools.tm_align`, returning the **inverse** transform so the query stays in its deposited frame, plus TM-score, RMSD and TM-align's own alignment (REQ-017) **Done.** orthogonality and determinant asserted; direction verified against the library's own reported RMSD
- [x] T046b Implement `build_superposition_view()`, passing the rotation to MolViewSpec **column-major** (REQ-017, REQ-V11) **Done.** transform node confirmed present in the emitted state, applied to the member only
- [x] T046c Implement `parse_secondary_structure()` from `_struct_conf`, collapsing helix subtypes to `H`, `STRN` to `E`, the rest to `T` (REQ-017c) **Done.** no external tool, no extra download
- [x] T046d Implement `plot_structural_alignment()` with `colour_by` of `scheme` or `plddt`, plus a secondary structure track per sequence, through pyMSAviz `set_custom_color_func` (REQ-017b, REQ-V12) **Done.** returning `None` from the custom func falls through to the scheme, so one code path serves both panels
- [x] T046e Refactor `fetch_model_plddt()` onto `fetch_verified_structure()`, so the superposition reuses **one** completeness-checked download for coordinates, pLDDT and secondary structure (NFR-012) **Done.**
- [x] T046f Append the superposition cells after section 7 **without editing any existing cell**, so the section can be deleted whole if it misbehaves (REQ-017) **Done.** all 34 pre-existing cells verified byte-identical after insertion
- [ ] T047 [P] Implement `fetch_msa()` with 403 detection by status and content type, never by parsing the body (REQ-012)
- [ ] T048 [P] Implement user-supplied a3m validation by gap-stripped sequence match (REQ-013)
- [ ] T049 [P] Implement `plot_msa_coverage()` in matplotlib, ColabFold style, stripping a3m lowercase insertions (REQ-V9)
- [x] T050 Write the bootstrap cell per D9: `git clone --depth 1` plus `sys.path.insert`, never `pip install` the package, guarded against duplicate path entries (NFR-001, NFR-007) **Done.** guarded against duplicate sys.path entries
- [x] T051 Write the install cell for molviewspec, biopython and pymsaviz, as a no-op when they already import (NFR-001) **Done.** no-op when packages present
- [x] T052 Write all markdown caveat cells, each placed before the output it bounds (REQ-015, REQ-E1 to REQ-E6) **Done.**
- [ ] T052a Write C025's prose: cascaded clustering and transitive linkage (REQ-E7), framed as a known phenomenon with citations to Steinegger and Soding 2018 and Barrio-Hernandez et al. 2023 (REQ-E9)
- [ ] T052b Write C012's disorder caveat: average pLDDT over a chain is substantially a function of disorder content, so a low family mean is not necessarily poor prediction (REQ-E8)
- [x] T053 Write the closing summary cell (REQ-016) **Done.**
- [x] T054 Implement `apply_plot_style()` as explicit, never applied on import (NFR-005) **Done.** not needed: no global state is mutated

## Phase E: Validation

- [ ] T060 Build the offline test harness at `tests/cluster_quality/`, using `pytest`, with fixtures located by a module-level path constant pointing at `specs/cluster_quality_diagnostic/fixtures/`. It replays pinned JSON against `src/insightfold/cluster_quality_utils.py`; the notebook keeps no fixture-loading path. Add `pytest` as a dev-only dependency, not a notebook one
- [ ] T061 Implement the snapshot assertions from T013 (REQ-003, REQ-004, REQ-009)
- [ ] T062 Assert FX-02's measured alignment against pinned sequences: **26.73% identity over 101 aligned columns, covering 21.3% of the query, residues 80-171 of 414**. Assert the aligned range, which is the stable signal; coverage is parameter-sensitive and must not be asserted as a threshold (REQ-011, REQ-E4)
- [ ] T063 Assert FX-05 (33.68%, 84.4% coverage) and FX-03 (25.09%, 74.6% coverage, residues 18-277 of 343) are reported as distant homology with no claim of an unrelated protein. FX-03's identity is **lower** than FX-02's, which is the point (REQ-011, REQ-E4)
- [ ] T064 [P] Test every refusal gate against its pinned response (REQ-014)
- [ ] T065 [P] Test collapse-rule consistency between the 3D and alignment sections (REQ-010, REQ-011)
- [ ] T066 [P] Test soft failure for every optional section: taxonomy, Mol\*, alignment figure, MSA (REQ-007, REQ-008, REQ-012)
- [ ] T067 [P] Test the missing-sequence path using a deleted-UniProt member as an extreme (REQ-011)
- [ ] T068 Instrument and assert per-endpoint request spacing (NFR-003)
- [ ] T069 Restart-and-run-all on every non-refusal fixture (NFR-002)
- [ ] T070 Run all twice in one kernel and confirm identical values and no duplicated `sys.path` entry (NFR-002)
- [ ] T071 Grep for prohibited imports, confirm lazy imports, and **check pyMSAviz's resolved transitive dependency tree** for anything prohibited (NFR-007)
- [ ] T071a [P] Test the collapse rule against `SYN-COLLAPSE`, both the index-0 and index-−1 cases (REQ-010, REQ-011)
- [ ] T071b [P] Test the MSA coverage plot against `SYN-A3M` so the plotting code is not shipped unexecuted (REQ-V9)
- [ ] T071c [P] Test that identity, aligned length and coverage still print with pyMSAviz uninstalled (REQ-011a)
- [ ] T071d [P] Test the alignment figure with `cleancolors.json` absent, confirming the built-in scheme is used (REQ-V7)
- [ ] T071e [P] Test partial degradation on `P00533`: sequence analysis completes, Foldseek section skipped (REQ-014a)
- [ ] T075a Assert the superposition transform on a pinned pair: orthogonal, determinant +1, and reproducing the library's reported RMSD when applied (REQ-017)
- [ ] T075b Assert `parse_secondary_structure()` against a pinned mmCIF, including a file with `STRN` records and one without (REQ-017c)
- [ ] T075c Test that the TM-score, RMSD and identity still print with `tmtools` present but pyMSAviz uninstalled (REQ-017a)
- [x] T076 Convert every user input to a Colab form field (REQ-018) **Done 2026-10-02.** accession textbox; model-source and colour-scheme dropdowns; PAE checkbox; MSA path textbox
- [ ] T076a Validate the `#@param` annotations as part of the offline harness: title first, widget type matching the assigned literal, dropdown containing its selected value, and the scheme list equal to `ALIGNMENT_SCHEME_CHOICES` (REQ-018, REQ-018a)
- [ ] T076b **Confirm the forms render as widgets on Colab.** The annotations are inert locally, so this is the one property that cannot be checked off Colab. Folds into T074
- [ ] T072 Confirm every plot function returns a bare Figure and the pyplot registry stays empty (NFR-005)
- [ ] T073 Grep for em dashes, `&mdash;` and `&#8212;`, including in generated HTML (NFR-009)
- [ ] T074 **Run the notebook on Colab free tier and record the evidence.** Not inferred; the repo has shipped an unverified Colab path before (NFR-001)
- [ ] T075 Confirm the offline harness passes with the network disabled

## Phase F: Documentation

- [ ] T080 Write the notebook's introductory tutorial cells per `docs-plan.md`
- [ ] T081 [P] Write the how-to section on adapting the notebook to another accession
- [ ] T082 [P] Write the reference section covering inputs, outputs, parameters and APIs
- [ ] T083 [P] Write the explanation section on scientific background and limitations, including why candidate recommendation is out of scope
- [ ] T084 Add a README entry for the notebook
- [ ] T085 Strip notebook outputs before commit (NFR-010)

## Phase G: Review

- [ ] T090 Run `$notebook-execution-validation`
- [ ] T091 Run `$notebook-review`
- [ ] T092 Computational structural biology review of the Foldseek description (OQ-2) and of whether cluster impurity has literature to cite
- [ ] T093 Independent reader review: do the three counter-intuitive explanations land?
- [ ] T094 **Grep for `TODO(merge)` and flip the bootstrap to `main` or a release tag before release**
