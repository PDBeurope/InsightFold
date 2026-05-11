# InsightFold Lifecycle Skill Roadmap

Source inputs:

- `PDBE-Lifecycle-050526-080147.pdf`
- existing InsightFold specs, PRD, agent, and skill files
- `/Users/mitsenkov/PycharmProjects/bioSkills`

## Lifecycle Assessment

The lifecycle document defines InsightFold as a notebook-driven development framework for AFDB and PDBe. The core path is:

1. Idea
2. Triage
3. Scope
4. Notebook prototype
5. Instrumented beta
6. Graduation review

The possible graduation outcomes are:

- integrate into AFDB or PDBe
- maintain as a standing notebook
- archive with a retrospective and preserved code state

Three cross-cutting workstreams span the lifecycle:

- development harness: stages 3 and 4
- feedback: stages 4 and 5
- documentation: stages 2, 3, 4, 5, and 6

The current repo is strongest around one stage-4 example: the homodimer diagnostic notebook. It already has notebook-oriented skills for AFDB fetching, ipSAE, pDockQ, LIS, PAE visualization, MolViewSpec rendering, notebook assembly, and reference validation. The repo also has lifecycle-level concept discovery, decision capture, and PRD normalization skills.

The next skills should not be homodimer-specific by default. They should support the lifecycle stages themselves: turning ideas into triage notes, triage notes into PRDs/specs, specs into notebook prototypes, prototypes into beta experiences, and beta evidence into graduation decisions.

## Stage-Based Skill Map

| Lifecycle stage | Existing support | New reusable skills to create | Purpose |
|---|---|---|---|
| 1. Idea | none | `idea-capture`, `idea-quality-check` | Capture one-liners consistently and test whether the idea is simple, novel, user-relevant, and worth triage |
| 2. Triage | partial via docs | `idea-to-triage-brief`, `evidence-scan`, `risk-and-assumption-log` | Build a short pitch with enough evidence, assumptions, risks, and stop/go criteria |
| 3. Scope | `.agents/skills/idea-scoping-interview/SKILL.md`, `.agents/skills/scoping-decision-capture/SKILL.md`, `.agents/skills/concept-to-prd/SKILL.md`, `.agents/skills/prd-to-notebook-spec/SKILL.md`, `.agents/skills/notebook-spec-review/SKILL.md`, `.agents/skills/fixture-selection/SKILL.md` | `acceptance-criteria-authoring`, `dependency-policy` | Turn approved ideas into concrete notebook specs with inputs, outputs, fixtures, dependency limits, review gates, and success criteria |
| 4. Notebook prototype | homodimer-specific legacy suite in `agent-skills/agents/skills/`, canonical skills in `skills/`, including `skills/molviewspec-rendering/SKILL.md`, `skills/notebook-from-spec/SKILL.md`, `skills/notebook-execution-validation/SKILL.md`, `skills/notebook-review/SKILL.md` | `notebook-architecture`, `data-contract-validation`, `structure-parsing-strategy`, `scientific-computation-patterns`, `notebook-result-interpretation` | Build first shippable notebooks that are runnable, scientifically defensible, educational, and validated |
| 5. Instrumented beta | none | `beta-feedback-instrumentation`, `usage-signal-summary`, `beta-release-packaging` | Move validated notebooks into user-facing beta surfaces and measure preference, return usage, failures, and comments |
| 6. Graduation review | none | `graduation-decision-brief`, `engineering-handoff-package`, `standing-notebook-maintenance-plan`, `archive-retrospective` | Decide integrate, standing notebook, or archive, and produce the correct handoff artifact |
| Cross-cutting documentation | `.agents/skills/concept-to-prd/SKILL.md` | `documentation-rollup`, `scientific-provenance-log` | Accumulate documentation throughout the lifecycle instead of writing it at the end |
| Cross-cutting development harness | existing project guidance and homodimer skills | `fixture-regression-harness`, `dependency-runtime-audit`, `agent-skill-suite-design` | Keep LLM-assisted notebook work testable, maintainable, and portable across tools |

## Agent Map

Agents should be organized around lifecycle responsibilities first, and domain responsibilities second. A useful project-wide suite would be:

| Agent | Primary skills | Lifecycle role |
|---|---|---|
| `lifecycle-orchestrator` | all stage-transition skills | Route work through the correct stage and prevent jumping straight from idea to code |
| `triage-reviewer` | `idea-to-triage-brief`, `evidence-scan`, `risk-and-assumption-log` | Decide whether an idea deserves scope/spec work |
| `prd-and-spec-architect` | `concept-to-prd`, `prd-to-notebook-spec`, `acceptance-criteria-authoring` | Convert scoped concepts into PRDs and executable notebook specs |
| `spec-reviewer` | `notebook-spec-review` | Review the spec pack before implementation and block vague or untestable specs |
| `fixture-curator` | `fixture-selection`, `data-contract-validation` | Choose pinned examples and expected outputs |
| `notebook-builder` | `notebook-from-spec`, domain skills, `scientific-computation-patterns` | Implement the notebook from the reviewed spec and fixture manifest |
| `notebook-validator` | `notebook-execution-validation`, `fixture-regression-harness`, `dependency-runtime-audit` | Check the prototype before it is considered shippable |
| `notebook-reviewer` | `notebook-review`, `notebook-result-interpretation`, `scientific-figure-style` | Review the final notebook as a scientific and user-facing artifact |
| `beta-analyst` | `beta-feedback-instrumentation`, `usage-signal-summary` | Summarize beta evidence |
| `graduation-reviewer` | `graduation-decision-brief`, `engineering-handoff-package`, `standing-notebook-maintenance-plan`, `archive-retrospective` | Recommend integrate, standing notebook, or archive |

## Existing InsightFold Skill Gaps

| Gap | Why it matters | Suggested skill |
|---|---|---|
| Idea quality gate | Prevents weak ideas from becoming premature notebooks | `idea-to-triage-brief` |
| Scope discipline | Keeps notebooks from discovering their own shape while being built | existing `prd-to-notebook-spec`, `notebook-spec-review` |
| PRD normalization | `.agents/skills/concept-to-prd/SKILL.md` exists as the general PRD-generation skill | improve `concept-to-prd` rather than turning PRDs into implementation specs |
| Fixture selection | Notebook prototypes need pinned examples and expected outputs | `fixture-selection` |
| Notebook implementation | Turns the reviewed spec into a notebook without losing requirements, fixtures, validation criteria, or documentation intent | `notebook-from-spec` |
| Data contract validation | AFDB and PDBe APIs can change; notebooks need explicit schema checks | `data-contract-validation` |
| Structure parsing strategy | Choose Biopython, manual parsing, or another parser based on runtime budget, dependency policy, and notebook scope | `structure-parsing-strategy` |
| Interface geometry QC | Contact definitions are high-risk and central to pDockQ/pDockQ2 | `interface-geometry-qc` |
| Metric reference locking | Formula constants, thresholds, and edge cases must not drift | `confidence-metric-reference-pack` |
| Result interpretation | Converts computed results into clear scientific interpretation for the notebook's target audience | `notebook-result-interpretation` |
| Notebook execution checks | Confirms sequential execution, no hidden state, no GPU-only dependency creep | `notebook-execution-validation` |
| Notebook quality review | Separates "runs successfully" from "is worth sharing, beta testing, or graduating" | `notebook-review` |
| Scientific figure style | Keeps PAE, pLDDT, contact, and summary figures consistent | `scientific-figure-style` |
| Beta telemetry | Measures preference, repeat usage, and failure points without overinstrumenting | `beta-feedback-instrumentation` |
| Graduation package | Produces the handoff package for AFDB/PDBe engineering or standing-notebook maintainers | `graduation-handoff-package` |
| Archive retrospective | Preserves the reason an idea was retired and when to revisit it | `archive-retrospective` |

## Recommended Skill Build Order

1. existing `.agents/skills/concept-to-prd/SKILL.md`
2. implemented `skills/prd-to-notebook-spec/SKILL.md`
3. implemented `skills/notebook-spec-review/SKILL.md`
4. implemented `skills/fixture-selection/SKILL.md`
5. implemented `skills/notebook-from-spec/SKILL.md`
6. implemented `skills/notebook-execution-validation/SKILL.md`
7. implemented `skills/notebook-review/SKILL.md`
8. next: `structure-parsing-strategy`
9. next: `data-contract-validation`
10. next: `beta-feedback-instrumentation`
11. next: `graduation-handoff-package`

This order improves the highest-risk part of the current lifecycle first: turning scoped ideas into notebooks that are scientifically correct, runnable, and reviewable.

## Proposed Skill Definitions

### `idea-to-triage-brief`

Use for stage 1 to stage 2. Convert a one-line idea into a triage pitch with problem, user, novelty, AFDB/PDBe fit, why now, known evidence, expected notebook artifact, and stop criteria.

Reusable resources:

- `references/triage-template.md`
- `references/scoring-rubric.md`

### `prd-to-notebook-spec`

Use for stage 3. Convert an approved idea or PRD into a concise notebook specification: inputs, outputs, in scope, out of scope, fixture accession, execution path, acceptance criteria, dependency budget, and scientific risks.

Reusable resources:

- `assets/notebook-spec-template.md`
- `references/nbdd-stage-gates.md`

### `fixture-selection`

Use when selecting pinned AFDB/PDB examples for a prototype. Require one happy-path fixture, one negative/edge fixture, expected score snapshots, cached metadata shape, and reason each fixture was chosen.

Reusable resources:

- `scripts/collect_afdb_fixture.py`
- `references/fixture-manifest-schema.md`

### `data-contract-validation`

Use when fetching or validating AFDB/PDBe metadata, mmCIF, PAE, pLDDT, BinaryCIF, or UniProt cross-references. Enforce explicit field contracts and graceful degradation when optional fields are absent.

Reusable resources:

- `references/afdb-api-contracts.md`
- `scripts/validate_afdb_response.py`

### `structure-parsing-strategy`

Use when deciding how a notebook or prototype should parse structural files. Select between Biopython, manual mmCIF parsing, BinaryCIF tooling, or a project-local parser based on execution target, install budget, dependency policy, and scientific requirements.

For Colab-first notebooks with a strict dependency budget, use a lightweight parser that extracts atom-site rows, CA/CB coordinates, residue names, residue numbers, chain IDs, B-factors, and model filtering using standard-library parsing plus NumPy. For development, validation, richer structure handling, or less constrained runtimes, Biopython is a reasonable default and should not be rejected automatically.

Reusable resources:

- `scripts/parse_mmcif_atom_site.py`
- `references/mmcif-field-map.md`
- `references/parser-selection-rubric.md`

### `interface-geometry-qc`

Use when detecting interfaces, contacts, residue coverage, clashes, or comparing CA-CA versus CB-CB definitions. Preserve GLY to CA fallback and report contact-count provenance.

Reusable resources:

- `scripts/contact_matrix.py`
- `references/interface-definition-table.md`

### `confidence-metric-reference-pack`

Use when implementing or reviewing ipTM, ipSAE variants, pDockQ, pDockQ2, LIS, pLDDT summaries, thresholds, and reference tolerances. This should supersede scattered formula copies.

Reusable resources:

- `references/confidence-metrics.md`
- `scripts/check_metric_regression.py`
- pinned fixture expected values

### `notebook-result-interpretation`

Use when converting notebook outputs into plain-language scientific interpretation. Start from the notebook's stated use case, target audience, computed results, uncertainty, limitations, and any validation checks. Produce a concise interpretation that explains what the result supports, what it does not support, and what a user should inspect next.

For the homodimer diagnostic notebook, one reference file can include score-agreement and score-disagreement patterns. For other notebooks, add domain-specific interpretation patterns as separate references rather than making the skill itself homodimer-specific.

Reusable resources:

- `references/interpretation-patterns.md`
- `references/homodimer-score-patterns.md`
- `assets/diagnostic-summary-template.md`

### `notebook-execution-validation`

Use before a notebook is considered a shippable prototype. Execute top-to-bottom, check hidden-state hazards, dependency budget, Colab compatibility, output presence, fixture match, and optional-section isolation.

Reusable resources:

- `scripts/validate_notebook.py`
- `references/notebook-acceptance-checklist.md`

### `beta-feedback-instrumentation`

Use for stage 5. Add low-friction telemetry and feedback prompts for hosted notebook/beta pages: preference, return usage, failed input, runtime failures, and qualitative comments.

Reusable resources:

- `references/feedback-event-schema.md`
- `assets/feedback-copy.md`

### `graduation-handoff-package`

Use for stage 6. Build the integration, standing-notebook, or archive package: final PRD/spec, fixture manifest, known limitations, user feedback summary, runtime/dependency profile, engineering risks, and owner recommendation.

Reusable resources:

- `assets/handoff-template.md`
- `assets/archive-retrospective-template.md`

## Extracted `bioSkills` Candidates

The most relevant `bioSkills` should not be copied wholesale into the runtime notebook path. Many use Biopython, RDKit, BLAST, GATK, or external command-line tools. For InsightFold, treat them as source material for domain logic, review checklists, examples, optional development-only workflows, and runtime workflows where the dependency budget allows them.

### Directly Relevant to Protein Structures

| Source skill | Extract | Use in InsightFold |
|---|---|---|
| `/Users/mitsenkov/PycharmProjects/bioSkills/structural-biology/alphafold-predictions/SKILL.md` | AFDB URL patterns, metadata lookup, pLDDT from B-factors, PAE loading and plotting | Merge into `afdb-pdbe-data-contracts`; adapt to AFDB homodimer API shape already used in InsightFold |
| `/Users/mitsenkov/PycharmProjects/bioSkills/structural-biology/structure-io/SKILL.md` | PDB/mmCIF format awareness, BinaryCIF mention, metadata extraction concepts | Use as reference for `structure-parsing-strategy`; Biopython is appropriate when the notebook/runtime allows it |
| `/Users/mitsenkov/PycharmProjects/bioSkills/structural-biology/structure-navigation/SKILL.md` | Structure-model-chain-residue-atom hierarchy, residue IDs, atom properties, three-letter to one-letter mapping | Use as conceptual reference for manual parser output model |
| `/Users/mitsenkov/PycharmProjects/bioSkills/structural-biology/geometric-analysis/SKILL.md` | Distance matrices, residue contacts, superimposition, RMSD, neighbor search | Adapt NumPy-only distance/contact parts into `interface-geometry-qc`; keep RMSD/superposition as future optional skill |
| `/Users/mitsenkov/PycharmProjects/bioSkills/structural-biology/structure-modification/SKILL.md` | B-factor/property modification and coordinate transforms | Useful for generated visualization states or local review artifacts, not core notebook scoring |
| `/Users/mitsenkov/PycharmProjects/bioSkills/structural-biology/modern-structure-prediction/SKILL.md` | Comparing AlphaFold3, ESMFold, Chai-1, Boltz-1 predictions | Future skill for model comparison notebooks, outside current homodimer diagnostic core |

### Protein Annotation and Context

| Source skill | Extract | Use in InsightFold |
|---|---|---|
| `/Users/mitsenkov/PycharmProjects/bioSkills/database-access/uniprot-access/SKILL.md` | UniProt sequences, names, GO terms, domains, annotations | Add optional metadata enrichment skill for notebook summaries |
| `/Users/mitsenkov/PycharmProjects/bioSkills/database-access/interaction-databases/SKILL.md` | STRING/BioGRID/IntAct interaction confidence and network context | Add optional biological context panel for whether an interface has supporting interaction evidence |
| `/Users/mitsenkov/PycharmProjects/bioSkills/database-access/sequence-similarity/SKILL.md` | orthologs, distant homologs, protein family search | Future comparative confidence notebook: inspect whether interface residues are conserved |
| `/Users/mitsenkov/PycharmProjects/bioSkills/database-access/blast-searches/SKILL.md` and `/local-blast/SKILL.md` | sequence similarity workflows | Development-only; too heavy for lightweight notebooks unless explicitly scoped |
| `/Users/mitsenkov/PycharmProjects/bioSkills/sequence-manipulation/sequence-properties/SKILL.md` | protein molecular weight, isoelectric point, composition | Metadata enrichment; implement with Biopython or a small local function depending on runtime constraints |
| `/Users/mitsenkov/PycharmProjects/bioSkills/sequence-io/read-sequences/SKILL.md` and `/format-conversion/SKILL.md` | FASTA parsing and conversion patterns | Useful for future local upload/offline mode |

### Variants and Clinical/Population Annotation

These are useful because InsightFold already has a `src/insightfold/variants/` module and variant notebooks.

| Source skill | Extract | Use in InsightFold |
|---|---|---|
| `/Users/mitsenkov/PycharmProjects/bioSkills/variant-calling/vcf-basics/SKILL.md` | VCF fields, querying, field extraction | Future `variant-input-normalization` skill |
| `/Users/mitsenkov/PycharmProjects/bioSkills/variant-calling/variant-normalization/SKILL.md` | normalized variant representation, multiallelic decomposition | Improve variant mapping robustness before structure projection |
| `/Users/mitsenkov/PycharmProjects/bioSkills/variant-calling/variant-annotation/SKILL.md` | VEP/SnpEff/ANNOVAR consequence concepts | Optional annotation context, not notebook core |
| `/Users/mitsenkov/PycharmProjects/bioSkills/clinical-databases/clinvar-lookup/SKILL.md` | ClinVar pathogenicity and review status | Optional RUO annotation overlay for variant notebooks |
| `/Users/mitsenkov/PycharmProjects/bioSkills/clinical-databases/gnomad-frequencies/SKILL.md` | population frequencies | Optional rarity context for variant notebooks |
| `/Users/mitsenkov/PycharmProjects/bioSkills/clinical-databases/myvariant-queries/SKILL.md` | aggregate variant annotation API | Convenient optional lookup if network access is acceptable |
| `/Users/mitsenkov/PycharmProjects/bioSkills/clinical-databases/variant-prioritization/SKILL.md` | pathogenicity/frequency/evidence ranking | Future prioritization workflow; keep RUO guardrails explicit |

### Notebook Quality, Reporting, and Visualization

| Source skill | Extract | Use in InsightFold |
|---|---|---|
| `/Users/mitsenkov/PycharmProjects/bioSkills/reporting/jupyter-reports/SKILL.md` | parameterized notebooks, papermill-style reproducibility | Seed `notebook-execution-validation` and batch fixture checks |
| `/Users/mitsenkov/PycharmProjects/bioSkills/reporting/figure-export/SKILL.md` | figure sizing, resolution, export hygiene | Seed `scientific-figure-style` |
| `/Users/mitsenkov/PycharmProjects/bioSkills/reporting/automated-qc-reports/SKILL.md` | standardized QC report thinking | Adapt into notebook acceptance summaries, not MultiQC itself |
| `/Users/mitsenkov/PycharmProjects/bioSkills/data-visualization/color-palettes/SKILL.md` | colorblind-friendly palettes | Use in all PAE/contact/pLDDT figure skills |
| `/Users/mitsenkov/PycharmProjects/bioSkills/data-visualization/heatmaps-clustering/SKILL.md` | heatmap annotation conventions | Adapt for PAE heatmaps without clustering unless explicitly useful |
| `/Users/mitsenkov/PycharmProjects/bioSkills/data-visualization/interactive-visualization/SKILL.md` | hoverable exploratory plots | Future beta experience, not required for first Colab notebook |
| `/Users/mitsenkov/PycharmProjects/bioSkills/data-visualization/multipanel-figures/SKILL.md` | multi-panel layout rules | Improve score decomposition and summary figures |
| `/Users/mitsenkov/PycharmProjects/bioSkills/data-visualization/network-visualization/SKILL.md` | biological network visualization | Optional interaction-evidence panels |

### Workflow and ML Skills Worth Keeping Adjacent

| Source skill | Extract | Use in InsightFold |
|---|---|---|
| `/Users/mitsenkov/PycharmProjects/bioSkills/workflow-management/snakemake-workflows/SKILL.md` | make-like dependency thinking | Could inform fixture regression workflows; probably overkill for notebook runtime |
| `/Users/mitsenkov/PycharmProjects/bioSkills/workflow-management/nextflow-pipelines/SKILL.md` | containerized production pipelines | Relevant only if a prototype graduates to production-scale batch processing |
| `/Users/mitsenkov/PycharmProjects/bioSkills/machine-learning/model-validation/SKILL.md` | leakage prevention, validation discipline | Useful if future InsightFold notebooks train classifiers or thresholds from data |
| `/Users/mitsenkov/PycharmProjects/bioSkills/machine-learning/prediction-explanation/SKILL.md` | feature attribution language | Future explainability skill if score synthesis becomes model-assisted |
| `/Users/mitsenkov/PycharmProjects/bioSkills/experimental-design/multiple-testing/SKILL.md` | FDR/statistical testing discipline | Relevant for batch analyses or benchmarking, not single-accession diagnostics |

## Skills to Avoid Copying Directly

Avoid installing or copying these into the core InsightFold skill path unless a future PRD explicitly needs them:

- proteomics mass-spectrometry skills: valuable biologically, but unrelated to AFDB structure-confidence notebook lifecycle
- chemoinformatics and docking skills: potentially useful for ligand-centric PDBe work, not current homodimer scoring
- heavy workflow engines: Nextflow, WDL, CWL are production/batch concerns, not stage 4 Colab prototype concerns
- clinical interpretation skills: only use in RUO variant notebooks with explicit disclaimers and no clinical decision framing

## Agent Implications

The current agent suite is implementation-heavy. Add lifecycle agents only after the corresponding skills exist:

| Agent | Skills it would use | Lifecycle role |
|---|---|---|
| idea-reviewer | `insightfold-idea-triage` | stage 2 gate |
| scope-architect | `insightfold-scope-spec`, `afdb-pdbe-data-contracts` | stage 3 |
| fixture-curator | `insightfold-fixture-curation`, `afdb-pdbe-data-contracts` | stages 3 and 4 |
| notebook-builder | existing AFDB/scoring/visualization/notebook skills | stage 4 |
| scientific-validator | `confidence-metric-reference-pack`, `interface-geometry-qc`, `notebook-execution-validation` | stage 4 review |
| beta-analyst | `beta-feedback-instrumentation`, `notebook-result-interpretation` | stage 5 |
| graduation-reviewer | `graduation-handoff-package`, `archive-retrospective` | stage 6 |

## Immediate Next Step

Start by turning `confidence-metric-reference-pack`, `structure-parsing-strategy`, and `notebook-execution-validation` into first-class skills. They directly raise quality for the existing homodimer notebook and will also generalize to future InsightFold notebook prototypes.
