## mapping Ca distances Specification

## Overview

## Purpose
The Mapping Cα Distances system is an ensemble-aware structural perturbation analysis notebook designed to compare experimentally resolved wild-type (WT) and pathogenic mutant protein structures using Cα distance mapping and contact perturbation analysis.

The system addresses the scientific question:

Can experimentally observed structural perturbations between WT and mutant protein structures be transformed into interpretable mechanistic hypotheses about disease-associated molecular effects?

The notebook focuses on:

 - mechanistic structural interpretation
 - biologically contextualized perturbation analysis
 - explainable evidence aggregation
 - ensemble-aware contact analysis for NMR structures

The system explicitly does not attempt:
 - pathogenicity prediction
 - ML-based classification
 - clinical diagnosis
 - energetic simulation

Intended Users
Primary users:
 - PhDs
 - postdoctoral researchers
 - structural bioinformaticians
 - structural biologists
 - computational biologists

Users are expected to understand:

 - protein structure
 - molecular interactions
 - mutation biology
 - structural terminology

Users are not assumed to be:

expert programmers

experts in structural comparison tooling

Primary Value Delivered
The notebook enables users to:

compare WT and mutant structures systematically

identify structural perturbations

quantify contact rewiring

inspect ensemble-aware contact persistence

rank biologically relevant perturbations

generate mechanistic structural hypotheses

The notebook functions as:

an evidence-driven structural reasoning assistant

rather than:

a predictive pathogenicity classifier.

Ecosystem Context
The system fits into:

structural bioinformatics workflows

pathogenic mutation interpretation pipelines

mechanistic hypothesis generation workflows

experimental structure analysis workflows

The strongest differentiator is:

ensemble-aware mechanistic structural interpretation using experimentally resolved structures.

Unlike many existing systems, this notebook prioritizes:

transparency

explainability

uncertainty communication

interpretable evidence aggregation

mechanistic reasoning over prediction accuracy

Domain / Business / Scientific Context
Scientific Background
Protein mutations may alter:

folding stability

packing interactions

conformational flexibility

interface interactions

allosteric communication

secondary structure stability

contact persistence

Many mutation effects are difficult to interpret from sequence alone.

Structural comparison enables:

direct observation of geometric perturbations

interaction-network analysis

mechanistic hypothesis generation

Core Structural Concept
The notebook uses:

Cα distance matrices

contact maps

ensemble contact persistence

to characterize structural differences between:

WT structures

mutant structures

Contact Maps
A contact map is a matrix representation describing which residues are spatially proximal.

Primary representation:

continuous Cα distance matrices

Secondary representation:

binary threshold contact maps

Continuous distance matrices are prioritized because:

binary contacts introduce threshold artifacts

small geometric differences can produce misleading binary changes

continuous distances preserve perturbation magnitude

NMR Ensemble Context
NMR structures are treated as:

ensembles
not:

single static structures

The notebook calculates:

contact persistence frequencies

ensemble variability metrics

perturbation consistency across ensemble members

Key metric:

P(contact)= 
total models
models with contact
​
 
This allows:

dynamic sensitivity

flexibility-aware interpretation

confidence-aware perturbation scoring

Mechanistic Interpretation Philosophy
The system performs:

mechanistic interpretation
not:

pathogenicity prediction

Outputs are:

evidence-backed hypotheses

not deterministic conclusions

Interpretations must use:

probabilistic language

uncertainty-aware phrasing

Examples:

“consistent with”

“may indicate”

“suggests”

Avoid:

“proves”

“demonstrates”

“causes”

Critical Scientific Caveats
Structural differences may not be mutation-driven
Observed differences may result from:

crystallization conditions

ligand states

oligomeric differences

experimental method differences

construct truncation

pH differences

rather than mutation effects.

The notebook must explicitly communicate:

structural comparability limitations

experimental-context mismatch warnings

Contact maps are abstractions
Cα distances capture:

topology
not:

detailed chemistry

The notebook does not directly model:

hydrogen bonds

electrostatics

side-chain energetics

free energy

Proteins are dynamic
Single structures do not fully represent:

biological conformational landscapes

Absence of perturbation:

does not imply absence of biological effect

Input Specification
Accepted Inputs
Required Inputs
Field	Type	Required	Description
wt_pdb_id	string	Yes	Experimental WT PDB identifier
mutant_pdb_id	string	Yes	Experimental mutant PDB identifier
Optional Inputs
Field	Type	Required	Description
wt_chain_id	string	No	WT chain identifier
mutant_chain_id	string	No	Mutant chain identifier
contact_threshold	float	No	Binary contact threshold in Å
ensemble_contact_cutoff	float	No	Minimum ensemble persistence threshold
alignment_method	string	No	Structural alignment method
distance_metric	string	No	Distance metric type
residue_subset	list[int]	No	Optional residue subset
output_directory	string	No	Output path
Supported Structure Types
Supported:

X-ray crystallography

NMR ensembles

cryo-EM structures with resolved coordinates

Unsupported:

AlphaFold predictions (v1)

MD trajectories (v1)

hybrid models

External APIs / Services
RCSB PDB API
Used for:

structure retrieval

metadata retrieval

Endpoints:

structure download

metadata lookup

Input Validation Rules
Hard Rejection Conditions
Reject if:

PDB ID invalid

structure unavailable

no Cα coordinates

chain absent

structure unreadable

zero sequence overlap

Warning Conditions
Warn if:

experimental methods differ

resolution differs substantially

oligomeric states differ

ligand states differ

missing residues exceed threshold

sequence overlap < 80%

ensemble size mismatch

chain mapping ambiguous

Warnings do not block execution.

Preprocessing Requirements
Structure Parsing
Required:

atom extraction

residue indexing

chain extraction

ensemble extraction

coordinate normalization

Residue Mapping
Required:

residue-number reconciliation

sequence alignment

insertion-code handling

Structural Alignment
Structures must be aligned before:

distance comparison

perturbation calculation

Default:

backbone-based superposition

Error Handling
Retry Logic
Retry:

network failures

transient API errors

Retry count:

3 attempts

Backoff:

exponential

Fatal Errors
Fatal errors terminate workflow with:

descriptive error message

suggested corrective action

System Execution / Processing Logic
Execution Flow
Input →
Validation →
Structure Parsing →
Alignment →
Distance Matrix Generation →
Contact Persistence Calculation →
Perturbation Scoring →
Biological Annotation →
Mechanistic Interpretation →
Visualization →
Export
Installation Requirements
Python Version
Python >= 3.10
Core Dependencies
biopython
numpy
scipy
pandas
matplotlib
seaborn
plotly
py3Dmol
scikit-learn
networkx
requests
Core Algorithms
Distance Matrix Calculation
For each residue pair:

d ij =∣∣Cαi−Cαj∣∣

Computed:

per ensemble member

per structure

Ensemble Distance Aggregation
Default aggregation:

mean distance

standard deviation

persistence frequency

Binary Contact Definition
Default threshold:

8.0 Å

Contact:

Cij ={ 1 dij <threshold
       0 otherwise  
     }
​

​
 
Threshold configurable.

Contact Persistence
For NMR ensembles:

Pij = models with contact/total models

​

 
Distance Difference Matrix
Primary perturbation metric:

ΔD = D
m
u
t
a
n
t
−
D
W
T
ΔD=D 
mutant
​
 −D 
WT
​
 
Per-Residue Perturbation Score
Initial v1 formula:

S
i
=
∑
j
∣
Δ
d
i
j
∣
S 
i
​
 = 
j
∑
​
 ∣Δd 
ij
​
 ∣
Used for:

ranking

perturbation prioritization

Intermediate Artifacts
Generated artifacts:

Artifact	Format
aligned_structures	PDB
wt_distance_matrix	NumPy array
mutant_distance_matrix	NumPy array
difference_matrix	NumPy array
contact_frequency_matrix	NumPy array
perturbation_scores	CSV
annotations	JSON
Output Artifacts
Generated outputs:

heatmaps

perturbation rankings

mechanistic summaries

validation reports

ensemble statistics

annotated residue tables

Performance Expectations
Expected scale:

proteins ≤ 1500 residues

ensembles ≤ 50 models

Expected runtime:

< 2 minutes for typical structures

Output Parsing and Interpretation
Output Schema
Structural Validation Report
{
  "sequence_overlap": 0.92,
  "missing_residue_fraction": 0.08,
  "experimental_method_match": true,
  "warnings": []
}
Perturbation Matrix
2D matrix:

Residue_i	Residue_j	DeltaDistance
Residue Ranking Table
Residue	PerturbationScore	PersistenceShift	Interpretation
Ranking Logic
Primary ranking criteria:

perturbation magnitude

ensemble consistency

contact persistence disruption

biological annotation overlap

buried residue involvement

Confidence Interpretation
Confidence categories:

High

Moderate

Low

Based on:

ensemble consistency

annotation support

perturbation robustness

structural comparability

Empty Result Behaviour
If no major perturbations detected:

Display:

“No substantial structural perturbations were detected under the current analysis parameters. This does not exclude biologically relevant effects outside the structural features captured by Cα distance analysis.”

Ambiguity Handling
Ambiguous results must:

expose uncertainty

avoid deterministic language

surface competing interpretations

Visualisations / UI Components
1. WT Distance Heatmap
Purpose
Visualize WT structural topology.

Rendering
Plotly heatmap

Interaction
hover residue lookup

zoom

export

2. Mutant Distance Heatmap
Purpose:

visualize mutant topology

Same interaction model as WT.

3. Difference Heatmap
Purpose
Primary perturbation visualization.

Data
Δ
D
ΔD
Color semantics
red: increased distance

blue: decreased distance

4. Ensemble Persistence Heatmap
Purpose
Visualize contact stability.

Interpretation Guidance
Display tooltip:

“High persistence indicates contacts consistently observed across ensemble members.”

5. Residue Perturbation Ranking Table
Purpose
Rank structurally perturbed residues.

Features
sortable

filterable

annotation overlays

6. 3D Structural Overlay
Purpose
Map perturbations onto structure.

Rendering
py3Dmol

Highlights
highly perturbed residues

persistent contact disruptions

Workflow / User Journey
Step 1 — User Input
User enters:

WT PDB ID

mutant PDB ID

System retrieves structures.

UI displays:

metadata

experimental method

chain information

Step 2 — Validation
System performs:

sequence comparison

chain mapping

missing residue analysis

comparability assessment

UI displays:

warnings

reliability concerns

comparability summary

Step 3 — Structural Computation
System:

aligns structures

computes matrices

calculates perturbations

UI displays:

loading state

progress indicators

Step 4 — Visualization
System renders:

heatmaps

perturbation maps

ensemble statistics

User explores:

residue-level perturbations

contact changes

Step 5 — Interpretation
System generates:

ranked perturbations

mechanistic hypotheses

Examples:

possible core destabilization

interface perturbation

flexibility increase

Step 6 — Evidence Inspection
User inspects:

supporting contacts

ensemble persistence

structural neighborhoods

System exposes:

raw evidence

uncertainty indicators

Component / Module Breakdown
Structure Retrieval Module
Responsibility
Retrieve PDB structures and metadata.

Dependencies
requests

RCSB API

Structure Parser Module
Responsibility
Extract:

coordinates

residues

ensembles

chains

Dependencies
BioPython

Alignment Module
Responsibility
Align WT and mutant structures.

Distance Matrix Module
Responsibility
Compute:

Cα matrices

difference matrices

persistence matrices

Perturbation Scoring Module
Responsibility
Compute:

residue scores

contact perturbation scores

Annotation Module
Responsibility
Attach:

conservation

secondary structure

accessibility

interface annotations

Mechanistic Interpretation Module
Responsibility
Generate:

evidence-backed hypotheses

confidence summaries

Visualization Module
Responsibility
Generate:

heatmaps

overlays

interactive visualizations

Reporting Module
Responsibility
Export:

CSV

JSON

images

summary reports

User-Facing Explanations and Copy
Onboarding Explanation
“This notebook compares experimentally resolved wild-type and mutant protein structures using Cα distance mapping and ensemble-aware contact analysis. The goal is to identify structural perturbations that may help explain mutation-associated molecular mechanisms.”

Structural Caveat Warning
“Structural differences may arise from experimental conditions, ligand states, or construct differences rather than mutation effects alone. Interpret perturbations within the broader experimental context.”

Contact Interpretation Guidance
“Cα distance analysis captures topological structural changes rather than detailed physicochemical interactions.”

Confidence Guidance
“Confidence estimates reflect structural consistency and ensemble support, not clinical certainty.”

Low Confidence Warning
“Interpretation confidence is reduced due to low ensemble consistency or limited structural comparability.”

Empty State Message
“No strong perturbation signatures were detected. Subtle functional effects may still exist outside the sensitivity of Cα distance analysis.”

Data, Logging, and Observability
Stored Data
Store:

input metadata

validation outputs

matrices

perturbation scores

interpretation outputs

Logging Requirements
Log:

API requests

validation failures

parsing errors

alignment statistics

runtime metrics

Metrics
Track:

runtime

matrix generation time

ensemble size

perturbation distribution

warning frequency

Debugging Requirements
Support:

verbose logging mode

matrix export

alignment inspection

Privacy Considerations
No sensitive clinical data stored.

Integration Points
Upstream Dependencies
RCSB PDB

BioPython

Export Formats
Supported:

CSV

JSON

PNG

HTML

Reusable Outputs
Outputs reusable for:

downstream structural analysis

variant prioritization workflows

reporting pipelines

Non-Functional Requirements
Performance
Target:

<2 minute runtime for standard structures

Reliability
System must:

fail gracefully

preserve intermediate artifacts

expose warnings clearly

Maintainability
Codebase must:

separate computation from interpretation

support modular scoring systems

Portability
Notebook must run:

locally

in Google Colab

Usability
Notebook must:

require minimal coding

expose parameters clearly

provide inline interpretation guidance

Out of Scope for v1
Feature	Rationale
Pathogenicity prediction	Requires major validation effort
ML classifiers	Reduces interpretability
Molecular dynamics	Computationally expensive
Free-energy calculations	Outside intended scope
AlphaFold integration	Experimental structures prioritized
Clinical interpretation	Insufficient validation
Side-chain energetic modeling	Too complex for v1
Open Questions Requiring Resolution Before Build
1. Perturbation Scoring Framework
Question
How should perturbation scores be weighted and aggregated?

Why It Matters
Scoring determines ranking quality.

Next Step
Benchmark against known pathogenic variants.

2. Contact Threshold Selection
Question
What binary threshold values should be default?

Why It Matters
Threshold sensitivity affects interpretation.

Next Step
Perform robustness analysis.

3. Confidence Scoring Model
Question
How should evidence confidence be quantified?

Why It Matters
Critical for trust calibration.

Next Step
Expert review and user testing.

4. Annotation Source Selection
Question
Which biological annotation systems should be integrated first?

Why It Matters
Affects interpretation quality and implementation complexity.

Next Step
Prioritize based on user needs and implementation effort.

5. Visualization Prioritization
Question
Which visualization modality should become primary?

Why It Matters
Large proteins create visualization scalability challenges.

Next Step
Conduct UX evaluation with expert users.