We have now completed our scoping conversation for the EnzyMM + AlphaFold 
notebook. I need you to produce a structured summary that captures 
everything we covered. This summary will be used directly to write the 
development spec, so it must be comprehensive and precise.

Please structure your summary as follows:

---

### 1. The biological question
What is the core scientific question this notebook answers? State it 
clearly in one or two sentences. Then describe why this is meaningful in 
the context of predicted protein structures at AFDB scale.

### 2. How EnzyMM works (as understood for this notebook)
Summarise the key technical aspects of EnzyMM relevant to this notebook: 
the geometric matching approach, the M-CSA template library, the key 
output metrics (RMSD, log E-value, orientation, completeness, 
predicted_correct), and the known limitations we discussed.

### 3. AlphaFold-specific caveats
List the specific considerations introduced by running EnzyMM on an 
AlphaFold model. Include how pLDDT relates to catalytic site confidence, 
side chain rotamer reliability, and when to flag results as unreliable.

### 4. Notebook scope decisions
For each of the following, state what we decided and the reasoning:
- Input method (UniProt accession / AFDB API vs user-uploaded structure)
- pLDDT filtering (yes/no, threshold, how explained to the user)
- Multimeric models (in scope for v1 or deferred)
- Which output fields to surface and interpret
- Hit prioritisation logic
- Visualisations to include

### 5. Narrative flow
Describe the intended user journey through the notebook, section by 
section. What does the user do, what do they see, what do they learn at 
each step?

### 6. Open questions and deferred decisions
List anything we explicitly did not resolve, with enough context that 
another session can pick it up cleanly.

### 7. Signal for graduation
What usage signal or biological validation outcome would justify 
graduating this notebook to a pre-computed AFDB pipeline?

---

Be precise. Where we made a specific decision, state it as a decision 
("We decided X because Y"), not as an open option. Where something was 
left open, say so explicitly.