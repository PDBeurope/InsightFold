---
name: molviewspec-rendering
description: This skill produces 3D structural views using MolViewSpec, including chain coloring, pLDDT overlays, and interface highlights for interactive educational displays.
---

# Goal

Generate MolViewSpec visualizations for the notebook.

# Required Views

1. Chain overview
2. pLDDT colouring
3. Interface highlighting
4. Diagnostic disagreement mapping

# Representation Rules

Prefer:
- cartoon representations
- uncluttered views
- clear chain colours

# Colour Rules

Chain A:
#009688

Chain B:
#FF7043

Use standard AlphaFold pLDDT colours.

# Technical Constraints

- prefer BinaryCIF
- generate .mvsj states
- support notebook embedding
- support Mol* URLs

# Dependencies

- molviewspec library
