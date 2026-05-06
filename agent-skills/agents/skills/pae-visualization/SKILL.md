---
name: pae-visualization
description: This skill creates heatmaps and overlays for PAE matrices, highlighting regions for different metrics with annotations and color conventions to teach score interpretations.
---

# Goal

Generate educational PAE visualizations.

# Required Plots

- full PAE heatmap
- inter-chain quadrants
- masked score regions
- contact overlays

# Visualization Rules

Use:
- seaborn heatmaps
- matplotlib
- visible chain boundaries
- colour bars
- annotations

# Dependencies

- seaborn
- matplotlib

# Required Masks

Visualize:
- ipTM regions
- ipSAE regions
- LIS regions
- pDockQ2 regions

# Educational Focus

Explain:
- what each mask means
- why score regions differ
- how cutoff logic changes interpretation
