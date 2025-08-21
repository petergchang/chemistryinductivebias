# Key Components

A scannable list; each item has a one‑sentence description.

- **Candidate generation**: enumerate acceptor A and donor B–H pairs, simulate A–H and B changes, sanitize with RDKit.
- **Property computation**: compute charge deltas, electronegativity differences, atomic radii, resonance flags, and an inductive stabilization score.
- **Object rules**: convert properties to qualitative labels and numeric scores via thresholds and mappings.
- **Aggregation**: weighted sum of component scores → final favorability.
- **Parameters**: weights, score map, EN/radius thresholds, inductive hyperparameters.
- **Diagnostics**: structured records for ablations and error analysis.
