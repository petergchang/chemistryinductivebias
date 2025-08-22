# Layers

A three-layer view that separates decision intent from implementation details.

## Higher-level (Intent)

- **Objective**: Predict thermodynamic favorability of proton transfer candidates.
- **Assumptions**: aprotic solvent, RDKit sanitization for chemical validity, ignore kinetics.
- **Scoring principle**: weighted sum of five object rules.

## Mid-level (Architecture)

- **Pipeline**:
  1. Generate candidates (acceptors, acidic hydrogens, simulate transfer, sanitize).
  2. Compute properties for each candidate (charge deltas, EN deltas, radii, resonance flags, inductive score).
  3. Score with object rules (formal charge, EN, resonance, radius, inductive) and aggregate.
- **Data model**: `Transformation` → `Components` → `MechanismProducts`.

## Lower-level (Mechanics)

- **Algorithms**: resonance detection, distance-attenuated inductive scoring, threshold discretization.
- **Parameters**: `weights`, `SCORE_MAP`, EN/Radius thresholds, inductive hyperparameters (`max_bonds`, `power`, `charge_boost`, `en_margin`).
- **Diagnostics**: per-candidate component scores and labels for debugging and ablations.


