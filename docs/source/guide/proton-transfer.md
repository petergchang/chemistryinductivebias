# Proton Transfer Mechanism

## Core Workflow

A three-step pipeline from candidates to a final score.

## 1) Apply Process Rule (enumerate candidates)

- Find acceptors A (features + optional exhaustive mode incl. non-H atoms and hydride).
- Find acidic hydrogens B–H (heavy atoms bonded to non-bridging hydrogens).
- Apply transfer: break B–H, form A–H, increment charge on A, decrement on B.
- Sanitize and discard chemically invalid products (via RDKit sanitization).

## 2) Compute Properties (per candidate "Transformation")

- Change in formal charge, electronegativity deltas, atomic radii, resonance flags.
- Inductive stabilization score using distance-attenuated relative EN and boosted charged centers.

## 3) Score and Aggregate

- Apply object rules: formal charge, electronegativity, resonance, atomic radius, inductive.
- Convert to labels/scores via thresholds and `SCORE_MAP` and take a weighted sum.

```python
# pseudo
cands = generate_candidates(mol)
props = [compute_transformation_properties(c) for c in cands]
comp = [score_transformation_components(p, params) for p in props]
final = [aggregate_components(cs, params.weights) for cs in comp]
```


## Modules

Left column is English; right is pseudocode-like steps you can map to code.

## Process pipeline

- English: Enumerate candidate transfers, compute properties, score via object rules, aggregate.
- Pseudocode:

```python
# pseudo
candidates = generate_candidates(mol)
records = [compute_properties(c) for c in candidates]
scores = [score_components(r, params) for r in records]
return aggregate(scores)
```

## Object rules

- English: Each rule maps properties → label → numeric score.
- Pseudocode:

```python
# pseudo
components = {
    'formal_charge': formal_charge_object_rule(record),
    'electronegativity': electronegativity_object_rule(record, thresholds),
    'resonance': resonance_object_rule(record),
    'atomic_radius': atomic_radius_object_rule(record, thresholds),
    'inductive': inductive_object_rule(record, hyperparams),
}
final_score = weighted_sum(components, weights)
```

## Parameters

- English: Tunables that control thresholds and weighting.
- Pseudocode:

```python
# pseudo
params = {
    'weights': (0.45, 0.25, 0.15, 0.10, 0.05),
    'SCORE_MAP': { '++': 2, '+': 1, '0': 0, '-': -1, '--': -2 },
    'en_thresholds': (0.05, 0.5),
    'radius_thresholds': (0.05, 0.15),
    'inductive': {
        'max_bonds': 4,
        'power': 2.0,
        'charge_boost': 0.5,
        'en_margin': 0.10,
    },
}
```



## Key Components

A scannable list; each item has a one‑sentence description.

- **Candidate generation**: enumerate acceptor A and donor B–H pairs, simulate A–H and B changes, sanitize with RDKit.
- **Property computation**: compute charge deltas, electronegativity differences, atomic radii, resonance flags, and an inductive stabilization score.
- **Object rules**: convert properties to qualitative labels and numeric scores via thresholds and mappings.
- **Aggregation**: weighted sum of component scores → final favorability.
- **Parameters**: weights, score map, EN/radius thresholds, inductive hyperparameters.
- **Diagnostics**: structured records for ablations and error analysis.



## Layers

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




## Parameters

Defaults listed here are the current working values.

## Scoring weights (weights)

- Five floats for: Formal Charge, Electronegativity, Resonance, Atomic Radius, Inductive.
- Default: `(0.45, 0.25, 0.15, 0.10, 0.05)`.

## Score-to-label mapping (SCORE_MAP)

- Numeric value assigned to qualitative labels: `{ '++': 2, '+': 1, '0': 0, '-': -1, '--': -2 }`.

## Thresholds

- Electronegativity ΔEN thresholds: neutrality boundary `0.05`, favorability boundary `0.5`.
- Atomic radius Δr thresholds: neutrality boundary `0.05 Å`, favorability boundary `0.15 Å`.

## Inductive effect

- Final-score label thresholds: `> 0.5` → FAVORABLE, `> 2.0` → VERY_FAVORABLE.
- Hyperparameters:
  - `max_bonds`: 4 (maximum bond distance considered)
  - `power`: 2.0 (distance attenuation exponent, 1/d^power)
  - `charge_boost`: 0.5 (extra contribution from formal positive centers)
  - `en_margin`: 0.10 (ignore negligible ΔEN)

```python
# pseudo snapshot of tunables
params = {
    'weights': (0.45, 0.25, 0.15, 0.10, 0.05),
    'SCORE_MAP': { '++': 2, '+': 1, '0': 0, '-': -1, '--': -2 },
    'en_thresholds': (0.05, 0.5),           # neutral, favorable
    'radius_thresholds': (0.05, 0.15),      # neutral, favorable (Å)
    'inductive': {
        'label_thresholds': (0.5, 2.0),     # FAVORABLE, VERY_FAVORABLE
        'max_bonds': 4,
        'power': 2.0,
        'charge_boost': 0.5,
        'en_margin': 0.10,
    },
}
```


## Significant Changes: V1 → V2

## Formal charge rule

- V1: Checked only if final atoms were neutral; ++ if both neutral, + if one neutral, − otherwise.
- V2: Prioritizes neutralization, heavily penalizes charge pileup (e.g., +1 → +2), distinguishes unfavorable intramolecular separation from plausible intermolecular charge creation.

## Atomic radius rule

- V1: Incorrectly favored larger proton acceptor (not donor).
- V2: Correctly favors larger proton donor.

## Resonance rule

- V1: A few SMARTS patterns to detect resonance.
- V2: Algorithmic check for adjacency to a π-system; excludes cases where the anion is directly part of the π system.

## Inductive rule

- V1: Sum of EN/distance² using a fixed list of withdrawing atoms; ignored local environment and formal charges.
- V2: Uses relative EN differences between neighbors, distance attenuation, and boosted contributions from positively charged centers.

## EN and radius discretization

- V1: Coarse thresholds.
- V2: Five-level scoring system for finer granularity.
