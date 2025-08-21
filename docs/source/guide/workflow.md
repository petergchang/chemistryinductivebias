# Workflow

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
