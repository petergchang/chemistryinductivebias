# Modules

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
