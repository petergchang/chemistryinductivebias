# Parameters

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
