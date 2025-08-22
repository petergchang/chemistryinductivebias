# Proton Transfer Mechanism
The proton transfer mechanism is implemented as {py:func}`proton_transfer_predict(.) <chemistryinductivebias.src.proton_transfer.proton_transfer_predict>`.

## Core Workflow: An Overview

We use a four-step pipeline from reactants to the final prediction.

### {ref}`1) Enumerate Proton Acceptors and Donors <step-1>`
- Find proton acceptors `A` (all non-H atoms and hydride).
- Find acidic hydrogens `B-H` (heavy atoms bonded to non-bridging hydrogens).

### {ref}`2) Apply Proton Transfer <step-2>`
- For each `(A, B, H)`, break `B-H` bond and form `A-H` bond.
- Increment charge on `A`, decrement charge on `B`.
- Sanitize and discard chemically invalid products.

### {ref}`3) Compute Transformation Properties <step-3>`
- Compute delta in formal charge, electronegativity, and atomic radii.
- Resonance flag computed by identifying {math}`\pi`-bonds adjacent to the new anion.
- Inductive stabilization score using distance-attenuated relative EN and boosted charged centers.

### {ref}`4) Score and Aggregate <step-4>`

- Apply object rules for formal charge, electronegativity, resonance, atomic radius, and inductive effects to generate favorability labels.
- Convert favorability labels to scores and take a weighted sum.

```python
## Pseudocode
# Step 1: Find potential proton acceptors and donors
acceptors = find_acceptors(reactants)
donors = find_acidic_hydrogens(reactants)

candidates = []
for A in acceptors:
    for B, H in donors:
        # Step 2: Apply proton transfer
        candidate = apply_proton_transfer(A, H, B)

        # Step 3: Compute transformation properties
        properties = compute_transformation_properties(candidate)
        candidates.append((candidate, properties))

# Step 4: Apply object rule to assess favorability
scores = score_transformation_components(candidates)
```

(step-1)=
## Step 1: Enumerate Proton Acceptors and Donors
Note: Steps 1 and 2 combined are implemented as {py:func}`proton_transfer_process_rule(.) <chemistryinductivebias.src.proton_transfer.proton_transfer_process_rule>`

Finding proton acceptors is implemented as {py:func}`find_acceptors(.) <chemistryinductivebias.src.proton_transfer.find_acceptors>`. It takes an `rdchem.Mol` instance and a `MechanismOptions` instance with the default `enumerate_all_acceptors=True` as arguments. Given this, it returns all non-hydrogen atoms in the molecule. Otherwise, it uses `RDKit`'s feature factory with family name `"Acceptor"` to return proton acceptors.

```py
## Pseudocode
def find_acceptors(mol):
    acceptors = []
    for a in mol.GetAtoms():
        if a != hydrogen or a == hydride:
            acceptors.append(a)
    return acceptors
```

Finding proton donors is implemented as {py:func}`find_acidic_hydrogens(.) <chemistryinductivebias.src.proton_transfer.find_acidic_hydrogens>`. It finds all non-bridging hydrogens and returns the donor atom and the hydrogen connected to it.

```py
## Pseudocode
def find_acidic_hydrogens(mol):
    donor_pairs = []
    for a in mol.GetAtoms():
        neighbors = a.neighbors
        if a == hydrogen and len(neighbors) == 1:
            donor_pairs.append((neighbors[0], a))
    return donor_pairs
```

(step-2)=
## Step 2: Apply Proton Transfer
Applying the proton transfer transformation for a given potential acceptor-donor pair is implemented as {py:func}`apply_proton_transfer(.) <chemistryinductivebias.src.proton_transfer.apply_proton_transfer>`. Given the acceptor `A`, and donor-hydrogen pair `B-H`, it first breaks the `B-H` bond, forms the new `A-H` bond, increments `A`'s formal charge, decreases `B`'s formal charge, and finally discard invalid structures by sanitizing using `rdkit.Chem.SanitizeMol(.)`.

```py
## Pseudocode
def apply_proton_transfer(mol, A, H, B):
    rwm = rdchem.RWMol(mol)  # Make an editable copy
    rwm.RemoveBond(B, H)  # Break B-H bond
    rwm.AddBond(A, H, BondType.SINGLE)  #Form A-H bond
    
    # Update formal charges
    A.SetFormalCharge(A.GetFormalCharge() + 1)
    B.SetFormalCharge(B.GetFormalCharge() - 1) 

    newmol = rwm.GetMol()
    try:
        Chem.SanitizeMol(newmol)
    except:
        return None  # Discard invalid structures.
    return newmol
```
### Examples of invalid structures that will be discarded

Mainly reactions that create hypervalent main-group atoms will be discarded, with RDKit raising `ValueError: Sanitization error: Explicit valence for atom # C 5 is greater than permitted`

- Attempting to protonate methane (`CH4`), the invalid product is `[CH5+]`.
- Attempting to protonate a hydronium ion (`[OH3+])`), the invalid product is `[H4O]2+`.
- Attempting to protonate hydrogen fluoride (`HF`), the invalid product is `[H2F+]`.

(step-3)=
## Step 3: Compute Transformation Properties

Computing the chemical properties of each valid transformation candidate is implemented as {py:func}`compute_transformation_properties(.) <chemistryinductivebias.src.proton_transfer.compute_transformation_properties>`. It gathers properties of the acceptor `A` and donor `B` atoms before and after the transfer, for the object rules to ultimately assess the transformation's favorability. This includes changes in formal charge, electronegativity, and atomic radius.

The two key stabilizing effects, resonance and induction, are also evaluated using helper functions.

```py
## Pseudocode
def compute_transformation_properties(old_mol, new_mol, A, B):
    a_old = old.GetAtomWithIdx(A)
    b_old = old.GetAtomWithIdx(B)
    a_new = new.GetAtomWithIdx(A)
    b_new = new.GetAtomWithIdx(B)
    a_old_charge = a_old.GetFormalCharge()
    b_old_charge = b_old.GetFormalCharge()
    a_new_charge = a_new.GetFormalCharge()
    b_new_charge = b_new.GetFormalCharge()

    return {
        "delta_charge_on_A": float(a_new_charge - a_old_charge),
        "delta_charge_on_B": float(b_new_charge - b_old_charge),
        # ... other properties ...
        "EN_A": get_en(a_old),  # lookup table
        "EN_B": get_en(b_old),
        "Radius_A": get_radius_angstrom(a_old),  # lookup table
        "Radius_B": get_radius_angstrom(b_old),
        "is_resonance_stabilized": is_resonance_stabilized(new, B),
        "inductive_score": calculate_inductive_score(new, B),
    }
```

### Resonance Stabilization
Evaluating resonance stabilization is implemented as
{py:func}`is_resonance_stabilized(.) <chemistryinductivebias.src.proton_transfer.is_resonance_stabilized>`. It checks whether the newly-formed anion (the deprotonated atom `B`) is adjacent to a {math}`\pi`-system, which allows the negative charge to be delocalized.

```py
## Pseudocode
def is_resonance_stabilized(mol, anion):
    # Rule 1: The anion itself must NOT be a part of pi system.
    # Its lone pair must be available for delocalization.
    for bond in anion.GetBonds():
        if bond.GetBondType() in (BondType.DOUBLE, BondType.TRIPLE):
            return False

    # Rule 2: The anion must be one bond away from a pi-system.
    for neighbor in anion.GetNeighbors():
        for bond in neighbor.GetBonds():
            if bond.GetBondType() in (BondType.DOUBLE, BondType.TRIPLE):
                return True
    
    return False
```
#### Why are anions *within* a {math}`\pi`-system not considered stabilized?
For a lone pair (negative charge) to be stabilized by resonance, its orbital must be able to overlap with the adjacent {math}`\pi`-system's p-orbitals.

When an anion is a *part* of a double bond, the deprotonated atom is {math}`sp^2`-hybridized. Its lone pair occupies an {math}`sp^2` orbital, which lies in the plane of the molecule and is **orthogonal** to the p-orbitals forming the {math}`\pi`-bond. Since they cannot overlap, the charge is localized and not stabilized by resonance.

### Inductive Effects
The inductive effect describes the stabilization of a charge through the {math}`\sigma`-bond framework of a molecule. Electronegative atoms can pull electron density through these bonds, delocalizing the negative charge on the newly-formed anion and making it more stable.

The stabilization is quantified by the function as {py:func}`calculate_inductive_score(.) <chemistryinductivebias.src.proton_transfer.calculate_inductive_score>`. It calculates a score by summing up the electron-withdrawing contributions from all relevant atoms in the molecule, attenuated by their distance from the anion. By default, it considers atoms up to 4 bonds away from the anion.

```py
## Pseudocode
def calculate_inductive_score(mol, anion_idx):
    total_score = 0
    for atom in mol.GetAtoms():
        if atom == anion: continue

        path = find_shortest_path(atom, anion)
        distance = len(path) - 1

        # 1. Path must be through sigma-bonds only to isolate inductive effects
        if not is_sigma_only(path):
            continue

        # 2. Calculate the atom's electron-withdrawing strength
        # Based on its electronegativity and formal charge
        ewg_strength = (
            get_electronegativity_pull(atom)
            + get_formal_charge_effect(atom)
        )

        # 3. Attenuate the effect by distance
        attenuation_factor = 1 / (distance ** 2)
        contribution = ewg_strength * attenuation_factor
        total_score += contribution

    return total_score
```
#### Key Principles of the Calculation:
- **{math}`\sigma`-Bond Framework Only:** Effects transmitted through {math}`\pi`-systems are considered resonance, not induction.
- **Distance Attenuation:** The score follows a power law to compute the weakining of the inductive effect with distance, but can also be modeled with exponential decay.
- **Local Electronegative Difference:** The score computes how much more electronegative an atom is than its closest neighbor along the path toward the anion.
- **Boost for Charged Centers:** Atoms with positive formal charge are considered more powerful electron-withdrawing groups, and vice versa.


(step-4)=
## Step 4: Score and Aggregate
After computing the chemical properties of each candidate transformation in Step 3, the final step is to evaluate and rank them. This is implemented by a series of **object rules**, or heuristics based on fundamental principles of acid-base chemistry. Each rule assigns a qualitative favorability label to the transformation, which is then converted into a numerical score and then finally combined in a weighted sum to produce a final ranking.

The entire scoring process for a single candidate is handled by {py:func}`score_transformation_components(.) <chemistryinductivebias.src.proton_transfer.score_transformation_components>`, which calls an individual object rule function for each component.

### Object Rule 1: Formal Charge Rule
The formal charge object rule, implemented as {py:func}`formal_charge_object_rule(.) <chemistryinductivebias.src.proton_transfer.formal_charge_object_rule>`, assesses the change in the overall charge distribution, following a strict hierarchy of preferences based on electrostatic principles:

- **Charge Neutralization is `VERY_FAVORABLE`**: A reaction that reduces the total number of charges (e.g., `(A+)+(B-) -> AB`) is strongly preferred.
- **Increasing Charge Magnitude is `VERY_UNFAVORABLE`**: Making an existing charge more extreme (e.g., `(A-)+B -> (A2-) + (B+)`) is strongly penalized.
- **Intermolecular Charge Separation is `NEUTRAL`**: Creating charges between two different neutral molecules (e.g., `A+B -> (A+) + (B-)`) is the basis of many reactions and is not penalized.
- **Intramolecular Charge Separation is `VERY_UNFAVORABLE`**: Creating opposite charges within the same molecule from a neutral precursor (e.g., `(ab) + C -> ((a+)(b-)) + C`) is strongly penalized.
- **Charge-Conserved Reaction is `NEUTRAL`**: Reactions whose total charge is unchanged (e.g., `A + (B-) -> (A-) + B`) is not penalized.

```py
## Pseudocode
def formal_charge_rule(props):
    old_abs_charge = abs(props["charge_A_old"]) + abs(props["charge_B_old"])
    new_abs_charge = abs(props["charge_A_new"]) + abs(props["charge_B_new"])

    # (a) Charge neutralization is VERY_FAVORABLE
    if new_abs_charge < old_abs_charge:
        return VERY_FAVORABLE

    # (b) Increasing magnitude of an existing charge is VERY_UNFAVORABLE
    is_pileup = (props["charge_A_new"] > props["charge_A_old"] > 0) or \
                (props["charge_A_new"] < props["charge_A_old"] < 0)
    if is_pileup:
        return VERY_UNFAVORABLE

    # (c, d) Creating new charge is handled differently based on context
    if new_abs_charge > old_abs_charge:
        if props["is_interfragment_transfer"]:
            return NEUTRAL  # Intermolecular charge separation
        else:
            return VERY_UNFAVORABLE # Intramolecular charge separation
    
    # (e) For charge-conserved reactions, producing a neutral atom is preferred
    if new_abs_charge == old_abs_charge:
        if props["charge_A_new"] == 0 or props["charge_B_new"] == 0:
            return NEUTRAL
        else:
            return UNFAVORABLE # Charge is just redistributed
```

### Object Rule 2: Electronegativity Rule
The electronegativity object rule, implemented as {py:func}`electronegativity_object_rule(.) <chemistryinductivebias.src.proton_transfer.electronegativity_object_rule>`, evaluates the difference in electronegativity between the new anion (`B-`) and the new cation (`A+`).

```py
## Pseudocode
def electronegativity_rule(props):
    delta_en = props["EN_B"] - props["EN_A"]
    if delta_en > 0.5: return VERY_FAVORABLE
    if delta_en > 0.05: return FAVORABLE
    if delta_en > -0.05: return NEUTRAL
    if delta_en > -0.5: return UNFAVORABLE
    return VERY_UNFAVORABLE
```

### Object Rule 3: Atomic Radius Rule
The atomic radius object rule, implemented as {py:func}`atomic_radius_object_rule(.) <chemistryinductivebias.src.proton_transfer.atomic_radius_object_rule>`, evaluates the difference in atomic radii between the new anion (`B-`) and the new cation (`A+`).

```py
## Pseudocode
def atomic_radius_rule(props):
    # A negative charge is more stable on a larger atom.
    # The new anion is B, so we check if B is larger than A.
    delta_radius = props["Radius_B"] - props["Radius_A"]
    if delta_radius > 0.15: return VERY_FAVORABLE
    if delta_radius > 0.05: return FAVORABLE
    if delta_radius > -0.05: return NEUTRAL
    if delta_radius > -0.15: return UNFAVORABLE
    return VERY_UNFAVORABLE
```

### Object Rule 4: Resonance Rule
The resonance object rule, implemented as {py:func}`resonance_object_rule(.) <chemistryinductivebias.src.proton_transfer.resonance_object_rule>`, directly uses the {py:func}`is_resonance_stabilized(.) <chemistryinductivebias.src.proton_transfer.is_resonance_stabilized>` flag computed in Step 3.

```py
## Pseudocode
def resonance_rule(props):
    if props["is_resonance_stabilized"]:
        return VERY_FAVORABLE
    else:
        return NEUTRAL
```

### Object Rule 5: Inductive Effects Rule
The inductive effects object rule, implemented as {py:func}`inductive_object_rule(.) <chemistryinductivebias.src.proton_transfer.inductive_object_rule>`, converts the numerical `inductive_score` computed in Step 3 to a qualitative label using fixed thresholds.

```py
## Pseudocode
def inductive_rule(props):
    score = props["inductive_score"]
    if score > 2.0: # Strong stabilization (e.g., multiple EWGs)
        return VERY_FAVORABLE
    if score > 0.5: # Moderate stabilization (e.g., one strong EWG)
        return FAVORABLE
    else:
        return NEUTRAL
```

### Aggregation to a Final Score
The qualitative labels from each object rule are first mapped to numerical scores, and then combined using a weighted sum to produce the final score for the transformation.
```py
## Pseudocode
def calculate_final_score(transformation):
    # 1. Apply all object rules to get labels
    fc_label = formal_charge_rule(transformation.properties)
    en_label = electronegativity_rule(transformation.properties)
    res_label = resonance_rule(transformation.properties)
    ar_label = atomic_radius_rule(transformation.properties)
    ind_label = inductive_rule(transformation.properties)

    # 2. Convert labels to numerical scores
    SCORE_MAP = {
        "++": 2, "+": 1, "0": 0, "-": -1, "--": -2
    }
    fc_score = SCORE_MAP[fc_label]
    en_score = SCORE_MAP[en_label]
    # ... and so on for all five rules

    # 3. Calculate the weighted sum
    WEIGHTS = {
        "FC": 0.30, "EN": 0.25, "RES": 0.20, "AR": 0.15, "IND": 0.10
    }
    total_score = (WEIGHTS["FC"] * fc_score +
                   WEIGHTS["EN"] * en_score +
                   WEIGHTS["RES"] * res_score +
                   WEIGHTS["AR"] * ar_score +
                   WEIGHTS["IND"] * ind_score)
    
    return total_score
```
This final score provides a quantitative ranking of all possible proton transfer products, allowing the {py:class}`MechanismProducts<chemistryinductivebias.src.proton_transfer.MechanismProducts>` class to identify the most plausible outcome.


## Tunable Parameters
The parameters arbitrarily chosen and further tunable are summarized as the following:
```python
params = {
    'weights': (0.30, 0.25, 0.20, 0.15, 0.10),
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

### Weights for object rule scores

- Five floats for: Formal Charge, Electronegativity, Resonance, Atomic Radius, Inductive.
- Default: `(0.30, 0.25, 0.20, 0.15, 0.10)`.

### Score-to-label mapping

- Numeric value assigned to qualitative labels: `{ '++': 2, '+': 1, '0': 0, '-': -1, '--': -2 }`.

### Thresholds

- Electronegativity ΔEN thresholds: neutrality boundary `0.05`, favorability boundary `0.5`.
- Atomic radius Δr thresholds: neutrality boundary `0.05 Å`, favorability boundary `0.15 Å`.

### Inductive effect

- Final-score label thresholds: `> 0.5` → `FAVORABLE`, `> 2.0` → `VERY_FAVORABLE`.
- Hyperparameters:
  - `max_bonds`: 4 (maximum bond distance considered)
  - `power`: 2.0 (distance attenuation exponent, 1/d^power)
  - `charge_boost`: 0.5 (extra contribution from formal positive centers)
  - `en_margin`: 0.10 (ignore negligible ΔEN)