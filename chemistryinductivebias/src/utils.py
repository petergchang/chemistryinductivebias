from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors
from IPython.display import display


def _label_all_atoms(mol):
    # Return a copy of the molecule with all atoms explicitly labeled (shows C's too)
    labeled = Chem.Mol(mol)
    for atom in labeled.GetAtoms():
        atom.SetProp("atomLabel", atom.GetSymbol())
    return labeled


def show_candidates(
    mp, top_k=-1, mols_per_row=3, size=(300, 300), legend_mode="canonical+formula"
):
    # Prepare and display reactants first (explicit H)
    reactant_mols = []
    reactants_explicit_smiles = []
    if hasattr(mp, "_reactants"):
        for m in getattr(mp, "_reactants"):
            mh = Chem.AddHs(m)
            reactant_mols.append(_label_all_atoms(mh))
            reactants_explicit_smiles.append(
                Chem.MolToSmiles(mh, canonical=True, allHsExplicit=True)
            )
    if reactant_mols:
        reactant_legends = list(reactants_explicit_smiles)
        r_img = Draw.MolsToGridImage(
            reactant_mols,
            molsPerRow=mols_per_row,
            subImgSize=size,
            legends=reactant_legends,
        )
        display(r_img)

    # Prepare candidate product depictions
    entries = mp.ranked_unique()
    if top_k > 0:
        entries = entries[:top_k]
    mols = []
    legends = []
    for s in entries:
        prods = [Chem.AddHs(m) for m in s.products()]
        combined = prods[0]
        for m in prods[1:]:
            combined = Chem.CombineMols(combined, m)
        mols.append(_label_all_atoms(combined))
        canonical = [Chem.MolToSmiles(m, canonical=True) for m in s.products()]
        explicit = s.products_explicit_smiles()
        formulas = [
            rdMolDescriptors.CalcMolFormula(Chem.AddHs(m)) for m in s.products()
        ]
        if legend_mode == "explicit":
            desc = " . ".join(explicit)
        elif legend_mode == "formula":
            desc = " + ".join(formulas)
        elif legend_mode == "canonical":
            desc = " . ".join(canonical)
        else:  # canonical+formula
            desc = f"{' . '.join(canonical)} | {' + '.join(formulas)}"
        legends.append(f"{desc}\n{mp.describe(s)}")
    if mols:
        img = Draw.MolsToGridImage(
            mols, molsPerRow=mols_per_row, subImgSize=size, legends=legends
        )
        display(img)

    # Build textual summary with explicit-H SMILES
    total_considered = len(entries)
    summary_lines = []
    summary_lines.append(f"{total_considered} protonation candidates considered")
    if reactants_explicit_smiles:
        summary_lines.append(f"Reactants: {' . '.join(reactants_explicit_smiles)}")
    for i, s in enumerate(entries, start=1):
        explicit = s.products_explicit_smiles()
        comps = s.meta().get("components", {})
        score = s.score
        # Force deterministic component order: FC, EN, RES, AR, IND
        order = ("FC", "EN", "RES", "AR", "IND")
        comp_str = "{" + ", ".join(f"'{k}': {comps.get(k)}" for k in order) + "}"
        summary_lines.append(
            f"{i}. {' . '.join(explicit):<50} | score: {score:<10.3f} | components: {comp_str}"
        )
    summary_text = "\n".join(summary_lines)
    return summary_text
