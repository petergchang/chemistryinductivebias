from dataclasses import dataclass
from contextlib import suppress
from enum import Enum
from functools import cache
import operator
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import BondType

if TYPE_CHECKING:
    # Only for static typing; won't run at import-time
    from rdkit.Chem.ChemicalFeatures import (  # type: ignore[import-not-found]
        ChemicalFeatureFactory as ChemFeatureFactory,
    )
else:
    ChemFeatureFactory = Any  # type: ignore[misc]


class RuleLabel(Enum):
    """Enum for qualitative labels assigned by object rules."""

    VERY_FAVORABLE = "++"
    FAVORABLE = "+"
    NEUTRAL = "0"
    UNFAVORABLE = "-"
    VERY_UNFAVORABLE = "--"


SCORE_MAP = {
    RuleLabel.VERY_FAVORABLE: 2,
    RuleLabel.FAVORABLE: 1,
    RuleLabel.NEUTRAL: 0,
    RuleLabel.UNFAVORABLE: -1,
    RuleLabel.VERY_UNFAVORABLE: -2,
}


@dataclass
class Transformation:
    products: rdchem.Mol
    properties: Dict[str, Any]
    A_idx: int
    B_idx: int
    H_idx: int


@dataclass
class ScoredTransformation:
    """Bundle a transformation with per-rule labels and component scores.

    The total score is computed as a weighted sum of components and exposed via
    the `score` property; weights can be overridden by MechanismProducts.
    """

    transformation: Transformation
    fc_label: RuleLabel
    en_label: RuleLabel
    res_label: RuleLabel
    ar_label: RuleLabel
    ind_label: RuleLabel
    fc_score: int
    en_score: int
    res_score: int
    ar_score: int
    ind_score: int
    _weights: Tuple[float, float, float, float, float] = (0.30, 0.25, 0.20, 0.15, 0.10)

    @property
    def labels(self) -> Dict[str, str]:
        return {
            "FC": self.fc_label.value,
            "EN": self.en_label.value,
            "RES": self.res_label.value,
            "AR": self.ar_label.value,
            "IND": self.ind_label.value,
        }

    @property
    def score(self) -> float:
        w_fc, w_en, w_res, w_ar, w_ind = self._weights
        return (
            w_fc * self.fc_score
            + w_en * self.en_score
            + w_res * self.res_score
            + w_ar * self.ar_score
            + w_ind * self.ind_score
        )

    def with_weights(
        self, weights: Tuple[float, float, float, float, float]
    ) -> "ScoredTransformation":
        self._weights = weights
        return self

    def products(self) -> List[rdchem.Mol]:
        frags = Chem.GetMolFrags(
            self.transformation.products,
            asMols=True,
            sanitizeFrags=False,
        )
        out: List[rdchem.Mol] = []
        for m in frags:
            # Skip RemoveHs for hydrogen-only fragments to avoid RDKit warnings
            has_heavy = any(a.GetAtomicNum() != 1 for a in m.GetAtoms())
            out.append(Chem.RemoveHs(m) if has_heavy else m)
        return out

    def meta(self) -> Dict[str, Any]:
        t = self.transformation
        return {
            "score": self.score,
            "components": {
                "FC": self.fc_score,
                "EN": self.en_score,
                "RES": self.res_score,
                "AR": self.ar_score,
                "IND": self.ind_score,
            },
            "labels": self.labels,
            "A_idx": t.A_idx,
            "B_idx": t.B_idx,
            "weights": {
                "FC": self._weights[0],
                "EN": self._weights[1],
                "RES": self._weights[2],
                "AR": self._weights[3],
                "IND": self._weights[4],
            },
        }

    def products_canonical_smiles(self) -> Tuple[str, ...]:
        """Canonical SMILES (without explicit Hs) for use as a dedup key.

        Fragments are converted to canonical SMILES and sorted to be
        order-invariant.
        """
        smiles = [Chem.MolToSmiles(m, canonical=True) for m in self.products()]
        return tuple(sorted(smiles))

    def products_explicit_smiles(self) -> List[str]:
        """SMILES with explicit hydrogens for display."""
        out: List[str] = []
        for m in self.products():
            mh = Chem.AddHs(m)
            out.append(Chem.MolToSmiles(mh, canonical=True, allHsExplicit=True))
        return out


@dataclass
class MechanismOptions:
    """Options to control enumeration/validation behavior.

    - enumerate_all_acceptors: if True, treat every non-hydrogen atom as a
      potential acceptor (no feature/charge gating).
    - include_bridging_hydrogens: if True and enumerate_all_donors is enabled,
      include hydrogens with degree > 1 by pairing that H with each heavy
      neighbor.
    """

    enumerate_all_acceptors: bool = True
    include_bridging_hydrogens: bool = False


class MechanismProducts:
    """Holds all proton-transfer candidates with scoring and ranking.

    This class manages a collection of scored proton transfer transformations,
    providing methods to access the best candidate, get ranked lists, and
    generate serializable representations.

    Attributes:
        _weights: Tuple of weights for scoring components (FC, EN, RES, AR, IND)
        _scored: List of ScoredTransformation objects with applied weights
        _reactants: Original reactant molecules

    Methods:
        best(): Returns the highest-scoring candidate as (products, metadata)
        ranked(): Returns all candidates sorted by descending score
        ranked_unique(): Returns ranked candidates with duplicates removed
        to_ranked_dicts(): Returns serializable dicts for all ranked candidates
        to_ranked_dicts_unique(): Returns serializable dicts for unique ranked
        describe(): Generates human-readable description of a transformation
    """

    def __init__(
        self,
        scored: List[ScoredTransformation],
        reactants: List[rdchem.Mol],
        weights: Tuple[float, float, float, float, float] = (
            0.30,
            0.25,
            0.20,
            0.15,
            0.10,
        ),
    ):
        # assign weights to scored item so that `score` reflects current policy
        self._weights = weights
        self._scored = [s.with_weights(weights) for s in scored]
        self._reactants = reactants.copy()

    def best(self) -> Tuple[List[rdchem.Mol], Dict[str, Any]]:
        if not self._scored:
            return self._reactants, {
                "score": 0,
                "labels": {},
                "note": "no valid proton transfers generated",
            }
        s = max(self._scored, key=lambda x: x.score)
        return s.products(), s.meta()

    def ranked(self) -> List[ScoredTransformation]:
        return sorted(self._scored, key=lambda x: x.score, reverse=True)

    def to_ranked_dicts(self) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = [
            {
                "products": s.products(),
                "score": s.score,
                "components": s.meta()["components"],
                "labels": s.labels,
                "A_idx": s.transformation.A_idx,
                "B_idx": s.transformation.B_idx,
                "description": self.describe(s),
                "weights": s.meta()["weights"],
            }
            for s in self.ranked()
        ]
        return out

    def ranked_unique(self) -> List[ScoredTransformation]:
        """Ranked list with duplicates removed by canonical product SMILES."""
        seen = set()
        unique: List[ScoredTransformation] = []
        for s in self.ranked():
            key = s.products_canonical_smiles()
            if key in seen:
                continue
            seen.add(key)
            unique.append(s)
        return unique

    def to_ranked_dicts_unique(self) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = [
            {
                "products": s.products(),
                "score": s.score,
                "components": s.meta()["components"],
                "labels": s.labels,
                "A_idx": s.transformation.A_idx,
                "B_idx": s.transformation.B_idx,
                "description": self.describe(s),
                "weights": s.meta()["weights"],
            }
            for s in self.ranked_unique()
        ]
        return out

    def _combined_with_offsets(self) -> Tuple[rdchem.Mol, List[int]]:
        # Rebuild the combined mol. (with explicit Hs), return offsets mapping
        mols_h = [Chem.AddHs(m) for m in self._reactants]
        offsets: List[int] = []
        n = 0
        for m in mols_h:
            offsets.append(n)
            n += m.GetNumAtoms()
        combined = mols_h[0]
        for m in mols_h[1:]:
            combined = Chem.CombineMols(combined, m)
        return combined, offsets

    def _map_global_to_reactant(self, idx: int, offsets: List[int]) -> Tuple[int, int]:
        # Return (reactant_index, local_atom_index)
        mols_h = [Chem.AddHs(m) for m in self._reactants]
        for i, m in enumerate(mols_h):
            start = offsets[i]
            end = start + m.GetNumAtoms()
            if start <= idx < end:
                return i, idx - start
        raise IndexError(idx)

    def describe(self, s: ScoredTransformation) -> str:
        """Return a human-readable sentence: 'Lone pair on A attacks H on B'."""
        combined, offsets = self._combined_with_offsets()
        A_idx = s.transformation.A_idx
        B_idx = s.transformation.B_idx
        atomA = combined.GetAtomWithIdx(A_idx)
        atomB = combined.GetAtomWithIdx(B_idx)
        # Map to reactants and create short names using SMILES
        a_r, _ = self._map_global_to_reactant(A_idx, offsets)
        b_r, _ = self._map_global_to_reactant(B_idx, offsets)
        reactant_names: List[str] = []
        for m in self._reactants:
            has_heavy = any(a.GetAtomicNum() != 1 for a in m.GetAtoms())
            mol_for_smiles = Chem.RemoveHs(m) if has_heavy else m
            reactant_names.append(Chem.MolToSmiles(mol_for_smiles, canonical=True))
        label_str = ", ".join([f"{k}: {v}" for k, v in s.labels.items()])
        return (
            f"Lone pair on {atomA.GetSymbol()} in {reactant_names[a_r]} "
            f"attacks H on {atomB.GetSymbol()} in {reactant_names[b_r]}\n"
            f"({label_str})"
        )


@cache
def periodic_en(atomic_number: int) -> float:
    """Return Pauling electronegativity for an atomic number.

    Uses pymatgen's periodic table to retrieve electronegativity values.
    Raises ValueError if pymatgen is unavailable or fails to return a value.

    Args:
        atomic_number (int): The atomic number of the element.

    Returns:
        float: The Pauling electronegativity value.

    Raises:
        ValueError: If pymatgen import fails or electronegativity is unavailable.
    """
    try:
        # Lazy import so pymatgen is an optional dependency
        from pymatgen.core.periodic_table import (  # pylint: disable=import-outside-toplevel
            Element as PmgElement,
        )

        x = PmgElement.from_Z(atomic_number).X
        if x is not None:
            return x
    except Exception as e:
        # Any import/runtime issue: fall back to local table
        raise ValueError(f"Failed to import pymatgen: {e}") from e
    return 0.0


def get_en(atom: rdchem.Atom) -> float:
    """Wrapper that accepts an RDKit Atom and returns Pauling EN."""
    return periodic_en(atom.GetAtomicNum())


def get_radius_angstrom(atom: rdchem.Atom) -> float:
    """
    Get the atomic radius in angstroms for a given atom.

    Prefers covalent radius; falls back to van der Waals radius if unavailable.

    Args:
        atom (rdchem.Atom): The atom whose radius is to be retrieved.

    Returns:
        float: The atomic radius in angstroms.
    """
    r = 0.0
    pt = rdchem.GetPeriodicTable()

    with suppress(Exception):
        r = float(pt.GetRcovalent(atom.GetAtomicNum()))
    if not r or r <= 0:
        with suppress(Exception):
            r = float(pt.GetRvdw(atom.GetAtomicNum()))
        if not r or r <= 0:
            r = 1.5
    return r


def _get_system_smiles_key(mol: rdchem.Mol) -> Tuple[str, ...]:
    """
    Generate a canonical, order-invariant SMILES key for a molecule.

    This function creates a unique identifier for molecules that may consist of
    multiple disconnected fragments by generating canonical SMILES for each
    fragment and sorting them to ensure order invariance.

    Args:
        mol (rdchem.Mol): The molecule to generate a SMILES key for.

    Returns:
        Tuple[str, ...]: A sorted tuple of canonical SMILES strings, one for
                        each fragment in the molecule.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    # Use canonical SMILES without explicit Hs for a robust key
    smiles_list = [Chem.MolToSmiles(f, canonical=True) for f in frags]
    return tuple(sorted(smiles_list))


def _build_feature_factory() -> Optional[ChemFeatureFactory]:
    try:
        # Local import so it only happens when needed
        from rdkit import RDConfig  # pylint: disable=import-outside-toplevel
        from rdkit.Chem import (  # pylint: disable=import-outside-toplevel
            ChemicalFeatures,
        )

        data_dir = RDConfig.RDDataDir
        # Make sure it's usable as a path
        try:
            data_dir_path = Path(data_dir)
        except TypeError:
            return None

        fdef = data_dir_path / "BaseFeatures.fdef"
        if fdef.exists():
            return ChemicalFeatures.BuildFeatureFactory(str(fdef))
    except Exception:  # pylint: disable=broad-exception-caught
        return None
    return None


FEATURE_FACTORY: Optional[ChemFeatureFactory] = None


def get_feature_factory() -> Optional[ChemFeatureFactory]:
    global FEATURE_FACTORY  # pylint: disable=global-statement
    if FEATURE_FACTORY is None:
        FEATURE_FACTORY = _build_feature_factory()
    return FEATURE_FACTORY


def _feature_atom_indices(mol: rdchem.Mol, family: str) -> List[int]:
    """
    Get atom indices for a specific chemical feature family using RDKit's
    feature factory.

    Args:
        mol (rdchem.Mol): The molecule to analyze for features.
        family (str): The feature family name (e.g., "Acceptor", "Donor").

    Returns:
        List[int]: Sorted list of unique atom indices that match the feature
        family.
    """
    # Ensure ring info/property cache are initialized for feature detection
    with suppress(Exception):
        mol.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.GetSymmSSSR(mol)
    factory = get_feature_factory()
    if factory is None:
        raise ValueError("Failed to get feature factory")
    feats = factory.GetFeaturesForMol(mol)
    indices: List[int] = []
    for f in feats:
        if f.GetFamily() == family:
            indices.extend(f.GetAtomIds())
    return sorted(set(indices))


def find_acceptors(mol: rdchem.Mol, options: MechanismOptions) -> List[int]:
    """Return indices of atoms that can act as proton acceptors.

    Uses RDKit's feature factory (Acceptor) and augments with
    any atoms carrying a negative formal charge.

    Args:
        mol (rdchem.Mol): The molecule to analyze for acceptor atoms.
        options (MechanismOptions): Configuration options.

    Returns:
        List[int]: Sorted list of atom indices that can act as proton acceptors.
    """
    # Exhaustive mode: all non-hydrogen atoms are potential acceptors,
    # plus allow hydride (H with -1 formal charge) to act as an acceptor
    if options.enumerate_all_acceptors:
        acceptors: List[int] = []
        for a in mol.GetAtoms():
            z = a.GetAtomicNum()
            q = a.GetFormalCharge()
            if z != 1:
                acceptors.append(a.GetIdx())
            else:
                # include hydride as acceptor so H2 formation is possible
                if q == -1:
                    acceptors.append(a.GetIdx())
        return acceptors

    acceptor_atoms = set(_feature_atom_indices(mol, "Acceptor"))
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            acceptor_atoms.add(atom.GetIdx())

    # # Also include neutral atoms with likely lone pairs (N/O/S) depending on
    # # charge policy
    # for atom in mol.GetAtoms():
    #     if atom.GetAtomicNum() in (7, 8, 16):
    #         acceptor_atoms.add(atom.GetIdx())
    return sorted(acceptor_atoms)


def find_acidic_hydrogens(
    mol: rdchem.Mol, options: MechanismOptions
) -> List[Tuple[int, int]]:
    """Find hydrogen atoms that can be donated in proton transfer reactions.

    Args:
        mol (rdchem.Mol): The molecule to analyze for acidic hydrogens.
        options (MechanismOptions): Configuration options.
    Returns:
        List[Tuple[int, int]]: List of (donor_atom_idx, hydrogen_idx) pairs
            representing potential proton donors. Each tuple contains the index
            of the heavy atom (donor) and the index of the hydrogen atom that
            can be transferred.

    Notes:
        By default, only considers hydrogens bonded to exactly one heavy atom
        (excludes bridging hydrogens). If options.include_bridging_hydrogens is
        True, includes all hydrogen-heavy atom pairs where the hydrogen is
        bonded to multiple heavy atoms.
    """
    pairs_set: set[Tuple[int, int]] = set()

    for h in mol.GetAtoms():
        if h.GetAtomicNum() != 1:
            continue
        neighbors = list(h.GetNeighbors())
        if not options.include_bridging_hydrogens:
            if len(neighbors) != 1:
                continue
            heavy = neighbors[0]
            if heavy.GetAtomicNum() != 1:
                pairs_set.add((heavy.GetIdx(), h.GetIdx()))
        else:
            for nb in neighbors:
                if nb.GetAtomicNum() != 1:
                    pairs_set.add((nb.GetIdx(), h.GetIdx()))
    return sorted(pairs_set)


def apply_proton_transfer(
    mol: rdchem.Mol, A_idx: int, H_idx: int, B_idx: int
) -> Optional[rdchem.Mol]:
    """
    Apply a proton transfer from atom B (donor) to atom A (acceptor) via
    hydrogen H.

    Args:
        mol (rdchem.Mol): The molecule to modify.
        A_idx (int): Index of the acceptor atom.
        H_idx (int): Index of the hydrogen atom to transfer.
        B_idx (int): Index of the donor atom.

    Returns:
        rdchem.Mol | None: The new molecule after proton transfer, or None if
        invalid.
    """
    # Make an editable copy
    rwm = rdchem.RWMol(mol)
    # Suppress RDKit warnings/errors during tentative edits and sanitize
    RDLogger.DisableLog("rdApp.error")
    RDLogger.DisableLog("rdApp.warning")
    try:
        if rwm.GetBondBetweenAtoms(B_idx, H_idx) is None:
            return None
        # Break B-H and form A-H
        try:
            rwm.RemoveBond(B_idx, H_idx)
        except Exception:  # pylint: disable=broad-exception-caught
            return None
        try:
            rwm.AddBond(A_idx, H_idx, BondType.SINGLE)
        except Exception:  # pylint: disable=broad-exception-caught
            return None

        atomA = rwm.GetAtomWithIdx(A_idx)
        atomB = rwm.GetAtomWithIdx(B_idx)
        atomA.SetFormalCharge(atomA.GetFormalCharge() + 1)
        atomB.SetFormalCharge(atomB.GetFormalCharge() - 1)

        newmol = rwm.GetMol()
        try:
            Chem.SanitizeMol(newmol)
        except Exception:  # pylint: disable=broad-exception-caught
            # Discard invalid structures (e.g., valence violations)
            return None
        return newmol
    finally:
        RDLogger.EnableLog("rdApp.error")
        RDLogger.EnableLog("rdApp.warning")


def _is_conjugated_to_pi(mol: rdchem.Mol, idx: int) -> bool:
    """
    Check if an atom is conjugated to a pi system within 2 bonds.

    This function determines if an atom is adjacent to a pi system
    by checking for aromatic character, conjugated bonds, or double/triple bonds
    within a 2-bond radius.

    Args:
        mol (rdchem.Mol): The molecule to analyze.
        idx (int): The index of the atom to check.

    Returns:
        bool: True if the atom is conjugated to a pi system, False otherwise.
    """
    a = mol.GetAtomWithIdx(idx)

    # If the anionic atom is itself part of an aromatic ring or a double/triple
    # bond, its lone pair is likely in an sp/sp2 hybrid orbital,
    # orthogonal to the adjacent pi system, making resonance impossible.
    if a.GetIsAromatic():
        return False

    for bond in a.GetBonds():
        if bond.GetBondType() in (BondType.DOUBLE, BondType.TRIPLE):
            return False

    for b in a.GetBonds():
        # 1 bond away
        nb = b.GetOtherAtom(a)
        if nb.GetIsAromatic():
            return True
        for b2 in nb.GetBonds():
            if b2.GetIdx() == b.GetIdx():
                continue
            if (
                b2.GetIsConjugated()
                or b2.GetIsAromatic()
                or b2.GetBondType() in (BondType.DOUBLE, BondType.TRIPLE)
            ):
                return True
    return False


def _is_resonance_stabilized(mol: rdchem.Mol, anion_idx: int) -> bool:
    """
    Check if an anion is resonance-stabilized by conjugation to a pi system.

    This function determines if a negatively charged atom can be stabilized
    through resonance by being conjugated to a pi system.

    Args:
        mol (rdchem.Mol): The molecule containing the anion.
        anion_idx (int): The index of the atom to check for anion stabilization.

    Returns:
        bool: True if the anion is resonance-stabilized, False otherwise.
    """
    anion_atom = mol.GetAtomWithIdx(anion_idx)
    if anion_atom.GetFormalCharge() >= 0:
        return False

    return _is_conjugated_to_pi(mol, anion_idx)


def _delta_en_above_baseline(
    mol, atom_idx: int, neighbor_toward_anion_idx: int | None
) -> float:
    """
    Calculate the electronegativity difference above a baseline for inductive
    effects.

    This function computes the inductive pull as the difference between an
    atom's electronegativity and a baseline value. The baseline is either the
    electronegativity of a neighboring atom (if provided) or carbon's
    electronegativity as a typical sigma framework reference.

    Args:
        mol: The molecule containing the atom.
        atom_idx (int): The index of the atom whose electronegativity is being
            compared.
        neighbor_toward_anion_idx (int | None): The index of a neighboring atom
            to use as baseline, or None to use carbon's electronegativity.

    Returns:
        float: The electronegativity difference above baseline, with a minimum
            of 0.0.
    """
    en_atom = get_en(mol.GetAtomWithIdx(atom_idx))
    if neighbor_toward_anion_idx is not None:
        en_base = get_en(mol.GetAtomWithIdx(neighbor_toward_anion_idx))
    else:
        # Pauling EN of carbon as a reasonable σ baseline
        en_base = get_en(Chem.MolFromSmiles("C"))
    return max(0.0, en_atom - en_base)


def _calculate_inductive_score(
    mol: rdchem.Mol,
    anion_idx: int,
    max_bonds: int = 4,
    power: float = 2.0,
    use_exponential: bool = False,
    rho: float = 0.6,  # if exponential, weight ~ rho**distance
    charge_boost: float = 0.5,
    en_margin: float = 0.10,  # ignore tiny ΔEN
) -> float:
    """
    Calculate distance-attenuated sigma-only inductive pull on an anion.

    This function computes the inductive stabilization of an anion by
    considering electron-withdrawing effects from nearby atoms through sigma
    bonds only. The contribution from each atom is attenuated by distance and
    filtered to exclude resonance pathways.

    Args:
        mol (rdchem.Mol): The molecule containing the anion.
        anion_idx (int): The index of the anion atom.
        max_bonds (int, optional): Maximum number of bonds to consider.
            Defaults to 4.
        power (float, optional): Power for distance attenuation when not using
            exponential. Defaults to 2.0.
        use_exponential (bool, optional): Whether to use exponential decay
            instead of power law. Defaults to False.
        rho (float, optional): Decay factor for exponential attenuation.
            Defaults to 0.6.
        charge_boost (float, optional): Multiplier for formal charge
            contributions. Defaults to 1.5.
        en_margin (float, optional): Minimum electronegativity difference to
            consider. Defaults to 0.10.

    Returns:
        float: The calculated inductive score representing stabilization
            strength.

    Notes:
        - Only considers atoms within max_bonds by shortest path
        - Discards paths that are not sigma-only to prevent resonance
          double-counting
        - Contribution scales as max(0, ΔEN above local baseline)/distance^power
          or ΔEN * rho^distance
        - Adds |formal charge| * charge_boost term for charged centers on
          sigma-only paths
    """
    score = 0.0
    for atom in mol.GetAtoms():
        ai = atom.GetIdx()
        if ai == anion_idx:
            continue

        path = Chem.GetShortestPath(mol, ai, anion_idx)
        d = len(path) - 1
        if d <= 0 or d > max_bonds:
            continue

        # Check if the path is sigma-only
        is_sigma_only = True
        for i in range(len(path) - 1):
            bond = mol.GetBondBetweenAtoms(path[i], path[i + 1])
            if bond is None or bond.GetBondType() != BondType.SINGLE:
                is_sigma_only = False
                break

        if not is_sigma_only:
            continue

        # neighbor along the path that is closer to the anion (for local EN
        # baseline)
        neighbor_toward_anion = path[-2] if d >= 1 else None
        dEN = _delta_en_above_baseline(mol, ai, neighbor_toward_anion)
        if dEN < en_margin and atom.GetFormalCharge() == 0:
            continue

        attenuation = (rho**d) if use_exponential else (d**power)
        contrib = (dEN / attenuation) if not use_exponential else (dEN * attenuation)

        # Positively charged centers withdraw strongly via field/inductive
        # effects. and vice versa for negatively charged centers.
        q = int(atom.GetFormalCharge())
        if use_exponential:
            contrib += charge_boost * q * (rho**d)
        else:
            contrib += charge_boost * q / (d**power)

        score += contrib

    return score


def compute_transformation_properties(
    old: rdchem.Mol, new: rdchem.Mol, A_idx: int, B_idx: int
) -> Dict[str, Any]:
    """
    Compute properties of a proton transfer transformation.

    Args:
        old (rdchem.Mol): The original molecule.
        new (rdchem.Mol): The molecule after transformation.
        A_idx (int): Index of the acceptor atom.
        B_idx (int): Index of the donor atom.

    Returns:
        Dict[str, Any]: Dictionary of computed properties.
    """
    a_old = old.GetAtomWithIdx(A_idx)
    b_old = old.GetAtomWithIdx(B_idx)
    a_new = new.GetAtomWithIdx(A_idx)
    b_new = new.GetAtomWithIdx(B_idx)
    a_old_charge = a_old.GetFormalCharge()
    b_old_charge = b_old.GetFormalCharge()
    a_new_charge = a_new.GetFormalCharge()
    b_new_charge = b_new.GetFormalCharge()

    # Determine whether A and B came from different initial fragments
    # (different reactants)
    frag_tuples = Chem.GetMolFrags(old, asMols=False, sanitizeFrags=False)
    atom_to_frag: Dict[int, int] = {}
    for frag_id, atoms in enumerate(frag_tuples):
        for ai in atoms:
            atom_to_frag[int(ai)] = frag_id
    is_interfragment = atom_to_frag.get(A_idx, -1) != atom_to_frag.get(B_idx, -1)

    props = {
        "delta_charge_on_A": float(a_new_charge - a_old_charge),
        "delta_charge_on_B": float(b_new_charge - b_old_charge),
        "charge_A_old": a_old_charge,
        "charge_B_old": b_old_charge,
        "charge_A_new": a_new_charge,
        "charge_B_new": b_new_charge,
        "is_interfragment_transfer": is_interfragment,
        "EN_A": get_en(a_old),
        "EN_B": get_en(b_old),
        "Radius_A": get_radius_angstrom(a_old),
        "Radius_B": get_radius_angstrom(b_old),
        "is_resonance_stabilized": _is_resonance_stabilized(new, B_idx),
        "inductive_score": _calculate_inductive_score(new, B_idx),
    }
    return props


def electronegativity_object_rule(transformation: Transformation) -> RuleLabel:
    """
    Assign a qualitative label to a transforma based on electronegativity diff.

    Args:
        transformation (Transformation): The transformation to evaluate.

    Returns:
        RuleLabel: One of VERY_FAVORABLE, FAVORABLE, UNFAVORABLE,
            VERY_UNFAVORABLE, or NEUTRAL.
    """
    prop = transformation.properties
    delta_en = float(prop["EN_A"] - prop["EN_B"])
    if prop["delta_charge_on_A"] > 0:
        if delta_en < -0.5:
            return RuleLabel.VERY_FAVORABLE
        if delta_en < -0.05:
            return RuleLabel.FAVORABLE
        if delta_en < 0.05:
            return RuleLabel.NEUTRAL
        if delta_en < 0.5:
            return RuleLabel.UNFAVORABLE
        return RuleLabel.VERY_UNFAVORABLE
    if prop["delta_charge_on_A"] < 0:
        if delta_en > 0.5:
            return RuleLabel.VERY_FAVORABLE
        if delta_en > 0.05:
            return RuleLabel.FAVORABLE
        if delta_en < -0.05:
            return RuleLabel.NEUTRAL
        if delta_en < -0.5:
            return RuleLabel.UNFAVORABLE
        return RuleLabel.VERY_UNFAVORABLE
    return RuleLabel.NEUTRAL


def atomic_radius_object_rule(transformation: Transformation) -> RuleLabel:
    """
    Assign a qualitative label to a transform based on atomic radius diff.

    Args:
        transformation (Transformation): The transformation to evaluate.

    Returns:
        RuleLabel: One of VERY_FAVORABLE, FAVORABLE, UNFAVORABLE,
            VERY_UNFAVORABLE, or NEUTRAL.
    """
    prop = transformation.properties
    if prop["delta_charge_on_A"] == 0:
        return RuleLabel.NEUTRAL
    delta_r = float(prop["Radius_B"] - prop["Radius_A"])  # in Angstrom
    # Thresholds adapted from 15 pm -> 0.15 Å
    if delta_r > 0.15:
        return RuleLabel.VERY_FAVORABLE
    if delta_r > 0.05:
        return RuleLabel.FAVORABLE
    if delta_r > -0.05:
        return RuleLabel.NEUTRAL
    if delta_r > -0.15:
        return RuleLabel.UNFAVORABLE
    return RuleLabel.VERY_UNFAVORABLE


def resonance_object_rule(transformation: Transformation) -> RuleLabel:
    """
    Assign a qualitative label based on the presence of resonance stabilization.

    This rule evaluates the resonance stabilization of the anion formed
    after proton transfer.

    Args:
        transformation (Transformation): The transformation to evaluate.

    Returns:
        RuleLabel: VERY_FAVORABLE if anion is resonance-stabilized, NEUTRAL
            otherwise.
    """
    if transformation.properties.get("is_resonance_stabilized", False):
        return RuleLabel.VERY_FAVORABLE
    return RuleLabel.NEUTRAL


def inductive_object_rule(transformation: Transformation) -> RuleLabel:
    """
    Assign a qualitative label based on the calculated inductive score.

    This rule evaluates the inductive stabilization of the conjugate base
    through electron-withdrawing effects from nearby atoms via sigma bonds.
    Higher scores indicate more stabilization.

    Args:
        transformation (Transformation): The transformation to evaluate.

    Returns:
        RuleLabel: VERY_FAVORABLE for strong inductive effects (score > 2.0),
                  FAVORABLE for moderate effects (score > 0.5), NEUTRAL
                  otherwise.
    """
    score = transformation.properties.get("inductive_score", 0.0)
    # These thresholds can be tuned based on desired sensitivity.
    if score > 2.0:  # e.g., Multiple strong EWGs like in trichloroacetic acid
        return RuleLabel.VERY_FAVORABLE
    if score > 0.5:  # e.g., A single strong EWG nearby like in chloroacetic acid
        return RuleLabel.FAVORABLE
    return RuleLabel.NEUTRAL


def formal_charge_object_rule(transformation: Transformation) -> RuleLabel:
    """
    Assign a qualitative label based on formal charge neutrality preference.

    This rule prefers transformations that result in neutral formal charges
    on both the acceptor and donor atoms after proton transfer.

    Args:
        transformation (Transformation): The transformation to evaluate.

    Returns:
        RuleLabel: VERY_FAVORABLE if both A and B are neutral in products,
                  FAVORABLE if exactly one is neutral, UNFAVORABLE if neither
                  is neutral.
    """
    props = transformation.properties
    a_q_old = int(props.get("charge_A_old", 0))
    b_q_old = int(props.get("charge_B_old", 0))
    a_q_new = int(props.get("charge_A_new", 0))
    b_q_new = int(props.get("charge_B_new", 0))

    old_abs_sum = abs(a_q_old) + abs(b_q_old)
    new_abs_sum = abs(a_q_new) + abs(b_q_new)

    # Rule 1: Neutralization is always highly favorable
    if new_abs_sum < old_abs_sum:
        return RuleLabel.VERY_FAVORABLE

    # Rule 2: Creation of more extreme charges is highly unfavorable
    charge_pileup_A = (a_q_new > a_q_old > 0) or (a_q_new < a_q_old < 0)
    charge_pileup_B = (b_q_new > b_q_old > 0) or (b_q_new < b_q_old < 0)
    if charge_pileup_A or charge_pileup_B:
        return RuleLabel.VERY_UNFAVORABLE

    # Rule 3: Creation of charge from neutral is neutral if intermolecular,
    # otherwise unfavorable
    if new_abs_sum > old_abs_sum:
        if props.get("is_interfragment_transfer", False):
            return RuleLabel.NEUTRAL
        return RuleLabel.VERY_UNFAVORABLE

    # Rule 4: If charge is conserved
    if 0 in (a_q_new, b_q_new):
        return RuleLabel.NEUTRAL

    return RuleLabel.UNFAVORABLE


def score_transformation_components(
    t: Transformation,
) -> Tuple[
    RuleLabel, int, RuleLabel, int, RuleLabel, int, RuleLabel, int, RuleLabel, int
]:
    """
    Return per-rule labels and scores for all five scoring components.

    Args:
        t (Transformation): The transformation to score.

    Returns:
        Tuple containing (FC_label, FC_score, EN_label, EN_score, RES_label,
            RES_score, AR_label, AR_score, IND_label, IND_score).
    """
    fc_label = formal_charge_object_rule(t)
    en_label = electronegativity_object_rule(t)
    res_label = resonance_object_rule(t)
    ar_label = atomic_radius_object_rule(t)
    ind_label = inductive_object_rule(t)
    fc_score = SCORE_MAP[fc_label]
    en_score = SCORE_MAP[en_label]
    res_score = SCORE_MAP[res_label]
    ar_score = SCORE_MAP[ar_label]
    ind_score = SCORE_MAP[ind_label]
    return (
        fc_label,
        fc_score,
        en_label,
        en_score,
        res_label,
        res_score,
        ar_label,
        ar_score,
        ind_label,
        ind_score,
    )


def combine_reactants(reactants: List[rdchem.Mol]) -> rdchem.Mol:
    """
    Combine a list of reactant molecules into a single molecule with explicit
    hydrogens.

    Args:
        reactants (List[rdchem.Mol]): List of reactant molecules.

    Returns:
        rdchem.Mol: The combined molecule.

    Raises:
        ValueError: If any reactant is None.
    """
    if any(m is None for m in reactants):
        bad_idx = [i for i, m in enumerate(reactants) if m is None]
        raise ValueError(
            f"Invalid reactant(s) at indices {bad_idx}: one or more SMILES "
            f"failed to parse"
        )
    mols_h = [Chem.AddHs(m) for m in reactants]
    combined = mols_h[0]
    for m in mols_h[1:]:
        combined = Chem.CombineMols(combined, m)
    # Initialize ring info on the combined molecule as features may require it
    with suppress(Exception):
        combined.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(combined)
    return combined


def proton_transfer_process_rule(
    reactants: List[rdchem.Mol], options: MechanismOptions
) -> List[Transformation]:
    """
    Generate all possible proton transfer transformations for a set of
    reactants.

    Args:
        reactants (List[rdchem.Mol]): List of reactant molecules.
        options (MechanismOptions): Configuration options controlling
            enumeration behavior.

    Returns:
        List[Transformation]: List of possible transformations.
    """
    combined = combine_reactants(reactants)
    reactant_key = _get_system_smiles_key(combined)
    acceptors = find_acceptors(combined, options)
    donors = find_acidic_hydrogens(combined, options)  # list of (B_idx, H_idx)

    candidates: List[Transformation] = []
    for A_idx in acceptors:
        for B_idx, H_idx in donors:
            if A_idx in (B_idx, H_idx):
                continue
            newmol = apply_proton_transfer(combined, A_idx, H_idx, B_idx)
            if newmol is None:
                continue

            # Skip non-productive transformations.
            if _get_system_smiles_key(newmol) == reactant_key:
                continue

            props = compute_transformation_properties(combined, newmol, A_idx, B_idx)
            candidates.append(
                Transformation(
                    products=newmol,
                    properties=props,
                    A_idx=A_idx,
                    B_idx=B_idx,
                    H_idx=H_idx,
                )
            )
    return candidates


def find_best_transformation(
    scored: List[Tuple[Transformation, int, Dict[str, str]]],
) -> Optional[Tuple[Transformation, int, Dict[str, str]]]:
    """
    Find the best transformation from a list of scored transformations.

    Args:
        scored (List[Tuple[Transformation, int, Dict[str, str]]]):
            List of (Transformation, score, labels) tuples.

    Returns:
        Tuple[Transformation, int, Dict[str, str]] | None:
            The best transformation tuple, or None if list is empty.
    """
    if not scored:
        return None
    return max(scored, key=operator.itemgetter(1))


def proton_transfer_predict(
    reactants: List[rdchem.Mol],
    options: MechanismOptions | None = None,
    weights: Tuple[float, float, float, float, float] = (0.30, 0.25, 0.20, 0.15, 0.10),
) -> MechanismProducts:
    """Generate all proton-transfer candidates and wrap them as
    MechanismProducts.

    This function enumerates all possible proton transfer transformations
    between the given reactant molecules, scores them using multiple criteria,
    and returns them wrapped in a MechanismProducts object for further analysis.

    Args:
        reactants (List[rdchem.Mol]): List of RDKit molecule objects
            representing the reactants for proton transfer reactions.
        options (MechanismOptions | None, optional): Configuration options
            controlling enumeration behavior such as which atoms to consider as
            acceptors/donors. If None, default options are used. Defaults to
            None.
        weights (Tuple[float, float, float, float, float], optional): Weights
            for combining the five scoring components (formal charge,
            electronegativity, resonance, atomic radius, inductive). Must sum
            to 1.0 or less. Defaults to (0.30, 0.25, 0.20, 0.15, 0.10).

    Returns:
        MechanismProducts: A container object holding all scored transformations
            with methods for ranking, filtering, and analysis.
    """
    if options is None:
        options = MechanismOptions()
    # Step 1: Apply process rule.
    candidates = proton_transfer_process_rule(reactants, options)

    # Step 2: Apply object rules.
    scored: List[ScoredTransformation] = []
    for t in candidates:
        (
            fc_label,
            fc_score,
            en_label,
            en_score,
            res_label,
            res_score,
            ar_label,
            ar_score,
            ind_label,
            ind_score,
        ) = score_transformation_components(t)

        scored.append(
            ScoredTransformation(
                transformation=t,
                fc_label=fc_label,
                en_label=en_label,
                res_label=res_label,
                ar_label=ar_label,
                ind_label=ind_label,
                fc_score=fc_score,
                en_score=en_score,
                res_score=res_score,
                ar_score=ar_score,
                ind_score=ind_score,
            )
        )

    # Step 3: Return the results.
    return MechanismProducts(scored=scored, reactants=reactants, weights=weights)
