import string
from collections import defaultdict
from functools import reduce

import hypothesis.strategies as st
import pkg_resources
from rdkit.Chem import (
    BRICS,
    CombineMols,
    Mol,
    MolFromSmiles,
    MolToMolBlock,
    MolToSmiles,
    RemoveHs,
    SanitizeMol,
)
from rdkit.Chem.BRICS import reverseReactions

__all__ = ["mols", "smiles", "mol_blocks"]


with pkg_resources.resource_stream(__name__, "fragments.smi") as f:
    fragments = [MolFromSmiles(line) for line in f]

possible_dummy_labels = set(
    [
        a.GetIsotope()
        for fragment in fragments
        for a in fragment.GetAtoms()
        if a.GetSymbol() == "*"
    ]
)


def has_dummy_atom(fragment, dummy_label):
    if isinstance(dummy_label, Mol):
        dummy_label = dummy_label.GetAtomWithIdx(0).GetIsotope()
    return any(
        (
            a.GetSymbol() == "*" and a.GetIsotope() == dummy_label
            for a in fragment.GetAtoms()
        )
    )


# all reactions and fragments compatible with each other
# reaction = fragment + X --> product
# reaction = X + fragment --> product
compatible_reactions = {
    dummy_label: [
        (r, True)
        for r in reverseReactions
        if has_dummy_atom(r._matchers[0], dummy_label)
    ]
    + [
        (r, False)
        for r in reverseReactions
        if has_dummy_atom(r._matchers[1], dummy_label)
    ]
    for dummy_label in possible_dummy_labels
}

# mapping dummy atom d to fragments containing d
compatible_fragments = defaultdict(list)
for fragment in fragments:
    for a in fragment.GetAtoms():
        if a.GetSymbol() == "*":
            dummy_label = a.GetIsotope()
            compatible_fragments[dummy_label].append(fragment)

# sort the lists in each key by number of dummy atoms in the fragment and fragment size
def num_dummy_atoms(mol):
    return sum([1 for a in mol.GetAtoms() if a.GetSymbol() == "*"])


for dummy_label in possible_dummy_labels:
    compatible_fragments[dummy_label] = sorted(
        compatible_fragments[dummy_label],
        key=lambda x: (num_dummy_atoms(x), x.GetNumAtoms()),
    )


dummy_pattern = MolFromSmiles("[*]")


@st.composite
def mols(
    draw,
    name=st.text(alphabet=list(string.ascii_lowercase), min_size=1),
    n_connected_components=st.just(1),
):
    """Strategy for generating random molecules.

    Parameters
    ----------
    name : str or hypothesis.strategies.SearchStrategy
        Name of the molecule.
    n_connected_components : int or hypothesis.strategies.SearchStrategy
        Number of connected components in the molecule.
    """
    if isinstance(name, str):
        name = st.just(name)
    if isinstance(n_connected_components, int):
        n_connected_components = st.just(n_connected_components)

    n_connected_components = draw(n_connected_components)

    # draw connected components independently
    if n_connected_components > 1:
        components = [
            draw(mols(n_connected_components=1)) for _ in range(n_connected_components)
        ]
        seed = reduce(CombineMols, components)
    else:
        # repeat this process until a valid molecule is drawn
        while True:
            # start by drawing a random fragment and use it as seed
            seed = draw(st.sampled_from(fragments))

            # extend the seed by running compatible reactions on the seed
            # stop when molecule has no dummy atoms anymore (i.e. is fully built)
            while seed.HasSubstructMatch(dummy_pattern):
                while True:
                    # find a free substitution site
                    for atom in seed.GetAtoms():
                        if atom.GetSymbol() == "*":
                            dummy_label = atom.GetIsotope()
                            break

                    # draw a random reaction and fragment compatible with the substitution site
                    # note: compatible_mapping[dummy_label] was sorted by number of dummy atoms and size
                    #       --> shrinking is done automatically
                    # run the reaction and collect the products
                    reaction, right = draw(
                        st.sampled_from(compatible_reactions[dummy_label])
                    )
                    other_label = (
                        reaction._matchers[int(right)].GetAtomWithIdx(0).GetIsotope()
                    )
                    fragment = draw(st.sampled_from(compatible_fragments[other_label]))
                    if right:
                        products = reaction.RunReactants((seed, fragment))
                    else:
                        products = reaction.RunReactants((fragment, seed))

                    # remove duplicate products
                    unique_smiles = set()
                    results = []
                    if products:
                        for product in products:
                            s = MolToSmiles(product[0])
                            if s not in unique_smiles:
                                unique_smiles.add(s)
                                results.append(product[0])

                    # if reactions were unsuccessful, try again (works only because of randomly drawing reactions)
                    if len(results) == 0:
                        continue

                    # if reactions were successful, pick a random product and use it as the new seed
                    # note: shrinking is not applicable here, because all products have the same size
                    seed = draw(st.sampled_from(results))

                    break

            try:
                seed.UpdatePropertyCache()
                SanitizeMol(seed)
                seed = RemoveHs(seed)
                break
            except:
                continue

    name = draw(name)
    if name is not None:
        seed.SetProp("_Name", name)

    return seed


@st.composite
def smiles(draw, **kwargs):
    mol = draw(mols(**kwargs))

    return f"{MolToSmiles(mol)} {mol.GetProp('_Name')}"


@st.composite
def mol_blocks(draw, **kwargs):
    mol = draw(mols(**kwargs))

    return MolToMolBlock(mol)
