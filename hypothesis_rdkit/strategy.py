import logging
import pickle
import string
from functools import reduce

import hypothesis.strategies as st
import pkg_resources
from rdkit.Chem import (
    BondType,
    CombineMols,
    MolFromSmiles,
    MolToMolBlock,
    MolToSmiles,
    RemoveHs,
    SanitizeMol,
)

from .fragment import precompute

__all__ = ["mols", "smiles", "mol_blocks"]

logger = logging.getLogger(__name__)

#
# A set of fragments was precomputed from the ChEMBL database
# Each fragment is annotated with its number of rotatable bonds
# Load this set of fragments
#
subset_size = 30_000

try:
    # try to load from the pickle file first
    # this might fail, because of wrong version of pickle, etc.
    # but it is faster than the backup solution (see except block)
    with pkg_resources.resource_stream(__name__, f"fragments_{subset_size}.pkl") as f:
        precomputed = pickle.load(f)

except:
    logger.info("Loading from pickle file failed. Parsing smi file instead.")
    # if loading from the pickle file fails, parse the smi file
    # this is slower than loading from the pickle file
    with pkg_resources.resource_stream(__name__, f"fragments_{subset_size}.smi") as f:
        fragments = [MolFromSmiles(line) for line in f]
        precomputed = precompute(fragments)

    # TODO: save precomputed data to resources(!) to speed up future runs
    # file_name_pkl = f"fragments_{subset_size}.pkl"
    # with open(file_name_pkl, "wb") as f:
    #     pickle.dump(precomputed, f)

fragments = precomputed["fragments"]
possible_dummy_labels = precomputed["possible_dummy_labels"]
compatible_fragments = precomputed["compatible_fragments"]
compatible_reactions = precomputed["compatible_reactions"]

dummy_pattern = MolFromSmiles("[*]")


@st.composite
def mols(
    draw,
    name=st.text(alphabet=list(string.ascii_lowercase), min_size=1),
    max_rotatable_bonds=st.just(10),
    n_connected_components=st.just(1),
):
    """Strategy for generating random molecules.

    Parameters
    ----------
    name : str or hypothesis.strategies.SearchStrategy
        Name of the molecule.
    max_rotatable_bonds : int or hypothesis.strategies.SearchStrategy
        Maximum number of rotatable bonds in the molecule. If
        n_connected_components > 1, each individual component will have less rotatable
        bonds than max_rotatable_bonds.
    n_connected_components : int or hypothesis.strategies.SearchStrategy
        Number of connected components in the molecule.
    """
    if isinstance(name, str):
        name = st.just(name)
    if isinstance(n_connected_components, int):
        n_connected_components = st.just(n_connected_components)
    if isinstance(max_rotatable_bonds, int):
        max_rotatable_bonds = st.just(max_rotatable_bonds)

    max_rotatable_bonds = draw(max_rotatable_bonds)
    n_connected_components = draw(n_connected_components)

    assert n_connected_components > 0, "n_connected_components must be > 0"
    assert max_rotatable_bonds >= 0, "max_rotatable_bonds must be >= 0"

    # draw connected components independently
    if n_connected_components > 1:
        components = [
            draw(
                mols(
                    max_rotatable_bonds=max_rotatable_bonds,
                    n_connected_components=1,
                )
            )
            for _ in range(n_connected_components)
        ]
        seed = reduce(CombineMols, components)
    else:
        # repeat this process until a valid molecule is drawn
        while True:
            # start by drawing a random fragment and use it as seed
            # make sure that the seed has less rotatable bonds than max_rotatable_bonds
            # TODO: expensive to compute this every time
            suitable_seeds = [
                (f, n) for (f, n) in fragments if n <= max_rotatable_bonds
            ]
            seed, n_rotatable_bonds_current = draw(st.sampled_from(suitable_seeds))

            # extend the seed by running compatible reactions on the seed
            # stop when molecule has no dummy atoms anymore (i.e. is fully built)
            while seed.HasSubstructMatch(dummy_pattern):
                while True:
                    # find a free substitution site
                    for atom in seed.GetAtoms():
                        if atom.GetSymbol() == "*":
                            dummy_label = atom.GetIsotope()
                            break

                    # draw a random reaction compatible with the chosen dummy atom
                    reaction, right = draw(
                        st.sampled_from(compatible_reactions[dummy_label])
                    )
                    other_label = (
                        reaction._matchers[int(right)].GetAtomWithIdx(0).GetIsotope()
                    )

                    # compute how many rotatable bonds are left to add
                    n_remaining_rotatable_bonds = (
                        max_rotatable_bonds - n_rotatable_bonds_current
                    )

                    # special case: some reactions connect fragments with a double bond
                    # --> no rotatable bond is added in the reaction
                    # --> add one from the remaining number of rotatable bonds
                    reaction_creates_nonrotatable_bond = (
                        reaction.GetProducts()[0].GetBondWithIdx(0).GetBondType()
                        != BondType.SINGLE
                    )
                    n_remaining_rotatable_bonds += reaction_creates_nonrotatable_bond

                    # draw a random fragment compatible with the drawn reaction, but
                    # filter all fragments that would exceed the maximum number of
                    # rotatable bonds
                    # note: compatible_mapping[dummy_label] was sorted by number of
                    #       dummy atoms and size --> shrinking is done automatically
                    # TODO: expensive to compute this every time
                    feasible_fragments = [
                        (f, n)
                        for (f, n) in compatible_fragments[other_label]
                        if n <= n_remaining_rotatable_bonds
                    ]

                    fragment, n = draw(st.sampled_from(feasible_fragments))
                    n -= reaction_creates_nonrotatable_bond

                    # run the reaction and collect the products
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

                    # if the reaction failed, try again (with another reaction)
                    if len(results) == 0:
                        continue

                    # if the reaction was successful, pick a random product and use it
                    # as the new seed
                    # note: shrinking is not applicable here, because all products have
                    # the same size
                    seed = draw(st.sampled_from(results))

                    # update the number of rotatable bonds
                    n_rotatable_bonds_current += n

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

    # TODO: create conformer

    return seed


@st.composite
def smiles(draw, **kwargs):
    mol = draw(mols(**kwargs))

    return f"{MolToSmiles(mol)} {mol.GetProp('_Name')}"


@st.composite
def mol_blocks(draw, **kwargs):
    mol = draw(mols(**kwargs))

    return MolToMolBlock(mol)
