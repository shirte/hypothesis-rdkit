import logging
import string
from functools import reduce

import hypothesis.strategies as st
from rdkit.Chem import (
    AddHs,
    AssignStereochemistry,
    BondType,
    CombineMols,
    GetMolFrags,
    MolFromInchi,
    MolFromMolBlock,
    MolFromSmiles,
    MolToInchi,
    MolToMolBlock,
    MolToSmiles,
    RemoveHs,
    SanitizeMol,
)
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

from .polyfills import BlockLogs
from .storage import load_fragment_data

__all__ = ["mols", "smiles", "mol_blocks", "inchis", "representations"]

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

dummy_pattern = MolFromSmiles("[*]")


@st.composite
def mols(
    draw,
    name=st.text(alphabet=list(string.ascii_lowercase), min_size=1),
    max_rotatable_bonds=st.just(10),
    n_connected_components=st.just(1),
    max_conformers=st.just(0),
    fragment_library_size="small",
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
    max_conformers : int or hypothesis.strategies.SearchStrategy
        Maximum number of conformers to generate. Be aware that generating conformers
        is slow.
    fragment_library_size : str or int
        Size of the fragment library to use. Can be 'small' (=30000), 'large' (=100000)
        or an integer specifying the number of fragments to use. If an integer is
        specified, the fragment library will be generated on the fly, which is slower.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        Random molecule.
    """
    if isinstance(name, str):
        name = st.just(name)
    if isinstance(n_connected_components, int):
        n_connected_components = st.just(n_connected_components)
    if isinstance(max_rotatable_bonds, int):
        max_rotatable_bonds = st.just(max_rotatable_bonds)
    if isinstance(max_conformers, int):
        max_conformers = st.just(max_conformers)
    if isinstance(fragment_library_size, str):
        assert fragment_library_size in [
            "small",
            "large",
            "full",
        ], "fragment_library_size must be 'small', 'large', 'full'"
        if fragment_library_size == "small":
            fragment_library_size = 30000
        elif fragment_library_size == "large":
            fragment_library_size = 100000
        elif fragment_library_size == "full":
            fragment_library_size = 365069

    max_rotatable_bonds = draw(max_rotatable_bonds)
    n_connected_components = draw(n_connected_components)
    max_conformers = draw(max_conformers)

    assert n_connected_components > 0, "n_connected_components must be > 0"
    assert max_rotatable_bonds >= 0, "max_rotatable_bonds must be >= 0"
    assert max_conformers >= 0, "max_conformers must be >= 0"
    assert isinstance(
        fragment_library_size, int
    ), "fragment_library_size must be an integer"

    # load precomputed fragment data
    fragment_data = load_fragment_data(fragment_library_size)

    # draw connected components independently
    if n_connected_components > 1:
        components = [
            draw(
                mols(
                    max_rotatable_bonds=max_rotatable_bonds,
                    n_connected_components=1,
                    max_conformers=0,
                )
            )
            for _ in range(n_connected_components)
        ]
        seed = reduce(CombineMols, components)
    else:
        # turn off RDKit logging
        with BlockLogs():
            # repeat this process until a valid molecule is drawn
            while True:
                # start by drawing a random fragment and use it as seed
                # make sure that the seed has less rotatable bonds than max_rotatable_bonds
                valid_seed_indices = fragment_data.get_fragment_indices(
                    max_rotatable_bonds=max_rotatable_bonds
                )
                seed_index = draw(st.sampled_from(valid_seed_indices))
                seed = fragment_data.get_fragment(seed_index)
                n_rotatable_bonds_current = fragment_data.get_rotatable_bonds(
                    seed_index
                )

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
                            st.sampled_from(
                                fragment_data.get_compatible_reactions(dummy_label)
                            )
                        )
                        other_label = (
                            reaction._matchers[int(right)]
                            .GetAtomWithIdx(0)
                            .GetIsotope()
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
                        n_remaining_rotatable_bonds += (
                            reaction_creates_nonrotatable_bond
                        )

                        # draw a random fragment compatible with the drawn reaction
                        # note: the result of get_fragment_indices was sorted by number of
                        #       dummy atoms and size --> shrinking is done automatically
                        feasible_fragment_indices = fragment_data.get_fragment_indices(
                            n_remaining_rotatable_bonds, other_label
                        )

                        fragment_index = draw(
                            st.sampled_from(feasible_fragment_indices)
                        )
                        fragment = fragment_data.get_fragment(fragment_index)
                        n = fragment_data.get_rotatable_bonds(fragment_index)
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
                    AssignStereochemistry(seed, force=True, cleanIt=True)
                    seed = RemoveHs(seed)

                    # AssignStereochemistry might add rotatable bonds
                    # --> check if the number of rotatable bonds is still valid
                    if CalcNumRotatableBonds(seed) > max_rotatable_bonds:
                        raise ValueError("too many rotatable bonds")

                    break
                except:
                    continue

    name = draw(name)
    if name is not None:
        seed.SetProp("_Name", name)

    # generate conformers
    if max_conformers > 0 and seed.GetNumAtoms() > 1:
        # turn off RDKit logging
        with BlockLogs():
            # we have to add hydrogens
            seed = AddHs(seed)

            random_seed = draw(st.integers(min_value=0, max_value=2**16 - 1))
            try:
                EmbedMultipleConfs(
                    seed,
                    numConfs=max_conformers,
                    maxAttempts=max_conformers * 2,
                    randomSeed=random_seed,
                )
            except:
                # doesn't work for some molecules
                pass

    return seed


@st.composite
def representations(draw, format, **kwargs):
    with BlockLogs():
        while True:
            try:
                mol = draw(mols(**kwargs))

                assert format in ["smiles", "mol_block", "inchi"]

                if format == "smiles":
                    result = f"{MolToSmiles(mol)} {mol.GetProp('_Name')}"
                    reconstructed_mol = MolFromSmiles(result)
                elif format == "mol_block":
                    result = MolToMolBlock(mol)
                    reconstructed_mol = MolFromMolBlock(result)
                elif format == "inchi":
                    result = MolToInchi(mol)
                    reconstructed_mol = MolFromInchi(result)

                # check that the reconstructed molecule does not have more rotatable bonds
                # or connected components than the original molecule (this can happen
                # during InChI conversion, i.e. when metals are disconnected)
                old_rotatable_bonds = CalcNumRotatableBonds(mol)
                new_rotatable_bonds = CalcNumRotatableBonds(reconstructed_mol)
                assert new_rotatable_bonds <= old_rotatable_bonds, (
                    f"reconstructed molecule has more rotatable bonds "
                    f"({new_rotatable_bonds}) than original molecule "
                    f"({old_rotatable_bonds})"
                )

                new_n_connected_components = len(
                    GetMolFrags(reconstructed_mol, asMols=True)
                )
                old_n_connected_components = len(GetMolFrags(mol, asMols=True))
                assert new_n_connected_components <= old_n_connected_components, (
                    f"reconstructed molecule has more connected components "
                    f"({new_n_connected_components}) than original molecule "
                    f"({old_n_connected_components})"
                )

                break
            except Exception as e:
                logger.debug(e)
                continue

        return result


@st.composite
def smiles(draw, **kwargs):
    return draw(representations(format="smiles", **kwargs))


@st.composite
def mol_blocks(draw, **kwargs):
    return draw(representations(format="mol_block", **kwargs))


@st.composite
def inchis(draw, **kwargs):
    return draw(representations(format="inchi", **kwargs))
