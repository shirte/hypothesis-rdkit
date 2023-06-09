from hypothesis import given, settings
from rdkit.Chem import GetMolFrags, Mol, MolFromSmiles, rdMolDescriptors

from hypothesis_rdkit import mols

dummy_pattern = MolFromSmiles("[*]")


@given(...)
@settings(max_examples=1000)
def test_mols(mol: Mol):
    # mols have atoms
    assert len(mol.GetAtoms()) > 0

    # mols have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1


@given(mols(n_connected_components=3))
@settings(max_examples=1000)
def test_connected_components(mol: Mol):
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 3


@given(mols(max_rotatable_bonds=10))
@settings(max_examples=1000)
def test_rotatable_bonds(mol: Mol):
    assert rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
