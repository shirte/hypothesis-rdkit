from hypothesis import given
from hypothesis_rdkit import mols
from rdkit.Chem import Mol, MolFromSmiles, GetMolFrags


dummy_pattern = MolFromSmiles("[*]")


@given(...)
def test_mols(mol: Mol):
    # mols have atoms
    assert len(mol.GetAtoms()) > 0

    # mols have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1


@given(mols(n_connected_components=3))
def test_connected_components(mol: Mol):
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 3
