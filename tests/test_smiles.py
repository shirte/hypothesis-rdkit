from hypothesis import given
from hypothesis_rdkit import smiles
from rdkit.Chem import MolFromSmiles, GetMolFrags


dummy_pattern = MolFromSmiles("[*]")


@given(smiles())
def test_smiles(smiles: str):
    # smiles are valid
    mol = MolFromSmiles(smiles)
    assert mol is not None

    # mols from smiles have atoms
    assert len(mol.GetAtoms()) > 0

    # mols from smiles have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1
