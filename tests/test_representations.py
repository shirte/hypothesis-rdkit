from hypothesis import given
from rdkit.Chem import GetMolFrags, MolFromInchi, MolFromMolBlock, MolFromSmiles

from hypothesis_rdkit import inchis, mol_blocks, smiles

dummy_pattern = MolFromSmiles("[*]")


@given(smiles())
def test_smiles(smiles: str):
    # smiles are valid
    mol = MolFromSmiles(smiles)
    assert mol is not None

    # mols have atoms
    assert len(mol.GetAtoms()) > 0

    # mols have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1


@given(mol_blocks())
def test_mol_blocks(mol_block: str):
    # mol blocks are valid
    mol = MolFromMolBlock(mol_block)
    assert mol is not None

    # mols have atoms
    assert len(mol.GetAtoms()) > 0

    # mols have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1


@given(inchis())
def test_inchis(inchi: str):
    # inchis are valid
    mol = MolFromInchi(inchi)
    assert mol is not None

    # mols have atoms
    assert len(mol.GetAtoms()) > 0

    # mols have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1
