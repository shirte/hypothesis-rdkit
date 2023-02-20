from hypothesis import given
from hypothesis_rdkit import mol_blocks
from rdkit.Chem import MolFromMolBlock, MolFromSmiles, GetMolFrags


dummy_pattern = MolFromSmiles("[*]")


@given(mol_blocks())
def test_mol_blocks(mol_block: str):
    # mol blocks are valid
    mol = MolFromMolBlock(mol_block)
    assert mol is not None

    # mols from smiles have atoms
    assert len(mol.GetAtoms()) > 0

    # mols from smiles have no dummy atoms
    assert not mol.HasSubstructMatch(dummy_pattern)

    # by default, mols consist of one single connected component
    frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    assert len(frags) == 1
