from hypothesis import HealthCheck, given, settings
from rdkit.Chem import GetMolFrags, Mol, MolFromSmiles, rdMolDescriptors

from hypothesis_rdkit import mols

dummy_pattern = MolFromSmiles("[*]")


@given(...)
@settings(max_examples=1000)
def test_mols(mol: Mol):
    # mols have atoms
    assert mol.GetNumAtoms() > 0

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


@given(mols(max_conformers=10))
@settings(
    max_examples=100,
    suppress_health_check=[
        HealthCheck.too_slow,
    ],
)
def test_conformers(mol: Mol):
    assert mol.GetNumConformers() <= 10


@given(mols(fragment_library_size="small"))
def test_small_fragmentation_library(mol):
    assert mol.GetNumAtoms() > 0


@given(mols(fragment_library_size="large"))
def test_large_fragmentation_library(mol):
    assert mol.GetNumAtoms() > 0


@given(mols(fragment_library_size="full"))
def test_full_fragmentation_library(mol):
    assert mol.GetNumAtoms() > 0
