# hypothesis-rdkit

A strategy to generate random molecules for the hypothesis testing framework. It uses 
a collection of fragments generated from the ChEMBL database to construct plausible 
molecular graphs. The fragments were mined using the [BRICS method](
https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cmdc.200800178).

![demo](demo.png)


## Installation

```bash
pip install -U hypothesis-rdkit
# or
conda install -c conda-forge hypothesis-rdkit
```

## Usage

The module ```hypothesis-rdkit``` provides a strategy for generating RDKit 
molecules. During the installation of the package, this strategy is linked to the 
```rdkit.Chem.Mol``` type:


```python
from hypothesis import given
from rdkit.Chem import Mol

@given(...)
def test_molecule_method(mol : Mol):
    # mol is a randomly generated molecule
    assert mol.GetNumAtoms() > 0
```

You can use the ```mols``` strategy directly for further customization:

```python
from hypothesis import given
from hypothesis_rdkit import mols
from rdkit.Chem import GetMolFrags, Mol
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

@given(mols(n_connected_components=2, max_rotatable_bonds=5, n_conformers=10))
def test_molecule_mixtures(mol : Mol):
    frags = GetMolFrags(mol, asMols=True)
    assert len(frags) == 2

    for frag in frags:
        assert CalcNumRotatableBonds(frag) <= 5

    assert mol.GetNumConformers() <= 10
```

There are also strategies to generate molecules in SMILES, mol block or InChI 
representation accepting the same parameters as ```mols```:

```python
from hypothesis import given
from hypothesis_rdkit import smiles, mol_blocks, inchis
from rdkit.Chem import MolFromSmiles, MolFromMolBlock, MolFromInchi

@given(smiles())
def test_smiles(smiles : str):
    mol = MolFromSmiles(smiles)
    assert mol is not None and mol.GetNumAtoms() > 0

@given(mol_blocks())
def test_mol_block(mol_block : str):
    mol = MolFromMolBlock(mol_block)
    assert mol is not None and mol.GetNumAtoms() > 0

@given(inchis())
def test_inchi(inchi : str):
    mol = MolFromInchi(inchi)
    assert mol is not None and mol.GetNumAtoms() > 0
```


## Development

All fragment files are generated during a test run (```pytest```) in the user data 
directory. On Linux, this is ~/.local/share/hypothesis_rdkit/{version}/.