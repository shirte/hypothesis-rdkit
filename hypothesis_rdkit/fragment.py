import gzip
import logging
import os
import signal
from contextlib import closing, contextmanager
from functools import partial
from multiprocessing import Pool, cpu_count
from urllib.request import urlretrieve

from rdkit.Chem import BondType, ForwardSDMolSupplier, GetMolFrags, MolFromSmiles
from rdkit.Chem.BRICS import BRICSDecompose, reverseReactions
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from tqdm import tqdm

from .util import get_possible_dummy_labels, num_dummy_atoms

__all__ = ["generate_fragments"]

logger = logging.getLogger(__name__)

# this script decomposes molecules from ChEMBL into fragments using BRICS
# Chembl 32 contains around 2.4 million molecules
# script takes approx 2.5 hours on 16 cores with 128 GB RAM

# idea:
# * iterate over molecules in ChEMBL data and decompose them into fragments
# * use a producer that reads molecules from disk
# * multiple consumers decompose the molecules into fragments
# * note 1: we use a timeout for the decomposition, because BRICS takes a long time for
#           some molecules (there are examples with a weight of 1600 Da)
# * note 2: BRICS leaks memory and so we need to terminate consumers after a batch of
#           molecules


def download_chembl(chembl_version):
    # TODO: download to user dir
    filename = f"chembl_{chembl_version:02d}.sdf.gz"

    # check if file already exists
    if os.path.exists(filename):
        return filename

    url = (
        f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
        f"chembl_{chembl_version:02d}/{filename}"
    )

    # download ChEMBL data
    path, status = urlretrieve(url, filename)

    return path


# a helper context manager to set a timeout
@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutError("Timed out!")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


# generate molecules by reading them from disk
def producer(chembl_version, n_subset, max_weight):
    # load the chembl database from the url
    path = download_chembl(chembl_version)

    with gzip.open(path, "rb") as f:
        supplier = ForwardSDMolSupplier(f)

        for mol, _ in zip(supplier, range(n_subset)):
            if mol is not None:
                # get connected components
                frags = GetMolFrags(mol, asMols=True)

                for frag in frags:
                    # remove isotope information
                    for a in frag.GetAtoms():
                        a.SetIsotope(0)

                    weight = CalcExactMolWt(frag)
                    if weight < max_weight:
                        yield frag


def consumer(batch, timeout):
    result = set()
    for mol in batch:
        try:
            with time_limit(timeout):
                # decompose mol into fragments using BRICS
                # we cancel a BRICS decomposition if it takes too long
                mol_fragments = BRICSDecompose(mol)
                result.update(mol_fragments)
        except TimeoutError:
            pass
        except:
            pass

    # TODO: check additional conditions here, especially normalize fragments with
    #       MolToSmiles(MolFromSmiles(...))

    return result


# produce batches from given stream of size batch_size
def batch(stream, batch_size):
    batch = []
    for mol in stream:
        batch.append(mol)
        if len(batch) == batch_size:
            yield batch
            batch = []

    if len(batch) > 0:
        yield batch


def generate_fragments(
    chembl_version=32,
    n_subset=100_000_000_000,
    n_processes=int(cpu_count() * 3 / 2),
    batch_size=100,
    timeout=5 * 60,
    max_weight=3000,
):
    """
    Decompose molecules from ChEMBL into fragments using BRICS. If this consumes too
    much memory, try to reduce the number of processes, the batch size, the timeout
    or the maximum weight of the molecules.

    Parameters
    ----------
    n_subset : int, optional
        Maximum number of molecules in ChEMBL to decompose. The default is
        100_000_000_000.

    n_processes : int, optional
        Number of processes to use for decomposition. The default is 1.5 times
        cpu_count().

    batch_size : int, optional
        Number of molecules to process until a consumer process is terminated. BRICS
        has memory leak problems, and to avoid accumulating occupied memory, consumer
        processes are killed after processing a batch of molecules. The default number
        of molecules in a batch is 100.

    timeout : int, optional
        Timeout in seconds for the decomposition of a single molecule. The default is
        5 * 60.

    max_weight : int, optional
        Maximum weight of a molecule to decompose. The default is 3000.

    Returns
    -------
    fragments : list of rdkit.Chem.rdchem.Mol
        List of fragments.
    """
    assert n_processes > 0, "n_processes must be > 0"

    logger.info(f"Generating fragments using {n_processes} processes...")

    molecules = producer(chembl_version, n_subset, max_weight)
    batches = batch(tqdm(molecules), batch_size)

    smiles_fragments = set()

    if n_processes == 1:
        # for debugging purposes
        for b in batches:
            smiles_fragments.update(consumer(b, timeout=timeout))
    else:
        # BRICS has memory leak problems and the used memory accumulates if
        # processes are recycled
        # --> kill processes after each batch
        # --> use maxtasksperchild=1
        with closing(Pool(n_processes, maxtasksperchild=1)) as pool:
            for result in pool.imap_unordered(
                partial(consumer, timeout=timeout), batches
            ):
                smiles_fragments.update(result)

    # check fragments
    fragments = []
    for f in smiles_fragments:
        # remove fragments that are not valid smiles
        mol = MolFromSmiles(f)
        if mol is None:
            continue

        fragments.append(mol)

    # Add a fragment [*:d][H] for each *possible* dummy label d *if* there is no reverse
    # reaction connecting fragments with anything other than a single bond (we cannot
    # connect a fragment with hydrogen using e.g. a double bond).
    possible_dummy_labels = get_possible_dummy_labels(fragments)

    non_single_bond_dummy_atoms = set()
    for reaction in reverseReactions:
        for p in reaction.GetProducts():
            if p.GetBondWithIdx(0).GetBondType() != BondType.SINGLE:
                for matcher in reaction._matchers:
                    non_single_bond_dummy_atoms.add(
                        matcher.GetAtomWithIdx(0).GetIsotope()
                    )

    dummy_fragments = [
        MolFromSmiles(f"[{d}*][H]")
        for d in possible_dummy_labels
        if d not in non_single_bond_dummy_atoms
    ]
    fragments = fragments + dummy_fragments

    # sort fragments by number of dummy atoms and size (--> shrinking)
    fragments = sorted(fragments, key=lambda x: (num_dummy_atoms(x), x.GetNumAtoms()))

    return fragments
