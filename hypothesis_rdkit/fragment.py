import gzip
import logging
import os
import pickle
import signal
from contextlib import contextmanager
from functools import partial
from multiprocessing import Pool, cpu_count
from urllib.request import urlretrieve

import numpy as np
from rdkit.Chem import (
    BondType,
    ForwardSDMolSupplier,
    GetMolFrags,
    MolFromSmiles,
    MolToSmiles,
)
from rdkit.Chem.BRICS import BRICSDecompose, reverseReactions
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from tqdm import tqdm

from .util import (
    get_possible_dummy_labels,
    get_rotatable_bonds,
    num_dummy_atoms,
    precompute,
)

__all__ = ["get_fragments", "main"]

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


def get_fragments(
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
        Maximum number of molecules in ChEMBL to decompose. The default is 100_000_000_000.

    n_processes : int, optional
        Number of processes to use for decomposition. The default is 1.5 times
        cpu_count().

    batch_size : int, optional
        Number of molecules to process until a consumer process is terminated. BRICS
        has memory leak problems, and to avoid accumulating occupied memory, consumer
        processes are killed after processing a batch of molecules. The default number
        of molecules in a batch is 100.

    timeout : int, optional
        Timeout in seconds for the decomposition of a single molecule. The default is 5 * 60.

    max_weight : int, optional
        Maximum weight of a molecule to decompose. The default is 3000.
    """
    fragments_cache = f"fragments.pkl"
    if not os.path.exists(fragments_cache):
        molecules = producer(chembl_version, n_subset, max_weight)
        batches = batch(tqdm(molecules), batch_size)

        # note: we use maxtasksperchild=1 to keep the memory consumption low
        # BRICS has memory leak problems and the used memory accumulates if processes are
        # recycled --> kill processes after each batch
        pool = Pool(n_processes, maxtasksperchild=1)
        fragments = set()
        try:
            for result in pool.imap_unordered(
                partial(consumer, timeout=timeout), batches
            ):
                fragments.update(result)
        finally:
            pool.close()
            pool.join()

        fragments = [s for f in fragments if (s := MolFromSmiles(f)) is not None]

        with open(fragments_cache, "wb") as f:
            pickle.dump(fragments, f)
    else:
        fragments = pickle.load(open(fragments_cache, "rb"))

    # add a fragment [*:d][H] for each *possible* dummy label d if
    # * the minimum number of rotatable bonds is > 0 and
    # * there is no reverse reaction connecting fragments with anything other than a
    #   single bond (we cannot connect a fragment with hydrogen using e.g. a double bond)
    possible_dummy_labels = get_possible_dummy_labels(fragments)

    # for each dummy label, find the fragment that creates the fewest rotatable bonds
    min_rotatable_bonds_per_dummy_label = {
        dummy_label: np.inf for dummy_label in possible_dummy_labels
    }

    for fragment in fragments:
        num_rotatable_bonds = get_rotatable_bonds(fragment)
        for a in fragment.GetAtoms():
            if a.GetSymbol() == "*":
                dummy_label = a.GetIsotope()
                min_rotatable_bonds_per_dummy_label[dummy_label] = min(
                    min_rotatable_bonds_per_dummy_label[dummy_label],
                    num_rotatable_bonds,
                )

    # find dummy labels that are connected via something that is not a single bond
    non_single_bond_dummy_atoms = set()
    for reaction in reverseReactions:
        for p in reaction.GetProducts():
            s = MolToSmiles(p)
            if p.GetBondWithIdx(0).GetBondType() != BondType.SINGLE:
                for matcher in reaction._matchers:
                    non_single_bond_dummy_atoms.add(
                        matcher.GetAtomWithIdx(0).GetIsotope()
                    )

    dummy_fragments = [
        MolFromSmiles(f"[{d}*][H]")
        for d in possible_dummy_labels
        if min_rotatable_bonds_per_dummy_label[d] > 0
        and d not in non_single_bond_dummy_atoms
    ]
    fragments = fragments + dummy_fragments

    # sort fragments by number of dummy atoms and size (--> shrinking)
    fragments = sorted(fragments, key=lambda x: (num_dummy_atoms(x), x.GetNumAtoms()))

    return fragments


def sample_subset(fragments, subset_size, temperature=10):
    # sample subsets of fragments using a boltzmann distribution on the weights

    # compute the weight for each fragment
    weights = np.array([CalcExactMolWt(fragment) for fragment in fragments])

    probabilities = np.exp(-weights / temperature)
    probabilities = probabilities / probabilities.sum()

    subset_indices = np.random.choice(
        len(fragments), size=subset_size, replace=False, p=probabilities
    )

    fragment_subset = [
        fragment for i, fragment in enumerate(fragments) if i in subset_indices
    ]

    return fragment_subset


def main():
    fragments = get_fragments()

    # sample subsets of fragments using a boltzmann distribution on the weights
    # compute the weight for each fragment
    subset_sizes = [30_000, 100_000, len(fragments)]
    temperature = 10

    for subset_size in subset_sizes:
        logger.info(f"Sampling subset of size {subset_size}")
        fragment_subset = sample_subset(fragments, subset_size, temperature)

        # precompute data
        logger.info("Precompute data")
        precomputed = precompute(fragment_subset)

        logger.info("Write subset to disk")
        file_name_pkl = f"fragments_{subset_size}.pkl"
        with open(file_name_pkl, "wb") as f:
            pickle.dump(precomputed, f)

        # write only the fragments as smiles into a plain file (as backup if the pickle
        # file cannot be loaded)
        file_name_smi = f"fragments_{subset_size}.smi"
        with open(file_name_smi, "w") as f:
            for fragment in fragment_subset:
                f.write(f"{MolToSmiles(fragment)}\n")


if __name__ == "__main__":
    main()
