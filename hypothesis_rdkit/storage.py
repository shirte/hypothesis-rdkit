import gzip
import logging
import os
import pickle
import re
import sys
from collections import defaultdict
from functools import lru_cache

import numpy as np
from platformdirs import user_data_dir
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.BRICS import reverseReactions
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, ReactionToSmarts
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

from .fragment import generate_fragments
from .util import get_rotatable_bonds, num_dummy_atoms
from .version import __version__

# import importlib.resources or fall back to the backport
if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files

__all__ = ["load_fragment_data"]

logger = logging.getLogger(__name__)


class FragmentData:
    def __init__(self, fragments):
        # compute the number of rotatable bonds for each fragment
        self.fragments = fragments
        self.reverse_reactions = reverseReactions
        self.rotatable_bonds = [get_rotatable_bonds(f) for f in fragments]
        self.compatible_reactions = self._get_compatible_reactions_mapping()
        self.compatible_fragments = self._get_compatible_fragments_mapping()

    def _get_compatible_reactions_mapping(self):
        # All reactions have the form [*:1] + [*:2] --> [*:1][*:2]. We would like to
        # quickly get all reactions that match a given dummy atom. Also, we want to know
        # if the dummy atom is in the first or the second reactant and which reactant
        # remains free to choose.
        # For that reason, we build a mapping of the form
        #   dummy_label --> [(reaction, right_side_free)]
        # where
        #   * dummy_label is the BRICS reaction label (e.g. 1, 2, 3, ...)
        #   * reaction is a BRICS reaction (e.g. [*:1] + [*:2] --> [*:1][*:2])
        #   * right_side_free is True if the dummy atom is in the second reactant
        compatible_reactions = defaultdict(list)
        for r in reverseReactions:
            for i, right_side_free in zip([0, 1], [True, False]):
                for a in r._matchers[i].GetAtoms():  # should be only one dummy atom
                    if a.GetSymbol() == "*":
                        dummy_label = a.GetIsotope()
                        compatible_reactions[dummy_label].append((r, right_side_free))

        return compatible_reactions

    def _get_compatible_fragments_mapping(self):
        compatible_fragments = dict()

        # mapping dummy atom d to fragments containing d
        for i, fragment, rotatable_bonds in zip(
            range(len(self.fragments)), self.fragments, self.rotatable_bonds
        ):
            dummy_labels = set(
                [a.GetIsotope() for a in fragment.GetAtoms() if a.GetSymbol() == "*"]
            )
            for dummy_label in dummy_labels:
                if dummy_label not in compatible_fragments:
                    compatible_fragments[dummy_label] = defaultdict(list)

                compatible_fragments[dummy_label][rotatable_bonds].append(i)

        # At the moment, each fragment is contained in the list for the corresponding
        # number of rotatable bonds. At the end, the lists should contain all fragments
        # with *less* or equal number of rotatable bonds. Therefore, we need to
        # accumulate the lists:
        for dummy_label in compatible_fragments.keys():
            max_rotatable_bonds = max(compatible_fragments[dummy_label].keys())
            stash = []
            for rotatable_bonds in range(max_rotatable_bonds + 1):
                stash += compatible_fragments[dummy_label][rotatable_bonds]
                compatible_fragments[dummy_label][rotatable_bonds] = list(set(stash))

        # Dummy label "-1" is a special label that matches all dummy labels.
        possible_rotatable_bonds = set(
            sum(
                [
                    list(compatible_fragments[dummy_label].keys())
                    for dummy_label in compatible_fragments.keys()
                ],
                [],
            )
        )
        compatible_fragments[-1] = defaultdict(list)
        for rotatable_bonds in possible_rotatable_bonds:
            compatible_fragments[-1][rotatable_bonds] = list(
                set(
                    sum(
                        [
                            compatible_fragments[dummy_label][rotatable_bonds]
                            for dummy_label in compatible_fragments.keys()
                            if rotatable_bonds in compatible_fragments[dummy_label]
                        ],
                        [],
                    )
                )
            )

        # precompute a ranking key for each fragment
        # --> the next sorting calls will be faster
        ranking_keys = [(num_dummy_atoms(f), f.GetNumAtoms()) for f in self.fragments]

        # sort the lists in each key by number of dummy atoms in the fragment and
        # fragment size --> shrinking
        for dummy_label in compatible_fragments.keys():
            for rotatable_bonds in compatible_fragments[dummy_label].keys():
                compatible_fragments[dummy_label][rotatable_bonds] = sorted(
                    compatible_fragments[dummy_label][rotatable_bonds],
                    key=lambda x: ranking_keys[x],
                )

        return compatible_fragments

    def get_compatible_reactions(self, dummy_label):
        return self.compatible_reactions[dummy_label]

    def get_fragment_indices(self, max_rotatable_bonds, dummy_label=None):
        assert max_rotatable_bonds >= 0
        assert dummy_label is None or dummy_label in self.compatible_fragments

        if dummy_label is None:
            dummy_label = -1

        lookup_table = self.compatible_fragments[dummy_label]
        max_value = max(lookup_table.keys())
        if max_rotatable_bonds > max_value:
            max_rotatable_bonds = max_value

        return self.compatible_fragments[dummy_label][max_rotatable_bonds]

    def get_fragment(self, idx):
        return self.fragments[idx]

    def get_rotatable_bonds(self, idx):
        # TODO: n_rotatable_bonds_current is overestimating the number of rotatable
        #       bonds, because the first reaction bond is already included
        return self.rotatable_bonds[idx]

    def __getstate__(self):
        # ChemicalReaction is not backwards / upwards compatible
        # --> store the reaction as SMARTS
        # --> replace all occurences of reactions with their SMARTS representation
        #     or indices
        compatible_reactions = defaultdict(list)
        for dummy_label in self.compatible_reactions.keys():
            old_values = self.compatible_reactions[dummy_label]
            compatible_reactions[dummy_label] = [
                (self.reverse_reactions.index(reaction), right_side_free)
                for (reaction, right_side_free) in old_values
            ]

        return {
            "fragments": self.fragments,
            "reverse_reactions": [ReactionToSmarts(r) for r in self.reverse_reactions],
            "matchers": [r._matchers for r in self.reverse_reactions],
            "rotatable_bonds": self.rotatable_bonds,
            "compatible_reactions": compatible_reactions,
            "compatible_fragments": self.compatible_fragments,
        }

    def __setstate__(self, state):
        self.fragments = state["fragments"]
        self.reverse_reactions = [
            ReactionFromSmarts(r) for r in state["reverse_reactions"]
        ]
        self.rotatable_bonds = state["rotatable_bonds"]
        self.compatible_reactions = state["compatible_reactions"]
        self.compatible_fragments = state["compatible_fragments"]

        for r, matchers in zip(self.reverse_reactions, state["matchers"]):
            r._matchers = matchers

        for dummy_label in self.compatible_reactions.keys():
            old_values = self.compatible_reactions[dummy_label]
            self.compatible_reactions[dummy_label] = [
                (self.reverse_reactions[rid], right_side_free)
                for (rid, right_side_free) in old_values
            ]


def sample_subset(fragments, subset_size, temperature=10):
    assert len(fragments) >= subset_size

    # make sure that all fragments [H][*:d] are contained in the fragment subset
    mandatory_fragment_indices = [
        i
        for i, fragment in enumerate(fragments)
        if set([a.GetSymbol() for a in fragment.GetAtoms()]) == set(["*", "H"])
    ]
    mandatory_fragments = [fragments[i] for i in mandatory_fragment_indices]

    # remove these fragments from the list of fragments
    fragments = [
        fragment
        for i, fragment in enumerate(fragments)
        if i not in mandatory_fragment_indices
    ]

    # sample fragments using a boltzmann distribution on the molecular weights
    # compute the weight of each fragment
    weights = np.array([CalcExactMolWt(fragment) for fragment in fragments])

    # probability is exp(-weight / temperature) / Z
    probabilities = np.exp(-weights / temperature)
    probabilities = probabilities / probabilities.sum()

    subset_indices = np.random.choice(
        len(fragments),
        size=subset_size - len(mandatory_fragments),
        replace=False,
        p=probabilities,
    )

    fragment_subset = mandatory_fragments + [
        fragment for i, fragment in enumerate(fragments) if i in subset_indices
    ]

    # sort fragments by number of dummy atoms and size (--> shrinking)
    fragment_subset = sorted(
        fragment_subset, key=lambda x: (num_dummy_atoms(x), x.GetNumAtoms())
    )

    return fragment_subset


def load_fragment_data_from_spec(subset_spec):
    subset_size, source, extension = subset_spec

    assert source in ["resources", "user_dir"]
    assert extension in ["pkl.gz", "smi"]

    filename = f"fragments_{subset_size}.{extension}"
    mode = "rb" if extension == "pkl.gz" else "r"

    if source == "resources":
        handler = (
            lambda: files(__package__)
            .joinpath("resources")
            .joinpath(filename)
            .open(mode)
        )
    elif source == "user_dir":
        user_dir = user_data_dir(appname=__package__, version=__version__)
        path_to_file = os.path.join(user_dir, filename)
        handler = lambda: open(path_to_file, mode)

    with handler() as f:
        if extension == "pkl.gz":
            with gzip.open(f, "rb") as g:
                return pickle.load(g)
        elif extension == "smi":
            fragments = [MolFromSmiles(line) for line in f]

            logger.info("Prepare fragment data...")
            precomputed = FragmentData(fragments)
            save_precomputed_fragments_to_user_dir(precomputed)

            return precomputed


def save_precomputed_fragments_to_user_dir(precomputed):
    fragments = precomputed.fragments
    subset_size = len(fragments)

    # save the precomputed data to a pickle file
    user_dir = user_data_dir(appname=__package__, version=__version__)
    file_name_pkl = f"fragments_{subset_size}.pkl.gz"
    path_pkl = os.path.join(user_dir, file_name_pkl)

    os.makedirs(user_dir, exist_ok=True)

    with gzip.open(path_pkl, "wb") as f:
        # use protocol 3 to be compatible with python < 3.8
        pickle.dump(precomputed, f, protocol=3)

    # write the fragments as smiles into a plain file (as backup if the pickle file
    # cannot be loaded)
    file_name_smi = f"fragments_{subset_size}.smi"
    path_smi = os.path.join(user_dir, file_name_smi)
    with open(path_smi, "w") as f:
        for fragment in fragments:
            f.write(f"{MolToSmiles(fragment)}\n")


def iter_fragment_subsets(requested_subset_size):
    # try precomputed files with the requested subset size
    yield (requested_subset_size, "resources", "pkl.gz")
    yield (requested_subset_size, "user_dir", "pkl.gz")

    # try raw fragment list files with the requested subset size
    yield (requested_subset_size, "resources", "smi")
    yield (requested_subset_size, "user_dir", "smi")

    # helper function to parse the subset size from a file name
    # e.g. fragments_30000.pkl.gz --> 30000
    def parse_subset_size(filename, source):
        # check if the file name matches the pattern
        match = re.match(
            r"fragments_(?P<subset_size>\d+).(?P<extension>smi|pkl\.gz)$", filename
        )
        if match:
            # get the subset size from the file name
            subset_size = int(match.group("subset_size"))
            extension = match.group("extension")

            return (subset_size, source, extension)
        else:
            return None

    # search for larger subsets in package resources
    subset_specs_resources = [
        parse_subset_size(path.name, "resources")
        for path in files(__package__).joinpath("resources").iterdir()
    ]

    # search for larger subsets in user directory
    user_dir = user_data_dir(appname=__package__, version=__version__)

    if not os.path.exists(user_dir):
        return

    subset_specs_user_dir = [
        parse_subset_size(path.name, "user_dir")
        for path in os.scandir(user_dir)
        if path.is_file()
    ]

    # concat
    subset_specs = subset_specs_resources + subset_specs_user_dir

    # filter out None values
    subset_specs = [s for s in subset_specs if s is not None]

    # sort by subset size and prefer pickle files over smiles files
    subset_specs = sorted(
        subset_specs, key=lambda x: x[0] * 2 + (x[2] == "pkl.gz"), reverse=True
    )

    for spec in subset_specs:
        if spec[0] > requested_subset_size:
            yield spec


@lru_cache(maxsize=1)
def load_fragment_data(subset_size=30_000):
    """Load precomputed fragment data or generate it if necessary.

    Parameters
    ----------
    subset_size : int
        Size of the fragment subset to load.

    Returns
    -------
    dict
        Dictionary containing precomputed fragment data.
    """

    # iterate through available fragment files
    for subset_spec in iter_fragment_subsets(subset_size):
        subset_size_, source, extension = subset_spec
        try:
            logger.info(
                f"Loading data containing {subset_size_} fragments from {extension} "
                f"file in {source}..."
            )
            precomputed = load_fragment_data_from_spec(subset_spec)
        except FileNotFoundError:
            logger.info(f"Could not find fragment data in {source}")
            continue
        except Exception as e:
            import traceback

            logger.warning(f"Could not load fragment data: {e}")
            traceback.print_exc()
            continue

        # sample if subset_size is too large
        fragments = precomputed.fragments
        if len(fragments) > subset_size:
            logger.info(f"Sampling subset of size {subset_size}")
            fragment_subset = sample_subset(fragments, subset_size, temperature=10)
            precomputed = FragmentData(fragment_subset)
            save_precomputed_fragments_to_user_dir(precomputed)

        return precomputed

    # Compute the fragments from scratch (slowest option, takes hours).
    logger.info("Generating fragments from scratch...")
    fragments = generate_fragments()
    precomputed = FragmentData(fragments)
    save_precomputed_fragments_to_user_dir(precomputed)

    # Try again. This should work now since we computed a larger set of fragments than
    # required. The next call will sample a subset from it (see step above).
    return load_fragment_data(subset_size)
