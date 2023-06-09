from collections import defaultdict

from rdkit.Chem.BRICS import reverseReactions
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds


def get_possible_dummy_labels(fragments):
    possible_dummy_labels = set(
        [
            a.GetIsotope()
            for fragment in fragments
            for a in fragment.GetAtoms()
            if a.GetSymbol() == "*"
        ]
    )
    return possible_dummy_labels


# compute the number of rotatable bonds for each fragment
def get_rotatable_bonds(fragment):
    # for performance reasons, we precompute the number of rotatable bonds
    # --> each fragments is annotated with its number of rotatable bonds
    if fragment.GetNumAtoms() <= 2:
        # if the fragment has only one heavy atom (e.g. *-CH3), it will lead to a
        # product with no additional rotatable bonds (the bond created by the
        # reaction is not rotatable)
        rotatable_bonds = 0
    else:
        # otherwise, calculate the number of rotatable bonds and add one to account
        # for the bond created by the reaction
        # note: this is an approximation for two reasons:
        # 1. the bond created by the reaction might be a double bond (not rotatable)
        # 2. the bond created by the reaction might not be rotatable because
        #    this or the other fragment is symmetric (e.g. CC(Cl)(Cl)(Cl))
        rotatable_bonds = CalcNumRotatableBonds(fragment) + 1
    return rotatable_bonds


def num_dummy_atoms(mol):
    return sum([1 for a in mol.GetAtoms() if a.GetSymbol() == "*"])


def get_compatible_reactions_mapping():
    #
    # for each dummy label, collect all reactions that contain a dummy atom with that label
    # if the dummy atom is in the first reactant, the second reactant is free to choose
    # additionally annotate, whether the the right_side is free to choose
    # obtain a map compatible_reactions: dummy_label --> [(reaction, right_side_free)]
    #
    compatible_reactions = defaultdict(list)
    for r in reverseReactions:
        for i, right_side_free in zip([0, 1], [True, False]):
            for a in r._matchers[i].GetAtoms():  # should be only one dummy atom
                if a.GetSymbol() == "*":
                    dummy_label = a.GetIsotope()
                    compatible_reactions[dummy_label].append((r, right_side_free))

    return compatible_reactions


def get_compatible_fragments_mapping(fragments):
    # mapping dummy atom d to fragments containing d
    compatible_fragments = defaultdict(list)
    for fragment, n in fragments:
        for a in fragment.GetAtoms():
            if a.GetSymbol() == "*":
                dummy_label = a.GetIsotope()
                compatible_fragments[dummy_label].append((fragment, n))

    # sort the lists in each key by number of dummy atoms in the fragment and
    # fragment size --> shrinking
    for dummy_label in compatible_fragments.keys():
        compatible_fragments[dummy_label] = sorted(
            compatible_fragments[dummy_label],
            key=lambda x: (num_dummy_atoms(x[0]), x[0].GetNumAtoms()),
        )

    return compatible_fragments


def precompute(fragments):
    possible_dummy_labels = get_possible_dummy_labels(fragments)

    # compute the number of rotatable bonds for each fragment
    fragments = [(f, get_rotatable_bonds(f)) for f in fragments]

    compatible_reactions = get_compatible_reactions_mapping()
    compatible_fragments = get_compatible_fragments_mapping(fragments)

    return dict(
        fragments=fragments,
        possible_dummy_labels=possible_dummy_labels,
        compatible_reactions=compatible_reactions,
        compatible_fragments=compatible_fragments,
    )
