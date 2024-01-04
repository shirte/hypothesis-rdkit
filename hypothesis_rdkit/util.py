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
