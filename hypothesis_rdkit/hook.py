def _hypothesis_setup_hook():
    from rdkit.Chem import Mol
    import hypothesis.strategies as st
    from .strategy import mols

    st.register_type_strategy(Mol, mols(n_connected_components=1))
