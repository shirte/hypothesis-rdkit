def _hypothesis_setup_hook():
    import hypothesis.strategies as st
    from rdkit.Chem import Mol

    # Here, we would like to run the following code to register the mols strategy for
    # the rdkit.Chem.Mol type:
    #
    #   from .strategy import mols
    #   st.register_type_strategy(Mol, mols())
    #
    # However, this leads to a circular import when importing hypothesis_rdkit:
    # * import hypothesis_rdkit
    # * --> import hypothesis_rdkit/__init__.py
    # * --> import hypothesis_rdkit/strategy.py
    # * --> import hypothesis.strategies as st
    # * --> import hypothesis
    # * --> hypothesis runs all hooks from plugins (including this function)
    # * --> call hypothesis_rdkit.hook._hypothesis_setup_hook()  (this function)
    # * --> import hypothesis_rdkit/strategy.py (--> circular import!)
    # For this reason, we use a closure to delay the import of
    # hypothesis_rdkit.strategy.mols:

    def _strategy(t):
        from .strategy import mols

        if t is Mol:
            return mols()
        else:
            raise NotImplemented()

    st.register_type_strategy(Mol, _strategy)
