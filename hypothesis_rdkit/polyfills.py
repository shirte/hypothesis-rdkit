from rdkit import rdBase

__all__ = ["BlockLogs"]

# polyfill BlockLogs (is valid contextmanager only since RDKit 2022)
if not hasattr(rdBase, "BlockLogs") or "__enter__" not in rdBase.BlockLogs.__dict__:
    from contextlib import contextmanager

    @contextmanager
    def nullcontext(enter_result=None):
        yield enter_result

    BlockLogs = nullcontext
else:
    BlockLogs = rdBase.BlockLogs
