import sys

if sys.version_info < (3, 10):
    from importlib_metadata import version
else:
    from importlib.metadata import version

__all__ = ["__version__"]

__version__ = version(__package__)
