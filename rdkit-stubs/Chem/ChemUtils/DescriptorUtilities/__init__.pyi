"""
Collection of utilities to be used with descriptors

"""
from __future__ import annotations
import rdkit.Chem.ChemUtils.DescriptorUtilities
import typing
import math

__all__ = [
    "VectorDescriptorNamespace",
    "VectorDescriptorWrapper",
    "math",
    "setDescriptorVersion"
]


class VectorDescriptorNamespace(dict):
    pass
class VectorDescriptorWrapper():
    """
    Wrap a function that returns a vector and make it seem like there
        is one function for each entry.  These functions are added to the global
        namespace with the names provided
    """
    pass
