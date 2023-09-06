from __future__ import annotations
import rdkit.ML.InfoTheory.BitClusterer
import typing
import rdkit.DataStructs
import rdkit.SimDivFilters.rdSimDivPickers

__all__ = [
    "BitClusterer",
    "DataStructs",
    "rdsimdiv"
]


class BitClusterer():
    """
    Class to cluster a set of bits based on their correllation

       The correlation matrix is first built using by reading the fingerprints
       from a database or a list of fingerprints
       
    """
    pass
