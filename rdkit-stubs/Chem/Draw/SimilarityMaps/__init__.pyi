from __future__ import annotations
import rdkit.Chem.Draw.SimilarityMaps
import typing
from matplotlib.colors import LinearSegmentedColormap
import copy
import math
import matplotlib.cm
import numpy
import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.Draw.rdMolDraw2D
import rdkit.Chem.rdDepictor
import rdkit.Chem.rdMolDescriptors
import rdkit.DataStructs
import rdkit.Geometry
_Shape = typing.Tuple[int, ...]

__all__ = [
    "Chem",
    "DataStructs",
    "Draw",
    "Geometry",
    "GetAPFingerprint",
    "GetAtomicWeightsForFingerprint",
    "GetAtomicWeightsForModel",
    "GetMorganFingerprint",
    "GetRDKFingerprint",
    "GetSimilarityMapForFingerprint",
    "GetSimilarityMapForModel",
    "GetSimilarityMapFromWeights",
    "GetStandardizedWeights",
    "GetTTFingerprint",
    "LinearSegmentedColormap",
    "apDict",
    "cm",
    "copy",
    "math",
    "numpy",
    "rdDepictor",
    "rdMD",
    "rdMolDraw2D",
    "ttDict"
]


apDict: dict # value = {'normal': <function <lambda>>, 'hashed': <function <lambda>>, 'bv': <function <lambda>>}
ttDict: dict # value = {'normal': <function <lambda>>, 'hashed': <function <lambda>>, 'bv': <function <lambda>>}
