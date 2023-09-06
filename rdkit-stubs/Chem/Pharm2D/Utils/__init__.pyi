""" utility functionality for the 2D pharmacophores code

  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
from __future__ import annotations
import rdkit.Chem.Pharm2D.Utils
import typing
import itertools

__all__ = [
    "BinsTriangleInequality",
    "CountUpTo",
    "GetAllCombinations",
    "GetIndexCombinations",
    "GetPossibleScaffolds",
    "GetTriangles",
    "GetUniqueCombinations",
    "GetUniqueCombinations_new",
    "NumCombinations",
    "OrderTriangle",
    "ScaffoldPasses",
    "UniquifyCombinations",
    "comb",
    "itertools",
    "nDistPointDict",
    "nPointDistDict"
]


_countCache = {}
_indexCombinations = {}
_numCombDict = {(8, 2): 36, (8, 3): 120}
_trianglesInPharmacophore = {3: ((0, 1, 2),)}
_verbose = 0
nDistPointDict = {1: 2, 3: 3, 5: 4, 7: 5, 9: 6, 11: 7, 13: 8, 15: 9, 17: 10}
nPointDistDict = {2: ((0, 1),), 3: ((0, 1), (0, 2), (1, 2)), 4: ((0, 1), (0, 2), (0, 3), (1, 2), (2, 3)), 5: ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)), 6: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)), 7: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)), 8: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)), 9: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)), 10: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9))}
