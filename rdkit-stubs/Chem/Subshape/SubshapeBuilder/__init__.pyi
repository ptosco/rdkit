from __future__ import annotations
import rdkit.Chem.Subshape.SubshapeBuilder
import typing
import copy
import pickle
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.Subshape.BuilderUtils
import rdkit.Chem.Subshape.SubshapeObjects
import rdkit.Geometry
import time

__all__ = [
    "AllChem",
    "BuilderUtils",
    "Chem",
    "Geometry",
    "SubshapeBuilder",
    "SubshapeCombineOperations",
    "SubshapeObjects",
    "copy",
    "pickle",
    "time"
]


class SubshapeBuilder():
    featFactory = None
    fraction = 0.25
    gridDims = (20, 15, 10)
    gridSpacing = 0.5
    nbrCount = 7
    stepSize = 1.0
    terminalPtRadScale = 0.75
    winRad = 3.0
    pass
class SubshapeCombineOperations():
    INTERSECT = 2
    SUM = 1
    UNION = 0
    pass
