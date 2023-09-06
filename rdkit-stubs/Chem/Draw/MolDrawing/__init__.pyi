from __future__ import annotations
import rdkit.Chem.Draw.MolDrawing
import typing
import copy
import functools
import math
import numpy
import rdkit.Chem
import rdkit.Chem.rdchem
_Shape = typing.Tuple[int, ...]

__all__ = [
    "Chem",
    "DrawingOptions",
    "Font",
    "MolDrawing",
    "cmp",
    "copy",
    "functools",
    "math",
    "numpy",
    "periodicTable"
]


class DrawingOptions():
    atomLabelDeuteriumTritium = False
    atomLabelFontFace = 'sans'
    atomLabelFontSize = 12
    atomLabelMinFontSize = 7
    atomNumberOffset = 0
    bgColor = (1, 1, 1)
    bondLineWidth = 1.2
    colorBonds = True
    coordScale = 1.0
    dash = (4, 4)
    dblBondLengthFrac = 0.8
    dblBondOffset = 0.25
    defaultColor = (1, 0, 0)
    dotsPerAngstrom = 30
    elemDict = {1: (0.55, 0.55, 0.55), 7: (0, 0, 1), 8: (1, 0, 0), 9: (0.2, 0.8, 0.8), 15: (1, 0.5, 0), 16: (0.8, 0.8, 0), 17: (0, 0.8, 0), 35: (0.5, 0.3, 0.1), 53: (0.63, 0.12, 0.94), 0: (0.5, 0.5, 0.5)}
    includeAtomNumbers = False
    noCarbonSymbols = True
    radicalSymbol = '∙'
    selectColor = (1, 0, 0)
    showUnknownDoubleBonds = True
    useFraction = 0.85
    wedgeDashedBonds = True
    pass
class Font():
    pass
class MolDrawing():
    pass
periodicTable: rdkit.Chem.rdchem.PeriodicTable
