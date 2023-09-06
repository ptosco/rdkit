from __future__ import annotations
import rdkit.Chem.Draw.spingCanvas
import typing
from rdkit.Chem.Draw.canvasbase import CanvasBase
import rdkit.Chem.Draw.canvasbase
import rdkit.sping.pid
import re

__all__ = [
    "Canvas",
    "CanvasBase",
    "convertColor",
    "faceMap",
    "pid",
    "re"
]


class Canvas(rdkit.Chem.Draw.canvasbase.CanvasBase):
    pass
faceMap = {'sans': 'helvetica', 'serif': 'times'}
