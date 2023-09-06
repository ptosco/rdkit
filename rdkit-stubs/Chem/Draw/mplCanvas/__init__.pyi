from __future__ import annotations
import rdkit.Chem.Draw.mplCanvas
import typing
from rdkit.Chem.Draw.canvasbase import CanvasBase
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import rdkit.Chem.Draw.canvasbase

__all__ = [
    "Canvas",
    "CanvasBase",
    "Line2D",
    "Polygon",
    "figure"
]


class Canvas(rdkit.Chem.Draw.canvasbase.CanvasBase):
    pass
