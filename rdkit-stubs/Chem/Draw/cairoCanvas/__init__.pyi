from __future__ import annotations
import rdkit.Chem.Draw.cairoCanvas
import typing
from rdkit.Chem.Draw.canvasbase import CanvasBase
import PIL.Image
import array
import cairo
import math
import os
import rdkit.Chem.Draw.canvasbase
import re

__all__ = [
    "Canvas",
    "CanvasBase",
    "Image",
    "array",
    "cairo",
    "have_cairocffi",
    "have_pango",
    "libType",
    "math",
    "os",
    "pango",
    "pangocairo",
    "re",
    "scriptPattern"
]


class Canvas(rdkit.Chem.Draw.canvasbase.CanvasBase):
    pass
have_cairocffi = False
have_pango = None
libType = 'pangocairo'
pango = None
pangocairo = None
scriptPattern: re.Pattern # value = re.compile('\\<.+?\\>')
