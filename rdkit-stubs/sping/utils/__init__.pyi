from __future__ import annotations
import rdkit.sping.utils
import typing

__all__ = [
    "availableCanvases"
]


availableCanvases = {'pdf': ('PDF.pidPDF', 'PDFCanvas', 'PDF'), 'ps': ('PS.pidPS', 'PSCanvas', 'PS'), 'svg': ('SVG.pidSVG', 'SVGCanvas', 'SVG'), 'jpg': ('PIL.pidPIL', 'PILCanvas', 'JPEG'), 'png': ('PIL.pidPIL', 'PILCanvas', 'PNG')}
