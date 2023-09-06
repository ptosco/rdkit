from __future__ import annotations
import rdkit.Chem.Draw.IPythonConsole
import typing
from _io import BytesIO
from IPython.core.display import HTML
from PIL.PngImagePlugin import PngInfo
from IPython.core.display import SVG
import IPython
import IPython.display
import PIL.Image
import base64
import copy
import html
import py3Dmol
import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.Draw.InteractiveRenderer
import rdkit.Chem.Draw.rdMolDraw2D
import rdkit.Chem.rdChemReactions
import rdkit.Chem.rdchem
import warnings

__all__ = [
    "BytesIO",
    "Chem",
    "DisableSubstructMatchRendering",
    "Draw",
    "DrawMorganBit",
    "DrawMorganBits",
    "DrawRDKitBit",
    "DrawRDKitBits",
    "EnableSubstructMatchRendering",
    "HTML",
    "IPython",
    "Image",
    "InstallIPythonRenderer",
    "InteractiveRenderer",
    "PngInfo",
    "SVG",
    "ShowMols",
    "UninstallIPythonRenderer",
    "addMolToView",
    "base64",
    "bgcolor_3d",
    "copy",
    "display",
    "display_pil_image",
    "drawMol3D",
    "drawOptions",
    "drawing_type_3d",
    "highlightByReactant",
    "highlightSubstructs",
    "html",
    "ipython_3d",
    "ipython_maxProperties",
    "ipython_showProperties",
    "ipython_useSVG",
    "kekulizeStructures",
    "molSize",
    "molSize_3d",
    "py3Dmol",
    "rdChemReactions",
    "rdMolDraw2D",
    "rdchem",
    "warnings"
]


_canUse3D = True
_methodsToDelete: list # value = [(<class 'rdkit.Chem.rdchem.Mol'>, '_repr_png_'), (<class 'rdkit.Chem.rdchem.Mol'>, '_repr_svg_'), (<class 'rdkit.Chem.rdchem.Mol'>, '_repr_html_'), (<class 'rdkit.Chem.rdChemReactions.ChemicalReaction'>, '_repr_png_'), (<class 'rdkit.Chem.rdChemReactions.ChemicalReaction'>, '_repr_svg_'), (<class 'rdkit.Chem.rdchem.MolBundle'>, '_repr_png_'), (<class 'rdkit.Chem.rdchem.MolBundle'>, '_repr_svg_'), (<class 'PIL.Image.Image'>, '_repr_png_')]
_rendererInstalled = True
bgcolor_3d = '0xeeeeee'
drawOptions: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
drawing_type_3d = 'stick'
highlightByReactant = False
highlightSubstructs = True
ipython_3d = False
ipython_maxProperties = 10
ipython_showProperties = True
ipython_useSVG = False
kekulizeStructures = True
molSize = (450, 150)
molSize_3d = (400, 400)
_DrawMorganBitSaved = rdkit.Chem.Draw.DrawMorganBit
_DrawMorganBitsSaved = rdkit.Chem.Draw.DrawMorganBits
_DrawRDKitBitSaved = rdkit.Chem.Draw.DrawRDKitBit
_DrawRDKitBitsSaved = rdkit.Chem.Draw.DrawRDKitBits
_MolsToGridImageSaved = rdkit.Chem.Draw.MolsToGridImage
