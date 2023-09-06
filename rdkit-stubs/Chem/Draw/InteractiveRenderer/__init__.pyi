""" Interactive molecule rendering through rdkit-structure-renderer.js """
from __future__ import annotations
import rdkit.Chem.Draw.InteractiveRenderer
import typing
from IPython.core.display import HTML
import base64
import json
import logging
import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.Draw.rdMolDraw2D
import re
import uuid
import xml.dom.minidom

__all__ = [
    "Chem",
    "Draw",
    "HTML",
    "MolsToHTMLTable",
    "base64",
    "camelCaseOptToDataTag",
    "clearOpt",
    "clearOpts",
    "display",
    "filterDefaultDrawOpts",
    "generateHTMLBody",
    "generateHTMLFooter",
    "getOpts",
    "injectHTMLFooterAfterTable",
    "isEnabled",
    "isNoInteractive",
    "json",
    "log",
    "logging",
    "minidom",
    "minimalLibJsUrl",
    "newlineToXml",
    "parentNodeQuery",
    "rdMolDraw2D",
    "rdkitStructureRendererJsUrl",
    "re",
    "setEnabled",
    "setNoInteractive",
    "setOpt",
    "setOpts",
    "toDataMol",
    "uuid",
    "xmlToNewline"
]


_camelCaseOptToTagRe: re.Pattern # value = re.compile('[A-Z]')
_defaultDrawOptions: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
_defaultDrawOptionsDict = None
_disabled = '_disabled'
_enabled_div_uuid = None
_opts = '__rnrOpts'
log: logging.Logger # value = <Logger rdkit.Chem.Draw.InteractiveRenderer (WARNING)>
minimalLibJsUrl = 'https://unpkg.com/rdkit-structure-renderer/public/RDKit_minimal.js'
parentNodeQuery = 'div[class*=jp-NotebookPanel-notebook]'
rdkitStructureRendererJsUrl = 'https://unpkg.com/rdkit-structure-renderer/dist/rdkit-structure-renderer-module.js'
