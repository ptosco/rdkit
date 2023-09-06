from __future__ import annotations
import rdkit.Chem.PandasPatcher
import typing
from xml.parsers.expat import ExpatError
from rdkit.Chem.rdchem import Mol
from _io import StringIO
import importlib
import logging
import pandas
import pandas.io.formats
import pandas.io.formats.format
import pandas.io.formats.html
import re
import xml.dom.minidom

__all__ = [
    "ExpatError",
    "InteractiveRenderer",
    "Mol",
    "MolFormatter",
    "PrintAsImageString",
    "RDK_MOLS_AS_IMAGE_ATTR",
    "StringIO",
    "adj",
    "changeMoleculeRendering",
    "check_rdk_attr",
    "dataframeformatter_class",
    "get_adjustment_name",
    "html_formatter_class",
    "html_formatter_module",
    "html_formatter_module_name",
    "importlib",
    "is_molecule_image",
    "log",
    "logging",
    "minidom",
    "molJustify",
    "orig_get_adjustment",
    "orig_get_formatter",
    "orig_to_html",
    "orig_write_cell",
    "pandas_formats",
    "pandas_formats_name",
    "pandas_frame",
    "patchPandas",
    "patched_get_adjustment",
    "patched_get_formatter",
    "patched_to_html",
    "patched_write_cell",
    "pd",
    "pprint_thing",
    "re",
    "renderImagesInAllDataFrames",
    "set_rdk_attr",
    "styleRegex",
    "to_html_class",
    "to_html_class_name",
    "unpatchPandas"
]


class MolFormatter():
    """
    Format molecules as images
    """
    pass
InteractiveRenderer = None
PrintAsImageString = None
RDK_MOLS_AS_IMAGE_ATTR = '__rdkitMolAsImage'
adj: pandas.io.formats.format.TextAdjustment
get_adjustment_name = 'get_adjustment'
html_formatter_module_name = 'html'
log: logging.Logger # value = <Logger rdkit.Chem.PandasPatcher (WARNING)>
molJustify = None
pandas_formats_name = 'pandas.io'
pandas_frame = None
styleRegex: re.Pattern # value = re.compile('^(.*style=["\'][^"^\']*)(["\'].*)$')
to_html_class_name = 'DataFrameRenderer'
dataframeformatter_class = pandas.io.formats.format.DataFrameFormatter
html_formatter_class = pandas.io.formats.html.HTMLFormatter
orig_get_adjustment = pandas.io.formats.format.get_adjustment
orig_get_formatter = pandas.io.formats.format.DataFrameFormatter._get_formatter
orig_to_html = pandas.io.formats.format.DataFrameRenderer.to_html
orig_write_cell = pandas.io.formats.html.HTMLFormatter._write_cell
to_html_class = pandas.io.formats.format.DataFrameRenderer
