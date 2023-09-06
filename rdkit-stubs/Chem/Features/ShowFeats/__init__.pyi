from __future__ import annotations
import rdkit.Chem.Features.ShowFeats
import typing
from optparse import OptionParser
import math
import optparse
import os
import rdkit.Chem.Features.FeatDirUtilsRD
import rdkit.Geometry
import rdkit.RDConfig
import rdkit.RDLogger
import sys

__all__ = [
    "ALPHA",
    "BEGIN",
    "COLOR",
    "CYLINDER",
    "END",
    "FeatDirUtils",
    "Geometry",
    "NORMAL",
    "OptionParser",
    "RDConfig",
    "SPHERE",
    "ShowArrow",
    "ShowMolFeats",
    "TRIANGLE_FAN",
    "VERTEX",
    "logger",
    "logging",
    "math",
    "os",
    "parser",
    "sys"
]


ALPHA = 25
BEGIN = 2
COLOR = 6
CYLINDER = 9
END = 3
NORMAL = 5
SPHERE = 7
TRIANGLE_FAN = 6
VERTEX = 4
_canonArrowhead = None
_featColors = {'Donor': (0, 1, 1), 'Acceptor': (1, 0, 1), 'NegIonizable': (1, 0, 0), 'PosIonizable': (0, 0, 1), 'ZnBinder': (1, 0.5, 0.5), 'Aromatic': (1, 0.8, 0.2), 'LumpedHydrophobe': (0.5, 0.25, 0), 'Hydrophobe': (0.5, 0.25, 0)}
_globalArrowCGO = []
_globalSphereCGO = []
_usage = '\n   ShowFeats [optional args] <filenames>\n\n  if "-" is provided as a filename, data will be read from stdin (the console)\n'
_version = '0.3.2'
_welcomeMessage = 'This is ShowFeats version 0.3.2'
logger: rdkit.RDLogger.logger
parser: optparse.OptionParser
