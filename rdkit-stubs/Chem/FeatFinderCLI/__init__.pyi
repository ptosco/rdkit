from __future__ import annotations
import rdkit.Chem.FeatFinderCLI
import typing
import argparse
import os
import rdkit.Chem
import rdkit.Chem.ChemicalFeatures
import rdkit.RDLogger
import re

__all__ = [
    "Chem",
    "ChemicalFeatures",
    "GetAtomFeatInfo",
    "RDLogger",
    "argparse",
    "existingFile",
    "initParser",
    "logger",
    "main",
    "os",
    "processArgs",
    "re",
    "splitExpr"
]


_splashMessage = '\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n  FeatFinderCLI\n  Part of the RDKit (http://www.rdkit.org)\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n'
logger: rdkit.RDLogger.logger
splitExpr: re.Pattern # value = re.compile('[ \\t,]')
