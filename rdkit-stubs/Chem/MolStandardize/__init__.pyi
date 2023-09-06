"""
MolVS - Molecule Validation and Standardization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolVS is a python tool built on top of RDKit that performs validation and standardization of chemical structures.

Note that the C++ reimplementation of this is available in the module rdkit.Chem.MolStandardize.rdMolStandardize

:copyright: (c) 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""
from __future__ import annotations
import rdkit.Chem.MolStandardize
import typing
from rdkit.Chem.MolStandardize.errors import MolVSError
from rdkit.Chem.MolStandardize.errors import StandardizeError
from rdkit.Chem.MolStandardize.standardize import Standardizer
from rdkit.Chem.MolStandardize.errors import ValidateError
from rdkit.Chem.MolStandardize.validate import Validator
import logging
import rdkit.Chem

__all__ = [
    "Chem",
    "MolVSError",
    "ReorderTautomers",
    "StandardizeError",
    "Standardizer",
    "ValidateError",
    "Validator",
    "canonicalize_tautomer_smiles",
    "charge",
    "enumerate_tautomers_smiles",
    "errors",
    "fragment",
    "log",
    "logging",
    "metal",
    "normalize",
    "rdMolStandardize",
    "standardize",
    "standardize_smiles",
    "tautomer",
    "utils",
    "validate",
    "validate_smiles",
    "validations"
]


__author__ = 'Matt Swain'
__copyright__ = 'Copyright 2016 Matt Swain'
__email__ = 'm.swain@me.com'
__license__ = 'MIT'
__title__ = 'MolVS'
__version__ = '0.1.1'
__warningregistry__: dict # value = {'version': 74, ('The module rdkit.Chem.MolStandardize.standardize is deprecated and will be removed in the next release.', <class 'DeprecationWarning'>, 20): True, ('The module rdkit.Chem.MolStandardize.validate is deprecated and will be removed in the next release.', <class 'DeprecationWarning'>, 22): True}
log: logging.Logger # value = <Logger rdkit.Chem.MolStandardize (WARNING)>
