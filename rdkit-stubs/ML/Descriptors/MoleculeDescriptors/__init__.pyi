""" Various bits and pieces for calculating Molecular descriptors

"""
from __future__ import annotations
import rdkit.ML.Descriptors.MoleculeDescriptors
import typing
import pickle
import rdkit.Chem.Descriptors
import rdkit.ML.Descriptors.Descriptors
import rdkit.RDLogger
import re

__all__ = [
    "Descriptors",
    "DescriptorsMod",
    "MolecularDescriptorCalculator",
    "logger",
    "pickle",
    "re"
]


class MolecularDescriptorCalculator(rdkit.ML.Descriptors.Descriptors.DescriptorCalculator):
    """
    used for calculating descriptors for molecules

     
    """
    pass
logger: rdkit.RDLogger.logger
