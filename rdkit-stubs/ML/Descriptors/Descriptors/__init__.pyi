""" Various bits and pieces for calculating descriptors

"""
from __future__ import annotations
import rdkit.ML.Descriptors.Descriptors
import typing
import pickle

__all__ = [
    "DescriptorCalculator",
    "pickle"
]


class DescriptorCalculator():
    """
    abstract base class for descriptor calculators

     
    """
    pass
