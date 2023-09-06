""" handles doing cross validation with neural nets

This is, perhaps, a little misleading.  For the purposes of this module,
cross validation == evaluating the accuracy of a net.

"""
from __future__ import annotations
import rdkit.ML.Neural.CrossValidate
import typing
import math
import rdkit.ML.Data.SplitData
import rdkit.ML.Neural.Network
import rdkit.ML.Neural.Trainers

__all__ = [
    "CrossValidate",
    "CrossValidationDriver",
    "Network",
    "SplitData",
    "Trainers",
    "math"
]


