"""
$Id$

Scoring - Calculate rank statistics

Created by Sereina Riniker, October 2012
after a file from Peter Gedeck, Greg Landrum

\param scores: ordered list with descending similarity containing
               active/inactive information
\param col: column index in scores where active/inactive information is stored
\param fractions: list of fractions at which the value shall be calculated
\param alpha: exponential weight
"""
from __future__ import annotations
import rdkit.ML.Scoring.Scoring
import typing
import math

__all__ = [
    "CalcAUC",
    "CalcBEDROC",
    "CalcEnrichment",
    "CalcRIE",
    "CalcROC",
    "math",
    "namedtuple"
]


