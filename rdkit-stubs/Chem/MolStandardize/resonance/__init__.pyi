"""
molvs.resonance
~~~~~~~~~~~~~~~

Resonance (mesomeric) transformations.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""
from __future__ import annotations
import rdkit.Chem.MolStandardize.resonance
import typing
import logging
import rdkit.Chem

__all__ = [
    "Chem",
    "MAX_STRUCTURES",
    "ResonanceEnumerator",
    "enumerate_resonance_smiles",
    "log",
    "logging",
    "warn"
]


class ResonanceEnumerator():
    """
    Simple wrapper around RDKit ResonanceMolSupplier.

        
    """
    pass
MAX_STRUCTURES = 1000
log: logging.Logger # value = <Logger rdkit.Chem.MolStandardize.resonance (WARNING)>
