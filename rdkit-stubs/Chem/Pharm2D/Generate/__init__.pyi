""" generation of 2D pharmacophores

**Notes**

  - The terminology for this gets a bit rocky, so here's a glossary of what
    terms used here mean:

      1) *N-point pharmacophore* a combination of N features along with
         distances between them.

      2) *N-point proto-pharmacophore*: a combination of N feature
         definitions without distances.  Each N-point
         proto-pharmacophore defines a manifold of potential N-point
         pharmacophores.

      3) *N-point scaffold*: a collection of the distances defining
         an N-point pharmacophore without feature identities.

  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
from __future__ import annotations
import rdkit.Chem.Pharm2D.Generate
import typing
import rdkit.Chem.Pharm2D.SigFactory
import rdkit.Chem.Pharm2D.Utils
import rdkit.RDLogger

__all__ = [
    "Gen2DFingerprint",
    "SigFactory",
    "Utils",
    "logger"
]


_verbose = 0
logger: rdkit.RDLogger.logger
