""" Contains an implementation of Physicochemical property fingerprints, as
described in:
Kearsley, S. K. et al.
"Chemical Similarity Using Physiochemical Property Descriptors."
J. Chem.Inf. Model. 36, 118-127 (1996)

The fingerprints can be accessed through the following functions:
- GetBPFingerprint
- GetBTFingerprint

"""
from __future__ import annotations
import rdkit.Chem.AtomPairs.Sheridan
import typing
import os
import rdkit.Chem
import rdkit.Chem.rdMolDescriptors
import rdkit.RDConfig
import re

__all__ = [
    "AssignPattyTypes",
    "Chem",
    "GetAtomPairFingerprint",
    "GetBPFingerprint",
    "GetBTFingerprint",
    "GetTopologicalTorsionFingerprint",
    "RDConfig",
    "fpLen",
    "numFpBits",
    "numPathBits",
    "os",
    "rdMolDescriptors",
    "re",
    "typMap"
]


def GetAtomPairFingerprint( mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
    GetAtomPairFingerprint( mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect
        Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<int>* GetAtomPairFingerprint(RDKit::ROMol [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]])
    """
def GetTopologicalTorsionFingerprint( mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
    GetTopologicalTorsionFingerprint( mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect
        Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<long>* GetTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]])
    """
_maxPathLen = 31
_pattyDefs = None
fpLen = 8388608
numFpBits = 23
numPathBits = 5
typMap = {'CAT': 1, 'ANI': 2, 'POL': 3, 'DON': 4, 'ACC': 5, 'HYD': 6, 'OTH': 7}
