""" Contains an implementation of Atom-pair fingerprints, as
described in:

R.E. Carhart, D.H. Smith, R. Venkataraghavan;
"Atom Pairs as Molecular Features in Structure-Activity Studies:
Definition and Applications" JCICS 25, 64-73 (1985).

The fingerprints can be accessed through the following functions:
- GetAtomPairFingerprint
- GetHashedAtomPairFingerprint (identical to GetAtomPairFingerprint)
- GetAtomPairFingerprintAsIntVect
- GetAtomPairFingerprintAsBitVect

"""
from __future__ import annotations
import rdkit.Chem.AtomPairs.Pairs
import typing
import rdkit.Chem.AtomPairs.Utils
import rdkit.Chem.rdMolDescriptors
import rdkit.DataStructs

__all__ = [
    "DataStructs",
    "ExplainPairScore",
    "GetAtomPairFingerprint",
    "GetAtomPairFingerprintAsBitVect",
    "GetAtomPairFingerprintAsIntVect",
    "GetHashedAtomPairFingerprint",
    "Utils",
    "fpLen",
    "numFpBits",
    "numPathBits",
    "pyScorePair",
    "rdMolDescriptors"
]


def GetAtomPairFingerprint( mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
    GetAtomPairFingerprint( mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect
        Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<int>* GetAtomPairFingerprint(RDKit::ROMol [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]])
    """
def GetHashedAtomPairFingerprint( mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
    GetHashedAtomPairFingerprint( mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect
        Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<int>* GetHashedAtomPairFingerprint(RDKit::ROMol [,unsigned int=2048 [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]]])
    """
_maxPathLen = 31
fpLen = 8388608
numFpBits = 23
numPathBits = 5
