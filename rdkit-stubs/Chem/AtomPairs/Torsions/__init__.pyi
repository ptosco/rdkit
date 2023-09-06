"""
Contains an implementation of Topological-torsion fingerprints, as
described in:

R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
"Topological Torsion: A New Molecular Descriptor for SAR Applications.
Comparison with Other Descriptors" JCICS 27, 82-85 (1987).

The fingerprints can be accessed through the following functions:
- GetTopologicalTorsionFingerprint
- GetHashedTopologicalTorsionFingerprint
- GetTopologicalTorsionFingerprintAsIntVect (identical to GetTopologicalTorsionFingerprint)
- GetTopologicalTorsionFingerprintAsIds

"""
from __future__ import annotations
import rdkit.Chem.AtomPairs.Torsions
import typing
import rdkit.Chem.AtomPairs.Utils
import rdkit.Chem.rdMolDescriptors

__all__ = [
    "ExplainPathScore",
    "GetHashedTopologicalTorsionFingerprint",
    "GetTopologicalTorsionFingerprint",
    "GetTopologicalTorsionFingerprintAsIds",
    "GetTopologicalTorsionFingerprintAsIntVect",
    "Utils",
    "pyScorePath",
    "rdMolDescriptors"
]


def GetHashedTopologicalTorsionFingerprint( mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
    GetHashedTopologicalTorsionFingerprint( mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect
        Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<long>* GetHashedTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=2048 [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]]])
    """
def GetTopologicalTorsionFingerprint( mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
    GetTopologicalTorsionFingerprint( mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect
        Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<long>* GetTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]])
    """
