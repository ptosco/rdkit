from __future__ import annotations
import rdkit.Chem.ReducedGraphs
import typing
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "GenerateErGFingerprintForReducedGraph",
    "GenerateMolExtendedReducedGraph",
    "GetErGFingerprint",
    "TanimotoSimilarity",
    "numpy"
]


def GenerateErGFingerprintForReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object:
    """
    GenerateErGFingerprintForReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object
        Returns the ErG fingerprint vector for a reduced graph

        C++ signature :
            _object* GenerateErGFingerprintForReducedGraph(RDKit::ROMol [,boost::python::api::object=0 [,double=0.3 [,int=1 [,int=15]]]])
    """
def GenerateMolExtendedReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0) -> Mol:
    """
    GenerateMolExtendedReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0) -> Mol
        Returns the reduced graph for a molecule

        C++ signature :
            RDKit::ROMol* GenerateMolExtendedReducedGraph(RDKit::ROMol [,boost::python::api::object=0])
    """
def GetErGFingerprint( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object:
    """
    GetErGFingerprint( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object
        Returns the ErG fingerprint vector for a molecule

        C++ signature :
            _object* GetErGFingerprint(RDKit::ROMol [,boost::python::api::object=0 [,double=0.3 [,int=1 [,int=15]]]])
    """
