"""Module containing functions to generate hashes for molecules"""
from __future__ import annotations
import rdkit.Chem.rdMolHash
import typing
import Boost.Python

__all__ = [
    "HashFunction",
    "MolHash"
]


class HashFunction(Boost.Python.enum, int):
    AnonymousGraph = rdkit.Chem.rdMolHash.HashFunction.AnonymousGraph
    ArthorSubstructureOrder = rdkit.Chem.rdMolHash.HashFunction.ArthorSubstructureOrder
    AtomBondCounts = rdkit.Chem.rdMolHash.HashFunction.AtomBondCounts
    CanonicalSmiles = rdkit.Chem.rdMolHash.HashFunction.CanonicalSmiles
    DegreeVector = rdkit.Chem.rdMolHash.HashFunction.DegreeVector
    ElementGraph = rdkit.Chem.rdMolHash.HashFunction.ElementGraph
    ExtendedMurcko = rdkit.Chem.rdMolHash.HashFunction.ExtendedMurcko
    HetAtomProtomer = rdkit.Chem.rdMolHash.HashFunction.HetAtomProtomer
    HetAtomTautomer = rdkit.Chem.rdMolHash.HashFunction.HetAtomTautomer
    HetAtomTautomerv2 = rdkit.Chem.rdMolHash.HashFunction.HetAtomTautomerv2
    Mesomer = rdkit.Chem.rdMolHash.HashFunction.Mesomer
    MolFormula = rdkit.Chem.rdMolHash.HashFunction.MolFormula
    MurckoScaffold = rdkit.Chem.rdMolHash.HashFunction.MurckoScaffold
    NetCharge = rdkit.Chem.rdMolHash.HashFunction.NetCharge
    RedoxPair = rdkit.Chem.rdMolHash.HashFunction.RedoxPair
    Regioisomer = rdkit.Chem.rdMolHash.HashFunction.Regioisomer
    SmallWorldIndexBR = rdkit.Chem.rdMolHash.HashFunction.SmallWorldIndexBR
    SmallWorldIndexBRL = rdkit.Chem.rdMolHash.HashFunction.SmallWorldIndexBRL
    __slots__ = ()
    names = {'AnonymousGraph': rdkit.Chem.rdMolHash.HashFunction.AnonymousGraph, 'ElementGraph': rdkit.Chem.rdMolHash.HashFunction.ElementGraph, 'CanonicalSmiles': rdkit.Chem.rdMolHash.HashFunction.CanonicalSmiles, 'MurckoScaffold': rdkit.Chem.rdMolHash.HashFunction.MurckoScaffold, 'ExtendedMurcko': rdkit.Chem.rdMolHash.HashFunction.ExtendedMurcko, 'MolFormula': rdkit.Chem.rdMolHash.HashFunction.MolFormula, 'AtomBondCounts': rdkit.Chem.rdMolHash.HashFunction.AtomBondCounts, 'DegreeVector': rdkit.Chem.rdMolHash.HashFunction.DegreeVector, 'Mesomer': rdkit.Chem.rdMolHash.HashFunction.Mesomer, 'HetAtomTautomer': rdkit.Chem.rdMolHash.HashFunction.HetAtomTautomer, 'HetAtomProtomer': rdkit.Chem.rdMolHash.HashFunction.HetAtomProtomer, 'RedoxPair': rdkit.Chem.rdMolHash.HashFunction.RedoxPair, 'Regioisomer': rdkit.Chem.rdMolHash.HashFunction.Regioisomer, 'NetCharge': rdkit.Chem.rdMolHash.HashFunction.NetCharge, 'SmallWorldIndexBR': rdkit.Chem.rdMolHash.HashFunction.SmallWorldIndexBR, 'SmallWorldIndexBRL': rdkit.Chem.rdMolHash.HashFunction.SmallWorldIndexBRL, 'ArthorSubstructureOrder': rdkit.Chem.rdMolHash.HashFunction.ArthorSubstructureOrder, 'HetAtomTautomerv2': rdkit.Chem.rdMolHash.HashFunction.HetAtomTautomerv2}
    values = {1: rdkit.Chem.rdMolHash.HashFunction.AnonymousGraph, 2: rdkit.Chem.rdMolHash.HashFunction.ElementGraph, 3: rdkit.Chem.rdMolHash.HashFunction.CanonicalSmiles, 4: rdkit.Chem.rdMolHash.HashFunction.MurckoScaffold, 5: rdkit.Chem.rdMolHash.HashFunction.ExtendedMurcko, 6: rdkit.Chem.rdMolHash.HashFunction.MolFormula, 7: rdkit.Chem.rdMolHash.HashFunction.AtomBondCounts, 8: rdkit.Chem.rdMolHash.HashFunction.DegreeVector, 9: rdkit.Chem.rdMolHash.HashFunction.Mesomer, 10: rdkit.Chem.rdMolHash.HashFunction.HetAtomTautomer, 11: rdkit.Chem.rdMolHash.HashFunction.HetAtomProtomer, 12: rdkit.Chem.rdMolHash.HashFunction.RedoxPair, 13: rdkit.Chem.rdMolHash.HashFunction.Regioisomer, 14: rdkit.Chem.rdMolHash.HashFunction.NetCharge, 15: rdkit.Chem.rdMolHash.HashFunction.SmallWorldIndexBR, 16: rdkit.Chem.rdMolHash.HashFunction.SmallWorldIndexBRL, 17: rdkit.Chem.rdMolHash.HashFunction.ArthorSubstructureOrder, 18: rdkit.Chem.rdMolHash.HashFunction.HetAtomTautomerv2}
    pass
def MolHash( mol: Mol, func: HashFunction, useCxSmiles: bool = False, cxFlagsToSkip: int = 0) -> str:
    """
    MolHash( mol: Mol, func: HashFunction, useCxSmiles: bool = False, cxFlagsToSkip: int = 0) -> str
        Generate a hash for a molecule. The func argument determines which hash is generated.

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolHash(RDKit::ROMol,RDKit::MolHash::HashFunction [,bool=False [,unsigned int=0]])
    """
