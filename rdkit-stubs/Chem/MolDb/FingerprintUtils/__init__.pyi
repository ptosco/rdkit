from __future__ import annotations
import rdkit.Chem.MolDb.FingerprintUtils
import typing
import pickle
import rdkit.Chem
import rdkit.DataStructs

__all__ = [
    "BuildAtomPairFP",
    "BuildAvalonFP",
    "BuildMorganFP",
    "BuildPharm2DFP",
    "BuildRDKitFP",
    "BuildSigFactory",
    "BuildTorsionsFP",
    "Chem",
    "DataStructs",
    "DepickleFP",
    "LayeredOptions",
    "pickle",
    "similarityMethods",
    "supportedSimilarityMethods"
]


class LayeredOptions():
    fpSize = 1024
    loadLayerFlags = 4294967295
    maxPath = 6
    minPath = 1
    nWords = 32
    searchLayerFlags = 7
    wordSize = 32
    pass
similarityMethods: dict # value = {'RDK': <class 'rdkit.DataStructs.cDataStructs.ExplicitBitVect'>, 'AtomPairs': <class 'rdkit.DataStructs.cDataStructs.IntSparseIntVect'>, 'TopologicalTorsions': <class 'rdkit.DataStructs.cDataStructs.LongSparseIntVect'>, 'Pharm2D': <class 'rdkit.DataStructs.cDataStructs.SparseBitVect'>, 'Gobbi2D': <class 'rdkit.DataStructs.cDataStructs.SparseBitVect'>, 'Morgan': <class 'rdkit.DataStructs.cDataStructs.UIntSparseIntVect'>, 'Avalon': <class 'rdkit.DataStructs.cDataStructs.ExplicitBitVect'>}
supportedSimilarityMethods = ['RDK', 'AtomPairs', 'TopologicalTorsions', 'Pharm2D', 'Gobbi2D', 'Morgan', 'Avalon']
