"""Module containing an assortment of functionality for basic data structures.

At the moment the data structures defined are:
  Bit Vector classes (for storing signatures, fingerprints and the like:
    - ExplicitBitVect: class for relatively small (10s of thousands of bits) or
                       dense bit vectors.
    - SparseBitVect:   class for large, sparse bit vectors
  DiscreteValueVect:   class for storing vectors of integers
  SparseIntVect:       class for storing sparse vectors of integers
"""
from __future__ import annotations
import rdkit.DataStructs
import typing
from rdkit.DataStructs.cDataStructs import DiscreteValueType
from rdkit.DataStructs.cDataStructs import DiscreteValueVect
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.DataStructs.cDataStructs import FPBReader
from rdkit.DataStructs.cDataStructs import IntSparseIntVect
from rdkit.DataStructs.cDataStructs import LongSparseIntVect
from rdkit.DataStructs.cDataStructs import MultiFPBReader
from rdkit.DataStructs.cDataStructs import SparseBitVect
from rdkit.DataStructs.cDataStructs import UIntSparseIntVect
from rdkit.DataStructs.cDataStructs import ULongSparseIntVect
import rdkit.DataStructs.cDataStructs
import rdkit.rdBase

__all__ = [
    "AllBitSimilarity",
    "AllProbeBitsMatch",
    "AsymmetricSimilarity",
    "AsymmetricSimilarityNeighbors",
    "AsymmetricSimilarityNeighbors_sparse",
    "BitVectToBinaryText",
    "BitVectToFPSText",
    "BitVectToText",
    "BraunBlanquetSimilarity",
    "BraunBlanquetSimilarityNeighbors",
    "BraunBlanquetSimilarityNeighbors_sparse",
    "BulkAllBitSimilarity",
    "BulkAsymmetricSimilarity",
    "BulkBraunBlanquetSimilarity",
    "BulkCosineSimilarity",
    "BulkDiceSimilarity",
    "BulkKulczynskiSimilarity",
    "BulkMcConnaugheySimilarity",
    "BulkOnBitSimilarity",
    "BulkRogotGoldbergSimilarity",
    "BulkRusselSimilarity",
    "BulkSokalSimilarity",
    "BulkTanimotoSimilarity",
    "BulkTverskySimilarity",
    "ComputeL1Norm",
    "ConvertToExplicit",
    "ConvertToNumpyArray",
    "CosineSimilarity",
    "CosineSimilarityNeighbors",
    "CosineSimilarityNeighbors_sparse",
    "CreateFromBinaryText",
    "CreateFromBitString",
    "CreateFromFPSText",
    "DiceSimilarity",
    "DiceSimilarityNeighbors",
    "DiceSimilarityNeighbors_sparse",
    "DiscreteValueType",
    "DiscreteValueVect",
    "EIGHTBITVALUE",
    "ExplicitBitVect",
    "FOURBITVALUE",
    "FPBReader",
    "FingerprintSimilarity",
    "FoldFingerprint",
    "FoldToTargetDensity",
    "InitFromDaylightString",
    "IntSparseIntVect",
    "KulczynskiSimilarity",
    "KulczynskiSimilarityNeighbors",
    "KulczynskiSimilarityNeighbors_sparse",
    "LongSparseIntVect",
    "McConnaugheySimilarity",
    "McConnaugheySimilarityNeighbors",
    "McConnaugheySimilarityNeighbors_sparse",
    "MultiFPBReader",
    "NumBitsInCommon",
    "ONEBITVALUE",
    "OffBitProjSimilarity",
    "OffBitsInCommon",
    "OnBitProjSimilarity",
    "OnBitSimilarity",
    "OnBitsInCommon",
    "RogotGoldbergSimilarity",
    "RogotGoldbergSimilarityNeighbors",
    "RogotGoldbergSimilarityNeighbors_sparse",
    "RusselSimilarity",
    "RusselSimilarityNeighbors",
    "RusselSimilarityNeighbors_sparse",
    "SIXTEENBITVALUE",
    "SokalSimilarity",
    "SokalSimilarityNeighbors",
    "SokalSimilarityNeighbors_sparse",
    "SparseBitVect",
    "TWOBITVALUE",
    "TanimotoSimilarity",
    "TanimotoSimilarityNeighbors",
    "TanimotoSimilarityNeighbors_sparse",
    "TopNContainer",
    "TverskySimilarity",
    "UIntSparseIntVect",
    "ULongSparseIntVect",
    "cDataStructs",
    "rdBase",
    "similarityFunctions"
]


@typing.overload
def AllBitSimilarity( v1: SparseBitVect, v2: SparseBitVect) -> float:
    """
    AllBitSimilarity( v1: SparseBitVect, v2: SparseBitVect) -> float

        C++ signature :
            double AllBitSimilarity(SparseBitVect,SparseBitVect)

        C++ signature :
            double AllBitSimilarity(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def AllBitSimilarity( v1: ExplicitBitVect, v2: ExplicitBitVect) -> float:
    pass
@typing.overload
def AllProbeBitsMatch( arg1: SparseBitVect, arg2: SparseBitVect) -> bool:
    """
    AllProbeBitsMatch( arg1: SparseBitVect, arg2: SparseBitVect) -> bool

        C++ signature :
            bool AllProbeBitsMatch(SparseBitVect,SparseBitVect)

        C++ signature :
            bool AllProbeBitsMatch(ExplicitBitVect,ExplicitBitVect)

        C++ signature :
            bool AllProbeBitsMatch(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

        C++ signature :
            bool AllProbeBitsMatch(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
@typing.overload
def AllProbeBitsMatch( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> bool:
    pass
@typing.overload
def AllProbeBitsMatch( arg1: SparseBitVect, arg2: str) -> bool:
    pass
@typing.overload
def AllProbeBitsMatch( arg1: ExplicitBitVect, arg2: str) -> bool:
    pass
@typing.overload
def AsymmetricSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    AsymmetricSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double AsymmetricSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double AsymmetricSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double AsymmetricSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double AsymmetricSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def AsymmetricSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def AsymmetricSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def AsymmetricSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def AsymmetricSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    AsymmetricSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / min(B(bv1),B(bv2))

        C++ signature :
            boost::python::list AsymmetricSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def AsymmetricSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    AsymmetricSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / min(B(bv1),B(bv2))

        C++ signature :
            boost::python::list AsymmetricSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def BitVectToBinaryText( arg1: SparseBitVect) -> object:
    """
    BitVectToBinaryText( arg1: SparseBitVect) -> object

        C++ signature :
            boost::python::api::object BitVectToBinaryText(SparseBitVect)

        C++ signature :
            boost::python::api::object BitVectToBinaryText(ExplicitBitVect)
    """
@typing.overload
def BitVectToBinaryText( arg1: ExplicitBitVect) -> object:
    pass
@typing.overload
def BitVectToFPSText( arg1: SparseBitVect) -> str:
    """
    BitVectToFPSText( arg1: SparseBitVect) -> str

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToFPSText(SparseBitVect)

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToFPSText(ExplicitBitVect)
    """
@typing.overload
def BitVectToFPSText( arg1: ExplicitBitVect) -> str:
    pass
@typing.overload
def BitVectToText( arg1: SparseBitVect) -> str:
    """
    BitVectToText( arg1: SparseBitVect) -> str

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(SparseBitVect)

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(ExplicitBitVect)
    """
@typing.overload
def BitVectToText( arg1: ExplicitBitVect) -> str:
    pass
@typing.overload
def BraunBlanquetSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    BraunBlanquetSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double BraunBlanquetSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double BraunBlanquetSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double BraunBlanquetSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double BraunBlanquetSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def BraunBlanquetSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def BraunBlanquetSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def BraunBlanquetSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def BraunBlanquetSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    BraunBlanquetSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / max(B(bv1),B(bv2))

        C++ signature :
            boost::python::list BraunBlanquetSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def BraunBlanquetSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    BraunBlanquetSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / max(B(bv1),B(bv2))

        C++ signature :
            boost::python::list BraunBlanquetSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
def BulkAllBitSimilarity( v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkAllBitSimilarity( v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkAllBitSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkAllBitSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkAsymmetricSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkAsymmetricSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkAsymmetricSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkAsymmetricSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkAsymmetricSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkBraunBlanquetSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkBraunBlanquetSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkBraunBlanquetSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkBraunBlanquetSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkBraunBlanquetSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkCosineSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkCosineSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkCosineSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkCosineSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkCosineSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkDiceSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkDiceSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkDiceSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkDiceSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<int>,boost::python::list [,bool=False])

        C++ signature :
            boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

        C++ signature :
            boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

        C++ signature :
            boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
@typing.overload
def BulkDiceSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkDiceSimilarity( v1: IntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkDiceSimilarity( v1: LongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkDiceSimilarity( v1: UIntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkDiceSimilarity( v1: ULongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkKulczynskiSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkKulczynskiSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkKulczynskiSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkKulczynskiSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkKulczynskiSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkMcConnaugheySimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkMcConnaugheySimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkMcConnaugheySimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkMcConnaugheySimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkMcConnaugheySimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
def BulkOnBitSimilarity( v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkOnBitSimilarity( v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkOnBitSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkOnBitSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkRogotGoldbergSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkRogotGoldbergSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkRogotGoldbergSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkRogotGoldbergSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkRogotGoldbergSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkRusselSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkRusselSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkRusselSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkRusselSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkRusselSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkSokalSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkSokalSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkSokalSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkSokalSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
@typing.overload
def BulkSokalSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkTanimotoSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    """
    BulkTanimotoSimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkTanimotoSimilarity(SparseBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkTanimotoSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])

        C++ signature :
            boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<int>,boost::python::list [,bool=False])

        C++ signature :
            boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

        C++ signature :
            boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

        C++ signature :
            boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
@typing.overload
def BulkTanimotoSimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkTanimotoSimilarity( v1: IntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTanimotoSimilarity( v1: LongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTanimotoSimilarity( v1: UIntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTanimotoSimilarity( v1: ULongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTverskySimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, a: float, b: float, returnDistance: bool = 0) -> list:
    """
    BulkTverskySimilarity( bv1: SparseBitVect, bvList: AtomPairsParameters, a: float, b: float, returnDistance: bool = 0) -> list

        C++ signature :
            boost::python::list BulkTverskySimilarity(SparseBitVect const*,boost::python::api::object,double,double [,bool=0])

        C++ signature :
            boost::python::list BulkTverskySimilarity(ExplicitBitVect const*,boost::python::api::object,double,double [,bool=0])

        C++ signature :
            boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<int>,boost::python::list,double,double [,bool=False])

        C++ signature :
            boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<long>,boost::python::list,double,double [,bool=False])

        C++ signature :
            boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list,double,double [,bool=False])

        C++ signature :
            boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list,double,double [,bool=False])
    """
@typing.overload
def BulkTverskySimilarity( bv1: ExplicitBitVect, bvList: AtomPairsParameters, a: float, b: float, returnDistance: bool = 0) -> list:
    pass
@typing.overload
def BulkTverskySimilarity( v1: IntSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTverskySimilarity( v1: LongSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTverskySimilarity( v1: UIntSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    pass
@typing.overload
def BulkTverskySimilarity( v1: ULongSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    pass
def ComputeL1Norm( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> int:
    """
    ComputeL1Norm( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> int
        Compute the distance between two discrete vector values
        

        C++ signature :
            unsigned int ComputeL1Norm(RDKit::DiscreteValueVect,RDKit::DiscreteValueVect)
    """
def ConvertToExplicit( arg1: SparseBitVect) -> ExplicitBitVect:
    """
    ConvertToExplicit( arg1: SparseBitVect) -> ExplicitBitVect
        Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector

        C++ signature :
            ExplicitBitVect* ConvertToExplicit(SparseBitVect const*)
    """
@typing.overload
def ConvertToNumpyArray( bv: ExplicitBitVect, destArray: AtomPairsParameters) -> None:
    """
    ConvertToNumpyArray( bv: ExplicitBitVect, destArray: AtomPairsParameters) -> None

        C++ signature :
            void ConvertToNumpyArray(ExplicitBitVect,boost::python::api::object)

        C++ signature :
            void ConvertToNumpyArray(RDKit::DiscreteValueVect,boost::python::api::object)

        C++ signature :
            void ConvertToNumpyArray(RDKit::SparseIntVect<int>,boost::python::api::object)

        C++ signature :
            void ConvertToNumpyArray(RDKit::SparseIntVect<long>,boost::python::api::object)

        C++ signature :
            void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned int>,boost::python::api::object)

        C++ signature :
            void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned long>,boost::python::api::object)
    """
@typing.overload
def ConvertToNumpyArray( bv: DiscreteValueVect, destArray: AtomPairsParameters) -> None:
    pass
@typing.overload
def ConvertToNumpyArray( bv: IntSparseIntVect, destArray: AtomPairsParameters) -> None:
    pass
@typing.overload
def ConvertToNumpyArray( bv: LongSparseIntVect, destArray: AtomPairsParameters) -> None:
    pass
@typing.overload
def ConvertToNumpyArray( bv: UIntSparseIntVect, destArray: AtomPairsParameters) -> None:
    pass
@typing.overload
def ConvertToNumpyArray( bv: ULongSparseIntVect, destArray: AtomPairsParameters) -> None:
    pass
@typing.overload
def CosineSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    CosineSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double CosineSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double CosineSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double CosineSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double CosineSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def CosineSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def CosineSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def CosineSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def CosineSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    CosineSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

        C++ signature :
            boost::python::list CosineSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def CosineSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    CosineSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

        C++ signature :
            boost::python::list CosineSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
def CreateFromBinaryText( arg1: str) -> ExplicitBitVect:
    """
    CreateFromBinaryText( arg1: str) -> ExplicitBitVect
        Creates an ExplicitBitVect from a binary string (byte array).

        C++ signature :
            ExplicitBitVect* CreateFromBinaryText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def CreateFromBitString( arg1: str) -> ExplicitBitVect:
    """
    CreateFromBitString( arg1: str) -> ExplicitBitVect
        Creates an ExplicitBitVect from a bit string (string of 0s and 1s).

        C++ signature :
            ExplicitBitVect* CreateFromBitString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def CreateFromFPSText( arg1: str) -> ExplicitBitVect:
    """
    CreateFromFPSText( arg1: str) -> ExplicitBitVect
        Creates an ExplicitBitVect from an FPS string.

        C++ signature :
            ExplicitBitVect* CreateFromFPSText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
@typing.overload
def DiceSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    DiceSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double DiceSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double DiceSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double DiceSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double DiceSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double DiceSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

        C++ signature :
            double DiceSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

        C++ signature :
            double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

        C++ signature :
            double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
@typing.overload
def DiceSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def DiceSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def DiceSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def DiceSimilarity( siv1: IntSparseIntVect, siv2: IntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def DiceSimilarity( siv1: LongSparseIntVect, siv2: LongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def DiceSimilarity( siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def DiceSimilarity( siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
def DiceSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    DiceSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        2*B(bv1&bv2) / (B(bv1) + B(bv2))

        C++ signature :
            boost::python::list DiceSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def DiceSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    DiceSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        2*B(bv1&bv2) / (B(bv1) + B(bv2))

        C++ signature :
            boost::python::list DiceSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def FoldFingerprint( bv: SparseBitVect, foldFactor: int = 2) -> SparseBitVect:
    """
    FoldFingerprint( bv: SparseBitVect, foldFactor: int = 2) -> SparseBitVect

        C++ signature :
            SparseBitVect* FoldFingerprint(SparseBitVect [,unsigned int=2])

        C++ signature :
            ExplicitBitVect* FoldFingerprint(ExplicitBitVect [,unsigned int=2])
    """
@typing.overload
def FoldFingerprint( bv: ExplicitBitVect, foldFactor: int = 2) -> ExplicitBitVect:
    pass
@typing.overload
def InitFromDaylightString( arg1: SparseBitVect, arg2: str) -> None:
    """
    InitFromDaylightString( arg1: SparseBitVect, arg2: str) -> None

        C++ signature :
            void InitFromDaylightString(SparseBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

        C++ signature :
            void InitFromDaylightString(ExplicitBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
@typing.overload
def InitFromDaylightString( arg1: ExplicitBitVect, arg2: str) -> None:
    pass
@typing.overload
def KulczynskiSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    KulczynskiSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double KulczynskiSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double KulczynskiSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double KulczynskiSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double KulczynskiSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def KulczynskiSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def KulczynskiSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def KulczynskiSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def KulczynskiSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    KulczynskiSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

        C++ signature :
            boost::python::list KulczynskiSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def KulczynskiSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    KulczynskiSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

        C++ signature :
            boost::python::list KulczynskiSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def McConnaugheySimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    McConnaugheySimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double McConnaugheySimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double McConnaugheySimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double McConnaugheySimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double McConnaugheySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def McConnaugheySimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def McConnaugheySimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def McConnaugheySimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def McConnaugheySimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    McConnaugheySimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

        C++ signature :
            boost::python::list McConnaugheySimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def McConnaugheySimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    McConnaugheySimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

        C++ signature :
            boost::python::list McConnaugheySimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def NumBitsInCommon( arg1: SparseBitVect, arg2: SparseBitVect) -> int:
    """
    NumBitsInCommon( arg1: SparseBitVect, arg2: SparseBitVect) -> int

        C++ signature :
            int NumBitsInCommon(SparseBitVect,SparseBitVect)

        C++ signature :
            int NumBitsInCommon(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def NumBitsInCommon( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> int:
    pass
@typing.overload
def OffBitProjSimilarity( arg1: SparseBitVect, arg2: SparseBitVect) -> _vectd:
    """
    OffBitProjSimilarity( arg1: SparseBitVect, arg2: SparseBitVect) -> _vectd

        C++ signature :
            std::vector<double, std::allocator<double> > OffBitProjSimilarity(SparseBitVect,SparseBitVect)

        C++ signature :
            std::vector<double, std::allocator<double> > OffBitProjSimilarity(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def OffBitProjSimilarity( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> _vectd:
    pass
@typing.overload
def OffBitsInCommon( arg1: SparseBitVect, arg2: SparseBitVect) -> _vecti:
    """
    OffBitsInCommon( arg1: SparseBitVect, arg2: SparseBitVect) -> _vecti

        C++ signature :
            std::vector<int, std::allocator<int> > OffBitsInCommon(SparseBitVect,SparseBitVect)

        C++ signature :
            std::vector<int, std::allocator<int> > OffBitsInCommon(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def OffBitsInCommon( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> _vecti:
    pass
@typing.overload
def OnBitProjSimilarity( arg1: SparseBitVect, arg2: SparseBitVect) -> _vectd:
    """
    OnBitProjSimilarity( arg1: SparseBitVect, arg2: SparseBitVect) -> _vectd

        C++ signature :
            std::vector<double, std::allocator<double> > OnBitProjSimilarity(SparseBitVect,SparseBitVect)

        C++ signature :
            std::vector<double, std::allocator<double> > OnBitProjSimilarity(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def OnBitProjSimilarity( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> _vectd:
    pass
@typing.overload
def OnBitSimilarity( v1: SparseBitVect, v2: SparseBitVect) -> float:
    """
    OnBitSimilarity( v1: SparseBitVect, v2: SparseBitVect) -> float

        C++ signature :
            double OnBitSimilarity(SparseBitVect,SparseBitVect)

        C++ signature :
            double OnBitSimilarity(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def OnBitSimilarity( v1: ExplicitBitVect, v2: ExplicitBitVect) -> float:
    pass
@typing.overload
def OnBitsInCommon( arg1: SparseBitVect, arg2: SparseBitVect) -> _vecti:
    """
    OnBitsInCommon( arg1: SparseBitVect, arg2: SparseBitVect) -> _vecti

        C++ signature :
            std::vector<int, std::allocator<int> > OnBitsInCommon(SparseBitVect,SparseBitVect)

        C++ signature :
            std::vector<int, std::allocator<int> > OnBitsInCommon(ExplicitBitVect,ExplicitBitVect)
    """
@typing.overload
def OnBitsInCommon( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> _vecti:
    pass
@typing.overload
def RogotGoldbergSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    RogotGoldbergSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double RogotGoldbergSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double RogotGoldbergSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double RogotGoldbergSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double RogotGoldbergSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def RogotGoldbergSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def RogotGoldbergSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def RogotGoldbergSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def RogotGoldbergSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    RogotGoldbergSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / B(bv1)

        C++ signature :
            boost::python::list RogotGoldbergSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def RogotGoldbergSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    RogotGoldbergSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / B(bv1)

        C++ signature :
            boost::python::list RogotGoldbergSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def RusselSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    RusselSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double RusselSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double RusselSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double RusselSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double RusselSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def RusselSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def RusselSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def RusselSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def RusselSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    RusselSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / B(bv1)

        C++ signature :
            boost::python::list RusselSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def RusselSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    RusselSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / B(bv1)

        C++ signature :
            boost::python::list RusselSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def SokalSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    SokalSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double SokalSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double SokalSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double SokalSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double SokalSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
@typing.overload
def SokalSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def SokalSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def SokalSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
def SokalSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    SokalSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

        C++ signature :
            boost::python::list SokalSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def SokalSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    SokalSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

        C++ signature :
            boost::python::list SokalSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def TanimotoSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    TanimotoSimilarity( bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float

        C++ signature :
            double TanimotoSimilarity(SparseBitVect,SparseBitVect [,bool=0])

        C++ signature :
            double TanimotoSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

        C++ signature :
            double TanimotoSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double TanimotoSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

        C++ signature :
            double TanimotoSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

        C++ signature :
            double TanimotoSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

        C++ signature :
            double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

        C++ signature :
            double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
@typing.overload
def TanimotoSimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def TanimotoSimilarity( bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def TanimotoSimilarity( bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def TanimotoSimilarity( siv1: IntSparseIntVect, siv2: IntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def TanimotoSimilarity( siv1: LongSparseIntVect, siv2: LongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def TanimotoSimilarity( siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def TanimotoSimilarity( siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
def TanimotoSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    TanimotoSimilarityNeighbors( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

        C++ signature :
            boost::python::list TanimotoSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
def TanimotoSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list:
    """
    TanimotoSimilarityNeighbors_sparse( bvqueries: AtomPairsParameters, bvList: AtomPairsParameters) -> list
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

        C++ signature :
            boost::python::list TanimotoSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
@typing.overload
def TverskySimilarity( bv1: SparseBitVect, bv2: SparseBitVect, a: float, b: float, returnDistance: bool = 0) -> float:
    """
    TverskySimilarity( bv1: SparseBitVect, bv2: SparseBitVect, a: float, b: float, returnDistance: bool = 0) -> float

        C++ signature :
            double TverskySimilarity(SparseBitVect,SparseBitVect,double,double [,bool=0])

        C++ signature :
            double TverskySimilarity(ExplicitBitVect,ExplicitBitVect,double,double [,bool=0])

        C++ signature :
            double TverskySimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

        C++ signature :
            double TverskySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

        C++ signature :
            double TverskySimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int>,double,double [,bool=False [,double=0.0]])

        C++ signature :
            double TverskySimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long>,double,double [,bool=False [,double=0.0]])

        C++ signature :
            double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

        C++ signature :
            double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
@typing.overload
def TverskySimilarity( bv1: ExplicitBitVect, bv2: ExplicitBitVect, a: float, b: float, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def TverskySimilarity( bv1: SparseBitVect, pkl: str, a: float, b: float, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def TverskySimilarity( bv1: ExplicitBitVect, pkl: str, a: float, b: float, returnDistance: bool = 0) -> float:
    pass
@typing.overload
def TverskySimilarity( siv1: IntSparseIntVect, siv2: IntSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def TverskySimilarity( siv1: LongSparseIntVect, siv2: LongSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def TverskySimilarity( siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
@typing.overload
def TverskySimilarity( siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    pass
EIGHTBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE
FOURBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE
ONEBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE
SIXTEENBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE
TWOBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE
similarityFunctions: list # value = [('Tanimoto', <Boost.Python.function object>, ''), ('Dice', <Boost.Python.function object>, ''), ('Cosine', <Boost.Python.function object>, ''), ('Sokal', <Boost.Python.function object>, ''), ('Russel', <Boost.Python.function object>, ''), ('RogotGoldberg', <Boost.Python.function object>, ''), ('AllBit', <Boost.Python.function object>, ''), ('Kulczynski', <Boost.Python.function object>, ''), ('McConnaughey', <Boost.Python.function object>, ''), ('Asymmetric', <Boost.Python.function object>, ''), ('BraunBlanquet', <Boost.Python.function object>, '')]
