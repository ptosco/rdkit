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
import rdkit.DataStructs.cDataStructs
import typing
import Boost.Python

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
    "FoldFingerprint",
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
    "TverskySimilarity",
    "UIntSparseIntVect",
    "ULongSparseIntVect"
]


class DiscreteValueType(Boost.Python.enum, int):
    EIGHTBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE
    FOURBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE
    ONEBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE
    SIXTEENBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE
    TWOBITVALUE = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE
    __slots__ = ()
    names = {'ONEBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE, 'TWOBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, 'FOURBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE, 'EIGHTBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE, 'SIXTEENBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE}
    values = {0: rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE, 1: rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, 2: rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE, 3: rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE, 4: rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE}
    pass
class DiscreteValueVect(Boost.Python.instance):
    """
    A container class for storing unsigned integer
    values within a particular range.

    The length of the vector and type of its elements (determines the maximum value
    that can be stored) are both set at construction time.

    As you would expect, _DiscreteValueVects_ support a set of binary operations
    so you can do things like:
      dvv3 = dvv1 & dvv2  the result contains the smallest value in each entry
      dvv3 = dvv1 | dvv2  the result contains the largest value in each entry
      dvv1 += dvv2     values are truncated when necessary
      dvv3 = dvv1 + dvv2    values are truncated when necessary
      dvv1 -= dvv3    would-be negative values are set to zero
      dvv3 = dvv1 - dvv2    would-be negative values are set to zero

    Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])
    """
    @staticmethod
    def GetTotalVal( arg1: DiscreteValueVect) -> int: 
        """
        GetTotalVal( arg1: DiscreteValueVect) -> int
            Get the sum of the values in the vector, basically L1 norm

            C++ signature :
                unsigned int GetTotalVal(RDKit::DiscreteValueVect {lvalue})
        """
    @staticmethod
    def GetValueType( arg1: DiscreteValueVect) -> DiscreteValueType: 
        """
        GetValueType( arg1: DiscreteValueVect) -> DiscreteValueType
            Get the type of value stored in the vector

            C++ signature :
                RDKit::DiscreteValueVect::DiscreteValueType GetValueType(RDKit::DiscreteValueVect {lvalue})
        """
    @staticmethod
    def __add__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object: 
        """
        __add__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object

            C++ signature :
                _object* __add__(RDKit::DiscreteValueVect {lvalue},RDKit::DiscreteValueVect)
        """
    @staticmethod
    def __and__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object: 
        """
        __and__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object

            C++ signature :
                _object* __and__(RDKit::DiscreteValueVect {lvalue},RDKit::DiscreteValueVect)
        """
    @staticmethod
    def __getinitargs__( arg1: DiscreteValueVect) -> tuple: 
        """
        __getinitargs__( arg1: DiscreteValueVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::DiscreteValueVect)
        """
    @staticmethod
    def __getitem__( arg1: DiscreteValueVect, arg2: int) -> int: 
        """
        __getitem__( arg1: DiscreteValueVect, arg2: int) -> int
            Get the value at a specified location

            C++ signature :
                unsigned int __getitem__(RDKit::DiscreteValueVect {lvalue},unsigned int)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    def __iadd__( arg1: object, arg2: DiscreteValueVect) -> object: 
        """
        __iadd__( arg1: object, arg2: DiscreteValueVect) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::DiscreteValueVect&>,RDKit::DiscreteValueVect)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: DiscreteValueType, arg3: int) -> None: 
        """
        __init__( arg1: object, arg2: DiscreteValueType, arg3: int) -> None
            Constructor

            C++ signature :
                void __init__(_object*,RDKit::DiscreteValueVect::DiscreteValueType,unsigned int)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    def __isub__( arg1: object, arg2: DiscreteValueVect) -> object: 
        """
        __isub__( arg1: object, arg2: DiscreteValueVect) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::DiscreteValueVect&>,RDKit::DiscreteValueVect)
        """
    @staticmethod
    def __len__( arg1: DiscreteValueVect) -> int: 
        """
        __len__( arg1: DiscreteValueVect) -> int
            Get the number of entries in the vector

            C++ signature :
                unsigned int __len__(RDKit::DiscreteValueVect {lvalue})
        """
    @staticmethod
    def __or__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object: 
        """
        __or__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object

            C++ signature :
                _object* __or__(RDKit::DiscreteValueVect {lvalue},RDKit::DiscreteValueVect)
        """
    @staticmethod
    def __setitem__( arg1: DiscreteValueVect, arg2: int, arg3: int) -> None: 
        """
        __setitem__( arg1: DiscreteValueVect, arg2: int, arg3: int) -> None
            Set the value at a specified location

            C++ signature :
                void __setitem__(RDKit::DiscreteValueVect {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object: 
        """
        __sub__( arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> object

            C++ signature :
                _object* __sub__(RDKit::DiscreteValueVect {lvalue},RDKit::DiscreteValueVect)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 64
    __safe_for_unpickling__ = True
    pass
class ExplicitBitVect(Boost.Python.instance):
    """
    A class to store explicit bit vectors.

    This class is most useful for situations where the size of the vector
    is relatively small (tens of thousands or smaller).

    For larger vectors, use the _SparseBitVect_ class instead.

    As you would expect, _ExplicitBitVects_ support a set of binary operations
    so you can do things like:
      bv3 = bv1 & bv2  (bitwise and)
      bv3 = bv1 | bv2  (bitwise or)
      bv3 = bv1 ^ bv2  (bitwise xor)
      bv3 = ~bv1       (bitwise negation)

    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).
    """
    @staticmethod
    def FromBase64( arg1: ExplicitBitVect, arg2: str) -> None: 
        """
        FromBase64( arg1: ExplicitBitVect, arg2: str) -> None
            Initializes the vector from a base64 encoded binary string.
            

            C++ signature :
                void FromBase64(ExplicitBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetBit( arg1: ExplicitBitVect, arg2: int) -> bool: 
        """
        GetBit( arg1: ExplicitBitVect, arg2: int) -> bool
            Returns the value of a bit.
            

            C++ signature :
                bool GetBit(ExplicitBitVect {lvalue},unsigned int)
        """
    @staticmethod
    def GetNumBits( arg1: ExplicitBitVect) -> int: 
        """
        GetNumBits( arg1: ExplicitBitVect) -> int
            Returns the number of bits in the vector (the vector's size).
            

            C++ signature :
                unsigned int GetNumBits(ExplicitBitVect {lvalue})
        """
    @staticmethod
    def GetNumOffBits( arg1: ExplicitBitVect) -> int: 
        """
        GetNumOffBits( arg1: ExplicitBitVect) -> int
            Returns the number of off bits.
            

            C++ signature :
                unsigned int GetNumOffBits(ExplicitBitVect {lvalue})
        """
    @staticmethod
    def GetNumOnBits( arg1: ExplicitBitVect) -> int: 
        """
        GetNumOnBits( arg1: ExplicitBitVect) -> int
            Returns the number of on bits.
            

            C++ signature :
                unsigned int GetNumOnBits(ExplicitBitVect {lvalue})
        """
    @staticmethod
    def GetOnBits( arg1: ExplicitBitVect) -> _vecti: 
        """
        GetOnBits( arg1: ExplicitBitVect) -> _vecti
            Returns a tuple containing IDs of the on bits.
            

            C++ signature :
                std::vector<int, std::allocator<int> > GetOnBits(ExplicitBitVect)
        """
    @staticmethod
    def SetBit( arg1: ExplicitBitVect, arg2: int) -> bool: 
        """
        SetBit( arg1: ExplicitBitVect, arg2: int) -> bool
            Turns on a particular bit.  Returns the original state of the bit.
            

            C++ signature :
                bool SetBit(ExplicitBitVect {lvalue},unsigned int)
        """
    @staticmethod
    def SetBitsFromList( arg1: ExplicitBitVect, arg2: AtomPairsParameters) -> None: 
        """
        SetBitsFromList( arg1: ExplicitBitVect, arg2: AtomPairsParameters) -> None
            Turns on a set of bits.  The argument should be a tuple or list of bit ids.
            

            C++ signature :
                void SetBitsFromList(ExplicitBitVect*,boost::python::api::object)
        """
    @staticmethod
    def ToBase64( arg1: ExplicitBitVect) -> str: 
        """
        ToBase64( arg1: ExplicitBitVect) -> str
            Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ToBase64(ExplicitBitVect {lvalue})
        """
    @staticmethod
    def ToBinary( arg1: ExplicitBitVect) -> object: 
        """
        ToBinary( arg1: ExplicitBitVect) -> object
            Returns an internal binary representation of the vector.
            

            C++ signature :
                boost::python::api::object ToBinary(ExplicitBitVect)
        """
    @staticmethod
    def ToList( arg1: ExplicitBitVect) -> list: 
        """
        ToList( arg1: ExplicitBitVect) -> list
            Return the Bitvector as a python list (faster than list(vect))

            C++ signature :
                boost::python::list ToList(ExplicitBitVect)
        """
    @staticmethod
    def UnSetBit( arg1: ExplicitBitVect, arg2: int) -> bool: 
        """
        UnSetBit( arg1: ExplicitBitVect, arg2: int) -> bool
            Turns off a particular bit.  Returns the original state of the bit.
            

            C++ signature :
                bool UnSetBit(ExplicitBitVect {lvalue},unsigned int)
        """
    @staticmethod
    def UnSetBitsFromList( arg1: ExplicitBitVect, arg2: AtomPairsParameters) -> None: 
        """
        UnSetBitsFromList( arg1: ExplicitBitVect, arg2: AtomPairsParameters) -> None
            Turns off a set of bits.  The argument should be a tuple or list of bit ids.
            

            C++ signature :
                void UnSetBitsFromList(ExplicitBitVect*,boost::python::api::object)
        """
    @staticmethod
    def __add__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object: 
        """
        __add__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __add__(ExplicitBitVect {lvalue},ExplicitBitVect)
        """
    @staticmethod
    def __and__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object: 
        """
        __and__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __and__(ExplicitBitVect {lvalue},ExplicitBitVect)
        """
    @staticmethod
    def __eq__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object: 
        """
        __eq__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __eq__(ExplicitBitVect {lvalue},ExplicitBitVect)
        """
    @staticmethod
    def __getinitargs__( arg1: ExplicitBitVect) -> tuple: 
        """
        __getinitargs__( arg1: ExplicitBitVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(ExplicitBitVect)
        """
    @staticmethod
    def __getitem__( arg1: ExplicitBitVect, arg2: int) -> int: 
        """
        __getitem__( arg1: ExplicitBitVect, arg2: int) -> int

            C++ signature :
                int __getitem__(ExplicitBitVect,int)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    def __iadd__( arg1: object, arg2: ExplicitBitVect) -> object: 
        """
        __iadd__( arg1: object, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<ExplicitBitVect&>,ExplicitBitVect)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None

            C++ signature :
                void __init__(_object*,unsigned int)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,unsigned int,bool)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, size: int, bitsSet: bool) -> None: ...
    @staticmethod
    def __invert__( arg1: ExplicitBitVect) -> object: 
        """
        __invert__( arg1: ExplicitBitVect) -> object

            C++ signature :
                _object* __invert__(ExplicitBitVect {lvalue})
        """
    @staticmethod
    def __len__( arg1: ExplicitBitVect) -> int: 
        """
        __len__( arg1: ExplicitBitVect) -> int

            C++ signature :
                unsigned int __len__(ExplicitBitVect {lvalue})
        """
    @staticmethod
    def __ne__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object: 
        """
        __ne__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __ne__(ExplicitBitVect {lvalue},ExplicitBitVect)
        """
    @staticmethod
    def __or__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object: 
        """
        __or__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __or__(ExplicitBitVect {lvalue},ExplicitBitVect)
        """
    @staticmethod
    def __setitem__( arg1: ExplicitBitVect, arg2: int, arg3: int) -> int: 
        """
        __setitem__( arg1: ExplicitBitVect, arg2: int, arg3: int) -> int

            C++ signature :
                int __setitem__(ExplicitBitVect {lvalue},int,int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __xor__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object: 
        """
        __xor__( arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> object

            C++ signature :
                _object* __xor__(ExplicitBitVect {lvalue},ExplicitBitVect)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class FPBReader(Boost.Python.instance):
    """
    A class for reading and searching FPB files from Andrew Dalke's chemfp.
        Note that this functionality is still experimental and the API may
        change in future releases.
    """
    @staticmethod
    def GetBytes( arg1: FPBReader, arg2: int) -> object: 
        """
        GetBytes( arg1: FPBReader, arg2: int) -> object
            returns a particular fingerprint as bytes

            C++ signature :
                boost::python::api::object GetBytes(RDKit::FPBReader const*,unsigned int)
        """
    @staticmethod
    def GetContainingNeighbors( arg1: FPBReader, bv: str) -> tuple: 
        """
        GetContainingNeighbors( arg1: FPBReader, bv: str) -> tuple
            returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)

            C++ signature :
                boost::python::tuple GetContainingNeighbors(RDKit::FPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetFP( arg1: FPBReader, arg2: int) -> ExplicitBitVect: 
        """
        GetFP( arg1: FPBReader, arg2: int) -> ExplicitBitVect
            returns a particular fingerprint as an ExplicitBitVect

            C++ signature :
                boost::shared_ptr<ExplicitBitVect> GetFP(RDKit::FPBReader {lvalue},unsigned int)
        """
    @staticmethod
    def GetId( arg1: FPBReader, arg2: int) -> str: 
        """
        GetId( arg1: FPBReader, arg2: int) -> str
            returns the id of a particular fingerprint

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetId(RDKit::FPBReader {lvalue},unsigned int)
        """
    @staticmethod
    def GetNumBits( arg1: FPBReader) -> int: 
        """
        GetNumBits( arg1: FPBReader) -> int
            returns the number of bits in a fingerprint

            C++ signature :
                unsigned int GetNumBits(RDKit::FPBReader {lvalue})
        """
    @staticmethod
    def GetTanimoto( arg1: FPBReader, arg2: int, arg3: str) -> float: 
        """
        GetTanimoto( arg1: FPBReader, arg2: int, arg3: str) -> float
            return the tanimoto similarity of a particular fingerprint to the bytes provided

            C++ signature :
                double GetTanimoto(RDKit::FPBReader const*,unsigned int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetTanimotoNeighbors( arg1: FPBReader, bv: str, threshold: float = 0.7) -> tuple: 
        """
        GetTanimotoNeighbors( arg1: FPBReader, bv: str, threshold: float = 0.7) -> tuple
            returns tanimoto similarities to and indices of all neighbors above the specified threshold

            C++ signature :
                boost::python::tuple GetTanimotoNeighbors(RDKit::FPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,double=0.7])
        """
    @staticmethod
    def GetTversky( arg1: FPBReader, arg2: int, arg3: str, arg4: float, arg5: float) -> float: 
        """
        GetTversky( arg1: FPBReader, arg2: int, arg3: str, arg4: float, arg5: float) -> float
            return the Tverksy similarity of a particular fingerprint to the bytes provided

            C++ signature :
                double GetTversky(RDKit::FPBReader const*,unsigned int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double)
        """
    @staticmethod
    def GetTverskyNeighbors( arg1: FPBReader, bv: str, ca: float, cb: float, threshold: float = 0.7) -> tuple: 
        """
        GetTverskyNeighbors( arg1: FPBReader, bv: str, ca: float, cb: float, threshold: float = 0.7) -> tuple
            returns Tversky similarities to and indices of all neighbors above the specified threshold

            C++ signature :
                boost::python::tuple GetTverskyNeighbors(RDKit::FPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,double=0.7])
        """
    @staticmethod
    def Init( arg1: FPBReader) -> None: 
        """
        Init( arg1: FPBReader) -> None
            Read the fingerprints from the file. This can take a while.
            

            C++ signature :
                void Init(RDKit::FPBReader {lvalue})
        """
    @staticmethod
    def __getitem__( arg1: FPBReader, arg2: int) -> tuple: 
        """
        __getitem__( arg1: FPBReader, arg2: int) -> tuple

            C++ signature :
                boost::python::tuple __getitem__(RDKit::FPBReader const*,unsigned int)
        """
    @staticmethod
    def __init__( arg1: object, filename: str, lazy: bool = False) -> None: 
        """
        __init__( arg1: object, filename: str, lazy: bool = False) -> None
            docstring

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    @staticmethod
    def __len__( arg1: FPBReader) -> int: 
        """
        __len__( arg1: FPBReader) -> int

            C++ signature :
                unsigned int __len__(RDKit::FPBReader {lvalue})
        """
    __instance_size__ = 48
    pass
class IntSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """
    @staticmethod
    def GetLength( arg1: IntSparseIntVect) -> int: 
        """
        GetLength( arg1: IntSparseIntVect) -> int
            Returns the length of the vector

            C++ signature :
                int GetLength(RDKit::SparseIntVect<int> {lvalue})
        """
    @staticmethod
    def GetNonzeroElements( arg1: IntSparseIntVect) -> dict: 
        """
        GetNonzeroElements( arg1: IntSparseIntVect) -> dict
            returns a dictionary of the nonzero elements

            C++ signature :
                boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<int> {lvalue})
        """
    @staticmethod
    def GetTotalVal( arg1: IntSparseIntVect, useAbs: bool = False) -> int: 
        """
        GetTotalVal( arg1: IntSparseIntVect, useAbs: bool = False) -> int
            Get the sum of the values in the vector, basically L1 norm

            C++ signature :
                int GetTotalVal(RDKit::SparseIntVect<int> {lvalue} [,bool=False])
        """
    @staticmethod
    def ToBinary( arg1: IntSparseIntVect) -> object: 
        """
        ToBinary( arg1: IntSparseIntVect) -> object
            returns a binary (pickle) representation of the vector

            C++ signature :
                boost::python::api::object ToBinary(RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def ToList( arg1: IntSparseIntVect) -> list: 
        """
        ToList( arg1: IntSparseIntVect) -> list
            Return the SparseIntVect as a python list

            C++ signature :
                boost::python::list ToList(RDKit::SparseIntVect<int> {lvalue})
        """
    @staticmethod
    def UpdateFromSequence( arg1: IntSparseIntVect, arg2: AtomPairsParameters) -> None: 
        """
        UpdateFromSequence( arg1: IntSparseIntVect, arg2: AtomPairsParameters) -> None
            update the vector based on the values in the list or tuple

            C++ signature :
                void UpdateFromSequence(RDKit::SparseIntVect<int> {lvalue},boost::python::api::object {lvalue})
        """
    @staticmethod
    def __add__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object: 
        """
        __add__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __add__(RDKit::SparseIntVect<int> {lvalue},RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def __and__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object: 
        """
        __and__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __and__(RDKit::SparseIntVect<int> {lvalue},RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def __eq__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object: 
        """
        __eq__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __eq__(RDKit::SparseIntVect<int> {lvalue},RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def __getinitargs__( arg1: IntSparseIntVect) -> tuple: 
        """
        __getinitargs__( arg1: IntSparseIntVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def __getitem__( arg1: IntSparseIntVect, arg2: int) -> int: 
        """
        __getitem__( arg1: IntSparseIntVect, arg2: int) -> int
            Get the value at a specified location

            C++ signature :
                int __getitem__(RDKit::SparseIntVect<int> {lvalue},int)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: IntSparseIntVect) -> object: 
        """
        __iadd__( arg1: object, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<int>&>,RDKit::SparseIntVect<int>)

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<int>&>,int)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __idiv__( arg1: object, arg2: int) -> object: 
        """
        __idiv__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __idiv__(boost::python::back_reference<RDKit::SparseIntVect<int>&>,int)
        """
    @staticmethod
    def __imul__( arg1: object, arg2: int) -> object: 
        """
        __imul__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __imul__(boost::python::back_reference<RDKit::SparseIntVect<int>&>,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None
            Constructor

            C++ signature :
                void __init__(_object*,int)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: IntSparseIntVect) -> object: 
        """
        __isub__( arg1: object, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<int>&>,RDKit::SparseIntVect<int>)

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<int>&>,int)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __ne__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object: 
        """
        __ne__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __ne__(RDKit::SparseIntVect<int> {lvalue},RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def __or__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object: 
        """
        __or__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __or__(RDKit::SparseIntVect<int> {lvalue},RDKit::SparseIntVect<int>)
        """
    @staticmethod
    def __setitem__( arg1: IntSparseIntVect, arg2: int, arg3: int) -> None: 
        """
        __setitem__( arg1: IntSparseIntVect, arg2: int, arg3: int) -> None
            Set the value at a specified location

            C++ signature :
                void __setitem__(RDKit::SparseIntVect<int> {lvalue},int,int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object: 
        """
        __sub__( arg1: IntSparseIntVect, arg2: IntSparseIntVect) -> object

            C++ signature :
                _object* __sub__(RDKit::SparseIntVect<int> {lvalue},RDKit::SparseIntVect<int>)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class LongSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """
    @staticmethod
    def GetLength( arg1: LongSparseIntVect) -> int: 
        """
        GetLength( arg1: LongSparseIntVect) -> int
            Returns the length of the vector

            C++ signature :
                long GetLength(RDKit::SparseIntVect<long> {lvalue})
        """
    @staticmethod
    def GetNonzeroElements( arg1: LongSparseIntVect) -> dict: 
        """
        GetNonzeroElements( arg1: LongSparseIntVect) -> dict
            returns a dictionary of the nonzero elements

            C++ signature :
                boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<long> {lvalue})
        """
    @staticmethod
    def GetTotalVal( arg1: LongSparseIntVect, useAbs: bool = False) -> int: 
        """
        GetTotalVal( arg1: LongSparseIntVect, useAbs: bool = False) -> int
            Get the sum of the values in the vector, basically L1 norm

            C++ signature :
                int GetTotalVal(RDKit::SparseIntVect<long> {lvalue} [,bool=False])
        """
    @staticmethod
    def ToBinary( arg1: LongSparseIntVect) -> object: 
        """
        ToBinary( arg1: LongSparseIntVect) -> object
            returns a binary (pickle) representation of the vector

            C++ signature :
                boost::python::api::object ToBinary(RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def ToList( arg1: LongSparseIntVect) -> list: 
        """
        ToList( arg1: LongSparseIntVect) -> list
            Return the SparseIntVect as a python list

            C++ signature :
                boost::python::list ToList(RDKit::SparseIntVect<long> {lvalue})
        """
    @staticmethod
    def UpdateFromSequence( arg1: LongSparseIntVect, arg2: AtomPairsParameters) -> None: 
        """
        UpdateFromSequence( arg1: LongSparseIntVect, arg2: AtomPairsParameters) -> None
            update the vector based on the values in the list or tuple

            C++ signature :
                void UpdateFromSequence(RDKit::SparseIntVect<long> {lvalue},boost::python::api::object {lvalue})
        """
    @staticmethod
    def __add__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object: 
        """
        __add__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __add__(RDKit::SparseIntVect<long> {lvalue},RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def __and__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object: 
        """
        __and__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __and__(RDKit::SparseIntVect<long> {lvalue},RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def __eq__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object: 
        """
        __eq__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __eq__(RDKit::SparseIntVect<long> {lvalue},RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def __getinitargs__( arg1: LongSparseIntVect) -> tuple: 
        """
        __getinitargs__( arg1: LongSparseIntVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def __getitem__( arg1: LongSparseIntVect, arg2: int) -> int: 
        """
        __getitem__( arg1: LongSparseIntVect, arg2: int) -> int
            Get the value at a specified location

            C++ signature :
                int __getitem__(RDKit::SparseIntVect<long> {lvalue},long)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: LongSparseIntVect) -> object: 
        """
        __iadd__( arg1: object, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<long>&>,RDKit::SparseIntVect<long>)

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<long>&>,int)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __idiv__( arg1: object, arg2: int) -> object: 
        """
        __idiv__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __idiv__(boost::python::back_reference<RDKit::SparseIntVect<long>&>,int)
        """
    @staticmethod
    def __imul__( arg1: object, arg2: int) -> object: 
        """
        __imul__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __imul__(boost::python::back_reference<RDKit::SparseIntVect<long>&>,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None
            Constructor

            C++ signature :
                void __init__(_object*,long)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: LongSparseIntVect) -> object: 
        """
        __isub__( arg1: object, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<long>&>,RDKit::SparseIntVect<long>)

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<long>&>,int)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __ne__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object: 
        """
        __ne__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __ne__(RDKit::SparseIntVect<long> {lvalue},RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def __or__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object: 
        """
        __or__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __or__(RDKit::SparseIntVect<long> {lvalue},RDKit::SparseIntVect<long>)
        """
    @staticmethod
    def __setitem__( arg1: LongSparseIntVect, arg2: int, arg3: int) -> None: 
        """
        __setitem__( arg1: LongSparseIntVect, arg2: int, arg3: int) -> None
            Set the value at a specified location

            C++ signature :
                void __setitem__(RDKit::SparseIntVect<long> {lvalue},long,int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object: 
        """
        __sub__( arg1: LongSparseIntVect, arg2: LongSparseIntVect) -> object

            C++ signature :
                _object* __sub__(RDKit::SparseIntVect<long> {lvalue},RDKit::SparseIntVect<long>)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class MultiFPBReader(Boost.Python.instance):
    """
    A class for reading and searching multiple FPB files from Andrew Dalke's chemfp.
        Note that this functionality is still experimental and the API may
        change in future releases.
    """
    @staticmethod
    def AddReader( arg1: MultiFPBReader, arg2: FPBReader) -> int: 
        """
        AddReader( arg1: MultiFPBReader, arg2: FPBReader) -> int
            adds an FPBReader to our set of readers

            C++ signature :
                unsigned int AddReader(RDKit::MultiFPBReader {lvalue},RDKit::FPBReader*)
        """
    @staticmethod
    def GetContainingNeighbors( arg1: MultiFPBReader, bv: str, numThreads: int = 1) -> tuple: 
        """
        GetContainingNeighbors( arg1: MultiFPBReader, bv: str, numThreads: int = 1) -> tuple
            returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)

            C++ signature :
                boost::python::tuple GetContainingNeighbors(RDKit::MultiFPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned int=1])
        """
    @staticmethod
    def GetNumBits( arg1: MultiFPBReader) -> int: 
        """
        GetNumBits( arg1: MultiFPBReader) -> int
            returns the number of bits in a fingerprint

            C++ signature :
                unsigned int GetNumBits(RDKit::MultiFPBReader {lvalue})
        """
    @staticmethod
    def GetReader( arg1: MultiFPBReader, arg2: int) -> FPBReader: 
        """
        GetReader( arg1: MultiFPBReader, arg2: int) -> FPBReader
            returns one of our readers

            C++ signature :
                RDKit::FPBReader* GetReader(RDKit::MultiFPBReader {lvalue},unsigned int)
        """
    @staticmethod
    def GetTanimotoNeighbors( arg1: MultiFPBReader, bv: str, threshold: float = 0.7, numThreads: int = 1) -> tuple: 
        """
        GetTanimotoNeighbors( arg1: MultiFPBReader, bv: str, threshold: float = 0.7, numThreads: int = 1) -> tuple
            returns tanimoto similarities to and indices of all neighbors above the specified threshold

            C++ signature :
                boost::python::tuple GetTanimotoNeighbors(RDKit::MultiFPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,double=0.7 [,unsigned int=1]])
        """
    @staticmethod
    def GetTverskyNeighbors( arg1: MultiFPBReader, bv: str, ca: float, cb: float, threshold: float = 0.7, numThreads: int = 1) -> tuple: 
        """
        GetTverskyNeighbors( arg1: MultiFPBReader, bv: str, ca: float, cb: float, threshold: float = 0.7, numThreads: int = 1) -> tuple
            returns Tversky similarities to and indices of all neighbors above the specified threshold

            C++ signature :
                boost::python::tuple GetTverskyNeighbors(RDKit::MultiFPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,double=0.7 [,unsigned int=1]])
        """
    @staticmethod
    def Init( arg1: MultiFPBReader) -> None: 
        """
        Init( arg1: MultiFPBReader) -> None
            Call Init() on each of our children. This can take a while.
            

            C++ signature :
                void Init(RDKit::MultiFPBReader {lvalue})
        """
    @staticmethod
    def __init__( arg1: object, initOnSearch: bool = False) -> None: 
        """
        __init__( arg1: object, initOnSearch: bool = False) -> None
            docstring

            C++ signature :
                void __init__(_object* [,bool=False])
        """
    @staticmethod
    def __len__( arg1: MultiFPBReader) -> int: 
        """
        __len__( arg1: MultiFPBReader) -> int

            C++ signature :
                unsigned int __len__(RDKit::MultiFPBReader {lvalue})
        """
    __instance_size__ = 56
    pass
class SparseBitVect(Boost.Python.instance):
    """
    A class to store sparse bit vectors.

    This class is most useful for situations where the size of the vector
    is large and relatively few bits are set

    For smaller or denser vectors, the _ExplicitBitVect_ class is much faster.

    As you would expect, _SparseBitVects_ support a set of binary operations
    so you can do things like:
      bv3 = bv1 & bv2  (bitwise and)
      bv3 = bv1 | bv2  (bitwise or)
      bv3 = bv1 ^ bv2  (bitwise xor)
      bv3 = ~bv1       (bitwise negation) NOTE: this operation is likely
                        to be VERY slow and inefficient.

    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).
    """
    @staticmethod
    def FromBase64( arg1: SparseBitVect, arg2: str) -> None: 
        """
        FromBase64( arg1: SparseBitVect, arg2: str) -> None
            Initializes the vector from a base64 encoded binary string.
            

            C++ signature :
                void FromBase64(SparseBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetBit( arg1: SparseBitVect, arg2: int) -> bool: 
        """
        GetBit( arg1: SparseBitVect, arg2: int) -> bool
            Returns the value of a bit.
            

            C++ signature :
                bool GetBit(SparseBitVect {lvalue},unsigned int)
        """
    @staticmethod
    def GetNumBits( arg1: SparseBitVect) -> int: 
        """
        GetNumBits( arg1: SparseBitVect) -> int
            Returns the number of bits in the vector (the vector's size).
            

            C++ signature :
                unsigned int GetNumBits(SparseBitVect {lvalue})
        """
    @staticmethod
    def GetNumOffBits( arg1: SparseBitVect) -> int: 
        """
        GetNumOffBits( arg1: SparseBitVect) -> int
            Returns the number of off bits.
            

            C++ signature :
                unsigned int GetNumOffBits(SparseBitVect {lvalue})
        """
    @staticmethod
    def GetNumOnBits( arg1: SparseBitVect) -> int: 
        """
        GetNumOnBits( arg1: SparseBitVect) -> int
            Returns the number of on bits.
            

            C++ signature :
                unsigned int GetNumOnBits(SparseBitVect {lvalue})
        """
    @staticmethod
    def GetOnBits( arg1: SparseBitVect) -> _vecti: 
        """
        GetOnBits( arg1: SparseBitVect) -> _vecti
            Returns a tuple containing IDs of the on bits.
            

            C++ signature :
                std::vector<int, std::allocator<int> > GetOnBits(SparseBitVect)
        """
    @staticmethod
    def SetBit( arg1: SparseBitVect, arg2: int) -> bool: 
        """
        SetBit( arg1: SparseBitVect, arg2: int) -> bool
            Turns on a particular bit.  Returns the original state of the bit.
            

            C++ signature :
                bool SetBit(SparseBitVect {lvalue},unsigned int)
        """
    @staticmethod
    def SetBitsFromList( arg1: SparseBitVect, arg2: AtomPairsParameters) -> None: 
        """
        SetBitsFromList( arg1: SparseBitVect, arg2: AtomPairsParameters) -> None
            Turns on a set of bits.  The argument should be a tuple or list of bit ids.
            

            C++ signature :
                void SetBitsFromList(SparseBitVect*,boost::python::api::object)
        """
    @staticmethod
    def ToBase64( arg1: SparseBitVect) -> str: 
        """
        ToBase64( arg1: SparseBitVect) -> str
            Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ToBase64(SparseBitVect {lvalue})
        """
    @staticmethod
    def ToBinary( arg1: SparseBitVect) -> object: 
        """
        ToBinary( arg1: SparseBitVect) -> object
            Returns an internal binary representation of the vector.
            

            C++ signature :
                boost::python::api::object ToBinary(SparseBitVect)
        """
    @staticmethod
    def ToList( arg1: SparseBitVect) -> list: 
        """
        ToList( arg1: SparseBitVect) -> list
            Return the BitVector as a python list

            C++ signature :
                boost::python::list ToList(SparseBitVect)
        """
    @staticmethod
    def UnSetBit( arg1: SparseBitVect, arg2: int) -> bool: 
        """
        UnSetBit( arg1: SparseBitVect, arg2: int) -> bool
            Turns off a particular bit.  Returns the original state of the bit.
            

            C++ signature :
                bool UnSetBit(SparseBitVect {lvalue},unsigned int)
        """
    @staticmethod
    def UnSetBitsFromList( arg1: SparseBitVect, arg2: AtomPairsParameters) -> None: 
        """
        UnSetBitsFromList( arg1: SparseBitVect, arg2: AtomPairsParameters) -> None
            Turns off a set of bits.  The argument should be a tuple or list of bit ids.
            

            C++ signature :
                void UnSetBitsFromList(SparseBitVect*,boost::python::api::object)
        """
    @staticmethod
    def __and__( arg1: SparseBitVect, arg2: SparseBitVect) -> object: 
        """
        __and__( arg1: SparseBitVect, arg2: SparseBitVect) -> object

            C++ signature :
                _object* __and__(SparseBitVect {lvalue},SparseBitVect)
        """
    @staticmethod
    def __eq__( arg1: SparseBitVect, arg2: SparseBitVect) -> object: 
        """
        __eq__( arg1: SparseBitVect, arg2: SparseBitVect) -> object

            C++ signature :
                _object* __eq__(SparseBitVect {lvalue},SparseBitVect)
        """
    @staticmethod
    def __getinitargs__( arg1: SparseBitVect) -> tuple: 
        """
        __getinitargs__( arg1: SparseBitVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(SparseBitVect)
        """
    @staticmethod
    def __getitem__( arg1: SparseBitVect, arg2: int) -> int: 
        """
        __getitem__( arg1: SparseBitVect, arg2: int) -> int

            C++ signature :
                int __getitem__(SparseBitVect,int)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None

            C++ signature :
                void __init__(_object*,unsigned int)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    def __invert__( arg1: SparseBitVect) -> object: 
        """
        __invert__( arg1: SparseBitVect) -> object

            C++ signature :
                _object* __invert__(SparseBitVect {lvalue})
        """
    @staticmethod
    def __len__( arg1: SparseBitVect) -> int: 
        """
        __len__( arg1: SparseBitVect) -> int

            C++ signature :
                unsigned int __len__(SparseBitVect {lvalue})
        """
    @staticmethod
    def __ne__( arg1: SparseBitVect, arg2: SparseBitVect) -> object: 
        """
        __ne__( arg1: SparseBitVect, arg2: SparseBitVect) -> object

            C++ signature :
                _object* __ne__(SparseBitVect {lvalue},SparseBitVect)
        """
    @staticmethod
    def __or__( arg1: SparseBitVect, arg2: SparseBitVect) -> object: 
        """
        __or__( arg1: SparseBitVect, arg2: SparseBitVect) -> object

            C++ signature :
                _object* __or__(SparseBitVect {lvalue},SparseBitVect)
        """
    @staticmethod
    def __setitem__( arg1: SparseBitVect, arg2: int, arg3: int) -> int: 
        """
        __setitem__( arg1: SparseBitVect, arg2: int, arg3: int) -> int

            C++ signature :
                int __setitem__(SparseBitVect {lvalue},int,int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __xor__( arg1: SparseBitVect, arg2: SparseBitVect) -> object: 
        """
        __xor__( arg1: SparseBitVect, arg2: SparseBitVect) -> object

            C++ signature :
                _object* __xor__(SparseBitVect {lvalue},SparseBitVect)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class UIntSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """
    @staticmethod
    def GetLength( arg1: UIntSparseIntVect) -> int: 
        """
        GetLength( arg1: UIntSparseIntVect) -> int
            Returns the length of the vector

            C++ signature :
                unsigned int GetLength(RDKit::SparseIntVect<unsigned int> {lvalue})
        """
    @staticmethod
    def GetNonzeroElements( arg1: UIntSparseIntVect) -> dict: 
        """
        GetNonzeroElements( arg1: UIntSparseIntVect) -> dict
            returns a dictionary of the nonzero elements

            C++ signature :
                boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<unsigned int> {lvalue})
        """
    @staticmethod
    def GetTotalVal( arg1: UIntSparseIntVect, useAbs: bool = False) -> int: 
        """
        GetTotalVal( arg1: UIntSparseIntVect, useAbs: bool = False) -> int
            Get the sum of the values in the vector, basically L1 norm

            C++ signature :
                int GetTotalVal(RDKit::SparseIntVect<unsigned int> {lvalue} [,bool=False])
        """
    @staticmethod
    def ToBinary( arg1: UIntSparseIntVect) -> object: 
        """
        ToBinary( arg1: UIntSparseIntVect) -> object
            returns a binary (pickle) representation of the vector

            C++ signature :
                boost::python::api::object ToBinary(RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def ToList( arg1: UIntSparseIntVect) -> list: 
        """
        ToList( arg1: UIntSparseIntVect) -> list
            Return the SparseIntVect as a python list

            C++ signature :
                boost::python::list ToList(RDKit::SparseIntVect<unsigned int> {lvalue})
        """
    @staticmethod
    def UpdateFromSequence( arg1: UIntSparseIntVect, arg2: AtomPairsParameters) -> None: 
        """
        UpdateFromSequence( arg1: UIntSparseIntVect, arg2: AtomPairsParameters) -> None
            update the vector based on the values in the list or tuple

            C++ signature :
                void UpdateFromSequence(RDKit::SparseIntVect<unsigned int> {lvalue},boost::python::api::object {lvalue})
        """
    @staticmethod
    def __add__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object: 
        """
        __add__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __add__(RDKit::SparseIntVect<unsigned int> {lvalue},RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def __and__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object: 
        """
        __and__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __and__(RDKit::SparseIntVect<unsigned int> {lvalue},RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def __eq__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object: 
        """
        __eq__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __eq__(RDKit::SparseIntVect<unsigned int> {lvalue},RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def __getinitargs__( arg1: UIntSparseIntVect) -> tuple: 
        """
        __getinitargs__( arg1: UIntSparseIntVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def __getitem__( arg1: UIntSparseIntVect, arg2: int) -> int: 
        """
        __getitem__( arg1: UIntSparseIntVect, arg2: int) -> int
            Get the value at a specified location

            C++ signature :
                int __getitem__(RDKit::SparseIntVect<unsigned int> {lvalue},unsigned int)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: UIntSparseIntVect) -> object: 
        """
        __iadd__( arg1: object, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<unsigned int>&>,RDKit::SparseIntVect<unsigned int>)

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<unsigned int>&>,int)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __idiv__( arg1: object, arg2: int) -> object: 
        """
        __idiv__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __idiv__(boost::python::back_reference<RDKit::SparseIntVect<unsigned int>&>,int)
        """
    @staticmethod
    def __imul__( arg1: object, arg2: int) -> object: 
        """
        __imul__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __imul__(boost::python::back_reference<RDKit::SparseIntVect<unsigned int>&>,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None
            Constructor

            C++ signature :
                void __init__(_object*,unsigned int)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: UIntSparseIntVect) -> object: 
        """
        __isub__( arg1: object, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<unsigned int>&>,RDKit::SparseIntVect<unsigned int>)

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<unsigned int>&>,int)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __ne__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object: 
        """
        __ne__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __ne__(RDKit::SparseIntVect<unsigned int> {lvalue},RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def __or__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object: 
        """
        __or__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __or__(RDKit::SparseIntVect<unsigned int> {lvalue},RDKit::SparseIntVect<unsigned int>)
        """
    @staticmethod
    def __setitem__( arg1: UIntSparseIntVect, arg2: int, arg3: int) -> None: 
        """
        __setitem__( arg1: UIntSparseIntVect, arg2: int, arg3: int) -> None
            Set the value at a specified location

            C++ signature :
                void __setitem__(RDKit::SparseIntVect<unsigned int> {lvalue},unsigned int,int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object: 
        """
        __sub__( arg1: UIntSparseIntVect, arg2: UIntSparseIntVect) -> object

            C++ signature :
                _object* __sub__(RDKit::SparseIntVect<unsigned int> {lvalue},RDKit::SparseIntVect<unsigned int>)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class ULongSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """
    @staticmethod
    def GetLength( arg1: ULongSparseIntVect) -> int: 
        """
        GetLength( arg1: ULongSparseIntVect) -> int
            Returns the length of the vector

            C++ signature :
                unsigned long GetLength(RDKit::SparseIntVect<unsigned long> {lvalue})
        """
    @staticmethod
    def GetNonzeroElements( arg1: ULongSparseIntVect) -> dict: 
        """
        GetNonzeroElements( arg1: ULongSparseIntVect) -> dict
            returns a dictionary of the nonzero elements

            C++ signature :
                boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<unsigned long> {lvalue})
        """
    @staticmethod
    def GetTotalVal( arg1: ULongSparseIntVect, useAbs: bool = False) -> int: 
        """
        GetTotalVal( arg1: ULongSparseIntVect, useAbs: bool = False) -> int
            Get the sum of the values in the vector, basically L1 norm

            C++ signature :
                int GetTotalVal(RDKit::SparseIntVect<unsigned long> {lvalue} [,bool=False])
        """
    @staticmethod
    def ToBinary( arg1: ULongSparseIntVect) -> object: 
        """
        ToBinary( arg1: ULongSparseIntVect) -> object
            returns a binary (pickle) representation of the vector

            C++ signature :
                boost::python::api::object ToBinary(RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def ToList( arg1: ULongSparseIntVect) -> list: 
        """
        ToList( arg1: ULongSparseIntVect) -> list
            Return the SparseIntVect as a python list

            C++ signature :
                boost::python::list ToList(RDKit::SparseIntVect<unsigned long> {lvalue})
        """
    @staticmethod
    def UpdateFromSequence( arg1: ULongSparseIntVect, arg2: AtomPairsParameters) -> None: 
        """
        UpdateFromSequence( arg1: ULongSparseIntVect, arg2: AtomPairsParameters) -> None
            update the vector based on the values in the list or tuple

            C++ signature :
                void UpdateFromSequence(RDKit::SparseIntVect<unsigned long> {lvalue},boost::python::api::object {lvalue})
        """
    @staticmethod
    def __add__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object: 
        """
        __add__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __add__(RDKit::SparseIntVect<unsigned long> {lvalue},RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def __and__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object: 
        """
        __and__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __and__(RDKit::SparseIntVect<unsigned long> {lvalue},RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def __eq__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object: 
        """
        __eq__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __eq__(RDKit::SparseIntVect<unsigned long> {lvalue},RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def __getinitargs__( arg1: ULongSparseIntVect) -> tuple: 
        """
        __getinitargs__( arg1: ULongSparseIntVect) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def __getitem__( arg1: ULongSparseIntVect, arg2: int) -> int: 
        """
        __getitem__( arg1: ULongSparseIntVect, arg2: int) -> int
            Get the value at a specified location

            C++ signature :
                int __getitem__(RDKit::SparseIntVect<unsigned long> {lvalue},unsigned long)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: ULongSparseIntVect) -> object: 
        """
        __iadd__( arg1: object, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<unsigned long>&>,RDKit::SparseIntVect<unsigned long>)

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDKit::SparseIntVect<unsigned long>&>,int)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __idiv__( arg1: object, arg2: int) -> object: 
        """
        __idiv__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __idiv__(boost::python::back_reference<RDKit::SparseIntVect<unsigned long>&>,int)
        """
    @staticmethod
    def __imul__( arg1: object, arg2: int) -> object: 
        """
        __imul__( arg1: object, arg2: int) -> object

            C++ signature :
                _object* __imul__(boost::python::back_reference<RDKit::SparseIntVect<unsigned long>&>,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None
            Constructor

            C++ signature :
                void __init__(_object*,unsigned long)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: ULongSparseIntVect) -> object: 
        """
        __isub__( arg1: object, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<unsigned long>&>,RDKit::SparseIntVect<unsigned long>)

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDKit::SparseIntVect<unsigned long>&>,int)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: int) -> object: ...
    @staticmethod
    def __ne__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object: 
        """
        __ne__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __ne__(RDKit::SparseIntVect<unsigned long> {lvalue},RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def __or__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object: 
        """
        __or__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __or__(RDKit::SparseIntVect<unsigned long> {lvalue},RDKit::SparseIntVect<unsigned long>)
        """
    @staticmethod
    def __setitem__( arg1: ULongSparseIntVect, arg2: int, arg3: int) -> None: 
        """
        __setitem__( arg1: ULongSparseIntVect, arg2: int, arg3: int) -> None
            Set the value at a specified location

            C++ signature :
                void __setitem__(RDKit::SparseIntVect<unsigned long> {lvalue},unsigned long,int)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object: 
        """
        __sub__( arg1: ULongSparseIntVect, arg2: ULongSparseIntVect) -> object

            C++ signature :
                _object* __sub__(RDKit::SparseIntVect<unsigned long> {lvalue},RDKit::SparseIntVect<unsigned long>)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
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
