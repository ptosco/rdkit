from __future__ import annotations
import rdkit.Chem.rdTautomerQuery
import typing
import Boost.Python

__all__ = [
    "PatternFingerprintTautomerTarget",
    "TautomerQuery",
    "TautomerQueryCanSerialize"
]


class TautomerQuery(Boost.Python.instance):
    """
    The Tautomer Query Class.
      Creates a query that enables structure search accounting for matching of
      Tautomeric forms
    """
    @staticmethod
    def GetModifiedAtoms( arg1: TautomerQuery) -> VectSizeT: 
        """
        GetModifiedAtoms( arg1: TautomerQuery) -> VectSizeT

            C++ signature :
                std::vector<unsigned long, std::allocator<unsigned long> > GetModifiedAtoms(RDKit::TautomerQuery {lvalue})
        """
    @staticmethod
    def GetModifiedBonds( arg1: TautomerQuery) -> VectSizeT: 
        """
        GetModifiedBonds( arg1: TautomerQuery) -> VectSizeT

            C++ signature :
                std::vector<unsigned long, std::allocator<unsigned long> > GetModifiedBonds(RDKit::TautomerQuery {lvalue})
        """
    @typing.overload
    def GetSubstructMatch(self, target: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object: 
        """
        GetSubstructMatch( self: TautomerQuery, target: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object

            C++ signature :
                _object* GetSubstructMatch(RDKit::TautomerQuery,RDKit::ROMol [,bool=False [,bool=False]])

            C++ signature :
                _object* GetSubstructMatch(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatch(self, target: Mol, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object: 
        """
        GetSubstructMatches( self: TautomerQuery, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object

            C++ signature :
                _object* GetSubstructMatches(RDKit::TautomerQuery,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

            C++ signature :
                _object* GetSubstructMatches(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatches(self, target: Mol, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatchesWithTautomers(self, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object: 
        """
        GetSubstructMatchesWithTautomers( self: TautomerQuery, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object

            C++ signature :
                _object* GetSubstructMatchesWithTautomers(RDKit::TautomerQuery,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

            C++ signature :
                _object* GetSubstructMatchesWithTautomers(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatchesWithTautomers(self, target: Mol, params: SubstructMatchParameters) -> object: ...
    @staticmethod
    def GetTautomers( arg1: TautomerQuery) -> object: 
        """
        GetTautomers( arg1: TautomerQuery) -> object

            C++ signature :
                _object* GetTautomers(RDKit::TautomerQuery)
        """
    @staticmethod
    def GetTemplateMolecule( arg1: TautomerQuery) -> Mol: 
        """
        GetTemplateMolecule( arg1: TautomerQuery) -> Mol

            C++ signature :
                RDKit::ROMol GetTemplateMolecule(RDKit::TautomerQuery {lvalue})
        """
    @typing.overload
    def IsSubstructOf(self, target: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool: 
        """
        IsSubstructOf( self: TautomerQuery, target: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool

            C++ signature :
                bool IsSubstructOf(RDKit::TautomerQuery,RDKit::ROMol [,bool=True [,bool=False [,bool=False]]])

            C++ signature :
                bool IsSubstructOf(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def IsSubstructOf(self, target: Mol, params: SubstructMatchParameters) -> bool: ...
    @staticmethod
    def PatternFingerprintTemplate( arg1: TautomerQuery, fingerprintSize: int = 2048) -> ExplicitBitVect: 
        """
        PatternFingerprintTemplate( arg1: TautomerQuery, fingerprintSize: int = 2048) -> ExplicitBitVect

            C++ signature :
                ExplicitBitVect* PatternFingerprintTemplate(RDKit::TautomerQuery {lvalue} [,unsigned int=2048])
        """
    @staticmethod
    def ToBinary( arg1: TautomerQuery) -> object: 
        """
        ToBinary( arg1: TautomerQuery) -> object

            C++ signature :
                boost::python::api::object ToBinary(RDKit::TautomerQuery)
        """
    @staticmethod
    def __getinitargs__( arg1: TautomerQuery) -> tuple: 
        """
        __getinitargs__( arg1: TautomerQuery) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::TautomerQuery)
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
    def __init__( arg1: AtomPairsParameters, arg2: Mol) -> object: 
        """
        __init__( arg1: AtomPairsParameters, arg2: Mol) -> object

            C++ signature :
                void* __init__(boost::python::api::object,RDKit::ROMol)

            C++ signature :
                void* __init__(boost::python::api::object,RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, arg2: Mol, arg3: str) -> object: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __safe_for_unpickling__ = True
    pass
def PatternFingerprintTautomerTarget( target: Mol, fingerprintSize: int = 2048) -> ExplicitBitVect:
    """
    PatternFingerprintTautomerTarget( target: Mol, fingerprintSize: int = 2048) -> ExplicitBitVect

        C++ signature :
            ExplicitBitVect* PatternFingerprintTautomerTarget(RDKit::ROMol [,unsigned int=2048])
    """
def TautomerQueryCanSerialize() -> bool:
    """
    TautomerQueryCanSerialize() -> bool
        Returns True if the TautomerQuery is serializable (requires that the RDKit was built with boost::serialization)

        C++ signature :
            bool TautomerQueryCanSerialize()
    """
