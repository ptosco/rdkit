from __future__ import annotations
import rdkit.Chem.rdfragcatalog
import typing
import Boost.Python

__all__ = [
    "FragCatGenerator",
    "FragCatParams",
    "FragCatalog",
    "FragFPGenerator"
]


class FragCatGenerator(Boost.Python.instance):
    @staticmethod
    def AddFragsFromMol( arg1: FragCatGenerator, arg2: Mol, arg3: FragCatalog) -> int: 
        """
        AddFragsFromMol( arg1: FragCatGenerator, arg2: Mol, arg3: FragCatalog) -> int

            C++ signature :
                unsigned int AddFragsFromMol(RDKit::FragCatGenerator {lvalue},RDKit::ROMol,RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int>*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 32
    pass
class FragCatParams(Boost.Python.instance):
    @staticmethod
    def GetFuncGroup( arg1: FragCatParams, arg2: int) -> Mol: 
        """
        GetFuncGroup( arg1: FragCatParams, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol const* GetFuncGroup(RDKit::FragCatParams {lvalue},int)
        """
    @staticmethod
    def GetLowerFragLength( arg1: FragCatParams) -> int: 
        """
        GetLowerFragLength( arg1: FragCatParams) -> int

            C++ signature :
                unsigned int GetLowerFragLength(RDKit::FragCatParams {lvalue})
        """
    @staticmethod
    def GetNumFuncGroups( arg1: FragCatParams) -> int: 
        """
        GetNumFuncGroups( arg1: FragCatParams) -> int

            C++ signature :
                unsigned int GetNumFuncGroups(RDKit::FragCatParams {lvalue})
        """
    @staticmethod
    def GetTolerance( arg1: FragCatParams) -> float: 
        """
        GetTolerance( arg1: FragCatParams) -> float

            C++ signature :
                double GetTolerance(RDKit::FragCatParams {lvalue})
        """
    @staticmethod
    def GetTypeString( arg1: FragCatParams) -> str: 
        """
        GetTypeString( arg1: FragCatParams) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetTypeString(RDKit::FragCatParams {lvalue})
        """
    @staticmethod
    def GetUpperFragLength( arg1: FragCatParams) -> int: 
        """
        GetUpperFragLength( arg1: FragCatParams) -> int

            C++ signature :
                unsigned int GetUpperFragLength(RDKit::FragCatParams {lvalue})
        """
    @staticmethod
    def Serialize( arg1: FragCatParams) -> str: 
        """
        Serialize( arg1: FragCatParams) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Serialize(RDKit::FragCatParams {lvalue})
        """
    @staticmethod
    def __init__( arg1: object, lLen: int, uLen: int, fgroupFilename: str, tol: float = 1e-08) -> None: 
        """
        __init__( arg1: object, lLen: int, uLen: int, fgroupFilename: str, tol: float = 1e-08) -> None

            C++ signature :
                void __init__(_object*,int,int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,double=1e-08])
        """
    __instance_size__ = 104
    pass
class FragCatalog(Boost.Python.instance):
    @staticmethod
    def GetBitDescription( arg1: FragCatalog, arg2: int) -> str: 
        """
        GetBitDescription( arg1: FragCatalog, arg2: int) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetBitDescription(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetBitDiscrims( arg1: FragCatalog, arg2: int) -> _vectd: 
        """
        GetBitDiscrims( arg1: FragCatalog, arg2: int) -> _vectd

            C++ signature :
                std::vector<double, std::allocator<double> > GetBitDiscrims(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetBitEntryId( arg1: FragCatalog, arg2: int) -> int: 
        """
        GetBitEntryId( arg1: FragCatalog, arg2: int) -> int

            C++ signature :
                unsigned int GetBitEntryId(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetBitFuncGroupIds( arg1: FragCatalog, arg2: int) -> _vecti: 
        """
        GetBitFuncGroupIds( arg1: FragCatalog, arg2: int) -> _vecti

            C++ signature :
                std::vector<int, std::allocator<int> > GetBitFuncGroupIds(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetBitOrder( arg1: FragCatalog, arg2: int) -> int: 
        """
        GetBitOrder( arg1: FragCatalog, arg2: int) -> int

            C++ signature :
                unsigned int GetBitOrder(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetCatalogParams( arg1: FragCatalog) -> FragCatParams: 
        """
        GetCatalogParams( arg1: FragCatalog) -> FragCatParams

            C++ signature :
                RDKit::FragCatParams* GetCatalogParams(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
    @staticmethod
    def GetEntryBitId( arg1: FragCatalog, arg2: int) -> int: 
        """
        GetEntryBitId( arg1: FragCatalog, arg2: int) -> int

            C++ signature :
                unsigned int GetEntryBitId(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryDescription( arg1: FragCatalog, arg2: int) -> str: 
        """
        GetEntryDescription( arg1: FragCatalog, arg2: int) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetEntryDescription(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryDownIds( arg1: FragCatalog, arg2: int) -> _vecti: 
        """
        GetEntryDownIds( arg1: FragCatalog, arg2: int) -> _vecti

            C++ signature :
                std::vector<int, std::allocator<int> > GetEntryDownIds(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryFuncGroupIds( arg1: FragCatalog, arg2: int) -> _vecti: 
        """
        GetEntryFuncGroupIds( arg1: FragCatalog, arg2: int) -> _vecti

            C++ signature :
                std::vector<int, std::allocator<int> > GetEntryFuncGroupIds(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryOrder( arg1: FragCatalog, arg2: int) -> int: 
        """
        GetEntryOrder( arg1: FragCatalog, arg2: int) -> int

            C++ signature :
                unsigned int GetEntryOrder(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetFPLength( arg1: FragCatalog) -> int: 
        """
        GetFPLength( arg1: FragCatalog) -> int

            C++ signature :
                unsigned int GetFPLength(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
    @staticmethod
    def GetNumEntries( arg1: FragCatalog) -> int: 
        """
        GetNumEntries( arg1: FragCatalog) -> int

            C++ signature :
                unsigned int GetNumEntries(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
    @staticmethod
    def Serialize( arg1: FragCatalog) -> str: 
        """
        Serialize( arg1: FragCatalog) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Serialize(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
    @staticmethod
    def __getinitargs__( arg1: FragCatalog) -> tuple: 
        """
        __getinitargs__( arg1: FragCatalog) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int>)
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
    def __init__( arg1: object, arg2: FragCatParams) -> None: 
        """
        __init__( arg1: object, arg2: FragCatParams) -> None

            C++ signature :
                void __init__(_object*,RDKit::FragCatParams*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
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
    __instance_size__ = 152
    __safe_for_unpickling__ = True
    pass
class FragFPGenerator(Boost.Python.instance):
    @staticmethod
    def GetFPForMol( arg1: FragFPGenerator, arg2: Mol, arg3: FragCatalog) -> ExplicitBitVect: 
        """
        GetFPForMol( arg1: FragFPGenerator, arg2: Mol, arg3: FragCatalog) -> ExplicitBitVect

            C++ signature :
                ExplicitBitVect* GetFPForMol(RDKit::FragFPGenerator {lvalue},RDKit::ROMol,RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int>)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 32
    pass
