from __future__ import annotations
import rdkit.Chem.rdMolCatalog
import typing
import Boost.Python

__all__ = [
    "CreateMolCatalog",
    "MolCatalog",
    "MolCatalogEntry"
]


class MolCatalog(Boost.Python.instance):
    @staticmethod
    def AddEdge( arg1: MolCatalog, arg2: int, arg3: int) -> None: 
        """
        AddEdge( arg1: MolCatalog, arg2: int, arg3: int) -> None

            C++ signature :
                void AddEdge(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def AddEntry( arg1: MolCatalog, arg2: MolCatalogEntry) -> int: 
        """
        AddEntry( arg1: MolCatalog, arg2: MolCatalogEntry) -> int

            C++ signature :
                unsigned int AddEntry(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int>*,RDKit::MolCatalogEntry*)
        """
    @staticmethod
    def GetBitDescription( arg1: MolCatalog, arg2: int) -> str: 
        """
        GetBitDescription( arg1: MolCatalog, arg2: int) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetBitDescription(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetBitEntryId( arg1: MolCatalog, arg2: int) -> int: 
        """
        GetBitEntryId( arg1: MolCatalog, arg2: int) -> int

            C++ signature :
                unsigned int GetBitEntryId(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryBitId( arg1: MolCatalog, arg2: int) -> int: 
        """
        GetEntryBitId( arg1: MolCatalog, arg2: int) -> int

            C++ signature :
                unsigned int GetEntryBitId(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryDescription( arg1: MolCatalog, arg2: int) -> str: 
        """
        GetEntryDescription( arg1: MolCatalog, arg2: int) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetEntryDescription(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetEntryDownIds( arg1: MolCatalog, arg2: int) -> _vecti: 
        """
        GetEntryDownIds( arg1: MolCatalog, arg2: int) -> _vecti

            C++ signature :
                std::vector<int, std::allocator<int> > GetEntryDownIds(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
    @staticmethod
    def GetFPLength( arg1: MolCatalog) -> int: 
        """
        GetFPLength( arg1: MolCatalog) -> int

            C++ signature :
                unsigned int GetFPLength(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue})
        """
    @staticmethod
    def GetNumEntries( arg1: MolCatalog) -> int: 
        """
        GetNumEntries( arg1: MolCatalog) -> int

            C++ signature :
                unsigned int GetNumEntries(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue})
        """
    @staticmethod
    def Serialize( arg1: MolCatalog) -> str: 
        """
        Serialize( arg1: MolCatalog) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Serialize(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue})
        """
    @staticmethod
    def __getinitargs__( arg1: MolCatalog) -> tuple: 
        """
        __getinitargs__( arg1: MolCatalog) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int>)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    def __init__( arg1: object, arg2: str) -> None: 
        """
        __init__( arg1: object, arg2: str) -> None

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
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
class MolCatalogEntry(Boost.Python.instance):
    @staticmethod
    def GetDescription( arg1: MolCatalogEntry) -> str: 
        """
        GetDescription( arg1: MolCatalogEntry) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetDescription(RDKit::MolCatalogEntry {lvalue})
        """
    @staticmethod
    def GetMol( arg1: MolCatalogEntry) -> Mol: 
        """
        GetMol( arg1: MolCatalogEntry) -> Mol

            C++ signature :
                RDKit::ROMol GetMol(RDKit::MolCatalogEntry {lvalue})
        """
    @staticmethod
    def GetOrder( arg1: MolCatalogEntry) -> int: 
        """
        GetOrder( arg1: MolCatalogEntry) -> int

            C++ signature :
                unsigned int GetOrder(RDKit::MolCatalogEntry {lvalue})
        """
    @staticmethod
    def SetDescription( arg1: MolCatalogEntry, arg2: str) -> None: 
        """
        SetDescription( arg1: MolCatalogEntry, arg2: str) -> None

            C++ signature :
                void SetDescription(RDKit::MolCatalogEntry {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetMol( arg1: MolCatalogEntry, arg2: Mol) -> None: 
        """
        SetMol( arg1: MolCatalogEntry, arg2: Mol) -> None

            C++ signature :
                void SetMol(RDKit::MolCatalogEntry*,RDKit::ROMol const*)
        """
    @staticmethod
    def SetOrder( arg1: MolCatalogEntry, arg2: int) -> None: 
        """
        SetOrder( arg1: MolCatalogEntry, arg2: int) -> None

            C++ signature :
                void SetOrder(RDKit::MolCatalogEntry {lvalue},unsigned int)
        """
    @staticmethod
    def __getinitargs__( arg1: MolCatalogEntry) -> tuple: 
        """
        __getinitargs__( arg1: MolCatalogEntry) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::MolCatalogEntry)
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
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

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
    __instance_size__ = 96
    __safe_for_unpickling__ = True
    pass
def CreateMolCatalog() -> MolCatalog:
    """
    CreateMolCatalog() -> MolCatalog

        C++ signature :
            RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int>* CreateMolCatalog()
    """
