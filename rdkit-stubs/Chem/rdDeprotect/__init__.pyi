from __future__ import annotations
import rdkit.Chem.rdDeprotect
import typing
import Boost.Python

__all__ = [
    "Deprotect",
    "DeprotectData",
    "DeprotectDataVect",
    "DeprotectInPlace",
    "GetDeprotections"
]


class DeprotectData(Boost.Python.instance):
    """
    DeprotectData class, contains a single deprotection reaction and information

     deprotectdata.deprotection_class - functional group being protected
     deprotectdata.reaction_smarts - reaction smarts used for deprotection
     deprotectdata.abbreviation - common abbreviation for the protecting group
     deprotectdata.full_name - full name for the protecting group
    """
    @staticmethod
    def __init__( arg1: object, deprotection_class: str, reaction_smarts: str, abbreviation: str, full_name: str) -> None: 
        """
        __init__( arg1: object, deprotection_class: str, reaction_smarts: str, abbreviation: str, full_name: str) -> None
            Construct a new DeprotectData instance.
              >>> reaction_class = "amine"
              >>> reaction_smarts = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]"
              >>> abbreviation = "Boc"
              >>> full_name = "tert-butyloxycarbonyl"
              >>> data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
              >>> assert data.isValid()
            
            

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def isValid( arg1: DeprotectData) -> bool: 
        """
        isValid( arg1: DeprotectData) -> bool
            Returns True if the DeprotectData has a valid reaction

            C++ signature :
                bool isValid(RDKit::Deprotect::DeprotectData {lvalue})
        """
    @property
    def abbreviation(self) -> None:
        """
        :type: None
        """
    @property
    def deprotection_class(self) -> None:
        """
        :type: None
        """
    @property
    def example(self) -> None:
        """
        :type: None
        """
    @property
    def full_name(self) -> None:
        """
        :type: None
        """
    @property
    def reaction_smarts(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 200
    pass
class DeprotectDataVect(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: DeprotectDataVect, arg2: object) -> bool: 
        """
        __contains__( arg1: DeprotectDataVect, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: DeprotectDataVect, arg2: object) -> None: 
        """
        __delitem__( arg1: DeprotectDataVect, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<RDKit::Deprotect::DeprotectData*, std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > > > __iter__(boost::python::back_reference<std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> >&>)
        """
    @staticmethod
    def __len__( arg1: DeprotectDataVect) -> int: 
        """
        __len__( arg1: DeprotectDataVect) -> int

            C++ signature :
                unsigned long __len__(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: DeprotectDataVect, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: DeprotectDataVect, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: DeprotectDataVect, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: DeprotectDataVect, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: DeprotectDataVect, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: DeprotectDataVect, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
def Deprotect( mol: Mol, deprotections: AtomPairsParameters = None) -> Mol:
    """
    Deprotect( mol: Mol, deprotections: AtomPairsParameters = None) -> Mol
        Return the deprotected version of the molecule.

        C++ signature :
            boost::shared_ptr<RDKit::ROMol> Deprotect(RDKit::ROMol [,boost::python::api::object=None])
    """
def DeprotectInPlace( mol: Mol, deprotections: AtomPairsParameters = None) -> bool:
    """
    DeprotectInPlace( mol: Mol, deprotections: AtomPairsParameters = None) -> bool
        Deprotects the molecule in place.

        C++ signature :
            bool DeprotectInPlace(RDKit::ROMol {lvalue} [,boost::python::api::object=None])
    """
def GetDeprotections() -> DeprotectDataVect:
    """
    GetDeprotections() -> DeprotectDataVect
        Return the default list of deprotections

        C++ signature :
            std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > GetDeprotections()
    """
