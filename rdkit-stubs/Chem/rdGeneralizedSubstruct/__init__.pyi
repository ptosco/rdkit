"""Module containing functions for generalized substructure searching"""
from __future__ import annotations
import rdkit.Chem.rdGeneralizedSubstruct
import typing
import Boost.Python

__all__ = [
    "CreateExtendedQueryMol",
    "ExtendedQueryMol",
    "MolGetSubstructMatch",
    "MolGetSubstructMatches",
    "MolHasSubstructMatch"
]


class ExtendedQueryMol(Boost.Python.instance):
    """
    Extended query molecule for use in generalized substructure searching.
    """
    @staticmethod
    def InitFromBinary( arg1: ExtendedQueryMol, arg2: str) -> None: 
        """
        InitFromBinary( arg1: ExtendedQueryMol, arg2: str) -> None

            C++ signature :
                void InitFromBinary(RDKit::GeneralizedSubstruct::ExtendedQueryMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def InitFromJSON( arg1: ExtendedQueryMol, arg2: str) -> None: 
        """
        InitFromJSON( arg1: ExtendedQueryMol, arg2: str) -> None

            C++ signature :
                void InitFromJSON(RDKit::GeneralizedSubstruct::ExtendedQueryMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def ToBinary( arg1: ExtendedQueryMol) -> object: 
        """
        ToBinary( arg1: ExtendedQueryMol) -> object

            C++ signature :
                boost::python::api::object ToBinary(RDKit::GeneralizedSubstruct::ExtendedQueryMol)
        """
    @staticmethod
    def ToJSON( arg1: ExtendedQueryMol) -> str: 
        """
        ToJSON( arg1: ExtendedQueryMol) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ToJSON(RDKit::GeneralizedSubstruct::ExtendedQueryMol {lvalue})
        """
    @staticmethod
    def __init__( arg1: object, text: str, isJSON: bool = False) -> None: 
        """
        __init__( arg1: object, text: str, isJSON: bool = False) -> None
            constructor from either a binary string (from ToBinary()) or a JSON string.

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    __instance_size__ = 40
    pass
def CreateExtendedQueryMol( mol: Mol, doEnumeration: bool = True, doTautomers: bool = True, adjustQueryProperties: bool = False, adjustQueryParameters: AdjustQueryParameters = None) -> ExtendedQueryMol:
    """
    CreateExtendedQueryMol( mol: Mol, doEnumeration: bool = True, doTautomers: bool = True, adjustQueryProperties: bool = False, adjustQueryParameters: AdjustQueryParameters = None) -> ExtendedQueryMol
        Creates an ExtendedQueryMol from the input molecule
        
          This takes a query molecule and, conceptually, performs the following steps to
          produce an ExtendedQueryMol:
        
            1. Enumerates features like Link Nodes and SRUs
            2. Converts everything into TautomerQueries
            3. Runs adjustQueryProperties()
        
          Each step is optional
        

        C++ signature :
            RDKit::GeneralizedSubstruct::ExtendedQueryMol* CreateExtendedQueryMol(RDKit::ROMol [,bool=True [,bool=True [,bool=False [,RDKit::MolOps::AdjustQueryParameters*=None]]]])
    """
def MolGetSubstructMatch( mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> object:
    """
    MolGetSubstructMatch( mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> object
        returns first match (if any) of a molecule to a generalized substructure query

        C++ signature :
            _object* MolGetSubstructMatch(RDKit::ROMol,RDKit::GeneralizedSubstruct::ExtendedQueryMol [,RDKit::SubstructMatchParameters*=None])
    """
def MolGetSubstructMatches( mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> object:
    """
    MolGetSubstructMatches( mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> object
        returns all matches (if any) of a molecule to a generalized substructure query

        C++ signature :
            _object* MolGetSubstructMatches(RDKit::ROMol,RDKit::GeneralizedSubstruct::ExtendedQueryMol [,RDKit::SubstructMatchParameters*=None])
    """
def MolHasSubstructMatch( mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> bool:
    """
    MolHasSubstructMatch( mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> bool
        determines whether or not a molecule is a match to a generalized substructure query

        C++ signature :
            bool MolHasSubstructMatch(RDKit::ROMol,RDKit::GeneralizedSubstruct::ExtendedQueryMol [,RDKit::SubstructMatchParameters*=None])
    """
