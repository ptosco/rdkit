"""Module containing functions for interchange of molecules.
Note that this should be considered beta and that the format
  and API will very likely change in future releases."""
from __future__ import annotations
import rdkit.Chem.rdMolInterchange
import typing
import Boost.Python

__all__ = [
    "JSONParseParameters",
    "JSONToMols",
    "JSONWriteParameters",
    "MolToJSON",
    "MolsToJSON"
]


class JSONParseParameters(Boost.Python.instance):
    """
    Parameters controlling the JSON parser
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def parseConformers(self) -> None:
        """
        parse conformers in the JSON

        :type: None
        """
    @property
    def parseProperties(self) -> None:
        """
        parse molecular properties in the JSON

        :type: None
        """
    @property
    def setAromaticBonds(self) -> None:
        """
        set bond types to aromatic for bonds flagged aromatic

        :type: None
        """
    @property
    def strictValenceCheck(self) -> None:
        """
        be strict when checking atom valences

        :type: None
        """
    @property
    def useHCounts(self) -> None:
        """
        use atomic H counts from the JSON. You may want to set this to False when parsing queries.

        :type: None
        """
    __instance_size__ = 32
    pass
class JSONWriteParameters(Boost.Python.instance):
    """
    Parameters controlling the JSON writer
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def useRDKitExtensions(self) -> None:
        """
        use RDKit extensions to the commonchem format

        :type: None
        """
    __instance_size__ = 32
    pass
def JSONToMols( jsonBlock: str, params: object = None) -> tuple:
    """
    JSONToMols( jsonBlock: str, params: object = None) -> tuple
        Convert JSON to a tuple of molecules
        
            ARGUMENTS:
              - jsonBlock: the molecule to work with
              - params: (optional) JSONParseParameters controlling the JSON parsing
            RETURNS:
              a tuple of Mols
        

        C++ signature :
            boost::python::tuple JSONToMols(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,boost::python::api::object=None])
    """
def MolToJSON( mol: Mol, params: object = None) -> str:
    """
    MolToJSON( mol: Mol, params: object = None) -> str
        Convert a single molecule to JSON
        
            ARGUMENTS:
              - mol: the molecule to work with
            RETURNS:
              a string
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToJSON(RDKit::ROMol [,boost::python::api::object=None])
    """
def MolsToJSON( mols: object, params: object = None) -> str:
    """
    MolsToJSON( mols: object, params: object = None) -> str
        Convert a set of molecules to JSON
        
            ARGUMENTS:
              - mols: the molecules to work with
            RETURNS:
              a string
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolsToJSON(boost::python::api::object [,boost::python::api::object=None])
    """
