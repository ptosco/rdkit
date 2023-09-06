"""Module containing functions for working with molecular abbreviations"""
from __future__ import annotations
import rdkit.Chem.rdAbbreviations
import typing
import Boost.Python

__all__ = [
    "AbbreviationDefinition",
    "CondenseAbbreviationSubstanceGroups",
    "CondenseMolAbbreviations",
    "GetDefaultAbbreviations",
    "GetDefaultLinkers",
    "LabelMolAbbreviations",
    "ParseAbbreviations",
    "ParseLinkers"
]


class AbbreviationDefinition(Boost.Python.instance):
    """
    Abbreviation Definition
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def displayLabel(self) -> None:
        """
        the label in a drawing when the bond comes from the right

        :type: None
        """
    @property
    def displayLabelW(self) -> None:
        """
        the label in a drawing when the bond comes from the west

        :type: None
        """
    @property
    def label(self) -> None:
        """
        the label

        :type: None
        """
    @property
    def mol(self) -> None:
        """
        the query molecule (should have a dummy as the first atom)

        :type: None
        """
    __instance_size__ = 192
    pass
class _vectN5RDKit13Abbreviations22AbbreviationDefinitionE(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: object) -> bool: 
        """
        __contains__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: object) -> None: 
        """
        __delitem__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> >&>,_object*)
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<RDKit::Abbreviations::AbbreviationDefinition*, std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > > > __iter__(boost::python::back_reference<std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> >&>)
        """
    @staticmethod
    def __len__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE) -> int: 
        """
        __len__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE) -> int

            C++ signature :
                unsigned long __len__(std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: _vectN5RDKit13Abbreviations22AbbreviationDefinitionE, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
def CondenseAbbreviationSubstanceGroups( mol: Mol) -> Mol:
    """
    CondenseAbbreviationSubstanceGroups( mol: Mol) -> Mol
        Finds and replaces abbreviation (i.e. "SUP") substance groups in a molecule. The result is not sanitized.

        C++ signature :
            RDKit::ROMol* CondenseAbbreviationSubstanceGroups(RDKit::ROMol const*)
    """
def CondenseMolAbbreviations( mol: Mol, abbrevs: AtomPairsParameters, maxCoverage: float = 0.4, sanitize: bool = True) -> Mol:
    """
    CondenseMolAbbreviations( mol: Mol, abbrevs: AtomPairsParameters, maxCoverage: float = 0.4, sanitize: bool = True) -> Mol
        Finds and replaces abbreviations in a molecule. The result is not sanitized.

        C++ signature :
            RDKit::ROMol* CondenseMolAbbreviations(RDKit::ROMol const*,boost::python::api::object [,double=0.4 [,bool=True]])
    """
def GetDefaultAbbreviations() -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    GetDefaultAbbreviations() -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE
        returns a list of the default abbreviation definitions

        C++ signature :
            std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > GetDefaultAbbreviations()
    """
def GetDefaultLinkers() -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    GetDefaultLinkers() -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE
        returns a list of the default linker definitions

        C++ signature :
            std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > GetDefaultLinkers()
    """
def LabelMolAbbreviations( mol: Mol, abbrevs: AtomPairsParameters, maxCoverage: float = 0.4) -> Mol:
    """
    LabelMolAbbreviations( mol: Mol, abbrevs: AtomPairsParameters, maxCoverage: float = 0.4) -> Mol
        Finds abbreviations and adds to them to a molecule as "SUP" SubstanceGroups

        C++ signature :
            RDKit::ROMol* LabelMolAbbreviations(RDKit::ROMol const*,boost::python::api::object [,double=0.4])
    """
def ParseAbbreviations( text: str, removeExtraDummies: bool = False, allowConnectionToDummies: bool = False) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    ParseAbbreviations( text: str, removeExtraDummies: bool = False, allowConnectionToDummies: bool = False) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE
        returns a set of abbreviation definitions from a string

        C++ signature :
            std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > ParseAbbreviations(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False [,bool=False]])
    """
def ParseLinkers( text: str) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    ParseLinkers( text: str) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE
        returns a set of linker definitions from a string

        C++ signature :
            std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > ParseLinkers(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
