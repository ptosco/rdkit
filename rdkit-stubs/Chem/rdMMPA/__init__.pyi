"""Module containing a C++ implementation of code for doing MMPA"""
from __future__ import annotations
import rdkit.Chem.rdMMPA
import typing

__all__ = [
    "FragmentMol"
]


def FragmentMol( mol: Mol, bondsToCut: AtomPairsParameters, minCuts: int = 1, maxCuts: int = 3, resultsAsMols: bool = True) -> tuple:
    """
    FragmentMol( mol: Mol, maxCuts: int = 3, maxCutBonds: int = 20, __CLOSE_SQUARE_BRACKET_TAG__!@!__EQUALS_TAG__!#__OPEN_SQUARE_BRACKET_TAG__*__CLOSE_SQUARE_BRACKET_TAG__': str)pattern='__OPEN_SQUARE_BRACKET_TAG__#6+0;!$(*__EQUALS_TAG__,#__OPEN_SQUARE_BRACKET_TAG__!#6__CLOSE_SQUARE_BRACKET_TAG__, resultsAsMols: bool = True) -> tuple
        Does the fragmentation necessary for an MMPA analysis

        C++ signature :
            boost::python::tuple FragmentMol(RDKit::ROMol [,unsigned int=3 [,unsigned int=20 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='[#6+0;!$(*=,#[!#6])]!@!=!#[*]' [,bool=True]]]])

        C++ signature :
            boost::python::tuple FragmentMol(RDKit::ROMol,unsigned int,unsigned int,unsigned int [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='[#6+0;!$(*=,#[!#6])]!@!=!#[*]' [,bool=True]])

        C++ signature :
            boost::python::tuple FragmentMol(RDKit::ROMol,boost::python::api::object [,unsigned int=1 [,unsigned int=3 [,bool=True]]])
    """
